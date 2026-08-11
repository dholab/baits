#!/usr/bin/env python3
"""Calibrate a Deacon absolute threshold from classified candidate reads."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, cast

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence

CLASSIFIED_FIELDS = (
    "metagenome_id", "fragment_id", "mate", "read_length", "distinct_bait_count",
    "fragment_distinct_bait_count", "representative_id", "classification",
    "best_target_bit_score", "best_non_target_bit_score",
)
PREPARATION_METRICS = (
    "design_id", "deacon_returned_read_count", "duplicate_fragment_count",
    "zero_bait_read_count", "candidate_read_count", "whole_read_blast_query_count",
)
SUMMARY_METRICS = (
    *PREPARATION_METRICS, "target_classified_read_count", "non_target_classified_read_count",
    "tied_read_count", "no_hit_read_count", "calibration_status",
    "recommended_deacon_absolute_threshold", "specificity_floor", "conclusion",
)
READ_COUNT_FIELDS = (
    "threshold", "target_read_count", "non_target_read_count", "tied_read_count",
    "no_hit_read_count",
)
CURVE_FIELDS = ("threshold", "retained_metagenome_count", "retained_fragment_count")
CLASSIFIED_SCHEMA = {
    "metagenome_id": pl.String,
    "fragment_id": pl.String,
    "mate": pl.String,
    "read_length": pl.Int64,
    "distinct_bait_count": pl.Int64,
    "fragment_distinct_bait_count": pl.Int64,
    "representative_id": pl.String,
    "classification": pl.String,
    "best_target_bit_score": pl.String,
    "best_non_target_bit_score": pl.String,
}


class DeaconThresholdCalibrationError(ValueError):
    """Raised when threshold calibration evidence is malformed."""

    @classmethod
    def malformed_classified_reads(cls, path: Path) -> DeaconThresholdCalibrationError:
        return cls(f"Classified-read evidence is malformed: {path}")

    @classmethod
    def duplicate_identity(cls) -> DeaconThresholdCalibrationError:
        return cls("Classified-read identity is duplicated")

    @classmethod
    def fragment_disagreement(cls) -> DeaconThresholdCalibrationError:
        return cls("Classified-read fragment-wide bait counts disagree")

    @classmethod
    def no_classified_reads(cls) -> DeaconThresholdCalibrationError:
        return cls("Classified-read evidence must contain at least one candidate read")

    @classmethod
    def malformed_preparation_summary(cls, path: Path) -> DeaconThresholdCalibrationError:
        return cls(f"Candidate-read preparation summary is malformed: {path}")

    @classmethod
    def design_mismatch(cls) -> DeaconThresholdCalibrationError:
        return cls("Candidate-read preparation summary design_id does not match")

    @classmethod
    def candidate_count_mismatch(cls) -> DeaconThresholdCalibrationError:
        return cls("Classified-read count disagrees with candidate-read preparation")


class CalibrationStatus(StrEnum):
    NO_CLASSIFIED_READS = "NO_CLASSIFIED_READS"
    RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD = "RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD"
    SPECIFICITY_FLOOR = "SPECIFICITY_FLOOR"


@dataclass(frozen=True)
class PreparationSummary:
    design_id: str
    deacon_returned_read_count: int
    duplicate_fragment_count: int
    zero_bait_read_count: int
    candidate_read_count: int
    whole_read_blast_query_count: int


@dataclass(frozen=True)
class CalibrationConclusion:
    status: CalibrationStatus
    recommended_threshold: int | None
    specificity_floor: int | None
    conclusion: str


@dataclass(frozen=True)
class ThresholdCalibration:
    read_counts: pl.LazyFrame
    curve: pl.LazyFrame
    summary: pl.LazyFrame
    conclusion: CalibrationConclusion


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calibrate a Deacon absolute threshold.")
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--classified-reads", type=Path, required=True)
    parser.add_argument("--preparation-summary", type=Path, required=True)
    parser.add_argument("--read-counts-out", type=Path, required=True)
    parser.add_argument("--curve-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    return parser.parse_args(argv)


def first_error(errors: pl.LazyFrame) -> str | None:
    error = errors.limit(1).collect()
    return None if error.is_empty() else error.item(0, "error")


def has_expected_header(path: Path, fields: tuple[str, ...]) -> bool:
    try:
        with path.open(newline="") as source:
            return next(csv.reader(source, delimiter="\t")) == list(fields)
    except (OSError, StopIteration, csv.Error):
        return False


def classified_read_errors(reads: pl.LazyFrame) -> pl.LazyFrame:
    valid_row = pl.all_horizontal(
        pl.col("metagenome_id").is_not_null(),
        pl.col("fragment_id").is_not_null(),
        pl.col("representative_id").is_not_null(),
        pl.col("classification").is_not_null(),
        pl.col("read_length").is_not_null(),
        pl.col("distinct_bait_count").is_not_null(),
        pl.col("fragment_distinct_bait_count").is_not_null(),
        pl.col("metagenome_id").str.strip_chars().str.len_chars() > 0,
        pl.col("fragment_id").str.strip_chars().str.len_chars() > 0,
        pl.col("representative_id").str.strip_chars().str.len_chars() > 0,
        pl.col("mate").is_in(("", "1", "2")),
        pl.col("read_length") > 0,
        pl.col("distinct_bait_count") > 0,
        pl.col("fragment_distinct_bait_count") > 0,
        pl.col("fragment_distinct_bait_count") >= pl.col("distinct_bait_count"),
        pl.col("classification").is_in(("TARGET", "NON_TARGET", "TIED", "NO_HIT")),
    )
    return pl.concat(
        [
            reads.filter(~valid_row.fill_null(value=False)).with_columns(
                pl.lit("malformed").alias("error"),
            ).select("error"),
            reads.group_by("metagenome_id", "fragment_id", "mate").len().filter(
                pl.col("len") > 1,
            ).with_columns(pl.lit("duplicate_identity").alias("error")).select("error"),
            reads.group_by("metagenome_id", "fragment_id").agg(
                pl.col("fragment_distinct_bait_count").n_unique().alias("count"),
            ).filter(pl.col("count") > 1).with_columns(
                pl.lit("fragment_disagreement").alias("error"),
            ).select("error"),
        ],
        how="vertical",
    )


def construct_classified_reads(path: Path) -> pl.LazyFrame:
    try:
        if not has_expected_header(path, CLASSIFIED_FIELDS):
            raise DeaconThresholdCalibrationError.malformed_classified_reads(path)
        reads = pl.scan_csv(
            path,
            separator="\t",
            schema=CLASSIFIED_SCHEMA,
            null_values=[],
            raise_if_empty=False,
        ).with_columns(pl.col("mate").fill_null(""))
        error = first_error(classified_read_errors(reads))
        count = reads.select(pl.len()).collect().item()
    except pl.exceptions.PolarsError as exception:
        raise DeaconThresholdCalibrationError.malformed_classified_reads(path) from exception
    if count == 0:
        raise DeaconThresholdCalibrationError.no_classified_reads()
    match error:
        case "duplicate_identity":
            raise DeaconThresholdCalibrationError.duplicate_identity()
        case "fragment_disagreement":
            raise DeaconThresholdCalibrationError.fragment_disagreement()
        case None:
            return reads
        case _:
            raise DeaconThresholdCalibrationError.malformed_classified_reads(path)


def parse_nonnegative_int64(value: str | None, path: Path) -> int:
    if value is None or not value.isascii() or not value.isdecimal():
        raise DeaconThresholdCalibrationError.malformed_preparation_summary(path)
    parsed = int(value)
    if parsed > 2**63 - 1:
        raise DeaconThresholdCalibrationError.malformed_preparation_summary(path)
    return parsed


def construct_preparation_summary(path: Path, design_id: str) -> PreparationSummary:
    try:
        if not has_expected_header(path, ("metric", "value")):
            raise DeaconThresholdCalibrationError.malformed_preparation_summary(path)
        summary = pl.scan_csv(
            path,
            separator="\t",
            schema={"metric": pl.String, "value": pl.String},
            null_values=[],
            raise_if_empty=False,
        ).collect()
    except pl.exceptions.PolarsError as exception:
        raise DeaconThresholdCalibrationError.malformed_preparation_summary(path) from exception
    if tuple(summary.get_column("metric")) != PREPARATION_METRICS:
        raise DeaconThresholdCalibrationError.malformed_preparation_summary(path)
    values = summary.get_column("value").to_list()
    if values[0] is None or values[0].strip() == "":
        raise DeaconThresholdCalibrationError.malformed_preparation_summary(path)
    if values[0] != design_id:
        raise DeaconThresholdCalibrationError.design_mismatch()
    counts = tuple(parse_nonnegative_int64(value, path) for value in values[1:])
    return PreparationSummary(str(values[0]), *counts)


def conclusion_from_counts(read_counts: pl.DataFrame) -> CalibrationConclusion:
    first = read_counts.row(0, named=True)
    if first["target_read_count"] + first["non_target_read_count"] + first["tied_read_count"] == 0:
        return CalibrationConclusion(
            CalibrationStatus.NO_CLASSIFIED_READS,
            None,
            None,
            "Every candidate read is a no-hit read; no threshold is supported.",
        )
    clean = read_counts.filter(
        (pl.col("non_target_read_count") == 0) & (pl.col("tied_read_count") == 0),
    ).row(0, named=True)
    threshold = clean["threshold"]
    if clean["target_read_count"] > 0:
        return CalibrationConclusion(
            CalibrationStatus.RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD,
            threshold,
            None,
            f"The recommended Deacon absolute threshold is {threshold} for this optimization read set.",
        )
    return CalibrationConclusion(
        CalibrationStatus.SPECIFICITY_FLOOR,
        None,
        threshold,
        f"The specificity floor is {threshold}; no target-classified read remains.",
    )


def summary_relation(preparation: PreparationSummary, totals: pl.DataFrame, conclusion: CalibrationConclusion) -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "metric": SUMMARY_METRICS,
            "value": (
                preparation.design_id,
                *(str(value) for value in (
                    preparation.deacon_returned_read_count,
                    preparation.duplicate_fragment_count,
                    preparation.zero_bait_read_count,
                    preparation.candidate_read_count,
                    preparation.whole_read_blast_query_count,
                )),
                *(str(totals.item(0, field)) for field in READ_COUNT_FIELDS[1:]),
                conclusion.status.value,
                "" if conclusion.recommended_threshold is None else str(conclusion.recommended_threshold),
                "" if conclusion.specificity_floor is None else str(conclusion.specificity_floor),
                conclusion.conclusion,
            ),
        },
    ).lazy()


def construct_threshold_calibration(
    classified_reads: pl.LazyFrame,
    preparation: PreparationSummary,
) -> ThresholdCalibration:
    reads = classified_reads.collect()
    if reads.height != preparation.candidate_read_count:
        raise DeaconThresholdCalibrationError.candidate_count_mismatch()
    maximum = cast("int", reads.get_column("fragment_distinct_bait_count").max())
    thresholds = pl.DataFrame({"threshold": range(1, maximum + 2)}).lazy()
    evidence = reads.lazy()
    retained = thresholds.join(evidence, how="cross").filter(
        pl.col("fragment_distinct_bait_count") >= pl.col("threshold"),
    )
    read_counts = thresholds.join(
        retained.group_by("threshold").agg(
            (pl.col("classification") == "TARGET").sum().alias("target_read_count"),
            (pl.col("classification") == "NON_TARGET").sum().alias("non_target_read_count"),
            (pl.col("classification") == "TIED").sum().alias("tied_read_count"),
            (pl.col("classification") == "NO_HIT").sum().alias("no_hit_read_count"),
        ),
        on="threshold",
        how="left",
    ).with_columns(pl.all().exclude("threshold").fill_null(0)).select(READ_COUNT_FIELDS).sort("threshold")
    fragments = evidence.select(
        "metagenome_id", "fragment_id", "fragment_distinct_bait_count",
    ).unique()
    curve = thresholds.join(
        thresholds.join(fragments, how="cross").filter(
            pl.col("fragment_distinct_bait_count") >= pl.col("threshold"),
        ).group_by("threshold").agg(
            pl.col("metagenome_id").n_unique().alias("retained_metagenome_count"),
            pl.len().alias("retained_fragment_count"),
        ),
        on="threshold",
        how="left",
    ).with_columns(pl.all().exclude("threshold").fill_null(0)).select(CURVE_FIELDS).sort("threshold")
    materialized_counts = read_counts.collect()
    materialized_curve = curve.collect()
    conclusion = conclusion_from_counts(materialized_counts)
    summary = summary_relation(preparation, materialized_counts.head(1), conclusion)
    return ThresholdCalibration(materialized_counts.lazy(), materialized_curve.lazy(), summary, conclusion)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    classified_reads = construct_classified_reads(args.classified_reads)
    preparation = construct_preparation_summary(args.preparation_summary, args.design_id)
    calibration = construct_threshold_calibration(classified_reads, preparation)
    calibration.read_counts.collect().write_csv(args.read_counts_out, separator="\t", quote_style="never")
    calibration.curve.collect().write_csv(args.curve_out, separator="\t", quote_style="never")
    calibration.summary.collect().write_csv(args.summary_out, separator="\t", quote_style="never")


if __name__ == "__main__":
    main()
