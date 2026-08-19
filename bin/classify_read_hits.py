#!/usr/bin/env python3
"""Classify candidate reads from read BLAST evidence."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from itertools import groupby
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

CANDIDATE_FIELDS = (
    "metagenome_id", "read_id", "read_length", "bait_count", "representative_id",
)
BLAST_FIELDS = (
    "qseqid", "qlen", "saccver", "staxids", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "stitle",
)
CLASSIFIED_FIELDS = (*CANDIDATE_FIELDS, "classification", "best_target_bit_score", "best_non_target_bit_score")
CANDIDATE_SCHEMA = {
    "metagenome_id": pl.String,
    "read_id": pl.String,
    "read_length": pl.Int64,
    "bait_count": pl.Int64,
    "representative_id": pl.String,
}
BLAST_SCHEMA = {
    "qseqid": pl.String,
    "qlen": pl.Int64,
    "saccver": pl.String,
    "staxids": pl.String,
    "pident": pl.Float64,
    "length": pl.Int64,
    "mismatch": pl.Int64,
    "gapopen": pl.Int64,
    "qstart": pl.Int64,
    "qend": pl.Int64,
    "sstart": pl.Int64,
    "send": pl.Int64,
    "evalue": pl.Float64,
    "bitscore": pl.String,
    "qcovhsp": pl.Float64,
    "stitle": pl.String,
}


class CandidateReadClassificationError(ValueError):
    """Raised when candidate-read classification evidence is malformed."""

    @classmethod
    def malformed_candidates(cls, path: Path) -> CandidateReadClassificationError:
        return cls(f"Candidate read counts are malformed: {path}")

    @classmethod
    def duplicate_identity(cls) -> CandidateReadClassificationError:
        return cls("Candidate read identity is duplicated")

    @classmethod
    def malformed_hits(cls, path: Path) -> CandidateReadClassificationError:
        return cls(f"Read BLAST hits are malformed: {path}")

    @classmethod
    def invalid_bitscore(cls) -> CandidateReadClassificationError:
        return cls("Read BLAST has an invalid bitscore")

    @classmethod
    def malformed_taxids(cls) -> CandidateReadClassificationError:
        return cls("Read BLAST has malformed staxids")

    @classmethod
    def unknown_query(cls) -> CandidateReadClassificationError:
        return cls("Read BLAST contains an unknown query")

    @classmethod
    def invalid_target_taxid(cls) -> CandidateReadClassificationError:
        return cls("target_taxid must be a positive integer")


def load_calibration_target_taxids(path: Path | None, target_taxid: str) -> frozenset[str]:
    """Load the calibration target scope, defaulting to the bait-set target taxid."""
    target_taxid = valid_target_taxid(target_taxid)
    if path is None:
        return frozenset({target_taxid})

    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["taxid"]:
            message = "Calibration target taxids must have exactly one header: taxid"
            raise CandidateReadClassificationError(message)

        taxids: list[str] = []
        for line_number, row in enumerate(reader, start=2):
            if len(row) != 1 or re.fullmatch(r"[1-9][0-9]*", row[0]) is None:
                message = f"Calibration target taxids row {line_number} must contain one canonical positive taxid"
                raise CandidateReadClassificationError(message)
            taxids.append(row[0])

    if len(set(taxids)) != len(taxids):
        message = "Calibration target taxids must be unique"
        raise CandidateReadClassificationError(message)
    if target_taxid not in taxids:
        message = f"Calibration target taxids must include target_taxid {target_taxid}"
        raise CandidateReadClassificationError(message)
    return frozenset(taxids)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify candidate reads from read BLAST evidence.")
    parser.add_argument("--candidate-read-counts", type=Path, required=True)
    parser.add_argument("--blast-hits", type=Path, required=True)
    parser.add_argument("--target-taxid", required=True)
    parser.add_argument("--calibration-target-taxids", type=Path)
    parser.add_argument("--output", type=Path, required=True)
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


def construct_candidate_read_counts(path: Path) -> pl.LazyFrame:
    try:
        if not has_expected_header(path, CANDIDATE_FIELDS):
            raise CandidateReadClassificationError.malformed_candidates(path)
        candidates = pl.scan_csv(
            path,
            separator="\t",
            schema=CANDIDATE_SCHEMA,
            null_values=[],
            raise_if_empty=False,
        ).with_row_index("_candidate_order")
        error = first_error(
            pl.concat(
                [
                    candidates.filter(
                        ~(
                            pl.all_horizontal(
                                pl.col("metagenome_id").is_not_null(),
                                pl.col("read_id").is_not_null(),
                                pl.col("representative_id").is_not_null(),
                                pl.col("read_length").is_not_null(),
                                pl.col("bait_count").is_not_null(),
                                pl.col("metagenome_id").str.len_chars() > 0,
                                pl.col("read_id").str.len_chars() > 0,
                                pl.col("representative_id").str.len_chars() > 0,
                                pl.col("read_length") > 0,
                                pl.col("bait_count") > 0,
                            )
                        ).fill_null(value=True),
                    ).with_columns(pl.lit("malformed").alias("error")).select("error"),
                    candidates.group_by("metagenome_id", "read_id").len().filter(
                        pl.col("len") > 1,
                    ).with_columns(pl.lit("duplicate_identity").alias("error")).select("error"),
                ],
                how="vertical",
            ),
        )
    except pl.exceptions.PolarsError as exception:
        raise CandidateReadClassificationError.malformed_candidates(path) from exception
    if error == "duplicate_identity":
        raise CandidateReadClassificationError.duplicate_identity()
    if error is not None:
        raise CandidateReadClassificationError.malformed_candidates(path)
    return candidates


def valid_target_taxid(target_taxid: str) -> str:
    if re.fullmatch(r"[1-9][0-9]*", target_taxid) is None:
        raise CandidateReadClassificationError.invalid_target_taxid()
    return target_taxid


def parse_bitscore(value: str | None) -> Decimal:
    if value is None:
        raise CandidateReadClassificationError.invalid_bitscore()
    try:
        bitscore = Decimal(value)
    except InvalidOperation as exception:
        raise CandidateReadClassificationError.invalid_bitscore() from exception
    if not bitscore.is_finite():
        raise CandidateReadClassificationError.invalid_bitscore()
    return bitscore


@dataclass(frozen=True)
class RepresentativeScores:
    representative_id: str
    best_target_bit_score: Decimal | None
    best_non_target_bit_score: Decimal | None


def construct_read_hits(
    path: Path,
    *,
    candidate_representatives: pl.LazyFrame,
    target_taxid: str,
    calibration_target_taxids: frozenset[str] | None = None,
) -> pl.LazyFrame:
    target_taxid = valid_target_taxid(target_taxid)
    target_scope = frozenset({target_taxid}) if calibration_target_taxids is None else calibration_target_taxids
    if target_taxid not in target_scope or any(valid_target_taxid(taxid) != taxid for taxid in target_scope):
        message = f"Calibration target taxids must be positive integers and include target_taxid {target_taxid}"
        raise CandidateReadClassificationError(message)
    try:
        if not has_expected_header(path, BLAST_FIELDS):
            raise CandidateReadClassificationError.malformed_hits(path)
        hits = pl.scan_csv(
            path,
            separator="\t",
            schema=BLAST_SCHEMA,
            null_values=[],
            raise_if_empty=False,
        ).with_columns(
            pl.col("staxids").str.replace_all(",", ";").str.split(";").alias("_taxids"),
        )
        malformed_taxids = (
            (pl.col("staxids") != "N/A")
            & (
                pl.col("staxids").is_null()
                | pl.col("_taxids").list.eval(
                    pl.element().str.strip_chars().str.contains(r"^[1-9][0-9]*$").not_(),
                ).list.any()
            )
        )
        error = first_error(
            pl.concat(
                [
                    hits.filter(pl.any_horizontal(*(pl.col(field).is_null() for field in BLAST_FIELDS if field != "bitscore"))).with_columns(
                        pl.lit("malformed").alias("error"),
                    ).select("error"),
                    hits.filter(malformed_taxids).with_columns(pl.lit("malformed_taxids").alias("error")).select("error"),
                    hits.join(
                        candidate_representatives.select("representative_id").unique(),
                        left_on="qseqid",
                        right_on="representative_id",
                        how="anti",
                    ).with_columns(pl.lit("unknown_query").alias("error")).select("error"),
                ],
                how="vertical",
            ),
        )
    except pl.exceptions.PolarsError as exception:
        raise CandidateReadClassificationError.malformed_hits(path) from exception
    try:
        for value in hits.select("bitscore").collect().get_column("bitscore"):
            parse_bitscore(value)
    except pl.exceptions.PolarsError as exception:
        raise CandidateReadClassificationError.malformed_hits(path) from exception
    match error:
        case "malformed_taxids":
            raise CandidateReadClassificationError.malformed_taxids()
        case "unknown_query":
            raise CandidateReadClassificationError.unknown_query()
        case None:
            return hits.with_columns(
                pl.col("_taxids").list.eval(pl.element().str.strip_chars()).alias("_taxids"),
            ).with_columns(
                pl.any_horizontal(
                    *(pl.col("_taxids").list.contains(taxid) for taxid in sorted(target_scope)),
                ).alias("_target"),
                pl.col("_taxids").list.eval(pl.element().is_in(sorted(target_scope)).not_()).list.any().alias("_non_target"),
            ).with_columns(
                pl.when(pl.col("staxids") == "N/A").then(pl.lit(value=False)).otherwise(pl.col("_target")).alias("_target"),
                pl.when(pl.col("staxids") == "N/A").then(pl.lit(value=True)).otherwise(pl.col("_non_target")).alias("_non_target"),
            )
        case _:
            raise CandidateReadClassificationError.malformed_hits(path)


def construct_representative_scores(hits: pl.LazyFrame) -> tuple[RepresentativeScores, ...]:
    hit_rows = hits.select(
        "qseqid",
        "bitscore",
        "_target",
        "_non_target",
    ).sort("qseqid").collect().iter_rows()

    def representative_scores(rows: Iterable[tuple[str, str, bool, bool]]) -> tuple[Decimal | None, Decimal | None]:
        best_target: Decimal | None = None
        best_non_target: Decimal | None = None
        for _, bitscore_text, target, non_target in rows:
            bitscore = parse_bitscore(bitscore_text)
            if target and (best_target is None or bitscore > best_target):
                best_target = bitscore
            if non_target and (best_non_target is None or bitscore > best_non_target):
                best_non_target = bitscore
        return best_target, best_non_target

    return tuple(
        RepresentativeScores(representative_id, *representative_scores(rows))
        for representative_id, rows in groupby(hit_rows, key=lambda row: row[0])
    )


def score_text(score: Decimal | None) -> str:
    if score is None:
        return ""
    text = format(score, "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def construct_representative_classifications(
    scores: tuple[RepresentativeScores, ...],
    tie_tolerance: Decimal,
) -> pl.LazyFrame:
    rows = tuple(
        (
            score.representative_id,
            "NO_HIT" if score.best_target_bit_score is None and score.best_non_target_bit_score is None
            else "NON_TARGET" if score.best_target_bit_score is None
            else "TARGET" if score.best_non_target_bit_score is None
            else "TIED" if abs(score.best_target_bit_score - score.best_non_target_bit_score) <= tie_tolerance
            else "TARGET" if score.best_target_bit_score > score.best_non_target_bit_score
            else "NON_TARGET",
            score_text(score.best_target_bit_score),
            score_text(score.best_non_target_bit_score),
        )
        for score in scores
    )
    return pl.LazyFrame(
        rows,
        schema={
            "representative_id": pl.String,
            "classification": pl.String,
            "best_target_bit_score": pl.String,
            "best_non_target_bit_score": pl.String,
        },
        orient="row",
    )


def construct_candidate_read_classifications(
    candidates: pl.LazyFrame,
    hits: pl.LazyFrame,
    *,
    tie_tolerance: Decimal = Decimal("0.1"),
) -> pl.LazyFrame:
    scores = construct_representative_classifications(
        construct_representative_scores(hits),
        tie_tolerance,
    )
    return (
        candidates.join(scores, on="representative_id", how="left")
        .with_columns(
            pl.col("classification").fill_null("NO_HIT"),
            pl.col("best_target_bit_score").fill_null(""),
            pl.col("best_non_target_bit_score").fill_null(""),
        )
        .sort("_candidate_order")
        .select(CLASSIFIED_FIELDS)
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    candidates = construct_candidate_read_counts(args.candidate_read_counts)
    calibration_target_taxids = load_calibration_target_taxids(
        args.calibration_target_taxids,
        args.target_taxid,
    )
    hits = construct_read_hits(
        args.blast_hits,
        candidate_representatives=candidates.select("representative_id").unique(),
        target_taxid=args.target_taxid,
        calibration_target_taxids=calibration_target_taxids,
    )
    classified = construct_candidate_read_classifications(candidates, hits).collect()
    classified.write_csv(args.output, separator="\t", null_value="", quote_style="never")


if __name__ == "__main__":
    main()
