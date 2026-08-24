#!/usr/bin/env python3
"""Project each calibration query's best target and other scores onto its reads."""

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

READ_FIELDS = (
    "metagenome_id",
    "read_id",
    "read_length",
    "bait_count",
    "query_id",
)
BLAST_FIELDS = (
    "qseqid",
    "qlen",
    "saccver",
    "staxids",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
    "stitle",
)
SCORED_READ_FIELDS = (*READ_FIELDS, "best_target_bit_score", "best_other_bit_score")
READ_SCHEMA = {
    "metagenome_id": pl.String,
    "read_id": pl.String,
    "read_length": pl.Int64,
    "bait_count": pl.Int64,
    "query_id": pl.String,
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


class CalibrationReadScoreError(ValueError):
    """Raised when calibration-read counts or matches are malformed."""

    @classmethod
    def malformed_reads(cls, path: Path) -> CalibrationReadScoreError:
        return cls(f"Read bait counts are malformed: {path}")

    @classmethod
    def duplicate_identity(cls) -> CalibrationReadScoreError:
        return cls("Calibration read identity is duplicated")

    @classmethod
    def malformed_hits(cls, path: Path) -> CalibrationReadScoreError:
        return cls(f"Calibration BLAST matches are malformed: {path}")

    @classmethod
    def invalid_bitscore(cls) -> CalibrationReadScoreError:
        return cls("Calibration BLAST has an invalid bitscore")

    @classmethod
    def malformed_taxids(cls) -> CalibrationReadScoreError:
        return cls("Calibration BLAST has malformed staxids")

    @classmethod
    def unknown_query(cls) -> CalibrationReadScoreError:
        return cls("Calibration BLAST contains an unknown query")

    @classmethod
    def invalid_target_taxid(cls) -> CalibrationReadScoreError:
        return cls("target_taxid must be a positive integer")


@dataclass(frozen=True)
class QueryScores:
    query_id: str
    best_target_bit_score: Decimal | None
    best_other_bit_score: Decimal | None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-bait-counts", type=Path, required=True)
    parser.add_argument("--blast-matches", type=Path, required=True)
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


def construct_read_bait_counts(path: Path) -> pl.LazyFrame:
    try:
        if not has_expected_header(path, READ_FIELDS):
            raise CalibrationReadScoreError.malformed_reads(path)
        reads = pl.scan_csv(
            path,
            separator="\t",
            schema=READ_SCHEMA,
            null_values=[],
            raise_if_empty=False,
        ).with_row_index("_read_order")
        error = first_error(
            pl.concat(
                [
                    reads.filter(
                        ~pl.all_horizontal(
                            pl.col("metagenome_id").is_not_null(),
                            pl.col("read_id").is_not_null(),
                            pl.col("query_id").is_not_null(),
                            pl.col("read_length").is_not_null(),
                            pl.col("bait_count").is_not_null(),
                            pl.col("metagenome_id").str.len_chars() > 0,
                            pl.col("read_id").str.len_chars() > 0,
                            pl.col("query_id").str.contains(r"^sequence_[0-9a-f]{64}$"),
                            pl.col("read_length") > 0,
                            pl.col("bait_count") > 0,
                        ).fill_null(value=False),
                    )
                    .with_columns(pl.lit("malformed").alias("error"))
                    .select("error"),
                    reads.group_by("metagenome_id", "read_id")
                    .len()
                    .filter(pl.col("len") > 1)
                    .with_columns(pl.lit("duplicate_identity").alias("error"))
                    .select("error"),
                ],
                how="vertical",
            ),
        )
    except pl.exceptions.PolarsError as error:
        raise CalibrationReadScoreError.malformed_reads(path) from error
    if error == "duplicate_identity":
        raise CalibrationReadScoreError.duplicate_identity()
    if error is not None:
        raise CalibrationReadScoreError.malformed_reads(path)
    return reads


def valid_target_taxid(target_taxid: str) -> str:
    if re.fullmatch(r"[1-9][0-9]*", target_taxid) is None:
        raise CalibrationReadScoreError.invalid_target_taxid()
    return target_taxid


def load_calibration_target_taxids(
    path: Path | None,
    target_taxid: str,
) -> frozenset[str]:
    target_taxid = valid_target_taxid(target_taxid)
    if path is None:
        return frozenset({target_taxid})
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["taxid"]:
            message = "Calibration target taxids must have exactly one header: taxid"
            raise CalibrationReadScoreError(message)
        taxids = tuple(row[0] if len(row) == 1 else "" for row in reader)
    if any(re.fullmatch(r"[1-9][0-9]*", taxid) is None for taxid in taxids):
        message = "Calibration target taxids must contain canonical positive taxids"
        raise CalibrationReadScoreError(message)
    if len(set(taxids)) != len(taxids):
        message = "Calibration target taxids must be unique"
        raise CalibrationReadScoreError(message)
    if target_taxid not in taxids:
        message = f"Calibration target taxids must include target_taxid {target_taxid}"
        raise CalibrationReadScoreError(message)
    return frozenset(taxids)


def parse_bitscore(value: str | None) -> Decimal:
    if value is None:
        raise CalibrationReadScoreError.invalid_bitscore()
    try:
        bitscore = Decimal(value)
    except InvalidOperation as error:
        raise CalibrationReadScoreError.invalid_bitscore() from error
    if not bitscore.is_finite():
        raise CalibrationReadScoreError.invalid_bitscore()
    return bitscore


def construct_read_hits(
    path: Path,
    *,
    query_ids: pl.LazyFrame,
    target_taxid: str,
    calibration_target_taxids: frozenset[str] | None = None,
) -> pl.LazyFrame:
    target_taxid = valid_target_taxid(target_taxid)
    target_scope = (
        frozenset({target_taxid})
        if calibration_target_taxids is None
        else calibration_target_taxids
    )
    if target_taxid not in target_scope or any(
        valid_target_taxid(taxid) != taxid for taxid in target_scope
    ):
        message = (
            f"Calibration target taxids must be positive integers and include target_taxid "
            f"{target_taxid}"
        )
        raise CalibrationReadScoreError(message)
    try:
        if not has_expected_header(path, BLAST_FIELDS):
            raise CalibrationReadScoreError.malformed_hits(path)
        hits = pl.scan_csv(
            path,
            separator="\t",
            schema=BLAST_SCHEMA,
            null_values=[],
            raise_if_empty=False,
        ).with_columns(
            pl.col("staxids").str.replace_all(",", ";").str.split(";").alias("_taxids"),
        )
        malformed_taxids = (pl.col("staxids") != "N/A") & (
            pl.col("staxids").is_null()
            | pl.col("_taxids")
            .list.eval(pl.element().str.strip_chars().str.contains(r"^[1-9][0-9]*$").not_())
            .list.any()
        )
        error = first_error(
            pl.concat(
                [
                    hits.filter(
                        pl.any_horizontal(
                            *(
                                pl.col(field).is_null()
                                for field in BLAST_FIELDS
                                if field != "bitscore"
                            ),
                        ),
                    )
                    .with_columns(pl.lit("malformed").alias("error"))
                    .select("error"),
                    hits.filter(malformed_taxids)
                    .with_columns(pl.lit("malformed_taxids").alias("error"))
                    .select("error"),
                    hits.join(
                        query_ids,
                        left_on="qseqid",
                        right_on="query_id",
                        how="anti",
                    )
                    .with_columns(pl.lit("unknown_query").alias("error"))
                    .select("error"),
                ],
                how="vertical",
            ),
        )
    except pl.exceptions.PolarsError as error:
        raise CalibrationReadScoreError.malformed_hits(path) from error
    try:
        for value in hits.select("bitscore").collect().get_column("bitscore"):
            parse_bitscore(value)
    except pl.exceptions.PolarsError as error:
        raise CalibrationReadScoreError.malformed_hits(path) from error
    match error:
        case "malformed_taxids":
            raise CalibrationReadScoreError.malformed_taxids()
        case "unknown_query":
            raise CalibrationReadScoreError.unknown_query()
        case None:
            return hits.with_columns(
                pl.col("_taxids").list.eval(pl.element().str.strip_chars()).alias("_taxids"),
            ).with_columns(
                pl.any_horizontal(
                    *(pl.col("_taxids").list.contains(taxid) for taxid in sorted(target_scope)),
                ).alias("_target"),
                pl.col("_taxids")
                .list.eval(pl.element().is_in(sorted(target_scope)).not_())
                .list.any()
                .alias("_other"),
            ).with_columns(
                pl.when(pl.col("staxids") == "N/A")
                .then(pl.lit(value=False))
                .otherwise(pl.col("_target"))
                .alias("_target"),
                pl.when(pl.col("staxids") == "N/A")
                .then(pl.lit(value=True))
                .otherwise(pl.col("_other"))
                .alias("_other"),
            )
        case _:
            raise CalibrationReadScoreError.malformed_hits(path)


def score_text(score: Decimal | None) -> str:
    if score is None:
        return ""
    text = format(score, "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def construct_query_scores(hits: pl.LazyFrame) -> tuple[QueryScores, ...]:
    hit_rows = hits.select("qseqid", "bitscore", "_target", "_other").sort("qseqid").collect().iter_rows()

    def best_scores(rows: Iterable[tuple[str, str, bool, bool]]) -> tuple[Decimal | None, Decimal | None]:
        best_target: Decimal | None = None
        best_other: Decimal | None = None
        for _, bitscore_text, target, other in rows:
            bitscore = parse_bitscore(bitscore_text)
            if target and (best_target is None or bitscore > best_target):
                best_target = bitscore
            if other and (best_other is None or bitscore > best_other):
                best_other = bitscore
        return best_target, best_other

    return tuple(
        QueryScores(query_id, *best_scores(rows))
        for query_id, rows in groupby(hit_rows, key=lambda row: row[0])
    )


def construct_scored_reads(reads: pl.LazyFrame, hits: pl.LazyFrame) -> pl.LazyFrame:
    scores = pl.LazyFrame(
        (
            (score.query_id, score_text(score.best_target_bit_score), score_text(score.best_other_bit_score))
            for score in construct_query_scores(hits)
        ),
        schema={
            "query_id": pl.String,
            "best_target_bit_score": pl.String,
            "best_other_bit_score": pl.String,
        },
        orient="row",
    )
    return (
        reads.join(scores, on="query_id", how="left")
        .with_columns(
            pl.col("best_target_bit_score").fill_null(""),
            pl.col("best_other_bit_score").fill_null(""),
        )
        .sort("_read_order")
        .select(SCORED_READ_FIELDS)
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    reads = construct_read_bait_counts(args.read_bait_counts)
    target_scope = load_calibration_target_taxids(
        args.calibration_target_taxids,
        args.target_taxid,
    )
    hits = construct_read_hits(
        args.blast_matches,
        query_ids=reads.select("query_id").unique(),
        target_taxid=args.target_taxid,
        calibration_target_taxids=target_scope,
    )
    construct_scored_reads(reads, hits).collect().write_csv(
        args.output,
        separator="\t",
        null_value="",
        quote_style="never",
    )


if __name__ == "__main__":
    main()
