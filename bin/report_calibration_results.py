#!/usr/bin/env python3
"""Publish run-wide calibration read rows and raw BLAST matches."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import BinaryIO

PRIVATE_READ_SCHEMA = {
    "metagenome_id": pl.String,
    "read_id": pl.String,
    "read_length": pl.UInt64,
    "bait_count": pl.UInt64,
    "query_id": pl.String,
    "best_target_bit_score": pl.String,
    "best_other_bit_score": pl.String,
}
PRIVATE_MATCH_SCHEMA = {
    "qseqid": pl.String,
    "qlen": pl.UInt64,
    "saccver": pl.String,
    "staxids": pl.String,
    "pident": pl.Float64,
    "length": pl.UInt64,
    "mismatch": pl.UInt64,
    "gapopen": pl.UInt64,
    "qstart": pl.UInt64,
    "qend": pl.UInt64,
    "sstart": pl.UInt64,
    "send": pl.UInt64,
    "evalue": pl.Float64,
    "bitscore": pl.Float64,
    "qcovhsp": pl.Float64,
    "stitle": pl.String,
}
PUBLIC_READ_COLUMNS = ["design_id", *PRIVATE_READ_SCHEMA]
PUBLIC_MATCH_COLUMNS = [
    "design_id",
    "query_id",
    "query_length",
    "subject_accession",
    "taxids",
    "percent_identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bit_score",
    "query_coverage",
    "subject_title",
]


class CalibrationReportError(ValueError):
    """Raised when private calibration tables cannot form the public report."""


@dataclass(frozen=True)
class DesignTables:
    design_id: str
    reads: Path
    matches: Path


@dataclass(frozen=True)
class CalibrationReport:
    reads: pl.LazyFrame
    matches: pl.LazyFrame


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-ids", type=Path, required=True)
    parser.add_argument("--calibration-reads-dir", type=Path, required=True)
    parser.add_argument("--calibration-matches-dir", type=Path, required=True)
    parser.add_argument("--reads-out", type=Path, required=True)
    parser.add_argument("--matches-out", type=Path, required=True)
    return parser.parse_args(argv)


def read_design_ids(path: Path) -> tuple[str, ...]:
    design_ids = tuple(path.read_text().splitlines())
    invalid = next(
        (
            design_id
            for design_id in design_ids
            if not design_id or design_id != design_id.strip() or "\t" in design_id
        ),
        None,
    )
    if invalid is not None:
        message = "Calibration report design IDs must be nonblank canonical TSV values"
        raise CalibrationReportError(message)
    if len(set(design_ids)) != len(design_ids):
        message = "Calibration report design IDs must be unique"
        raise CalibrationReportError(message)
    return design_ids


def discover_staged_tables(root: Path, filename: str) -> tuple[Path, ...]:
    if not root.exists():
        return ()
    directories = tuple(sorted(path for path in root.iterdir() if path.is_dir()))
    tables = tuple(directory / filename for directory in directories)
    missing = next((table for table in tables if not table.is_file()), None)
    if missing is not None:
        message = f"Calibration report input is missing {missing}"
        raise CalibrationReportError(message)
    return tables


def pair_design_tables(
    design_ids: Sequence[str],
    read_tables: Sequence[Path],
    match_tables: Sequence[Path],
) -> tuple[DesignTables, ...]:
    if len(design_ids) != len(read_tables) or len(design_ids) != len(match_tables):
        message = "Calibration report design and table counts disagree"
        raise CalibrationReportError(message)
    return tuple(
        DesignTables(design_id, reads, matches)
        for design_id, reads, matches in zip(
            design_ids,
            read_tables,
            match_tables,
            strict=True,
        )
    )


def scan_exact_tsv(
    path: Path,
    schema: Mapping[str, type[pl.DataType]],
    label: str,
) -> pl.LazyFrame:
    with path.open() as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    if header != list(schema):
        message = f"{label} has an unexpected schema: {path}"
        raise CalibrationReportError(message)
    return pl.scan_csv(
        path,
        separator="\t",
        schema=dict(schema),
        quote_char=None,
        null_values=[],
        raise_if_empty=False,
    )


def report_one_design(tables: DesignTables) -> CalibrationReport:
    reads = scan_exact_tsv(
        tables.reads,
        PRIVATE_READ_SCHEMA,
        "Private calibration read table",
    )
    matches = scan_exact_tsv(
        tables.matches,
        PRIVATE_MATCH_SCHEMA,
        "Private calibration match table",
    )
    unknown_query = (
        matches.select("qseqid")
        .unique()
        .join(
            reads.select(pl.col("query_id").alias("qseqid")).unique(),
            on="qseqid",
            how="anti",
        )
        .sort("qseqid")
        .limit(1)
        .collect()
    )
    if not unknown_query.is_empty():
        message = (
            f"Calibration matches for {tables.design_id} contain unknown query "
            f"{unknown_query.item(0, 'qseqid')}"
        )
        raise CalibrationReportError(message)
    return CalibrationReport(
        reads=reads.select(
            pl.lit(tables.design_id).alias("design_id"),
            *PRIVATE_READ_SCHEMA,
        ),
        matches=matches.select(
            pl.lit(tables.design_id).alias("design_id"),
            pl.col("qseqid").alias("query_id"),
            pl.col("qlen").alias("query_length"),
            pl.col("saccver").alias("subject_accession"),
            pl.col("staxids").alias("taxids"),
            pl.col("pident").alias("percent_identity"),
            pl.col("length").alias("alignment_length"),
            pl.col("mismatch").alias("mismatches"),
            pl.col("gapopen").alias("gap_opens"),
            pl.col("qstart").alias("query_start"),
            pl.col("qend").alias("query_end"),
            pl.col("sstart").alias("subject_start"),
            pl.col("send").alias("subject_end"),
            "evalue",
            pl.col("bitscore").alias("bit_score"),
            pl.col("qcovhsp").alias("query_coverage"),
            pl.col("stitle").alias("subject_title"),
        ),
    )


def construct_report(inputs: Sequence[DesignTables]) -> CalibrationReport:
    reports = tuple(report_one_design(tables) for tables in inputs)
    if not reports:
        return CalibrationReport(
            reads=pl.LazyFrame(schema=dict.fromkeys(PUBLIC_READ_COLUMNS, pl.String)),
            matches=pl.LazyFrame(schema=dict.fromkeys(PUBLIC_MATCH_COLUMNS, pl.String)),
        )
    return CalibrationReport(
        reads=pl.concat([report.reads for report in reports], how="vertical").sort(
            "design_id",
            "metagenome_id",
            "read_id",
        ),
        matches=pl.concat([report.matches for report in reports], how="vertical").sort(
            PUBLIC_MATCH_COLUMNS,
        ),
    )


def write_gzip_tsv(frame: pl.DataFrame, path: Path) -> None:
    with (
        path.open("wb") as raw_handle,
        gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as gzip_handle,
    ):
        frame.write_csv(
            cast("BinaryIO", gzip_handle),
            separator="\t",
            quote_style="never",
        )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    inputs = pair_design_tables(
        read_design_ids(args.design_ids),
        discover_staged_tables(args.calibration_reads_dir, "calibration_reads.tsv"),
        discover_staged_tables(args.calibration_matches_dir, "read_blast_hits.tsv"),
    )
    report = construct_report(inputs)
    reads, matches = pl.collect_all([report.reads, report.matches])
    reads.select(PUBLIC_READ_COLUMNS).write_csv(
        args.reads_out,
        separator="\t",
        quote_style="never",
        null_value="",
    )
    write_gzip_tsv(matches.select(PUBLIC_MATCH_COLUMNS), args.matches_out)


if __name__ == "__main__":
    main()
