#!/usr/bin/env python3
"""Publish run-wide source-discovery loci and matches."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import BinaryIO

PRIVATE_LOCUS_SCHEMA = {
    "candidate_locus_id": pl.String,
    "source_sequence_id": pl.String,
    "query_id": pl.String,
    "query_group": pl.String,
    "assembly_sequence_id": pl.String,
    "assembly_start": pl.UInt64,
    "assembly_end": pl.UInt64,
    "strand": pl.String,
    "percent_identity": pl.Float64,
    "query_coverage": pl.Float64,
}
PRIVATE_MATCH_SCHEMA = {
    "query_id": pl.String,
    "query_length": pl.UInt64,
    "query_start": pl.UInt64,
    "query_end": pl.UInt64,
    "assembly_sequence_id": pl.String,
    "assembly_start": pl.UInt64,
    "assembly_end": pl.UInt64,
    "alignment_length": pl.UInt64,
    "percent_identity": pl.Float64,
}
PUBLIC_LOCUS_SCHEMA = {"design_id": pl.String, **PRIVATE_LOCUS_SCHEMA}
PUBLIC_MATCH_SCHEMA = {"design_id": pl.String, **PRIVATE_MATCH_SCHEMA}


@dataclass(frozen=True)
class DiscoveryTables:
    """Private source-discovery tables for one design."""

    design_id: str
    candidate_loci: Path
    source_matches: Path


@dataclass(frozen=True)
class SourceDiscoveryReport:
    """Run-wide public source-discovery relations."""

    candidate_loci: pl.LazyFrame
    source_matches: pl.LazyFrame


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-ids", type=Path, required=True)
    parser.add_argument("--candidate-loci-dir", type=Path, required=True)
    parser.add_argument("--source-matches-dir", type=Path, required=True)
    parser.add_argument("--candidate-loci-out", type=Path, required=True)
    parser.add_argument("--source-matches-out", type=Path, required=True)
    return parser.parse_args(argv)


def discover_staged_tables(root: Path, filename: str) -> tuple[Path, ...]:
    if not root.exists():
        return ()
    return tuple(
        directory / filename
        for directory in sorted(path for path in root.iterdir() if path.is_dir())
    )


def pair_design_tables(
    design_ids: Sequence[str],
    locus_tables: Sequence[Path],
    match_tables: Sequence[Path],
) -> tuple[DiscoveryTables, ...]:
    if len(design_ids) != len(locus_tables) or len(design_ids) != len(match_tables):
        message = "Source-discovery design and table counts disagree"
        raise ValueError(message)
    return tuple(
        DiscoveryTables(design_id, candidate_loci, source_matches)
        for design_id, candidate_loci, source_matches in zip(
            design_ids,
            locus_tables,
            match_tables,
            strict=True,
        )
    )


def scan_matches(path: Path) -> pl.LazyFrame:
    if path.stat().st_size == 0:
        return pl.LazyFrame(schema=PRIVATE_MATCH_SCHEMA)
    return pl.scan_csv(
        path,
        separator="\t",
        has_header=False,
        quote_char=None,
        schema=PRIVATE_MATCH_SCHEMA,
    )


def report_one_design(tables: DiscoveryTables) -> SourceDiscoveryReport:
    candidate_loci = pl.scan_csv(
        tables.candidate_loci,
        separator="\t",
        quote_char=None,
        schema=PRIVATE_LOCUS_SCHEMA,
    )
    source_matches = scan_matches(tables.source_matches)
    return SourceDiscoveryReport(
        candidate_loci=candidate_loci.select(
            pl.lit(tables.design_id).alias("design_id"),
            *PRIVATE_LOCUS_SCHEMA,
        ),
        source_matches=source_matches.select(
            pl.lit(tables.design_id).alias("design_id"),
            *PRIVATE_MATCH_SCHEMA,
        ),
    )


def construct_report(inputs: Sequence[DiscoveryTables]) -> SourceDiscoveryReport:
    reports = tuple(report_one_design(tables) for tables in inputs)
    if not reports:
        return SourceDiscoveryReport(
            candidate_loci=pl.LazyFrame(schema=PUBLIC_LOCUS_SCHEMA),
            source_matches=pl.LazyFrame(schema=PUBLIC_MATCH_SCHEMA),
        )
    return SourceDiscoveryReport(
        candidate_loci=pl.concat(
            [report.candidate_loci for report in reports],
            how="vertical",
        ).sort("design_id", "candidate_locus_id"),
        source_matches=pl.concat(
            [report.source_matches for report in reports],
            how="vertical",
        ).sort(
            "design_id",
            "query_id",
            "assembly_sequence_id",
            "query_start",
            "query_end",
            "assembly_start",
            "assembly_end",
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
        args.design_ids.read_text().splitlines(),
        discover_staged_tables(args.candidate_loci_dir, "candidate_loci.tsv"),
        discover_staged_tables(args.source_matches_dir, "query_blast_hits.txt"),
    )
    report = construct_report(inputs)
    candidate_loci, source_matches = pl.collect_all(
        [report.candidate_loci, report.source_matches],
    )
    candidate_loci.write_csv(args.candidate_loci_out, separator="\t", quote_style="never")
    write_gzip_tsv(source_matches, args.source_matches_out)


if __name__ == "__main__":
    main()
