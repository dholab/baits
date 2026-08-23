#!/usr/bin/env python3
"""Publish run-wide taxonomic exact-match rows."""

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

PRIVATE_CANDIDATE_SCHEMA = {
    "candidate_kmer_id": pl.String,
    "kmer": pl.String,
    "source_copy_count": pl.UInt64,
    "background_occurrences": pl.UInt64,
    "outcome": pl.String,
}
PRIVATE_MATCH_SCHEMA = {
    "qseqid": pl.String,
    "saccver": pl.String,
    "staxids": pl.String,
    "sscinames": pl.String,
    "pident": pl.Float64,
    "length": pl.UInt64,
    "mismatch": pl.UInt64,
    "gapopen": pl.UInt64,
    "qstart": pl.UInt64,
    "qend": pl.UInt64,
    "sstart": pl.UInt64,
    "send": pl.UInt64,
    "stitle": pl.String,
}
PUBLIC_MATCH_COLUMNS = [
    "design_id",
    "kmer",
    "subject_accession",
    "taxids",
    "scientific_names",
    "percent_identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "subject_title",
]


class TaxonomicMatchReportError(ValueError):
    """Raised when private taxonomic tables cannot form the public report."""


@dataclass(frozen=True)
class DesignTables:
    """Private taxonomic tables for one design."""

    design_id: str
    candidates: Path
    matches: Path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-ids", type=Path, required=True)
    parser.add_argument("--candidate-manifests-dir", type=Path, required=True)
    parser.add_argument("--taxonomic-matches-dir", type=Path, required=True)
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
        message = "Taxonomic match report design IDs must be nonblank canonical TSV values"
        raise TaxonomicMatchReportError(message)
    if len(set(design_ids)) != len(design_ids):
        message = "Taxonomic match report design IDs must be unique"
        raise TaxonomicMatchReportError(message)
    return design_ids


def discover_staged_tables(root: Path, filename: str) -> tuple[Path, ...]:
    if not root.exists():
        return ()
    directories = tuple(sorted(path for path in root.iterdir() if path.is_dir()))
    tables = tuple(directory / filename for directory in directories)
    missing = next((table for table in tables if not table.is_file()), None)
    if missing is not None:
        message = f"Taxonomic match report input is missing {missing}"
        raise TaxonomicMatchReportError(message)
    return tables


def pair_design_tables(
    design_ids: Sequence[str],
    candidate_tables: Sequence[Path],
    match_tables: Sequence[Path],
) -> tuple[DesignTables, ...]:
    if len(design_ids) != len(candidate_tables) or len(design_ids) != len(match_tables):
        message = "Taxonomic match report design and table counts disagree"
        raise TaxonomicMatchReportError(message)
    return tuple(
        DesignTables(design_id, candidates, matches)
        for design_id, candidates, matches in zip(
            design_ids,
            candidate_tables,
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
        raise TaxonomicMatchReportError(message)
    return pl.scan_csv(
        path,
        separator="\t",
        schema=dict(schema),
        quote_char=None,
    )


def report_one_design(tables: DesignTables) -> pl.LazyFrame:
    candidates = scan_exact_tsv(
        tables.candidates,
        PRIVATE_CANDIDATE_SCHEMA,
        "Private candidate table",
    )
    matches = scan_exact_tsv(
        tables.matches,
        PRIVATE_MATCH_SCHEMA,
        "Private taxonomic match table",
    )
    unsupported_outcome = (
        candidates.filter(pl.col("outcome") == "off_target_match")
        .select("candidate_kmer_id")
        .join(
            matches.select(pl.col("qseqid").alias("candidate_kmer_id")).unique(),
            on="candidate_kmer_id",
            how="anti",
        )
        .sort("candidate_kmer_id")
        .limit(1)
        .collect()
    )
    if not unsupported_outcome.is_empty():
        message = (
            f"Taxonomic matches for {tables.design_id} do not support off-target candidate "
            f"{unsupported_outcome.item(0, 'candidate_kmer_id')}"
        )
        raise TaxonomicMatchReportError(message)
    joined = matches.join(
        candidates.select("candidate_kmer_id", "kmer"),
        left_on="qseqid",
        right_on="candidate_kmer_id",
        how="left",
    )
    unknown_candidate = (
        joined.filter(pl.col("kmer").is_null())
        .select("qseqid")
        .sort("qseqid")
        .limit(1)
        .collect()
    )
    if not unknown_candidate.is_empty():
        message = (
            f"Taxonomic matches for {tables.design_id} contain an unknown candidate "
            f"identifier {unknown_candidate.item(0, 'qseqid')}"
        )
        raise TaxonomicMatchReportError(message)
    return joined.select(
        pl.lit(tables.design_id).alias("design_id"),
        "kmer",
        pl.col("saccver").alias("subject_accession"),
        pl.col("staxids").alias("taxids"),
        pl.col("sscinames").alias("scientific_names"),
        pl.col("pident").alias("percent_identity"),
        pl.col("length").alias("alignment_length"),
        pl.col("mismatch").alias("mismatches"),
        pl.col("gapopen").alias("gap_opens"),
        pl.col("qstart").alias("query_start"),
        pl.col("qend").alias("query_end"),
        pl.col("sstart").alias("subject_start"),
        pl.col("send").alias("subject_end"),
        pl.col("stitle").alias("subject_title"),
    )


def construct_report(inputs: Sequence[DesignTables]) -> pl.LazyFrame:
    reports = tuple(report_one_design(tables) for tables in inputs)
    if not reports:
        return pl.LazyFrame(schema=dict.fromkeys(PUBLIC_MATCH_COLUMNS, pl.String))
    return pl.concat(reports, how="vertical").sort(PUBLIC_MATCH_COLUMNS)


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
        discover_staged_tables(args.candidate_manifests_dir, "candidate_kmers.tsv"),
        discover_staged_tables(args.taxonomic_matches_dir, "taxonomic_blast_hits.tsv"),
    )
    matches = construct_report(inputs).collect()
    write_gzip_tsv(matches.select(PUBLIC_MATCH_COLUMNS), args.matches_out)


if __name__ == "__main__":
    main()
