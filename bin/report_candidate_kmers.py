#!/usr/bin/env python3
"""Publish run-wide candidate k-mer and source-occurrence tables."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

PRIVATE_CANDIDATE_SCHEMA = {
    "candidate_kmer_id": pl.String,
    "bait_id": pl.String,
    "kmer": pl.String,
    "source_copy_count": pl.UInt64,
    "background_occurrences": pl.UInt64,
    "outcome": pl.String,
}
PRIVATE_OCCURRENCE_SCHEMA = {
    "candidate_kmer_id": pl.String,
    "source_sequence_id": pl.String,
    "start": pl.UInt64,
    "query_group": pl.String,
}
PUBLIC_CANDIDATE_COLUMNS = ["design_id", "kmer", "background_occurrences", "outcome"]
PUBLIC_OCCURRENCE_COLUMNS = [
    "design_id",
    "kmer",
    "source_sequence_id",
    "start",
    "query_group",
]
class CandidateReportError(ValueError):
    """Raised when private candidate tables cannot form the public report."""


@dataclass(frozen=True)
class DesignTables:
    """Private candidate tables for one design."""

    design_id: str
    candidates: Path
    occurrences: Path


@dataclass(frozen=True)
class CandidateReport:
    """Run-wide public candidate relations."""

    candidates: pl.LazyFrame
    occurrences: pl.LazyFrame


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-ids", type=Path, required=True)
    parser.add_argument("--candidate-manifests-dir", type=Path, required=True)
    parser.add_argument("--candidate-occurrences-dir", type=Path, required=True)
    parser.add_argument("--candidates-out", type=Path, required=True)
    parser.add_argument("--occurrences-out", type=Path, required=True)
    return parser.parse_args(argv)


def read_design_ids(path: Path) -> tuple[str, ...]:
    design_ids = tuple(path.read_text().splitlines())
    if not design_ids:
        message = "Candidate reporting requires at least one design"
        raise CandidateReportError(message)
    invalid = next(
        (
            design_id
            for design_id in design_ids
            if not design_id or design_id != design_id.strip() or "\t" in design_id
        ),
        None,
    )
    if invalid is not None:
        message = "Candidate report design IDs must be nonblank canonical TSV values"
        raise CandidateReportError(message)
    if len(set(design_ids)) != len(design_ids):
        message = "Candidate report design IDs must be unique"
        raise CandidateReportError(message)
    return design_ids


def discover_staged_tables(root: Path, filename: str) -> tuple[Path, ...]:
    directories = tuple(sorted(path for path in root.iterdir() if path.is_dir()))
    tables = tuple(directory / filename for directory in directories)
    missing = next((table for table in tables if not table.is_file()), None)
    if missing is not None:
        message = f"Candidate report input is missing {missing}"
        raise CandidateReportError(message)
    return tables


def pair_design_tables(
    design_ids: Sequence[str],
    candidate_tables: Sequence[Path],
    occurrence_tables: Sequence[Path],
) -> tuple[DesignTables, ...]:
    if len(design_ids) != len(candidate_tables) or len(design_ids) != len(occurrence_tables):
        message = "Candidate report design and table counts disagree"
        raise CandidateReportError(message)
    return tuple(
        DesignTables(design_id, candidates, occurrences)
        for design_id, candidates, occurrences in zip(
            design_ids,
            candidate_tables,
            occurrence_tables,
            strict=True,
        )
    )


def scan_exact_tsv(
    path: Path,
    schema: Mapping[str, type[pl.DataType]],
    label: str,
) -> pl.LazyFrame:
    frame = pl.scan_csv(
        path,
        separator="\t",
        schema=dict(schema),
        quote_char=None,
    )
    if frame.collect_schema().names() != list(schema):
        message = f"{label} has an unexpected schema: {path}"
        raise CandidateReportError(message)
    return frame


def report_one_design(tables: DesignTables) -> CandidateReport:
    candidates = scan_exact_tsv(
        tables.candidates,
        PRIVATE_CANDIDATE_SCHEMA,
        "Private candidate table",
    )
    occurrences = scan_exact_tsv(
        tables.occurrences,
        PRIVATE_OCCURRENCE_SCHEMA,
        "Private candidate occurrence table",
    )

    public_candidates = candidates.select(
        pl.lit(tables.design_id).alias("design_id"),
        "kmer",
        "background_occurrences",
        "outcome",
    )
    public_occurrences = occurrences.join(
        candidates.select("candidate_kmer_id", "kmer"),
        on="candidate_kmer_id",
        how="left",
    )
    unknown_candidate = (
        public_occurrences.filter(pl.col("kmer").is_null())
        .select("candidate_kmer_id")
        .sort("candidate_kmer_id")
        .limit(1)
        .collect()
    )
    if not unknown_candidate.is_empty():
        message = (
            f"Candidate occurrence table for {tables.design_id} contains an unknown candidate "
            f"identifier {unknown_candidate.item(0, 'candidate_kmer_id')}"
        )
        raise CandidateReportError(message)
    return CandidateReport(
        candidates=public_candidates,
        occurrences=public_occurrences.select(
            pl.lit(tables.design_id).alias("design_id"),
            "kmer",
            "source_sequence_id",
            "start",
            pl.col("query_group").fill_null(""),
        ),
    )


def construct_report(inputs: Sequence[DesignTables]) -> CandidateReport:
    reports = tuple(report_one_design(tables) for tables in inputs)
    return CandidateReport(
        candidates=pl.concat([report.candidates for report in reports], how="vertical").sort(
            "design_id",
            "kmer",
        ),
        occurrences=pl.concat([report.occurrences for report in reports], how="vertical").sort(
            "design_id",
            "kmer",
            "source_sequence_id",
            "start",
            "query_group",
        ),
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    inputs = pair_design_tables(
        read_design_ids(args.design_ids),
        discover_staged_tables(args.candidate_manifests_dir, "candidate_kmers.tsv"),
        discover_staged_tables(
            args.candidate_occurrences_dir,
            "candidate_kmer_occurrences.tsv",
        ),
    )
    report = construct_report(inputs)
    candidates, occurrences = pl.collect_all([report.candidates, report.occurrences])
    candidates.select(PUBLIC_CANDIDATE_COLUMNS).write_csv(
        args.candidates_out,
        separator="\t",
        quote_style="never",
    )
    occurrences.select(PUBLIC_OCCURRENCE_COLUMNS).write_csv(
        args.occurrences_out,
        separator="\t",
        quote_style="never",
    )


if __name__ == "__main__":
    main()
