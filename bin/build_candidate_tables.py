#!/usr/bin/env python3
"""Build canonical Candidate K-mer evidence from Source Sequences."""

from __future__ import annotations

import argparse
import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from collections.abc import Sequence

MANIFEST_COLUMNS = [
    "candidate_kmer_id",
    "bait_id",
    "kmer",
    "source_copy_count",
    "background_occurrences",
    "status",
    "rejection_reason",
    "taxonomic_screening_status",
    "on_target_hits",
    "off_target_hits",
]
OCCURRENCE_COLUMNS = ["candidate_kmer_id", "source_sequence_id", "start", "query_group"]
QUERY_GROUP_SCHEMA = {"source_sequence_id": pl.String, "query_group": pl.String}
MERYL_COUNT_SCHEMA = {"kmer": pl.String, "count": pl.Int64}
DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")


class SourceSequenceError(ValueError):
    """Raised when Source Sequences are not valid input."""


class QueryGroupError(ValueError):
    """Raised when Source Sequence Query Groups are not valid input."""


class MerylCountError(ValueError):
    """Raised when a Meryl printout is not valid count evidence."""


class CandidateEvidenceError(ValueError):
    """Raised when independent Candidate evidence disagrees."""


@dataclass(frozen=True)
class CandidateEvidence:
    """The manifest, Source occurrences, and candidates for complexity filtering."""

    manifest: pl.LazyFrame
    occurrences: pl.LazyFrame
    complexity_candidates: pl.LazyFrame


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sequences", type=Path, required=True)
    parser.add_argument("--source-sequence-query-groups", type=Path, required=True)
    parser.add_argument("--meryl-source-counts", type=Path, required=True)
    parser.add_argument("--background-intersection-counts", type=Path, required=True)
    parser.add_argument("--kmer-size", type=int, required=True)
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    parser.add_argument("--occurrences-out", type=Path, required=True)
    parser.add_argument("--complexity-candidates-out", type=Path, required=True)
    parser.add_argument("--filtering-status-out", type=Path, required=True)
    parser.add_argument("--terminal-manifest-out", type=Path, required=True)
    return parser.parse_args(argv)


def open_fasta(path: Path) -> TextIO:
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_source_sequences(path: Path) -> pl.LazyFrame:
    with open_fasta(path) as handle:
        records = tuple(SeqIO.parse(handle, "fasta"))
    if not records:
        message = "Source Sequence FASTA must contain at least one record"
        raise SourceSequenceError(message)

    identifiers = tuple(record.id for record in records)
    blank_identifier = next((identifier for identifier in identifiers if not identifier), None)
    if blank_identifier is not None:
        message = "Source Sequence FASTA record ID must not be blank"
        raise SourceSequenceError(message)
    duplicate_identifier = next(
        (
            identifier
            for identifier, count in sorted(Counter(identifiers).items())
            if count > 1
        ),
        None,
    )
    if duplicate_identifier is not None:
        message = f"Duplicate Source Sequence FASTA record ID: {duplicate_identifier}"
        raise SourceSequenceError(message)

    raw_sequences = tuple(str(record.seq) for record in records)
    non_normalized = next(
        (
            record.id
            for record, sequence in zip(records, raw_sequences, strict=True)
            if sequence != sequence.upper()
        ),
        None,
    )
    if non_normalized is not None:
        message = f"Source Sequence FASTA must contain uppercase bases: {non_normalized}"
        raise SourceSequenceError(message)
    sequences = raw_sequences
    invalid_sequence = next(
        (
            (record.id, sequence)
            for record, sequence in zip(records, sequences, strict=True)
            if not set(sequence) <= DNA_IUPAC
        ),
        None,
    )
    if invalid_sequence is not None:
        identifier, sequence = invalid_sequence
        invalid_bases = "".join(sorted(set(sequence) - DNA_IUPAC))
        message = (
            f"Source Sequence FASTA contains malformed DNA/IUPAC sequence for {identifier}: "
            f"{invalid_bases or sequence}"
        )
        raise SourceSequenceError(message)

    return pl.LazyFrame(
        {"source_sequence_id": identifiers, "sequence": sequences},
        schema={"source_sequence_id": pl.String, "sequence": pl.String},
    )


def scan_query_groups(path: Path, source_sequences: pl.LazyFrame) -> pl.LazyFrame:
    try:
        query_groups = pl.scan_csv(path, separator="\t", schema_overrides=QUERY_GROUP_SCHEMA)
        columns = query_groups.collect_schema().names()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Source Sequence Query Groups: {error}"
        raise QueryGroupError(message) from error
    if columns != list(QUERY_GROUP_SCHEMA):
        message = "Source Sequence Query Groups columns must be exactly: source_sequence_id, query_group"
        raise QueryGroupError(message)

    contextualized = query_groups.with_row_index("_row").with_columns(
        pl.when(pl.col("source_sequence_id").is_null() | (pl.col("source_sequence_id") == ""))
        .then(pl.lit("blank source_sequence_id in Source Sequence Query Groups"))
        .when(pl.len().over("source_sequence_id") > 1)
        .then(
            pl.concat_str(
                pl.lit("duplicate source_sequence_id in Source Sequence Query Groups: "),
                pl.col("source_sequence_id"),
            ),
        )
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias("_error"),
    )
    validation = contextualized.select(
        pl.col("_error").drop_nulls().first().alias("row_error"),
    ).collect()
    if row_error := validation.item(0, "row_error"):
        raise QueryGroupError(row_error)

    missing = (
        source_sequences.join(contextualized, on="source_sequence_id", how="anti")
        .select("source_sequence_id")
        .sort("source_sequence_id")
        .limit(1)
        .collect()
    )
    if missing.height:
        message = f"Missing Source Sequence Query Group ID: {missing.item(0, 'source_sequence_id')}"
        raise QueryGroupError(message)
    extra = (
        contextualized.join(source_sequences, on="source_sequence_id", how="anti")
        .select("source_sequence_id")
        .sort("source_sequence_id")
        .limit(1)
        .collect()
    )
    if extra.height:
        message = f"Unknown Source Sequence Query Group ID: {extra.item(0, 'source_sequence_id')}"
        raise QueryGroupError(message)
    return contextualized.with_columns(
        pl.col("query_group").fill_null(""),
    ).drop("_row", "_error")


def project_canonical(kmer: str) -> str:
    reverse_complement = str(Seq(kmer).reverse_complement())
    return min(kmer, reverse_complement)


def scan_meryl_counts(path: Path, kmer_size: int, evidence_name: str) -> pl.LazyFrame:
    if path.stat().st_size == 0:
        return pl.LazyFrame(schema=MERYL_COUNT_SCHEMA)
    try:
        counts = pl.scan_csv(
            path,
            separator="\t",
            has_header=False,
            quote_char=None,
            schema=MERYL_COUNT_SCHEMA,
        ).with_row_index("_row")
        malformed = counts.with_columns(
            pl.when(
                pl.col("kmer").is_null()
                | ~pl.col("kmer").str.contains(r"^[ACGT]+$")
                | (pl.col("kmer").str.len_chars() != kmer_size),
            )
            .then(pl.concat_str(pl.lit("malformed k-mer in "), pl.lit(evidence_name)))
            .when(pl.col("count").is_null() | (pl.col("count") <= 0))
            .then(pl.concat_str(pl.lit("nonpositive count in "), pl.lit(evidence_name)))
            .otherwise(pl.lit(None, dtype=pl.String))
            .alias("_error"),
        )
        validation = malformed.select(pl.col("_error").drop_nulls().first().alias("error")).collect()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse {evidence_name}: {error}"
        raise MerylCountError(message) from error
    if error := validation.item(0, "error"):
        raise MerylCountError(error)

    canonicalized = malformed.filter(pl.col("_error").is_null()).with_columns(
        pl.col("kmer").map_elements(project_canonical, return_dtype=pl.String).alias("kmer"),
    )
    duplicate = (
        canonicalized.group_by("kmer")
        .len()
        .filter(pl.col("len") > 1)
        .sort("kmer")
        .select("kmer")
        .limit(1)
        .collect()
    )
    if duplicate.height:
        message = f"Duplicate canonical k-mer in {evidence_name}: {duplicate.item(0, 'kmer')}"
        raise MerylCountError(message)
    return canonicalized.select("kmer", "count")


def construct_candidate_evidence(
    source_sequences: pl.LazyFrame,
    query_groups: pl.LazyFrame,
    meryl_source_counts: pl.LazyFrame,
    background_intersection_counts: pl.LazyFrame,
    kmer_size: int,
) -> CandidateEvidence:
    """Construct Candidate identity, Source occurrences, and cancellation evidence."""
    contextualized_sources = source_sequences.join(query_groups, on="source_sequence_id", how="inner")
    occurrences = (
        contextualized_sources.with_columns(
            pl.int_ranges(
                0,
                (
                    pl.col("sequence").str.len_chars().cast(pl.Int64)
                    - kmer_size
                    + 1
                ).clip(lower_bound=0),
            ).alias("start"),
        )
        .explode("start", empty_as_null=False)
        .with_columns(pl.col("sequence").str.slice(pl.col("start"), kmer_size).alias("_window"))
        .filter(pl.col("_window").str.contains(r"^[ACGT]+$"))
        .with_columns(
            pl.col("_window").map_elements(project_canonical, return_dtype=pl.String).alias("kmer"),
        )
        .select("source_sequence_id", "start", "query_group", "kmer")
    )
    candidates = occurrences.group_by("kmer").len().rename({"len": "source_copy_count"})

    mismatched_source_count = (
        candidates.join(meryl_source_counts, on="kmer", how="full", coalesce=True)
        .filter(
            pl.col("source_copy_count").is_null()
            | pl.col("count").is_null()
            | (pl.col("source_copy_count") != pl.col("count")),
        )
        .sort("kmer")
        .select("kmer")
        .limit(1)
        .collect()
    )
    if mismatched_source_count.height:
        message = (
            "Meryl Source counts must equal Python Candidate keys and positional occurrence "
            f"counts: {mismatched_source_count.item(0, 'kmer')}"
        )
        raise CandidateEvidenceError(message)

    unexpected_background = (
        background_intersection_counts.join(candidates, on="kmer", how="anti")
        .sort("kmer")
        .select("kmer")
        .limit(1)
        .collect()
    )
    if unexpected_background.height:
        message = (
            "Background intersection contains non-Candidate K-mer: "
            f"{unexpected_background.item(0, 'kmer')}"
        )
        raise CandidateEvidenceError(message)

    numbered_candidates = (
        candidates.sort("kmer")
        .with_row_index("_candidate_number", offset=1)
        .with_columns(
            pl.concat_str(
                pl.lit("candidate_kmer_"),
                pl.col("_candidate_number").cast(pl.String).str.pad_start(6, "0"),
            ).alias("candidate_kmer_id"),
        )
        .drop("_candidate_number")
    )
    manifest = (
        numbered_candidates.join(background_intersection_counts, on="kmer", how="left")
        .with_columns(pl.col("count").fill_null(0).alias("background_occurrences"))
        .with_columns(
            pl.lit("").alias("bait_id"),
            pl.when(pl.col("background_occurrences") > 0)
            .then(pl.lit("REJECT_INTERFERENCE_BACKGROUND"))
            .otherwise(pl.lit("PASS"))
            .alias("status"),
            pl.when(pl.col("background_occurrences") > 0)
            .then(pl.lit("background_occurrence"))
            .otherwise(pl.lit("none"))
            .alias("rejection_reason"),
            pl.when(pl.col("background_occurrences") > 0)
            .then(pl.lit("NOT_APPLICABLE"))
            .otherwise(pl.lit("NOT_RUN"))
            .alias("taxonomic_screening_status"),
            pl.lit("").alias("on_target_hits"),
            pl.lit("").alias("off_target_hits"),
        )
        .select(MANIFEST_COLUMNS)
        .sort("kmer")
    )
    return CandidateEvidence(
        manifest=manifest,
        occurrences=(
            occurrences.join(numbered_candidates.select("candidate_kmer_id", "kmer"), on="kmer")
            .select(OCCURRENCE_COLUMNS)
            .sort("candidate_kmer_id", "source_sequence_id", "start", "query_group")
        ),
        complexity_candidates=(
            manifest.filter(pl.col("status") == "PASS").select("candidate_kmer_id", "kmer").sort("kmer")
        ),
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.kmer_size <= 0:
        message = "kmer_size must be positive"
        raise ValueError(message)

    source_sequences = read_source_sequences(args.source_sequences)
    query_groups = scan_query_groups(args.source_sequence_query_groups, source_sequences)
    meryl_source_counts = scan_meryl_counts(args.meryl_source_counts, args.kmer_size, "Meryl Source counts")
    background_counts = scan_meryl_counts(
        args.background_intersection_counts,
        args.kmer_size,
        "Background intersection counts",
    )
    evidence = construct_candidate_evidence(
        source_sequences,
        query_groups,
        meryl_source_counts,
        background_counts,
        args.kmer_size,
    )

    manifest, occurrences, complexity_candidates = pl.collect_all(
        [
            evidence.manifest,
            evidence.occurrences,
            evidence.complexity_candidates,
        ],
    )
    manifest.write_csv(args.manifest_out, separator="\t", quote_style="never")
    occurrences.write_csv(args.occurrences_out, separator="\t", quote_style="never")
    if complexity_candidates.height:
        SeqIO.write(
            (
                SeqRecord(Seq(kmer), id=candidate_kmer_id, description="")
                for candidate_kmer_id, kmer in complexity_candidates.iter_rows()
            ),
            args.complexity_candidates_out,
            "fasta",
        )
        return

    terminal_stage = "candidate_kmer_enumeration" if manifest.is_empty() else "explicit_background_cancellation"
    pl.DataFrame(
        {
            "metric": [
                "design_id",
                "candidate_kmer_count",
                "rejected_interference_background_count",
                "rejected_low_complexity_count",
                "locally_filtered_bait_count",
                "terminal_stage",
            ],
            "value": [args.design_id, str(manifest.height), str(manifest.height), "0", "0", terminal_stage],
        },
        schema={"metric": pl.String, "value": pl.String},
    ).write_csv(args.filtering_status_out, separator="\t", quote_style="never")
    manifest.write_csv(args.terminal_manifest_out, separator="\t", quote_style="never")


if __name__ == "__main__":
    main()
