#!/usr/bin/env python3
"""Apply Deacon low-complexity results to Candidate K-mers."""

from __future__ import annotations

import argparse
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from collections.abc import Sequence

MANIFEST_FIELDS = [
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
MANIFEST_SCHEMA = dict.fromkeys(MANIFEST_FIELDS, pl.String)


class ComplexityFilterError(ValueError):
    """Raised when complexity-filter inputs do not satisfy their contract."""


class SourceSequenceOrigin(StrEnum):
    """Provenance of the sequences that yielded Candidate K-mers."""

    CURATED_INPUT = "curated_input"
    QUERY_GUIDED_DISCOVERY = "query_guided_discovery"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-id", required=True)
    parser.add_argument(
        "--source-sequence-origin",
        type=SourceSequenceOrigin,
        choices=list(SourceSequenceOrigin),
        required=True,
    )
    parser.add_argument("--manifest-in", type=Path, required=True)
    parser.add_argument("--passing-kmers", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    parser.add_argument("--baits-out", type=Path, required=True)
    parser.add_argument("--filtering-status-out", type=Path, required=True)
    parser.add_argument("--bait-set-status-out", type=Path, required=True)
    return parser.parse_args(argv)


def canonical_kmer(kmer: str) -> str:
    """Return the project representation shared by a k-mer and its reverse complement."""
    return min(kmer, str(Seq(kmer).reverse_complement()))


def scan_manifest(path: Path) -> pl.LazyFrame:
    """Scan and construct the valid pre-complexity Candidate K-mer manifest."""
    try:
        header = pl.scan_csv(
            path,
            separator="\t",
            infer_schema_length=0,
        ).collect_schema().names()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Candidate manifest: {error}"
        raise ComplexityFilterError(message) from error
    if header != MANIFEST_FIELDS:
        message = "Candidate manifest has an unexpected schema"
        raise ComplexityFilterError(message)

    try:
        manifest = pl.scan_csv(
            path,
            separator="\t",
            schema_overrides=MANIFEST_SCHEMA,
        ).with_row_index("_row", offset=1)
        manifest.collect_schema()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Candidate manifest: {error}"
        raise ComplexityFilterError(message) from error

    source_copy_count = pl.col("source_copy_count").cast(pl.Int64, strict=False)
    background_occurrences = pl.col("background_occurrences").cast(pl.Int64, strict=False)
    kmer = pl.col("kmer")
    expected_candidate_id = pl.concat_str(
        pl.lit("candidate_kmer_"),
        pl.col("_row").cast(pl.String).str.pad_start(6, "0"),
    )
    is_pass = (
        (pl.col("status") == "PASS")
        & (pl.col("rejection_reason") == "none")
        & (pl.col("taxonomic_screening_status") == "NOT_RUN")
        & (background_occurrences == 0)
    )
    is_background_rejection = (
        (pl.col("status") == "REJECT_INTERFERENCE_BACKGROUND")
        & (pl.col("rejection_reason") == "background_occurrence")
        & (pl.col("taxonomic_screening_status") == "NOT_APPLICABLE")
        & (background_occurrences > 0)
    )
    evidence_fields = ("bait_id", "on_target_hits", "off_target_hits")
    empty_evidence = pl.all_horizontal([pl.col(field).fill_null("") == "" for field in evidence_fields])
    invalid_row = pl.coalesce(
        pl.when(pl.any_horizontal([pl.col(field).is_null() for field in MANIFEST_FIELDS if field not in evidence_fields]))
        .then(pl.lit("Candidate manifest contains a missing field")),
        pl.when(pl.col("candidate_kmer_id") != expected_candidate_id).then(
            pl.concat_str(pl.lit("invalid or non-sequential Candidate K-mer ID: "), pl.col("candidate_kmer_id")),
        ),
        pl.when((kmer < kmer.shift(1)).fill_null(value=False)).then(
            pl.concat_str(pl.lit("Candidate manifest is not sorted by k-mer at "), pl.col("candidate_kmer_id")),
        ),
        pl.when(~kmer.str.contains(r"^[ACGT]+$")).then(
            pl.concat_str(pl.lit("invalid Candidate K-mer: "), kmer),
        ),
        pl.when(kmer != kmer.map_elements(canonical_kmer, return_dtype=pl.String)).then(
            pl.concat_str(pl.lit("noncanonical Candidate K-mer: "), kmer),
        ),
        pl.when(pl.len().over("kmer") > 1).then(
            pl.concat_str(pl.lit("duplicate Candidate K-mer: "), kmer),
        ),
        pl.when(~pl.col("source_copy_count").str.contains(r"^[0-9]+$") | (source_copy_count <= 0)).then(
            pl.concat_str(pl.lit("invalid source_copy_count for "), pl.col("candidate_kmer_id")),
        ),
        pl.when(
            ~pl.col("background_occurrences").str.contains(r"^[0-9]+$")
            | background_occurrences.is_null(),
        ).then(pl.concat_str(pl.lit("invalid background_occurrences for "), pl.col("candidate_kmer_id"))),
        pl.when(~empty_evidence).then(
            pl.concat_str(pl.lit("pre-complexity evidence must be empty for "), pl.col("candidate_kmer_id")),
        ),
        pl.when(~(is_pass | is_background_rejection)).then(
            pl.concat_str(pl.lit("inconsistent pre-complexity state for "), pl.col("candidate_kmer_id")),
        ),
    ).alias("_invalid_row")
    parsed = manifest.with_columns(
        *[pl.col(field).fill_null("").alias(field) for field in evidence_fields],
        source_copy_count.alias("source_copy_count"),
        background_occurrences.alias("background_occurrences"),
        invalid_row,
    )
    try:
        validation = parsed.select(
            pl.len().alias("row_count"),
            (pl.col("status") == "PASS").sum().alias("eligible_count"),
            pl.col("kmer").str.len_chars().n_unique().alias("kmer_lengths"),
            pl.col("_invalid_row").drop_nulls().first().alias("first_error"),
        ).collect()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Candidate manifest: {error}"
        raise ComplexityFilterError(message) from error
    if validation.item(0, "row_count") == 0:
        message = "Complexity filtering requires at least one eligible Candidate K-mer"
        raise ComplexityFilterError(message)
    if first_error := validation.item(0, "first_error"):
        raise ComplexityFilterError(first_error)
    if validation.item(0, "eligible_count") == 0:
        message = "Complexity filtering requires at least one eligible Candidate K-mer"
        raise ComplexityFilterError(message)
    if validation.item(0, "kmer_lengths") != 1:
        message = "Candidate K-mers must have one common length"
        raise ComplexityFilterError(message)
    return parsed.drop("_invalid_row")


def read_passing_kmers(path: Path, kmer_size: int) -> pl.LazyFrame:
    """Read Deacon's FASTA dump and construct its unique canonical passing keys."""
    text = path.read_text()
    first_content = next((line for line in text.splitlines() if line.strip()), "")
    if first_content and not first_content.startswith(">"):
        message = "Deacon passing FASTA contains sequence text before the first header"
        raise ComplexityFilterError(message)
    try:
        records = tuple(SeqIO.parse(path, "fasta"))
    except ValueError as error:
        message = f"Could not parse Deacon passing FASTA: {error}"
        raise ComplexityFilterError(message) from error
    if text.strip() and not records:
        message = "Deacon passing FASTA contains malformed nonempty text"
        raise ComplexityFilterError(message)

    sequences = tuple(str(record.seq) for record in records)
    if any(sequence == "" for sequence in sequences):
        message = "Deacon passing FASTA contains an empty record"
        raise ComplexityFilterError(message)
    if invalid := next((sequence for sequence in sequences if not sequence.isupper() or set(sequence) - set("ACGT")), None):
        message = f"Deacon passing FASTA contains malformed non-ACGT k-mer: {invalid}"
        raise ComplexityFilterError(message)
    if invalid := next((sequence for sequence in sequences if len(sequence) != kmer_size), None):
        message = f"Deacon passing FASTA contains wrong-length k-mer: {invalid}"
        raise ComplexityFilterError(message)

    canonical_sequences = tuple(map(canonical_kmer, sequences))
    if len(set(canonical_sequences)) != len(canonical_sequences):
        message = "Deacon passing FASTA contains duplicate canonical k-mer"
        raise ComplexityFilterError(message)
    return pl.LazyFrame({"kmer": canonical_sequences}, schema={"kmer": pl.String})


def apply_complexity_results(manifest: pl.LazyFrame, passing_kmers: pl.LazyFrame) -> pl.LazyFrame:
    """Join Deacon's passing keys onto valid rows and assign local-filter outcomes."""
    passing_marker = passing_kmers.with_columns(pl.lit(value=True).alias("_passes_complexity"))
    transitioned = manifest.join(passing_marker, on="kmer", how="left").with_columns(
        pl.when((pl.col("status") == "PASS") & pl.col("_passes_complexity").is_null())
        .then(pl.lit("REJECT_LOW_COMPLEXITY"))
        .otherwise(pl.col("status"))
        .alias("status"),
        pl.when((pl.col("status") == "PASS") & pl.col("_passes_complexity").is_null())
        .then(pl.lit("low_complexity"))
        .otherwise(pl.col("rejection_reason"))
        .alias("rejection_reason"),
        pl.when((pl.col("status") == "PASS") & pl.col("_passes_complexity").is_null())
        .then(pl.lit("NOT_APPLICABLE"))
        .otherwise(pl.col("taxonomic_screening_status"))
        .alias("taxonomic_screening_status"),
    )
    survivors = (
        transitioned.filter(pl.col("status") == "PASS")
        .sort("kmer")
        .with_row_index("_bait_number", offset=1)
        .select(
            "kmer",
            pl.concat_str(
                pl.lit("bait_"),
                pl.col("_bait_number").cast(pl.String).str.pad_start(6, "0"),
            ).alias("bait_id"),
        )
    )
    return (
        transitioned.drop("_passes_complexity")
        .join(survivors, on="kmer", how="left", suffix="_survivor")
        .with_columns(pl.coalesce("bait_id_survivor", "bait_id").alias("bait_id"))
        .drop("bait_id_survivor", "_row")
        .select(MANIFEST_FIELDS)
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    manifest = scan_manifest(args.manifest_in)
    kmer_size = manifest.select(pl.col("kmer").str.len_chars().first()).collect().item()
    passing_kmers = read_passing_kmers(args.passing_kmers, kmer_size)
    eligible = manifest.filter(pl.col("status") == "PASS").select("kmer")
    unexpected = passing_kmers.join(eligible, on="kmer", how="anti").limit(1).collect()
    if unexpected.height:
        message = "Deacon passing sequences are not eligible Candidate K-mers"
        raise ComplexityFilterError(message)

    results = apply_complexity_results(manifest, passing_kmers)
    manifest_output, survivors = pl.collect_all(
        [
            results.sort("candidate_kmer_id"),
            results.filter(pl.col("status") == "PASS").sort("kmer"),
        ],
    )
    manifest_output.write_csv(args.manifest_out, separator="\t", quote_style="never")

    counts = manifest_output.select(
        pl.len().alias("candidate_kmer_count"),
        (pl.col("status") == "REJECT_INTERFERENCE_BACKGROUND")
        .sum()
        .alias("rejected_interference_background_count"),
        (pl.col("status") == "REJECT_LOW_COMPLEXITY").sum().alias("rejected_low_complexity_count"),
        (pl.col("status") == "PASS").sum().alias("locally_filtered_bait_count"),
    )
    count_values = {column: str(counts.item(0, column)) for column in counts.columns}
    terminal_stage = "" if survivors.height else "low_complexity_filtering"
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
            "value": [
                args.design_id,
                count_values["candidate_kmer_count"],
                count_values["rejected_interference_background_count"],
                count_values["rejected_low_complexity_count"],
                count_values["locally_filtered_bait_count"],
                terminal_stage,
            ],
        },
        schema={"metric": pl.String, "value": pl.String},
    ).write_csv(args.filtering_status_out, separator="\t", quote_style="never")
    if survivors.height:
        SeqIO.write(
            (
                SeqRecord(Seq(kmer), id=bait_id, description="")
                for bait_id, kmer in survivors.select("bait_id", "kmer").iter_rows()
            ),
            args.baits_out,
            "fasta",
        )
        pl.DataFrame(
            {
                "metric": [
                    "design_id",
                    "source_sequence_origin",
                    "candidate_kmer_count",
                    "locally_filtered_bait_count",
                    "taxonomic_screening_status",
                    "taxonomically_screened_bait_count",
                    "deepest_bait_set",
                    "deacon_index_source",
                ],
                "value": [
                    args.design_id,
                    args.source_sequence_origin.value,
                    count_values["candidate_kmer_count"],
                    count_values["locally_filtered_bait_count"],
                    "NOT_RUN",
                    "",
                    "locally_filtered",
                    "",
                ],
            },
            schema={"metric": pl.String, "value": pl.String},
        ).write_csv(args.bait_set_status_out, separator="\t", quote_style="never")


if __name__ == "__main__":
    main()
