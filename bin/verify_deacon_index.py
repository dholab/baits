#!/usr/bin/env python3
"""Verify a Deacon index against its bait set and interference background."""

from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from Bio import SeqIO

if TYPE_CHECKING:
    from collections.abc import Sequence

    from Bio.SeqRecord import SeqRecord

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
BAIT_SET_METRICS = [
    "design_id",
    "source_sequence_origin",
    "candidate_kmer_count",
    "locally_filtered_bait_count",
    "taxonomic_screening_status",
    "taxonomically_screened_bait_count",
    "deepest_bait_set",
    "deacon_index_source",
]
FASTA_SCHEMA = {
    "record_id": pl.String,
    "sequence": pl.String,
}
CONCLUSION = (
    "The Deacon index reproduces the bait set and retains no interference background records."
)
MAX_COUNT = 2**63 - 1


class DeaconVerificationError(ValueError):
    """Raised when Deacon verification evidence violates its contract."""


class SourceSequenceOrigin(StrEnum):
    """Provenance of the sequences that yielded candidate k-mers."""

    CURATED_INPUT = "curated_input"
    QUERY_GUIDED_DISCOVERY = "query_guided_discovery"


class TaxonomicScreeningStatus(StrEnum):
    """Taxonomic screening state of the selected bait set."""

    NOT_RUN = "NOT_RUN"
    SCREENED = "PASS"


class BaitSetSource(StrEnum):
    """Bait set selected to build the Deacon index."""

    LOCALLY_FILTERED = "locally_filtered"
    TAXONOMICALLY_SCREENED = "taxonomically_screened"


@dataclass(frozen=True)
class BaitSetStatus:
    """Validated status before Deacon index verification."""

    design_id: str
    source_sequence_origin: SourceSequenceOrigin
    candidate_kmer_count: int
    locally_filtered_bait_count: int
    taxonomic_screening_status: TaxonomicScreeningStatus
    taxonomically_screened_bait_count: int | None
    deepest_bait_set: BaitSetSource
    deacon_index_source: BaitSetSource | None


@dataclass(frozen=True)
class VerificationResult:
    """Evidence from one successful Deacon index verification."""

    design_id: str
    bait_set_source: BaitSetSource
    kmer_size: int
    deacon_window: int
    bait_count: int
    roundtrip_count: int


@dataclass(frozen=True)
class VerificationRelations:
    """Typed rectangular evidence consumed by index verification."""

    baits: pl.LazyFrame
    roundtrip: pl.LazyFrame
    background: pl.LazyFrame
    manifest: pl.LazyFrame


@dataclass(frozen=True)
class RoundTripDigests:
    """Byte identities retained by the existing output contract."""

    bait: str
    background: str


def positive_integer(value: str) -> int:
    """Parse one canonical positive integer for argparse."""
    if re.fullmatch(r"[1-9][0-9]*", value) is None:
        message = "must be a positive integer"
        raise argparse.ArgumentTypeError(message)
    return int(value)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--kmer-size", type=positive_integer, required=True)
    parser.add_argument("--deacon-window", type=positive_integer, required=True)
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--bait-roundtrip", type=Path, required=True)
    parser.add_argument("--background-roundtrip", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--bait-set-status-in", type=Path, required=True)
    parser.add_argument("--bait-set-status-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    return parser.parse_args(argv)


def construct_fasta_relation(
    records: Sequence[SeqRecord],
    label: str,
    expected_length: int | None,
    *,
    allow_empty: bool = False,
) -> pl.LazyFrame:
    """Construct one valid ordered FASTA relation from Biopython records."""
    if not records:
        if allow_empty:
            return pl.LazyFrame(schema=FASTA_SCHEMA)
        message = f"{label} must not be empty"
        raise DeaconVerificationError(message)

    record_ids = tuple(record.description for record in records)
    sequences = tuple(str(record.seq) for record in records)
    if any(not record_id for record_id in record_ids):
        message = f"{label} has an empty record identifier"
        raise DeaconVerificationError(message)
    if len(set(record_ids)) != len(record_ids):
        message = f"{label} record identifiers must be unique"
        raise DeaconVerificationError(message)
    if len(set(sequences)) != len(sequences):
        message = f"{label} sequences must be unique"
        raise DeaconVerificationError(message)

    invalid_sequence = next(
        (
            (record_id, sequence)
            for record_id, sequence in zip(record_ids, sequences, strict=True)
            if not sequence
            or sequence != sequence.upper()
            or bool(set(sequence) - set("ACGT"))
            or (expected_length is not None and len(sequence) != expected_length)
        ),
        None,
    )
    if invalid_sequence is not None:
        record_id, _ = invalid_sequence
        expected = f" {expected_length}-mer" if expected_length is not None else " DNA sequence"
        message = f"{label} has an invalid A/C/G/T{expected}: {record_id}"
        raise DeaconVerificationError(message)
    return pl.LazyFrame(
        {
            "record_id": record_ids,
            "sequence": sequences,
        },
        schema=FASTA_SCHEMA,
    )


def read_fasta(
    path: Path,
    label: str,
    expected_length: int | None,
    *,
    allow_empty: bool = False,
) -> pl.LazyFrame:
    """Parse FASTA with Biopython and construct its valid relation."""
    text = path.read_text()
    first_content = next((line for line in text.splitlines() if line.strip()), "")
    if first_content and not first_content.startswith(">"):
        message = f"{label} contains sequence text before the first header"
        raise DeaconVerificationError(message)
    try:
        records = tuple(SeqIO.parse(path, "fasta"))
    except ValueError as error:
        message = f"Could not parse {label}: {error}"
        raise DeaconVerificationError(message) from error
    return construct_fasta_relation(records, label, expected_length, allow_empty=allow_empty)


def _table_header(path: Path, label: str) -> list[str]:
    try:
        return pl.scan_csv(path, separator="\t", infer_schema_length=0).collect_schema().names()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse {label}: {error}"
        raise DeaconVerificationError(message) from error


def scan_manifest(path: Path, kmer_size: int) -> pl.LazyFrame:
    """Scan and construct the valid candidate k-mer manifest."""
    label = "candidate k-mer manifest"
    if _table_header(path, label) != MANIFEST_FIELDS:
        message = f"{label.capitalize()} has an unexpected schema"
        raise DeaconVerificationError(message)
    try:
        manifest = pl.scan_csv(
            path,
            separator="\t",
            schema_overrides=dict.fromkeys(MANIFEST_FIELDS, pl.String),
        ).with_row_index("_manifest_order", offset=1)
        manifest.collect_schema()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse {label}: {error}"
        raise DeaconVerificationError(message) from error

    source_count = pl.col("source_copy_count").cast(pl.Int64, strict=False)
    background_count = pl.col("background_occurrences").cast(pl.Int64, strict=False)
    on_target = pl.col("on_target_hits").cast(pl.Int64, strict=False)
    off_target = pl.col("off_target_hits").cast(pl.Int64, strict=False)
    expected_candidate_id = pl.concat_str(
        pl.lit("candidate_kmer_"),
        pl.col("_manifest_order").cast(pl.String).str.pad_start(6, "0"),
    )
    kmer = pl.col("kmer")
    reverse_complement = (
        kmer.str.reverse()
        .str.replace_many(["A", "C", "G", "T"], ["T", "G", "C", "A"])
    )
    empty_hits = pl.all_horizontal(
        pl.col("on_target_hits").fill_null("") == "",
        pl.col("off_target_hits").fill_null("") == "",
    )
    numeric_hits = pl.all_horizontal(
        pl.col("on_target_hits").fill_null("") != "",
        pl.col("off_target_hits").fill_null("") != "",
        on_target.is_not_null(),
        off_target.is_not_null(),
        on_target >= 0,
        off_target >= 0,
    )
    locally_passing = (
        (pl.col("status") == "PASS")
        & (pl.col("rejection_reason") == "none")
        & (pl.col("taxonomic_screening_status") == "NOT_RUN")
        & (pl.col("bait_id").fill_null("") != "")
        & (background_count == 0)
        & empty_hits
    )
    screened_passing = (
        (pl.col("status") == "PASS")
        & (pl.col("rejection_reason") == "none")
        & (pl.col("taxonomic_screening_status") == "PASS")
        & (pl.col("bait_id").fill_null("") != "")
        & (background_count == 0)
        & numeric_hits
        & (off_target == 0)
    )
    background_rejected = (
        (pl.col("status") == "REJECT_INTERFERENCE_BACKGROUND")
        & (pl.col("rejection_reason") == "background_occurrence")
        & (pl.col("taxonomic_screening_status") == "NOT_APPLICABLE")
        & (pl.col("bait_id").fill_null("") == "")
        & (background_count > 0)
        & empty_hits
    )
    complexity_rejected = (
        (pl.col("status") == "REJECT_LOW_COMPLEXITY")
        & (pl.col("rejection_reason") == "low_complexity")
        & (pl.col("taxonomic_screening_status") == "NOT_APPLICABLE")
        & (pl.col("bait_id").fill_null("") == "")
        & (background_count == 0)
        & empty_hits
    )
    off_target_rejected = (
        (pl.col("status") == "REJECT_OFF_TARGET_HIT")
        & (pl.col("rejection_reason") == "off_target_exact_match")
        & (pl.col("taxonomic_screening_status") == "REJECT_OFF_TARGET_HIT")
        & (pl.col("bait_id").fill_null("") != "")
        & (background_count == 0)
        & numeric_hits
        & (off_target > 0)
    )
    state_is_valid = pl.any_horizontal(
        locally_passing,
        screened_passing,
        background_rejected,
        complexity_rejected,
        off_target_rejected,
    )
    invalid_row = pl.coalesce(
        pl.when(
            pl.any_horizontal(
                pl.col(field).is_null()
                for field in (
                    "candidate_kmer_id",
                    "kmer",
                    "status",
                    "rejection_reason",
                    "taxonomic_screening_status",
                )
            ),
        )
        .then(pl.lit("contains a null required field")),
        pl.when(
            (pl.col("candidate_kmer_id") == "")
            | (kmer == "")
            | ~kmer.str.contains(r"^[ACGT]+$")
            | (kmer.str.len_chars() != kmer_size),
        ).then(pl.lit("contains an invalid candidate identity or k-mer")),
        pl.when(pl.col("candidate_kmer_id") != expected_candidate_id).then(
            pl.lit("contains an invalid or non-sequential candidate k-mer identifier"),
        ),
        pl.when((kmer < kmer.shift(1)).fill_null(value=False)).then(
            pl.lit("is not sorted by k-mer"),
        ),
        pl.when(kmer > reverse_complement).then(pl.lit("contains a noncanonical k-mer")),
        pl.when(pl.len().over("kmer") > 1).then(pl.lit("contains a duplicate k-mer")),
        pl.when(
            (pl.col("bait_id").fill_null("") != "")
            & (pl.len().over("bait_id") > 1),
        ).then(pl.lit("contains a duplicate bait identifier")),
        pl.when(source_count.is_null() | (source_count <= 0)).then(
            pl.lit("contains an invalid source_copy_count"),
        ),
        pl.when(background_count.is_null() | (background_count < 0)).then(
            pl.lit("contains an invalid background_occurrences"),
        ),
        pl.when(~state_is_valid).then(pl.lit("contains an incoherent candidate state")),
    )
    typed = manifest.with_columns(
        pl.col("bait_id").fill_null(""),
        source_count.alias("source_copy_count"),
        background_count.alias("background_occurrences"),
        invalid_row.alias("_error"),
    )
    validation = typed.select(
        pl.col("_error").drop_nulls().first().alias("row_error"),
    ).collect()
    if row_error := validation.item(0, "row_error"):
        message = f"Candidate k-mer manifest {row_error}"
        raise DeaconVerificationError(message)
    return typed.select(
        "candidate_kmer_id",
        "bait_id",
        "kmer",
        "status",
        "rejection_reason",
        "taxonomic_screening_status",
    )


def _parse_count(value: str, metric: str, *, optional: bool = False) -> int | None:
    if optional and value == "":
        return None
    if re.fullmatch(r"(?:0|[1-9][0-9]*)", value) is None or int(value) > MAX_COUNT:
        message = f"Bait set status has an invalid {metric}"
        raise DeaconVerificationError(message)
    return int(value)


def read_bait_set_status(path: Path, design_id: str) -> BaitSetStatus:
    """Read and construct the valid pre-verification bait set status."""
    label = "bait set status"
    if _table_header(path, label) != ["metric", "value"]:
        message = "Bait set status has an unexpected schema"
        raise DeaconVerificationError(message)
    try:
        table = pl.read_csv(
            path,
            separator="\t",
            schema_overrides={"metric": pl.String, "value": pl.String},
        ).with_columns(pl.col("value").fill_null(""))
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse {label}: {error}"
        raise DeaconVerificationError(message) from error
    if table.get_column("metric").to_list() != BAIT_SET_METRICS:
        message = "Bait set status has unexpected metrics or metric order"
        raise DeaconVerificationError(message)
    values = dict(table.select("metric", "value").iter_rows())
    if values["design_id"] != design_id:
        message = "Bait set status design_id does not match the design"
        raise DeaconVerificationError(message)
    try:
        origin = SourceSequenceOrigin(values["source_sequence_origin"])
        screening = TaxonomicScreeningStatus(values["taxonomic_screening_status"])
        deepest = BaitSetSource(values["deepest_bait_set"])
        index_source = (
            BaitSetSource(values["deacon_index_source"])
            if values["deacon_index_source"]
            else None
        )
    except ValueError as error:
        message = "Bait set status contains an unsupported state"
        raise DeaconVerificationError(message) from error
    candidate_count = _parse_count(values["candidate_kmer_count"], "candidate_kmer_count")
    locally_filtered_count = _parse_count(
        values["locally_filtered_bait_count"],
        "locally_filtered_bait_count",
    )
    screened_count = _parse_count(
        values["taxonomically_screened_bait_count"],
        "taxonomically_screened_bait_count",
        optional=True,
    )
    if index_source is not None:
        message = "Bait set status claims a Deacon index before verification"
        raise DeaconVerificationError(message)
    local_state = (
        deepest is BaitSetSource.LOCALLY_FILTERED
        and screening is TaxonomicScreeningStatus.NOT_RUN
        and screened_count is None
    )
    screened_state = (
        deepest is BaitSetSource.TAXONOMICALLY_SCREENED
        and screening is TaxonomicScreeningStatus.SCREENED
        and screened_count is not None
        and screened_count > 0
    )
    if not (local_state or screened_state):
        message = "Bait set status is not coherent before Deacon verification"
        raise DeaconVerificationError(message)
    if candidate_count is None or locally_filtered_count is None:
        message = "Bait set status is missing a required count"
        raise DeaconVerificationError(message)
    return BaitSetStatus(
        design_id=design_id,
        source_sequence_origin=origin,
        candidate_kmer_count=candidate_count,
        locally_filtered_bait_count=locally_filtered_count,
        taxonomic_screening_status=screening,
        taxonomically_screened_bait_count=screened_count,
        deepest_bait_set=deepest,
        deacon_index_source=index_source,
    )


def construct_verification_result(
    relations: VerificationRelations,
    status: BaitSetStatus,
    *,
    kmer_size: int,
    deacon_window: int,
) -> VerificationResult:
    """Verify all relations and construct immutable Deacon evidence."""
    facts = relations.manifest.select(
        pl.len().alias("candidate_count"),
        pl.col("status").is_in(["PASS", "REJECT_OFF_TARGET_HIT"]).sum().alias(
            "locally_filtered_count",
        ),
        (pl.col("status") == "PASS").sum().alias("passing_count"),
    )
    expected_bait_count = (
        status.locally_filtered_bait_count
        if status.deepest_bait_set is BaitSetSource.LOCALLY_FILTERED
        else status.taxonomically_screened_bait_count
    )
    passing = relations.manifest.filter(pl.col("status") == "PASS")
    expected_screening = status.taxonomic_screening_status.value
    manifest_baits = passing.select(
        pl.col("bait_id").alias("record_id"),
        pl.col("kmer").alias("sequence"),
    )
    errors = pl.concat(
        [
            facts.filter(pl.col("candidate_count") != status.candidate_kmer_count)
            .with_columns(
                pl.lit("Candidate k-mer manifest count disagrees with bait set status").alias(
                    "error",
                ),
            )
            .select("error"),
            facts.filter(pl.col("locally_filtered_count") != status.locally_filtered_bait_count)
            .with_columns(
                pl.lit("Locally filtered bait count disagrees with bait set status").alias(
                    "error",
                ),
            )
            .select("error"),
            facts.filter(pl.col("passing_count") != expected_bait_count)
            .with_columns(
                pl.lit("Selected bait count disagrees with bait set status").alias("error"),
            )
            .select("error"),
            passing.filter(
                (pl.col("bait_id").fill_null("") == "")
                | (pl.col("rejection_reason") != "none")
                | (pl.col("taxonomic_screening_status") != expected_screening),
            )
            .with_columns(
                pl.lit("Passing candidate k-mer rows are not coherent with bait set status").alias(
                    "error",
                ),
            )
            .select("error"),
            manifest_baits.join(
                relations.baits.select("record_id", "sequence"),
                on=["record_id", "sequence"],
                how="anti",
            )
            .with_columns(pl.lit("Passing candidate k-mer is absent from the bait FASTA").alias("error"))
            .select("error"),
            relations.baits.select("record_id", "sequence")
            .join(manifest_baits, on=["record_id", "sequence"], how="anti")
            .with_columns(pl.lit("Bait FASTA contains no matching passing candidate k-mer").alias("error"))
            .select("error"),
            relations.baits.select("sequence")
            .join(relations.roundtrip.select("sequence"), on="sequence", how="anti")
            .with_columns(pl.lit("Bait sequence set differs after the Deacon round trip").alias("error"))
            .select("error"),
            relations.roundtrip.select("sequence")
            .join(relations.baits.select("sequence"), on="sequence", how="anti")
            .with_columns(pl.lit("Bait sequence set differs after the Deacon round trip").alias("error"))
            .select("error"),
            relations.background.with_columns(
                pl.lit(
                    "Deacon retained interference background records at the permissive threshold",
                ).alias("error"),
            ).select("error"),
        ],
        how="vertical",
    ).limit(1)
    error, bait_count, roundtrip_count = pl.collect_all(
        [
            errors,
            relations.baits.select(pl.len()),
            relations.roundtrip.select(pl.len()),
        ],
    )
    if error.height:
        raise DeaconVerificationError(error.item(0, "error"))
    return VerificationResult(
        design_id=status.design_id,
        bait_set_source=status.deepest_bait_set,
        kmer_size=kmer_size,
        deacon_window=deacon_window,
        bait_count=bait_count.item(),
        roundtrip_count=roundtrip_count.item(),
    )


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file's bytes."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def status_table(result: VerificationResult, status: BaitSetStatus) -> pl.DataFrame:
    """Render the final ordered bait set status table."""
    screened_count = (
        ""
        if status.taxonomically_screened_bait_count is None
        else str(status.taxonomically_screened_bait_count)
    )
    return pl.DataFrame(
        [
            ("design_id", status.design_id),
            ("source_sequence_origin", status.source_sequence_origin.value),
            ("candidate_kmer_count", str(status.candidate_kmer_count)),
            ("locally_filtered_bait_count", str(status.locally_filtered_bait_count)),
            ("taxonomic_screening_status", status.taxonomic_screening_status.value),
            ("taxonomically_screened_bait_count", screened_count),
            ("deepest_bait_set", status.deepest_bait_set.value),
            ("deacon_index_source", result.bait_set_source.value),
        ],
        schema={"metric": pl.String, "value": pl.String},
        orient="row",
    )


def summary_table(result: VerificationResult, digests: RoundTripDigests) -> pl.DataFrame:
    """Render the ordered machine-readable verification evidence."""
    return pl.DataFrame(
        [
            ("design_id", result.design_id),
            ("kmer_size", str(result.kmer_size)),
            ("deacon_window", str(result.deacon_window)),
            ("deacon_index_entropy_threshold", "0"),
            ("deacon_filter_absolute_threshold", "1"),
            ("deacon_filter_relative_threshold", "0"),
            ("bait_roundtrip_record_count", str(result.roundtrip_count)),
            ("bait_roundtrip_sha256", digests.bait),
            ("bait_sequence_sets_equal", "true"),
            ("interference_background_roundtrip_record_count", "0"),
            ("interference_background_roundtrip_sha256", digests.background),
            ("conclusion", CONCLUSION),
        ],
        schema={"metric": pl.String, "value": pl.String},
        orient="row",
    )


def render_report(result: VerificationResult, digests: RoundTripDigests) -> str:
    """Render the readable verification report."""
    bait_set_name = (
        "locally filtered bait set"
        if result.bait_set_source is BaitSetSource.LOCALLY_FILTERED
        else "taxonomically screened bait set"
    )
    return "\n".join(
        [
            "# Deacon index verification",
            "",
            "## Conclusion",
            "",
            CONCLUSION,
            "",
            "## Evidence",
            "",
            f"- Design: `{result.design_id}`",
            f"- Source: {bait_set_name}",
            f"- Deacon index parameters: `-k {result.kmer_size} -w {result.deacon_window} -e 0`",
            "- Deacon filter parameters: `-a 1 -r 0`",
            f"- Input bait records: {result.bait_count}",
            f"- Bait round-trip records: {result.roundtrip_count}",
            "- Bait sequence sets equal: yes",
            f"- Bait round-trip SHA-256: `{digests.bait}`",
            "- Retained interference background records: 0",
            f"- Interference background round-trip SHA-256: `{digests.background}`",
            "",
        ],
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    baits = read_fasta(args.baits, "bait FASTA", args.kmer_size)
    roundtrip = read_fasta(args.bait_roundtrip, "bait round-trip FASTA", args.kmer_size)
    background = read_fasta(
        args.background_roundtrip,
        "interference background round-trip FASTA",
        None,
        allow_empty=True,
    )
    manifest = scan_manifest(args.manifest, args.kmer_size)
    status = read_bait_set_status(args.bait_set_status_in, args.design_id)
    result = construct_verification_result(
        VerificationRelations(
            baits=baits,
            roundtrip=roundtrip,
            background=background,
            manifest=manifest,
        ),
        status,
        kmer_size=args.kmer_size,
        deacon_window=args.deacon_window,
    )
    digests = RoundTripDigests(
        bait=sha256_file(args.bait_roundtrip),
        background=sha256_file(args.background_roundtrip),
    )
    status_table(result, status).write_csv(
        args.bait_set_status_out,
        separator="\t",
        quote_style="never",
    )
    summary_table(result, digests).write_csv(
        args.summary_out,
        separator="\t",
        quote_style="never",
    )
    args.report_out.write_text(render_report(result, digests))


if __name__ == "__main__":
    main()
