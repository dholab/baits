#!/usr/bin/env python3
"""Apply exact-match taxonomic evidence to locally filtered Baits."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
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
BLAST_FIELDS = [
    "qseqid",
    "saccver",
    "staxids",
    "sscinames",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "stitle",
]
DECISION_FIELDS = [
    "bait_id",
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
EXACT_PERCENT_IDENTITY = 100.0


class TaxonomicScreeningError(ValueError):
    """Raised when taxonomic-screening evidence violates its contract."""


class SourceSequenceOrigin(StrEnum):
    """Provenance of the sequences that yielded Candidate K-mers."""

    CURATED_INPUT = "curated_input"
    QUERY_GUIDED_DISCOVERY = "query_guided_discovery"


@dataclass(frozen=True)
class BaitSetStatus:
    """Validated status of the locally filtered Bait Set."""

    design_id: str
    source_sequence_origin: SourceSequenceOrigin
    candidate_kmer_count: int
    locally_filtered_bait_count: int
    deacon_index_source: str


@dataclass(frozen=True)
class TaxonomicScreeningResult:
    """Lazy projections of one validated taxonomic-screening decision."""

    manifest: pl.LazyFrame
    decisions: pl.LazyFrame
    survivors: pl.LazyFrame


def positive_integer(value: str) -> int:
    """Parse one canonical positive integer for argparse."""
    if re.fullmatch(r"[1-9][0-9]*", value) is None:
        message = "must be a positive integer"
        raise argparse.ArgumentTypeError(message)
    return int(value)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--target-taxid", type=positive_integer, required=True)
    parser.add_argument("--kmer-size", type=positive_integer, required=True)
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--blast-hits", type=Path, required=True)
    parser.add_argument("--manifest-in", type=Path, required=True)
    parser.add_argument("--bait-set-status-in", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    parser.add_argument("--baits-out", type=Path, required=True)
    parser.add_argument("--decisions-out", type=Path, required=True)
    parser.add_argument("--screening-status-out", type=Path, required=True)
    parser.add_argument("--bait-set-status-out", type=Path, required=True)
    parser.add_argument("--terminal-bait-set-status-out", type=Path)
    return parser.parse_args(argv)


def canonical_kmer(kmer: str) -> str:
    """Return the project representation shared by a k-mer and its reverse complement."""
    return min(kmer, str(Seq(kmer).reverse_complement()))


def _header(path: Path, label: str) -> list[str]:
    try:
        return pl.scan_csv(path, separator="\t", infer_schema_length=0).collect_schema().names()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse {label}: {error}"
        raise TaxonomicScreeningError(message) from error


def construct_baits(records: Sequence[SeqRecord], kmer_size: int) -> pl.LazyFrame:
    """Construct the ordered, canonical Bait relation from parsed FASTA records."""
    if not records:
        message = "Locally filtered Bait FASTA must not be empty"
        raise TaxonomicScreeningError(message)

    identifiers = tuple(record.id for record in records)
    sequences = tuple(str(record.seq) for record in records)
    invalid_description = next(
        (
            record.description
            for record in records
            if record.id == "" or record.description != record.id
        ),
        None,
    )
    if invalid_description is not None:
        message = f"Bait FASTA has an invalid identifier or description: {invalid_description}"
        raise TaxonomicScreeningError(message)
    if len(set(identifiers)) != len(identifiers):
        message = "Bait FASTA record IDs must be unique"
        raise TaxonomicScreeningError(message)
    if len(set(sequences)) != len(sequences):
        message = "Bait FASTA sequences must be unique"
        raise TaxonomicScreeningError(message)
    invalid_sequence = next(
        (
            sequence
            for sequence in sequences
            if len(sequence) != kmer_size
            or not sequence.isupper()
            or set(sequence) - set("ACGT")
        ),
        None,
    )
    if invalid_sequence is not None:
        message = f"Bait FASTA contains a malformed {kmer_size}-mer: {invalid_sequence}"
        raise TaxonomicScreeningError(message)
    if invalid := next(
        (sequence for sequence in sequences if sequence != canonical_kmer(sequence)),
        None,
    ):
        message = f"Bait FASTA contains a noncanonical k-mer: {invalid}"
        raise TaxonomicScreeningError(message)

    expected_identifiers = tuple(f"bait_{number:06d}" for number in range(1, len(records) + 1))
    if identifiers != expected_identifiers:
        message = "Bait FASTA IDs must be sequential in k-mer order"
        raise TaxonomicScreeningError(message)
    if sequences != tuple(sorted(sequences)):
        message = "Bait FASTA is not sorted by k-mer"
        raise TaxonomicScreeningError(message)
    return pl.LazyFrame(
        {
            "_bait_order": range(1, len(records) + 1),
            "bait_id": identifiers,
            "kmer": sequences,
        },
        schema={"_bait_order": pl.UInt32, "bait_id": pl.String, "kmer": pl.String},
    )


def read_baits(path: Path, kmer_size: int) -> pl.LazyFrame:
    """Read the locally filtered Bait FASTA and construct its domain relation."""
    text = path.read_text()
    first_content = next((line for line in text.splitlines() if line.strip()), "")
    if first_content and not first_content.startswith(">"):
        message = "Bait FASTA contains sequence text before the first header"
        raise TaxonomicScreeningError(message)
    try:
        records = tuple(SeqIO.parse(path, "fasta"))
    except ValueError as error:
        message = f"Could not parse Bait FASTA: {error}"
        raise TaxonomicScreeningError(message) from error
    return construct_baits(records, kmer_size)


def scan_manifest(path: Path, kmer_size: int) -> pl.LazyFrame:
    """Scan and construct the valid pre-screening Candidate K-mer manifest."""
    if _header(path, "Candidate manifest") != MANIFEST_FIELDS:
        message = "Candidate manifest has an unexpected schema"
        raise TaxonomicScreeningError(message)
    try:
        manifest = pl.scan_csv(
            path,
            separator="\t",
            schema_overrides=dict.fromkeys(MANIFEST_FIELDS, pl.String),
        ).with_row_index("_manifest_order", offset=1)
        manifest.collect_schema()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Candidate manifest: {error}"
        raise TaxonomicScreeningError(message) from error

    source_count = pl.col("source_copy_count").cast(pl.Int64, strict=False)
    background_count = pl.col("background_occurrences").cast(pl.Int64, strict=False)
    expected_candidate_id = pl.concat_str(
        pl.lit("candidate_kmer_"),
        pl.col("_manifest_order").cast(pl.String).str.pad_start(6, "0"),
    )
    empty_hit_evidence = pl.all_horizontal(
        pl.col("on_target_hits").fill_null("") == "",
        pl.col("off_target_hits").fill_null("") == "",
    )
    pass_state = (
        (pl.col("status") == "PASS")
        & (pl.col("rejection_reason") == "none")
        & (pl.col("taxonomic_screening_status") == "NOT_RUN")
        & (pl.col("bait_id").fill_null("") != "")
        & (background_count == 0)
    )
    background_state = (
        (pl.col("status") == "REJECT_INTERFERENCE_BACKGROUND")
        & (pl.col("rejection_reason") == "background_occurrence")
        & (pl.col("taxonomic_screening_status") == "NOT_APPLICABLE")
        & (pl.col("bait_id").fill_null("") == "")
        & (background_count > 0)
    )
    complexity_state = (
        (pl.col("status") == "REJECT_LOW_COMPLEXITY")
        & (pl.col("rejection_reason") == "low_complexity")
        & (pl.col("taxonomic_screening_status") == "NOT_APPLICABLE")
        & (pl.col("bait_id").fill_null("") == "")
        & (background_count == 0)
    )
    kmer = pl.col("kmer")
    invalid_row = pl.coalesce(
        pl.when(
            pl.any_horizontal(
                pl.col(field).is_null()
                for field in MANIFEST_FIELDS
                if field not in {"bait_id", "on_target_hits", "off_target_hits"}
            ),
        ).then(pl.lit("Candidate manifest contains a missing field")),
        pl.when(pl.col("candidate_kmer_id") != expected_candidate_id).then(
            pl.concat_str(
                pl.lit("invalid or non-sequential Candidate K-mer ID: "),
                pl.col("candidate_kmer_id"),
            ),
        ),
        pl.when((kmer < kmer.shift(1)).fill_null(value=False)).then(
            pl.concat_str(
                pl.lit("Candidate manifest is not sorted by k-mer at "),
                pl.col("candidate_kmer_id"),
            ),
        ),
        pl.when(~kmer.str.contains(r"^[ACGT]+$") | (kmer.str.len_chars() != kmer_size)).then(
            pl.concat_str(pl.lit("invalid Candidate K-mer: "), kmer),
        ),
        pl.when(kmer != kmer.map_elements(canonical_kmer, return_dtype=pl.String)).then(
            pl.concat_str(pl.lit("noncanonical Candidate K-mer: "), kmer),
        ),
        pl.when(pl.len().over("kmer") > 1).then(
            pl.concat_str(pl.lit("duplicate Candidate K-mer: "), kmer),
        ),
        pl.when(
            (pl.col("bait_id").fill_null("") != "")
            & (pl.len().over("bait_id") > 1),
        ).then(
            pl.concat_str(pl.lit("duplicate Bait ID: "), pl.col("bait_id")),
        ),
        pl.when(
            ~pl.col("source_copy_count").str.contains(r"^[0-9]+$")
            | source_count.is_null()
            | (source_count <= 0),
        ).then(
            pl.concat_str(pl.lit("invalid source_copy_count for "), pl.col("candidate_kmer_id")),
        ),
        pl.when(
            ~pl.col("background_occurrences").str.contains(r"^[0-9]+$")
            | background_count.is_null(),
        ).then(
            pl.concat_str(
                pl.lit("invalid background_occurrences for "),
                pl.col("candidate_kmer_id"),
            ),
        ),
        pl.when(~empty_hit_evidence).then(
            pl.concat_str(
                pl.lit("pre-screening hit evidence must be empty for "),
                pl.col("candidate_kmer_id"),
            ),
        ),
        pl.when(~(pass_state | background_state | complexity_state)).then(
            pl.concat_str(
                pl.lit("inconsistent pre-screening state for "),
                pl.col("candidate_kmer_id"),
            ),
        ),
    ).alias("_invalid_row")
    parsed = manifest.with_columns(
        pl.col("bait_id").fill_null(""),
        pl.col("on_target_hits").fill_null(""),
        pl.col("off_target_hits").fill_null(""),
        source_count.alias("source_copy_count"),
        background_count.alias("background_occurrences"),
        invalid_row,
    )
    try:
        validation = parsed.select(
            pl.len().alias("row_count"),
            (pl.col("status") == "PASS").sum().alias("bait_count"),
            pl.col("_invalid_row").drop_nulls().first().alias("first_error"),
        ).collect()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Candidate manifest: {error}"
        raise TaxonomicScreeningError(message) from error
    if validation.item(0, "row_count") == 0 or validation.item(0, "bait_count") == 0:
        message = "Taxonomic screening requires at least one locally filtered Bait"
        raise TaxonomicScreeningError(message)
    if first_error := validation.item(0, "first_error"):
        raise TaxonomicScreeningError(first_error)
    return parsed.drop("_invalid_row")


def scan_blast_hits(path: Path, kmer_size: int, target_taxid: int) -> pl.LazyFrame:
    """Scan and construct valid exact-hit rows with their scope classification."""
    if _header(path, "BLAST TSV") != BLAST_FIELDS:
        message = "BLAST TSV has an unexpected schema"
        raise TaxonomicScreeningError(message)
    try:
        hits = pl.scan_csv(
            path,
            separator="\t",
            quote_char=None,
            schema_overrides=dict.fromkeys(BLAST_FIELDS, pl.String),
        ).with_row_index("_hit_order", offset=1)
        hits.collect_schema()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse BLAST TSV: {error}"
        raise TaxonomicScreeningError(message) from error

    pident = pl.col("pident").cast(pl.Float64, strict=False)
    integer_fields = ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send")
    parsed_integers = {
        field: pl.col(field).cast(pl.Int64, strict=False)
        for field in integer_fields
    }
    valid_integers = pl.all_horizontal(
        pl.col(field).str.contains(r"^[0-9]+$") & parsed_integers[field].is_not_null()
        for field in integer_fields
    )
    exact_geometry = (
        (pident == EXACT_PERCENT_IDENTITY)
        & (parsed_integers["length"] == kmer_size)
        & (parsed_integers["mismatch"] == 0)
        & (parsed_integers["gapopen"] == 0)
        & (parsed_integers["qstart"] == 1)
        & (parsed_integers["qend"] == kmer_size)
        & (parsed_integers["sstart"] > 0)
        & (parsed_integers["send"] > 0)
        & ((parsed_integers["send"] - parsed_integers["sstart"]).abs() + 1 == kmer_size)
    )
    valid_taxids = (pl.col("staxids") == "N/A") | pl.col("staxids").str.contains(
        r"^[1-9][0-9]*(;[1-9][0-9]*)*$",
    )
    taxids = pl.col("staxids").str.split(";").list.unique().list.sort()
    on_target = (
        (pl.col("staxids") != "N/A")
        & (taxids.list.len() == 1)
        & taxids.list.contains(str(target_taxid))
    )
    invalid_row = pl.coalesce(
        pl.when(pl.col("qseqid").is_null() | (pl.col("qseqid") == "")).then(
            pl.lit("BLAST TSV contains a missing query identifier"),
        ),
        pl.when(
            pl.col("pident").is_null() | pident.is_null() | ~pident.is_finite(),
        ).then(pl.lit("BLAST TSV contains an invalid pident")),
        pl.when(~valid_integers.fill_null(value=False)).then(
            pl.lit("BLAST TSV contains an invalid integer field"),
        ),
        pl.when(~exact_geometry.fill_null(value=False)).then(
            pl.lit("BLAST TSV contains a non-exact full-length match"),
        ),
        pl.when(pl.col("staxids").is_null() | ~valid_taxids).then(
            pl.lit("BLAST TSV contains malformed staxids"),
        ),
    ).alias("_invalid_row")
    parsed = hits.with_columns(
        pident.alias("pident"),
        *(expression.alias(field) for field, expression in parsed_integers.items()),
        on_target.alias("_on_target"),
        invalid_row,
    )
    try:
        first_error = parsed.select(pl.col("_invalid_row").drop_nulls().first()).collect().item()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse BLAST TSV: {error}"
        raise TaxonomicScreeningError(message) from error
    if first_error:
        raise TaxonomicScreeningError(first_error)
    return parsed.drop("_invalid_row")


def read_bait_set_status(
    path: Path,
    design_id: str,
    candidate_kmer_count: int,
    locally_filtered_bait_count: int,
) -> BaitSetStatus:
    """Read and construct the locally filtered Bait Set status."""
    if _header(path, "Bait Set status") != ["metric", "value"]:
        message = "Bait Set status has an unexpected schema"
        raise TaxonomicScreeningError(message)
    try:
        status = pl.read_csv(
            path,
            separator="\t",
            schema_overrides={"metric": pl.String, "value": pl.String},
        ).with_columns(pl.col("value").fill_null(""))
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Bait Set status: {error}"
        raise TaxonomicScreeningError(message) from error
    if status.get_column("metric").to_list() != BAIT_SET_METRICS:
        message = "Bait Set status has unexpected metrics or metric order"
        raise TaxonomicScreeningError(message)
    values = dict(status.select("metric", "value").iter_rows())
    try:
        origin = SourceSequenceOrigin(values["source_sequence_origin"])
    except ValueError as error:
        message = "Bait Set status has an invalid Source Sequence origin"
        raise TaxonomicScreeningError(message) from error
    expected = {
        "design_id": design_id,
        "candidate_kmer_count": str(candidate_kmer_count),
        "locally_filtered_bait_count": str(locally_filtered_bait_count),
        "taxonomic_screening_status": "NOT_RUN",
        "taxonomically_screened_bait_count": "",
        "deepest_bait_set": "locally_filtered",
        "deacon_index_source": "",
    }
    if any(values[key] != value for key, value in expected.items()):
        message = "Bait Set status is not coherent before taxonomic screening"
        raise TaxonomicScreeningError(message)
    return BaitSetStatus(
        design_id=design_id,
        source_sequence_origin=origin,
        candidate_kmer_count=candidate_kmer_count,
        locally_filtered_bait_count=locally_filtered_bait_count,
        deacon_index_source=values["deacon_index_source"],
    )


def construct_screening_result(
    baits: pl.LazyFrame,
    manifest: pl.LazyFrame,
    hits: pl.LazyFrame,
) -> TaxonomicScreeningResult:
    """Assign outcomes whose hit counts count accepted BLAST rows (HSPs)."""
    eligible = manifest.filter(pl.col("status") == "PASS").select("bait_id", "kmer")
    relation_errors = pl.concat(
        [
            eligible.join(baits.select("bait_id", "kmer"), on=["bait_id", "kmer"], how="anti")
            .with_columns(pl.lit("PASS Candidate K-mer is absent from the Bait FASTA").alias("error"))
            .select("error"),
            baits.select("bait_id", "kmer")
            .join(eligible, on=["bait_id", "kmer"], how="anti")
            .with_columns(pl.lit("Bait FASTA contains no matching PASS Candidate K-mer").alias("error"))
            .select("error"),
            hits.select("qseqid")
            .join(baits.select(pl.col("bait_id").alias("qseqid")), on="qseqid", how="anti")
            .with_columns(pl.lit("BLAST TSV contains an unknown Bait identifier").alias("error"))
            .select("error"),
        ],
        how="vertical",
    ).limit(1).collect()
    if relation_errors.height:
        raise TaxonomicScreeningError(relation_errors.item(0, "error"))

    hit_counts = hits.group_by(pl.col("qseqid").alias("bait_id")).agg(
        pl.col("_on_target").cast(pl.UInt64).sum().alias("on_target_hits"),
        (~pl.col("_on_target")).cast(pl.UInt64).sum().alias("off_target_hits"),
    )
    decisions = (
        baits.join(hit_counts, on="bait_id", how="left")
        .with_columns(
            pl.col("on_target_hits").fill_null(0),
            pl.col("off_target_hits").fill_null(0),
        )
        .with_columns(
            pl.when(pl.col("off_target_hits") > 0)
            .then(pl.lit("REJECT_OFF_TARGET_HIT"))
            .otherwise(pl.lit("PASS"))
            .alias("taxonomic_screening_status"),
        )
    )
    decision_evidence = decisions.select(
        "bait_id",
        pl.col("taxonomic_screening_status").alias("_screening_decision"),
        pl.col("on_target_hits").alias("_on_target_hits"),
        pl.col("off_target_hits").alias("_off_target_hits"),
    )
    transitioned_manifest = (
        manifest.join(decision_evidence, on="bait_id", how="left")
        .with_columns(
            pl.when(pl.col("_screening_decision") == "REJECT_OFF_TARGET_HIT")
            .then(pl.lit("REJECT_OFF_TARGET_HIT"))
            .otherwise(pl.col("status"))
            .alias("status"),
            pl.when(pl.col("_screening_decision") == "REJECT_OFF_TARGET_HIT")
            .then(pl.lit("off_target_exact_match"))
            .otherwise(pl.col("rejection_reason"))
            .alias("rejection_reason"),
            pl.coalesce("_screening_decision", "taxonomic_screening_status").alias(
                "taxonomic_screening_status",
            ),
            pl.when(pl.col("_screening_decision").is_not_null())
            .then(pl.col("_on_target_hits").cast(pl.String))
            .otherwise(pl.col("on_target_hits"))
            .alias("on_target_hits"),
            pl.when(pl.col("_screening_decision").is_not_null())
            .then(pl.col("_off_target_hits").cast(pl.String))
            .otherwise(pl.col("off_target_hits"))
            .alias("off_target_hits"),
        )
        .sort("_manifest_order")
        .select(MANIFEST_FIELDS)
    )
    decision_output = decisions.sort("_bait_order").select(DECISION_FIELDS)
    survivors = decisions.filter(pl.col("taxonomic_screening_status") == "PASS").sort(
        "_bait_order",
    )
    return TaxonomicScreeningResult(
        manifest=transitioned_manifest,
        decisions=decision_output,
        survivors=survivors,
    )


def status_table(metrics: Sequence[tuple[str, str]]) -> pl.DataFrame:
    """Construct one ordered status table from its domain metrics."""
    return pl.DataFrame(metrics, schema={"metric": pl.String, "value": pl.String}, orient="row")


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    baits = read_baits(args.baits, args.kmer_size)
    manifest = scan_manifest(args.manifest_in, args.kmer_size)
    hits = scan_blast_hits(args.blast_hits, args.kmer_size, args.target_taxid)
    counts = manifest.select(
        pl.len().alias("candidate_kmer_count"),
        (pl.col("status") == "PASS").sum().alias("locally_filtered_bait_count"),
    ).collect()
    candidate_kmer_count = counts.item(0, "candidate_kmer_count")
    locally_filtered_bait_count = counts.item(0, "locally_filtered_bait_count")
    bait_set = read_bait_set_status(
        args.bait_set_status_in,
        args.design_id,
        candidate_kmer_count,
        locally_filtered_bait_count,
    )
    result = construct_screening_result(baits, manifest, hits)
    manifest_output, decisions, survivors = pl.collect_all(
        [result.manifest, result.decisions, result.survivors],
    )

    screened_count = survivors.height
    screening_status = "PASS" if screened_count else "NO_BAITS"
    screening_status_output = status_table(
        [
            ("design_id", args.design_id),
            ("locally_filtered_bait_count", str(locally_filtered_bait_count)),
            ("taxonomic_screening_status", screening_status),
            ("taxonomically_screened_bait_count", str(screened_count)),
        ],
    )
    bait_set_status_output = status_table(
        [
            ("design_id", args.design_id),
            ("source_sequence_origin", bait_set.source_sequence_origin.value),
            ("candidate_kmer_count", str(candidate_kmer_count)),
            ("locally_filtered_bait_count", str(locally_filtered_bait_count)),
            ("taxonomic_screening_status", screening_status),
            ("taxonomically_screened_bait_count", str(screened_count)),
            ("deepest_bait_set", "taxonomically_screened" if screened_count else "locally_filtered"),
            ("deacon_index_source", bait_set.deacon_index_source),
        ],
    )

    manifest_output.write_csv(args.manifest_out, separator="\t", quote_style="never")
    decisions.write_csv(args.decisions_out, separator="\t", quote_style="never")
    screening_status_output.write_csv(
        args.screening_status_out,
        separator="\t",
        quote_style="never",
    )
    bait_set_status_output.write_csv(
        args.bait_set_status_out,
        separator="\t",
        quote_style="never",
    )
    if screened_count:
        SeqIO.write(
            (
                SeqRecord(Seq(kmer), id=bait_id, description="")
                for bait_id, kmer in survivors.select("bait_id", "kmer").iter_rows()
            ),
            args.baits_out,
            "fasta",
        )
    elif args.baits_out.exists():
        args.baits_out.unlink()
    if not screened_count and args.terminal_bait_set_status_out:
        bait_set_status_output.write_csv(
            args.terminal_bait_set_status_out,
            separator="\t",
            quote_style="never",
        )
    elif args.terminal_bait_set_status_out and args.terminal_bait_set_status_out.exists():
        args.terminal_bait_set_status_out.unlink()


if __name__ == "__main__":
    main()
