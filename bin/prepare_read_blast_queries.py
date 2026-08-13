#!/usr/bin/env python3
"""Prepare strand-deduplicated whole-read BLAST queries."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from collections.abc import Sequence

COUNT_SCHEMA = {
    "metagenome_id": pl.String,
    "fragment_id": pl.String,
    "mate": pl.String,
    "read_length": pl.Int64,
    "distinct_bait_count": pl.Int64,
    "fragment_distinct_bait_count": pl.Int64,
    "candidate_sequence_id": pl.String,
}
STATUS_METRICS = (
    "metagenome_id",
    "deacon_returned_read_count",
    "duplicate_fragment_count",
    "zero_bait_read_count",
    "candidate_read_count",
)
PUBLISHED_COUNT_FIELDS = (
    "metagenome_id", "fragment_id", "mate", "read_length", "distinct_bait_count",
    "fragment_distinct_bait_count", "representative_id",
)
SUMMARY_METRICS = ("design_id", *STATUS_METRICS[1:], "whole_read_blast_query_count")
THRESHOLD_SUMMARY_METRICS = (
    *SUMMARY_METRICS, "target_classified_read_count", "non_target_classified_read_count",
    "tied_read_count", "no_hit_read_count", "calibration_status",
    "recommended_deacon_absolute_threshold", "specificity_floor", "conclusion",
)
BLAST_FIELDS = (
    "qseqid", "qlen", "saccver", "staxids", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "stitle",
)
CLASSIFIED_FIELDS = (*PUBLISHED_COUNT_FIELDS, "classification", "best_target_bit_score", "best_non_target_bit_score")
THRESHOLD_COUNT_FIELDS = ("threshold", "target_read_count", "non_target_read_count", "tied_read_count", "no_hit_read_count")
THRESHOLD_CURVE_FIELDS = ("threshold", "retained_metagenome_count", "retained_fragment_count")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class BlastQueryPreparationError(ValueError):
    """Raised when candidate-read evidence cannot be prepared safely."""

    @classmethod
    def malformed_counts(cls, path: Path) -> BlastQueryPreparationError:
        return cls(f"Candidate read counts are malformed: {path}")

    @classmethod
    def malformed_status(cls, path: Path) -> BlastQueryPreparationError:
        return cls(f"Candidate read status is malformed: {path}")

    @classmethod
    def unexpected_status_metrics(cls, path: Path) -> BlastQueryPreparationError:
        return cls(f"Candidate read status has unexpected metrics: {path}")

    @classmethod
    def multiple_metagenomes(cls) -> BlastQueryPreparationError:
        return cls("Candidate read count relation identifies multiple metagenomes")

    @classmethod
    def metagenome_disagreement(cls) -> BlastQueryPreparationError:
        return cls("Candidate read counts and status identify different metagenomes")

    @classmethod
    def duplicate_identity(cls) -> BlastQueryPreparationError:
        return cls("Candidate read identity is duplicated")

    @classmethod
    def duplicate_sequence_id(cls) -> BlastQueryPreparationError:
        return cls("Candidate read sequence ID is duplicated")

    @classmethod
    def duplicate_fasta_id(cls) -> BlastQueryPreparationError:
        return cls("Candidate read FASTA record ID is duplicated")

    @classmethod
    def fasta_disagreement(cls) -> BlastQueryPreparationError:
        return cls("Candidate read count and FASTA records disagree")

    @classmethod
    def empty_fasta_disagreement(cls) -> BlastQueryPreparationError:
        return cls("An empty candidate-read count relation requires an empty FASTA")

    @classmethod
    def length_disagreement(cls) -> BlastQueryPreparationError:
        return cls("Candidate read sequence length disagrees with counts")

    @classmethod
    def status_count_disagreement(cls, metagenome_id: str) -> BlastQueryPreparationError:
        return cls(f"Candidate read status count disagrees for {metagenome_id}")

    @classmethod
    def duplicate_status_metagenome(cls) -> BlastQueryPreparationError:
        return cls("Candidate read statuses contain duplicate metagenome_id values")

    @classmethod
    def unequal_input_lists(cls) -> BlastQueryPreparationError:
        return cls("Candidate read count, FASTA, and status file lists have different lengths")

    @classmethod
    def missing_terminal_paths(cls) -> BlastQueryPreparationError:
        return cls("No-candidate evidence output paths are required")


@dataclass(frozen=True)
class CandidateReadStatus:
    metagenome_id: str
    deacon_returned_read_count: int
    duplicate_fragment_count: int
    zero_bait_read_count: int
    candidate_read_count: int


@dataclass(frozen=True)
class PreparationSummary:
    design_id: str
    deacon_returned_read_count: int
    duplicate_fragment_count: int
    zero_bait_read_count: int
    candidate_read_count: int
    whole_read_blast_query_count: int


@dataclass(frozen=True)
class TerminalCalibrationEvidence:
    blast_hits: pl.LazyFrame
    classified_reads: pl.LazyFrame
    threshold_read_counts: pl.LazyFrame
    threshold_curve: pl.LazyFrame
    summary: pl.LazyFrame


@dataclass(frozen=True)
class MetagenomeCandidateEvidence:
    metagenome_id: str
    counts: pl.LazyFrame
    sequences: pl.LazyFrame
    status: CandidateReadStatus


@dataclass(frozen=True)
class BlastQueryPreparation:
    candidate_read_counts: pl.LazyFrame
    query_records: tuple[SeqRecord, ...]
    summary: PreparationSummary
    terminal_evidence: TerminalCalibrationEvidence | None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare strand-deduplicated whole-read BLAST queries.")
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--counts", type=Path, nargs="+", required=True)
    parser.add_argument("--fastas", type=Path, nargs="+", required=True)
    parser.add_argument("--statuses", type=Path, nargs="+", required=True)
    parser.add_argument("--candidate-counts-out", type=Path, required=True)
    parser.add_argument("--query-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--terminal-blast-hits-out", type=Path)
    parser.add_argument("--terminal-classified-reads-out", type=Path)
    parser.add_argument("--terminal-read-counts-out", type=Path)
    parser.add_argument("--terminal-curve-out", type=Path)
    parser.add_argument("--terminal-summary-out", type=Path, required=True)
    return parser.parse_args(argv)


def canonical(sequence: str) -> str:
    uppercase_sequence = sequence.upper()
    return min(uppercase_sequence, uppercase_sequence.translate(COMPLEMENT)[::-1])


def has_expected_header(path: Path, fields: tuple[str, ...]) -> bool:
    try:
        with path.open(newline="") as source:
            return next(csv.reader(source, delimiter="\t")) == list(fields)
    except (OSError, StopIteration, csv.Error):
        return False


def scan_counts(path: Path) -> pl.LazyFrame:
    try:
        if not has_expected_header(path, tuple(COUNT_SCHEMA)):
            raise BlastQueryPreparationError.malformed_counts(path)
        return pl.scan_csv(path, separator="\t", schema=COUNT_SCHEMA).with_columns(
            pl.col("mate").fill_null(""),
        )
    except pl.exceptions.PolarsError as error:
        raise BlastQueryPreparationError.malformed_counts(path) from error


def read_status(path: Path) -> CandidateReadStatus:
    try:
        if not has_expected_header(path, ("metric", "value")):
            raise BlastQueryPreparationError.malformed_status(path)
        status = pl.scan_csv(path, separator="\t", schema={"metric": pl.String, "value": pl.String}, null_values=[]).collect()
    except pl.exceptions.PolarsError as error:
        raise BlastQueryPreparationError.malformed_status(path) from error
    if tuple(status.get_column("metric")) != STATUS_METRICS:
        raise BlastQueryPreparationError.unexpected_status_metrics(path)
    values = status.get_column("value")
    if values.null_count() or not values[0]:
        raise BlastQueryPreparationError.malformed_status(path)
    try:
        counts = tuple(int(value) for value in values[1:])
    except (TypeError, ValueError) as error:
        raise BlastQueryPreparationError.malformed_status(path) from error
    if any(count < 0 or count > 2**63 - 1 for count in counts):
        raise BlastQueryPreparationError.malformed_status(path)
    return CandidateReadStatus(str(values[0]), *counts)


def candidate_evidence_errors(
    counts: pl.LazyFrame,
    sequences: pl.LazyFrame,
    status: CandidateReadStatus,
) -> pl.LazyFrame:
    """Construct the ordered candidate-read evidence invariant violations."""
    status_relation = pl.LazyFrame(
        {"metagenome_id": [status.metagenome_id], "candidate_read_count": [status.candidate_read_count]},
    )
    valid_rows = (
        pl.col("metagenome_id").is_not_null()
        & pl.col("fragment_id").is_not_null()
        & pl.col("candidate_sequence_id").is_not_null()
        & (pl.col("metagenome_id").str.len_chars() > 0)
        & (pl.col("fragment_id").str.len_chars() > 0)
        & (pl.col("candidate_sequence_id").str.len_chars() > 0)
        & pl.col("mate").is_in(("", "1", "2"))
        & (pl.col("read_length") > 0)
        & (pl.col("distinct_bait_count") > 0)
        & (pl.col("fragment_distinct_bait_count") >= pl.col("distinct_bait_count"))
    )
    return pl.concat(
        [
            counts.filter(~valid_rows.fill_null(strategy="zero"))
            .with_columns(pl.lit("malformed_counts").alias("error"))
            .select("error"),
            counts.select(pl.col("metagenome_id").n_unique().alias("metagenome_count")).filter(
                pl.col("metagenome_count") > 1,
            ).with_columns(pl.lit("multiple_metagenomes").alias("error")).select("error"),
            counts.filter(pl.col("metagenome_id") != status.metagenome_id)
            .with_columns(pl.lit("metagenome_disagreement").alias("error"))
            .select("error"),
            counts.group_by("metagenome_id", "fragment_id", "mate").len().filter(
                pl.col("len") > 1,
            ).with_columns(pl.lit("duplicate_identity").alias("error")).select("error"),
            counts.group_by("candidate_sequence_id").len().filter(pl.col("len") > 1).select(
                "candidate_sequence_id",
            ).with_columns(pl.lit("duplicate_sequence_id").alias("error")).select("error"),
            sequences.group_by("candidate_sequence_id").len().filter(pl.col("len") > 1).select(
                "candidate_sequence_id",
            ).with_columns(pl.lit("duplicate_fasta_id").alias("error")).select("error"),
            counts.select("candidate_sequence_id").unique().join(
                sequences.select("candidate_sequence_id").unique(),
                on="candidate_sequence_id",
                how="anti",
            ).with_columns(pl.lit("fasta_disagreement").alias("error")).select("error"),
            sequences.select("candidate_sequence_id").unique().join(
                counts.select("candidate_sequence_id").unique(),
                on="candidate_sequence_id",
                how="anti",
            ).with_columns(pl.lit("fasta_disagreement").alias("error")).select("error"),
            counts.join(sequences, on="candidate_sequence_id", how="inner").filter(
                pl.col("read_length") != pl.col("sequence_length"),
            ).with_columns(pl.lit("length_disagreement").alias("error")).select("error"),
            counts.select(pl.len().alias("observed_count")).join(
                status_relation,
                how="cross",
            ).filter(pl.col("observed_count") != pl.col("candidate_read_count")).with_columns(
                pl.lit("status_count_disagreement").alias("error"),
            ).select("error"),
        ],
        how="vertical",
    )
def candidate_evidence_error(
    code: str,
    *,
    counts_path: Path,
    metagenome_id: str,
) -> BlastQueryPreparationError:
    errors = {
        "malformed_counts": BlastQueryPreparationError.malformed_counts(counts_path),
        "multiple_metagenomes": BlastQueryPreparationError.multiple_metagenomes(),
        "metagenome_disagreement": BlastQueryPreparationError.metagenome_disagreement(),
        "duplicate_identity": BlastQueryPreparationError.duplicate_identity(),
        "duplicate_sequence_id": BlastQueryPreparationError.duplicate_sequence_id(),
        "duplicate_fasta_id": BlastQueryPreparationError.duplicate_fasta_id(),
        "fasta_disagreement": BlastQueryPreparationError.fasta_disagreement(),
        "length_disagreement": BlastQueryPreparationError.length_disagreement(),
        "status_count_disagreement": BlastQueryPreparationError.status_count_disagreement(metagenome_id),
    }
    try:
        return errors[code]
    except KeyError as error:
        message = f"Unexpected candidate-read evidence error: {code}"
        raise AssertionError(message) from error


def construct_metagenome_candidate_evidence(
    counts_path: Path,
    fasta_path: Path,
    status_path: Path,
) -> MetagenomeCandidateEvidence:
    status = read_status(status_path)
    try:
        count_frame = scan_counts(counts_path).collect()
    except pl.exceptions.PolarsError as error:
        raise BlastQueryPreparationError.malformed_counts(counts_path) from error
    counts = count_frame.lazy()
    fasta_records = tuple(SeqIO.parse(fasta_path, "fasta"))
    sequences = pl.DataFrame(
        {
            "metagenome_id": [status.metagenome_id] * len(fasta_records),
            "candidate_sequence_id": [record.id for record in fasta_records],
            "canonical_sequence": [canonical(str(record.seq)) for record in fasta_records],
        },
        schema={
            "metagenome_id": pl.String,
            "candidate_sequence_id": pl.String,
            "canonical_sequence": pl.String,
        },
    ).lazy().with_columns(
        pl.col("canonical_sequence").str.len_chars().alias("sequence_length"),
    )
    error = candidate_evidence_errors(counts, sequences, status).limit(1).collect()
    if error.height:
        raise candidate_evidence_error(
            error.item(0, "error"),
            counts_path=counts_path,
            metagenome_id=status.metagenome_id,
        )
    return MetagenomeCandidateEvidence(status.metagenome_id, counts, sequences, status)


def status_relation(summary: PreparationSummary) -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "metric": SUMMARY_METRICS,
            "value": [
                str(summary.design_id),
                str(summary.deacon_returned_read_count),
                str(summary.duplicate_fragment_count),
                str(summary.zero_bait_read_count),
                str(summary.candidate_read_count),
                str(summary.whole_read_blast_query_count),
            ],
        },
    ).lazy()


def construct_blast_query_preparation(*, design_id: str, evidence: Sequence[MetagenomeCandidateEvidence]) -> BlastQueryPreparation:
    if len({item.metagenome_id for item in evidence}) != len(evidence):
        raise BlastQueryPreparationError.duplicate_status_metagenome()
    counts = pl.concat([item.counts for item in evidence], how="vertical") if evidence else pl.LazyFrame(schema=COUNT_SCHEMA)
    sequences = pl.concat([item.sequences for item in evidence], how="vertical") if evidence else pl.LazyFrame(schema={"metagenome_id": pl.String, "candidate_sequence_id": pl.String, "canonical_sequence": pl.String, "sequence_length": pl.Int64})
    joined = counts.join(sequences, on=["metagenome_id", "candidate_sequence_id"], how="inner")
    unique_sequences = tuple(joined.select("canonical_sequence").unique().sort("canonical_sequence").collect().get_column("canonical_sequence"))
    representatives = pl.DataFrame(
        {
            "canonical_sequence": unique_sequences,
            "representative_id": [
                f"representative_{number:06d}"
                for number in range(1, len(unique_sequences) + 1)
            ],
        },
        schema={"canonical_sequence": pl.String, "representative_id": pl.String},
    ).lazy()
    published = joined.join(representatives, on="canonical_sequence").select(PUBLISHED_COUNT_FIELDS).sort("metagenome_id", "fragment_id", "mate")
    published_count = published.select(pl.len()).collect().item()
    summary = PreparationSummary(
        design_id,
        sum(item.status.deacon_returned_read_count for item in evidence),
        sum(item.status.duplicate_fragment_count for item in evidence),
        sum(item.status.zero_bait_read_count for item in evidence),
        published_count,
        len(unique_sequences),
    )
    if unique_sequences:
        return BlastQueryPreparation(published, tuple(SeqRecord(Seq(sequence), id=f"representative_{number:06d}", description="") for number, sequence in enumerate(unique_sequences, 1)), summary, None)
    terminal_summary = pl.DataFrame(
        {
            "metric": THRESHOLD_SUMMARY_METRICS,
            "value": [
                str(summary.design_id),
                str(summary.deacon_returned_read_count),
                str(summary.duplicate_fragment_count),
                str(summary.zero_bait_read_count),
                str(summary.candidate_read_count),
                str(summary.whole_read_blast_query_count),
                "0",
                "0",
                "0",
                "0",
                "NO_CANDIDATE_READS",
                "",
                "",
                "The optimization read set contains no candidate reads.",
            ],
        },
    ).lazy()
    terminal = TerminalCalibrationEvidence(
        pl.LazyFrame(schema=dict.fromkeys(BLAST_FIELDS, pl.String)),
        pl.LazyFrame(schema=dict.fromkeys(CLASSIFIED_FIELDS, pl.String)),
        pl.LazyFrame({"threshold": [1], "target_read_count": [0], "non_target_read_count": [0], "tied_read_count": [0], "no_hit_read_count": [0]}),
        pl.LazyFrame({"threshold": [1], "retained_metagenome_count": [0], "retained_fragment_count": [0]}),
        terminal_summary,
    )
    return BlastQueryPreparation(published, (), summary, terminal)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if len({len(args.counts), len(args.fastas), len(args.statuses)}) != 1:
        raise BlastQueryPreparationError.unequal_input_lists()
    evidence = tuple(construct_metagenome_candidate_evidence(*triplet) for triplet in zip(args.counts, args.fastas, args.statuses, strict=True))
    preparation = construct_blast_query_preparation(design_id=args.design_id, evidence=evidence)
    preparation.candidate_read_counts.collect().write_csv(args.candidate_counts_out, separator="\t")
    status_relation(preparation.summary).collect().write_csv(args.summary_out, separator="\t")
    if preparation.terminal_evidence is None:
        SeqIO.write(preparation.query_records, args.query_out, "fasta")
        if args.terminal_summary_out.exists():
            args.terminal_summary_out.unlink()
        return
    terminal_paths = (args.terminal_blast_hits_out, args.terminal_classified_reads_out, args.terminal_read_counts_out, args.terminal_curve_out)
    if any(path is None for path in terminal_paths):
        raise BlastQueryPreparationError.missing_terminal_paths()
    terminal = preparation.terminal_evidence
    terminal.blast_hits.collect().write_csv(args.terminal_blast_hits_out, separator="\t")
    terminal.classified_reads.collect().write_csv(args.terminal_classified_reads_out, separator="\t")
    terminal.threshold_read_counts.collect().write_csv(args.terminal_read_counts_out, separator="\t")
    terminal.threshold_curve.collect().write_csv(args.terminal_curve_out, separator="\t")
    terminal.summary.collect().write_csv(args.terminal_summary_out, separator="\t")
    if args.query_out.exists():
        args.query_out.unlink()


if __name__ == "__main__":
    main()
