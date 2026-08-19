#!/usr/bin/env python3
"""Combine streamed read counts and prepare globally deduplicated BLAST queries."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from Bio import SeqIO

if TYPE_CHECKING:
    from collections.abc import Sequence

COUNT_FIELDS = (
    "metagenome_id",
    "read_id",
    "read_length",
    "bait_count",
    "representative_id",
)
STATUS_METRICS = (
    "metagenome_id",
    "deacon_returned_read_count",
    "candidate_read_count",
)
SUMMARY_METRICS = (
    "design_id",
    "deacon_returned_read_count",
    "candidate_read_count",
    "duplicate_sequence_count",
    "read_blast_query_count",
)
THRESHOLD_SUMMARY_METRICS = (
    *SUMMARY_METRICS,
    "target_classified_read_count",
    "non_target_classified_read_count",
    "tied_read_count",
    "no_hit_read_count",
    "calibration_status",
    "recommended_threshold",
    "conclusion",
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
CLASSIFIED_FIELDS = (
    *COUNT_FIELDS,
    "classification",
    "best_target_bit_score",
    "best_non_target_bit_score",
)
THRESHOLD_COUNT_FIELDS = (
    "threshold",
    "target_read_count",
    "non_target_read_count",
    "tied_read_count",
    "no_hit_read_count",
)
REPRESENTATIVE_PATTERN = re.compile(r"sequence_[0-9a-f]{64}")
RELATION_WIDTH = 2


class BlastQueryPreparationError(ValueError):
    """Raised when streamed candidate-read products disagree."""

    @classmethod
    def malformed_counts(cls, path: Path) -> BlastQueryPreparationError:
        return cls(f"Candidate read counts are malformed: {path}")

    @classmethod
    def malformed_status(cls, path: Path) -> BlastQueryPreparationError:
        return cls(f"Candidate read status is malformed: {path}")

    @classmethod
    def duplicate_status_metagenome(cls) -> BlastQueryPreparationError:
        return cls("Candidate read statuses contain duplicate metagenome_id values")

    @classmethod
    def status_count_disagreement(cls, metagenome_id: str) -> BlastQueryPreparationError:
        return cls(f"Candidate read status count disagrees for {metagenome_id}")

    @classmethod
    def malformed_query(cls, identifier: str) -> BlastQueryPreparationError:
        return cls(f"Deduplicated candidate query is malformed: {identifier}")

    @classmethod
    def candidate_query_disagreement(cls) -> BlastQueryPreparationError:
        return cls("Candidate reads and deduplicated BLAST queries disagree")

    @classmethod
    def unequal_input_lists(cls) -> BlastQueryPreparationError:
        return cls("Candidate read count and status file lists have different lengths")

    @classmethod
    def missing_terminal_paths(cls) -> BlastQueryPreparationError:
        return cls("No-candidate output paths are required")


@dataclass(frozen=True)
class CandidateReadStatus:
    metagenome_id: str
    deacon_returned_read_count: int
    candidate_read_count: int


@dataclass(frozen=True)
class PreparationSummary:
    design_id: str
    deacon_returned_read_count: int
    candidate_read_count: int
    duplicate_sequence_count: int
    read_blast_query_count: int


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine candidate-read rows and prepare unique BLAST queries.",
    )
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--counts", type=Path, nargs="+", required=True)
    parser.add_argument("--statuses", type=Path, nargs="+", required=True)
    parser.add_argument("--unique-fasta", type=Path, required=True)
    parser.add_argument("--candidate-counts-out", type=Path, required=True)
    parser.add_argument("--query-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--terminal-blast-hits-out", type=Path)
    parser.add_argument("--terminal-search-parameters-out", type=Path)
    parser.add_argument("--terminal-classified-reads-out", type=Path)
    parser.add_argument("--terminal-read-counts-out", type=Path)
    parser.add_argument("--terminal-summary-out", type=Path, required=True)
    return parser.parse_args(argv)


def read_nonnegative_integer(value: str, path: Path) -> int:
    if not value.isascii() or not value.isdecimal():
        raise BlastQueryPreparationError.malformed_status(path)
    return int(value)


def read_status(path: Path) -> CandidateReadStatus:
    try:
        with path.open(newline="") as handle:
            rows = tuple(csv.reader(handle, delimiter="\t"))
    except (OSError, csv.Error) as error:
        raise BlastQueryPreparationError.malformed_status(path) from error
    if not rows or rows[0] != ["metric", "value"] or tuple(row[0] for row in rows[1:]) != STATUS_METRICS:
        raise BlastQueryPreparationError.malformed_status(path)
    if any(len(row) != RELATION_WIDTH for row in rows[1:]) or not rows[1][1]:
        raise BlastQueryPreparationError.malformed_status(path)
    status = CandidateReadStatus(
        rows[1][1],
        read_nonnegative_integer(rows[2][1], path),
        read_nonnegative_integer(rows[3][1], path),
    )
    if status.deacon_returned_read_count != status.candidate_read_count:
        raise BlastQueryPreparationError.status_count_disagreement(status.metagenome_id)
    return status


def valid_count_row(row: dict[str, str | None], metagenome_id: str, ordinal: int) -> bool:
    return (
        row["metagenome_id"] == metagenome_id
        and row["read_id"] == f"read_{ordinal:012d}"
        and (row["read_length"] or "").isascii()
        and (row["read_length"] or "").isdecimal()
        and int(row["read_length"] or "0") > 0
        and (row["bait_count"] or "").isascii()
        and (row["bait_count"] or "").isdecimal()
        and int(row["bait_count"] or "0") > 0
        and REPRESENTATIVE_PATTERN.fullmatch(row["representative_id"] or "") is not None
    )


def combine_candidate_counts(
    count_paths: Sequence[Path],
    status_paths: Sequence[Path],
    output: Path,
) -> tuple[int, frozenset[str]]:
    status_by_source = {
        path.name.removesuffix(".candidate_read_status.tsv"): read_status(path)
        for path in status_paths
    }
    count_by_source = {
        path.name.removesuffix(".counted_reads.tsv"): path
        for path in count_paths
    }
    if status_by_source.keys() != count_by_source.keys():
        raise BlastQueryPreparationError.unequal_input_lists()
    total = 0
    representative_ids: set[str] = set()
    with output.open("w", newline="") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=COUNT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source in sorted(count_by_source):
            path = count_by_source[source]
            status = status_by_source[source]
            try:
                with path.open(newline="") as input_handle:
                    reader = csv.DictReader(input_handle, delimiter="\t")
                    if tuple(reader.fieldnames or ()) != COUNT_FIELDS:
                        raise BlastQueryPreparationError.malformed_counts(path)
                    rows = 0
                    for rows, row in enumerate(reader, start=1):
                        if not valid_count_row(row, status.metagenome_id, rows):
                            raise BlastQueryPreparationError.malformed_counts(path)
                        representative_ids.add(row["representative_id"] or "")
                        writer.writerow(row)
            except (OSError, csv.Error, KeyError, ValueError) as error:
                if isinstance(error, BlastQueryPreparationError):
                    raise
                raise BlastQueryPreparationError.malformed_counts(path) from error
            if rows != status.candidate_read_count:
                raise BlastQueryPreparationError.status_count_disagreement(status.metagenome_id)
            total += rows
    return total, frozenset(representative_ids)


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def write_unique_queries(source: Path, destination: Path) -> frozenset[str]:
    query_ids: set[str] = set()
    with open_text(source) as input_handle, destination.open("w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence = str(record.seq).upper()
            expected_id = f"sequence_{hashlib.sha256(sequence.encode()).hexdigest()}"
            if record.id != expected_id:
                raise BlastQueryPreparationError.malformed_query(record.id)
            query_ids.add(record.id)
            output_handle.write(f">{record.id}\n{sequence}\n")
    return frozenset(query_ids)


def write_relation(path: Path, fields: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def write_summary(path: Path, summary: PreparationSummary) -> None:
    write_relation(path, ("metric", "value"), tuple(zip(SUMMARY_METRICS, summary.__dict__.values(), strict=True)))


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if len(args.counts) != len(args.statuses):
        raise BlastQueryPreparationError.unequal_input_lists()
    statuses = tuple(read_status(path) for path in args.statuses)
    if len({status.metagenome_id for status in statuses}) != len(statuses):
        raise BlastQueryPreparationError.duplicate_status_metagenome()
    candidate_count, candidate_representatives = combine_candidate_counts(
        args.counts,
        args.statuses,
        args.candidate_counts_out,
    )
    query_ids = write_unique_queries(args.unique_fasta, args.query_out)
    if candidate_representatives != query_ids:
        raise BlastQueryPreparationError.candidate_query_disagreement()
    query_count = len(query_ids)
    summary = PreparationSummary(
        args.design_id,
        sum(status.deacon_returned_read_count for status in statuses),
        candidate_count,
        candidate_count - query_count,
        query_count,
    )
    write_summary(args.summary_out, summary)
    if candidate_count:
        args.terminal_summary_out.unlink(missing_ok=True)
        return
    args.query_out.unlink(missing_ok=True)
    blast_hits = args.terminal_blast_hits_out
    search_parameters = args.terminal_search_parameters_out
    classified_reads = args.terminal_classified_reads_out
    read_counts = args.terminal_read_counts_out
    if (
        blast_hits is None
        or search_parameters is None
        or classified_reads is None
        or read_counts is None
    ):
        raise BlastQueryPreparationError.missing_terminal_paths()
    write_relation(blast_hits, BLAST_FIELDS, ())
    write_relation(search_parameters, ("parameter", "value"), ())
    write_relation(classified_reads, CLASSIFIED_FIELDS, ())
    write_relation(read_counts, THRESHOLD_COUNT_FIELDS, ((1, 0, 0, 0, 0),))
    values = (
        *summary.__dict__.values(),
        0,
        0,
        0,
        0,
        "NO_CANDIDATE_READS",
        "",
        "The calibration reads contain no candidate reads.",
    )
    write_relation(
        args.terminal_summary_out,
        ("metric", "value"),
        tuple(zip(THRESHOLD_SUMMARY_METRICS, values, strict=True)),
    )


if __name__ == "__main__":
    main()
