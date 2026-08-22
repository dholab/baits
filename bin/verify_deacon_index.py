#!/usr/bin/env python3
"""Verify Deacon round trips against the selected candidate k-mers."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO

if TYPE_CHECKING:
    from collections.abc import Sequence

    from Bio.SeqRecord import SeqRecord

MANIFEST_FIELDS = (
    "candidate_kmer_id",
    "bait_id",
    "kmer",
    "source_copy_count",
    "background_occurrences",
    "outcome",
)
OUTCOMES = frozenset({"selected", "background_match", "low_complexity", "off_target_match"})
MAX_COUNT = 2**63 - 1


class DeaconVerificationError(ValueError):
    """Raised when Deacon verification inputs violate their contract."""


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """One validated FASTA record."""

    record_id: str
    sequence: str


@dataclass(frozen=True, slots=True)
class CandidateKmer:
    """One validated candidate k-mer manifest row."""

    candidate_kmer_id: str
    bait_id: str
    kmer: str
    source_copy_count: int
    background_occurrences: int
    outcome: str


def positive_integer(value: str) -> int:
    """Parse one canonical positive integer for argparse."""
    if re.fullmatch(r"[1-9][0-9]*", value) is None:
        message = "must be a positive integer"
        raise argparse.ArgumentTypeError(message)
    return int(value)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kmer-size", type=positive_integer, required=True)
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--bait-roundtrip", type=Path, required=True)
    parser.add_argument("--background-roundtrip", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args(argv)


def canonical_kmer(kmer: str) -> str:
    """Return the canonical representation shared with a reverse complement."""
    return min(kmer, kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1])


def construct_fasta_records(
    records: Sequence[SeqRecord],
    label: str,
    expected_length: int | None,
) -> tuple[FastaRecord, ...]:
    """Construct a unique A/C/G/T FASTA relation from parsed records."""
    relation = tuple(FastaRecord(record.id or "", str(record.seq)) for record in records)
    if len({record.record_id for record in relation}) != len(relation):
        message = f"{label} record identifiers must be unique"
        raise DeaconVerificationError(message)
    if len({record.sequence for record in relation}) != len(relation):
        message = f"{label} sequences must be unique"
        raise DeaconVerificationError(message)
    invalid = next(
        (
            record
            for record in relation
            if not record.record_id
            or not record.sequence
            or record.sequence != record.sequence.upper()
            or set(record.sequence) - set("ACGT")
            or (expected_length is not None and len(record.sequence) != expected_length)
        ),
        None,
    )
    if invalid is not None:
        message = f"{label} has an invalid record: {invalid.record_id}"
        raise DeaconVerificationError(message)
    return relation


def read_fasta(path: Path, label: str, expected_length: int | None) -> tuple[FastaRecord, ...]:
    """Read one FASTA file into a validated relation."""
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
    if text.strip() and not records:
        message = f"Could not parse {label}"
        raise DeaconVerificationError(message)
    return construct_fasta_records(records, label, expected_length)


def parse_count(value: str, field: str) -> int:
    """Parse one nonnegative manifest count."""
    if re.fullmatch(r"(?:0|[1-9][0-9]*)", value) is None or int(value) > MAX_COUNT:
        message = f"Candidate k-mer manifest has an invalid {field}"
        raise DeaconVerificationError(message)
    return int(value)


def parse_manifest_row(row: dict[str, str | None], kmer_size: int) -> CandidateKmer:
    """Validate one candidate manifest row independently of external artifacts."""
    if any(row[field] is None for field in MANIFEST_FIELDS):
        message = "Candidate k-mer manifest contains a missing field"
        raise DeaconVerificationError(message)
    candidate_kmer_id = row["candidate_kmer_id"] or ""
    bait_id = row["bait_id"] or ""
    kmer = row["kmer"] or ""
    outcome = row["outcome"] or ""
    if not candidate_kmer_id or not re.fullmatch(r"[ACGT]+", kmer) or len(kmer) != kmer_size:
        message = "Candidate k-mer manifest contains an invalid candidate identity or k-mer"
        raise DeaconVerificationError(message)
    if kmer != canonical_kmer(kmer):
        message = "Candidate k-mer manifest contains a noncanonical k-mer"
        raise DeaconVerificationError(message)
    source_copy_count = parse_count(row["source_copy_count"] or "", "source_copy_count")
    if source_copy_count == 0:
        message = "Candidate k-mer manifest has an invalid source_copy_count"
        raise DeaconVerificationError(message)
    background_occurrences = parse_count(row["background_occurrences"] or "", "background_occurrences")
    if outcome not in OUTCOMES:
        message = "Candidate k-mer manifest contains an unknown outcome"
        raise DeaconVerificationError(message)
    has_bait = bool(bait_id)
    valid_outcome = (
        (outcome == "selected" and has_bait and background_occurrences == 0)
        or (outcome == "off_target_match" and not has_bait and background_occurrences == 0)
        or (outcome == "background_match" and not has_bait and background_occurrences > 0)
        or (outcome == "low_complexity" and not has_bait and background_occurrences == 0)
    )
    if not valid_outcome:
        message = "Candidate k-mer manifest contains an incoherent candidate outcome"
        raise DeaconVerificationError(message)
    return CandidateKmer(
        candidate_kmer_id=candidate_kmer_id,
        bait_id=bait_id,
        kmer=kmer,
        source_copy_count=source_copy_count,
        background_occurrences=background_occurrences,
        outcome=outcome,
    )


def read_manifest(path: Path, kmer_size: int) -> tuple[CandidateKmer, ...]:
    """Read the candidate manifest and validate its mutually exclusive outcomes."""
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
                message = "Candidate k-mer manifest has an unexpected schema"
                raise DeaconVerificationError(message)
            candidates = tuple(parse_manifest_row(row, kmer_size) for row in reader)
    except csv.Error as error:
        message = f"Could not parse candidate k-mer manifest: {error}"
        raise DeaconVerificationError(message) from error
    if not candidates:
        message = "Candidate k-mer manifest must not be empty"
        raise DeaconVerificationError(message)
    if len({candidate.candidate_kmer_id for candidate in candidates}) != len(candidates):
        message = "Candidate k-mer manifest contains a duplicate candidate identifier"
        raise DeaconVerificationError(message)
    if len({candidate.kmer for candidate in candidates}) != len(candidates):
        message = "Candidate k-mer manifest contains a duplicate k-mer"
        raise DeaconVerificationError(message)
    bait_ids = tuple(candidate.bait_id for candidate in candidates if candidate.bait_id)
    if len(set(bait_ids)) != len(bait_ids):
        message = "Candidate k-mer manifest contains a duplicate bait identifier"
        raise DeaconVerificationError(message)
    return candidates


def verify_relations(
    manifest: Sequence[CandidateKmer],
    baits: Sequence[FastaRecord],
    bait_roundtrip: Sequence[FastaRecord],
    background_roundtrip: Sequence[FastaRecord],
) -> None:
    """Verify candidate selection and both Deacon round-trip relations."""
    selected = frozenset(
        (candidate.bait_id, candidate.kmer)
        for candidate in manifest
        if candidate.outcome == "selected"
    )
    bait_relation = frozenset((record.record_id, record.sequence) for record in baits)
    if selected != bait_relation:
        message = "Selected candidate k-mers do not equal the bait FASTA"
        raise DeaconVerificationError(message)
    if {record.sequence for record in bait_roundtrip} != {record.sequence for record in baits}:
        message = "Bait sequence set differs after the Deacon round trip"
        raise DeaconVerificationError(message)
    if background_roundtrip:
        message = "Background round-trip FASTA must be empty"
        raise DeaconVerificationError(message)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    manifest = read_manifest(args.manifest, args.kmer_size)
    baits = read_fasta(args.baits, "bait FASTA", args.kmer_size)
    bait_roundtrip = read_fasta(args.bait_roundtrip, "bait round-trip FASTA", args.kmer_size)
    background_roundtrip = read_fasta(args.background_roundtrip, "background round-trip FASTA", None)
    verify_relations(manifest, baits, bait_roundtrip, background_roundtrip)


if __name__ == "__main__":
    main()
