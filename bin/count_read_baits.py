#!/usr/bin/env python3
"""Stream candidate reads into count rows and sequence-deduplication input."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from Bio import SeqIO

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from Bio.SeqRecord import SeqRecord

COUNT_FIELDS = (
    "metagenome_id",
    "read_id",
    "read_length",
    "bait_count",
    "representative_id",
)
STATUS_FIELDS = (
    "metagenome_id",
    "deacon_returned_read_count",
    "candidate_read_count",
)
COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


class CandidateReadRecountError(ValueError):
    """Raised when candidate reads cannot be recounted safely."""

    @classmethod
    def nonpositive_kmer_size(cls) -> CandidateReadRecountError:
        return cls("kmer_size must be positive")

    @classmethod
    def empty_bait_set(cls) -> CandidateReadRecountError:
        return cls("Taxonomically screened bait set must not be empty")

    @classmethod
    def duplicate_baits(cls) -> CandidateReadRecountError:
        return cls("Taxonomically screened baits must be unique")

    @classmethod
    def unequal_bait_lengths(cls) -> CandidateReadRecountError:
        return cls("Taxonomically screened baits must have one k-mer size")

    @classmethod
    def disagreeing_kmer_size(cls) -> CandidateReadRecountError:
        return cls("Taxonomically screened bait length disagrees with kmer_size")

    @classmethod
    def zero_bait_candidate(cls) -> CandidateReadRecountError:
        return cls("Deacon returned a candidate read with zero taxonomically screened baits")


@dataclass(frozen=True)
class RecountStatus:
    metagenome_id: str
    deacon_returned_read_count: int
    candidate_read_count: int


@dataclass(frozen=True)
class CandidateWriters:
    counts: TextIO
    fasta: TextIO


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count baits while streaming independent candidate reads.",
    )
    parser.add_argument("--metagenome-id", required=True)
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--kmer-size", type=int, required=True)
    parser.add_argument("--read", type=Path, required=True)
    parser.add_argument("--counts-out", type=Path, required=True)
    parser.add_argument("--fasta-out", type=Path, required=True)
    parser.add_argument("--status-out", type=Path, required=True)
    return parser.parse_args(argv)


def canonical(sequence: str) -> str:
    uppercase_sequence = sequence.upper()
    return min(uppercase_sequence, uppercase_sequence.translate(COMPLEMENT)[::-1])


def load_bait_set(path: Path, kmer_size: int) -> frozenset[str]:
    if kmer_size <= 0:
        raise CandidateReadRecountError.nonpositive_kmer_size()
    with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as handle:
        baits: set[str] = set()
        bait_lengths: set[int] = set()
        duplicate = False
        for record in SeqIO.parse(handle, "fasta"):
            bait = canonical(str(record.seq))
            duplicate = duplicate or bait in baits
            baits.add(bait)
            bait_lengths.add(len(bait))
    if not baits:
        raise CandidateReadRecountError.empty_bait_set()
    if duplicate:
        raise CandidateReadRecountError.duplicate_baits()
    if len(bait_lengths) != 1:
        raise CandidateReadRecountError.unequal_bait_lengths()
    if bait_lengths.pop() != kmer_size:
        raise CandidateReadRecountError.disagreeing_kmer_size()
    return frozenset(baits)


def iter_reads(path: Path) -> Iterator[SeqRecord]:
    with path.open("rb") as raw_handle:
        compressed = raw_handle.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rt") if compressed else path.open() as handle:
        yield from SeqIO.parse(handle, "fasta")


def matched_baits(sequence: str, baits: frozenset[str], kmer_size: int) -> frozenset[str]:
    return frozenset(
        kmer
        for offset in range(max(0, len(sequence) - kmer_size + 1))
        if (kmer := canonical(sequence[offset : offset + kmer_size])) in baits
    )


def representative_id(sequence: str) -> str:
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    return f"sequence_{digest}"


def write_candidate_reads(
    *,
    metagenome_id: str,
    records: Iterator[SeqRecord],
    baits: frozenset[str],
    kmer_size: int,
    writers: CandidateWriters,
) -> RecountStatus:
    writer = csv.writer(writers.counts, delimiter="\t", lineterminator="\n")
    writer.writerow(COUNT_FIELDS)
    candidate_read_count = 0
    for candidate_read_count, record in enumerate(records, start=1):
        sequence = canonical(str(record.seq))
        bait_count = len(matched_baits(sequence, baits, kmer_size))
        if bait_count == 0:
            raise CandidateReadRecountError.zero_bait_candidate()
        query_id = representative_id(sequence)
        writer.writerow(
            (
                metagenome_id,
                f"read_{candidate_read_count:012d}",
                len(sequence),
                bait_count,
                query_id,
            ),
        )
        writers.fasta.write(f">{query_id}\n{sequence}\n")
    return RecountStatus(
        metagenome_id,
        candidate_read_count,
        candidate_read_count,
    )


def write_status(path: Path, status: RecountStatus) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(zip(STATUS_FIELDS, status.__dict__.values(), strict=True))


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    baits = load_bait_set(args.baits, args.kmer_size)
    with args.counts_out.open("w", newline="") as counts_handle, args.fasta_out.open("w") as fasta_handle:
        status = write_candidate_reads(
            metagenome_id=args.metagenome_id,
            records=iter_reads(args.read),
            baits=baits,
            kmer_size=args.kmer_size,
            writers=CandidateWriters(counts_handle, fasta_handle),
        )
    write_status(args.status_out, status)


if __name__ == "__main__":
    main()
