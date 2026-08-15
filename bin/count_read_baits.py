#!/usr/bin/env python3
"""Recount taxonomically screened baits on candidate reads."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

COUNT_SCHEMA = {
    "metagenome_id": pl.String,
    "fragment_id": pl.String,
    "mate": pl.String,
    "read_length": pl.Int64,
    "bait_count": pl.Int64,
    "candidate_sequence_id": pl.String,
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class CandidateReadRecountError(ValueError):
    """Raised when candidate-read recount input is malformed."""

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
    def unequal_paired_record_counts(cls) -> CandidateReadRecountError:
        return cls("Paired candidate-read files contain different record counts")

    @classmethod
    def invalid_read_file_count(cls) -> CandidateReadRecountError:
        return cls("Candidate reads must contain one single-end file or one paired file set")

    @classmethod
    def invalid_fragment_layout(cls) -> CandidateReadRecountError:
        return cls("Read fragment must contain one or two reads")

    @classmethod
    def mismatched_fragment_identities(cls) -> CandidateReadRecountError:
        return cls("Paired candidate-read records have different fragment identifiers")

    @classmethod
    def duplicate_fragment_identity(cls, identifier: str) -> CandidateReadRecountError:
        return cls(f"Candidate read fragment identifier is duplicated: {identifier}")


@dataclass(frozen=True)
class RecountStatus:
    metagenome_id: str
    deacon_returned_read_count: int
    duplicate_fragment_count: int
    zero_bait_read_count: int
    candidate_read_count: int


@dataclass(frozen=True)
class CandidateReadRecount:
    counts: pl.LazyFrame
    fasta_records: tuple[SeqRecord, ...]
    status: RecountStatus


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count baits on individual candidate reads.",
    )
    parser.add_argument("--metagenome-id", required=True)
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--kmer-size", type=int, required=True)
    parser.add_argument("--reads", type=Path, nargs="+", required=True)
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
        baits = tuple(canonical(str(record.seq)) for record in SeqIO.parse(handle, "fasta"))
    if not baits:
        raise CandidateReadRecountError.empty_bait_set()
    if len(baits) != len(set(baits)):
        raise CandidateReadRecountError.duplicate_baits()
    bait_lengths = {len(bait) for bait in baits}
    if len(bait_lengths) != 1:
        raise CandidateReadRecountError.unequal_bait_lengths()
    if bait_lengths.pop() != kmer_size:
        raise CandidateReadRecountError.disagreeing_kmer_size()
    return frozenset(baits)


def iter_read_fragments(reads: Sequence[Path]) -> Iterator[tuple[SeqRecord, ...]]:
    match reads:
        case (read,):
            with gzip.open(read, "rt") if read.suffix == ".gz" else read.open() as handle:
                yield from ((record,) for record in SeqIO.parse(handle, "fastq"))
        case (read_1, read_2):
            with (
                gzip.open(read_1, "rt") if read_1.suffix == ".gz" else read_1.open() as handle_1,
                gzip.open(read_2, "rt") if read_2.suffix == ".gz" else read_2.open() as handle_2,
            ):
                for first_record, second_record in zip_longest(
                    SeqIO.parse(handle_1, "fastq"),
                    SeqIO.parse(handle_2, "fastq"),
                ):
                    if first_record is None or second_record is None:
                        raise CandidateReadRecountError.unequal_paired_record_counts()
                    yield first_record, second_record
        case _:
            raise CandidateReadRecountError.invalid_read_file_count()


def fragment_id(identifier: str) -> str:
    return identifier.removesuffix("/1").removesuffix("/2")


def matched_baits(sequence: str, baits: frozenset[str], kmer_size: int) -> frozenset[str]:
    return frozenset(
        kmer
        for offset in range(max(0, len(sequence) - kmer_size + 1))
        if (kmer := canonical(sequence[offset : offset + kmer_size])) in baits
    )


def construct_candidate_read_recount(
    *,
    metagenome_id: str,
    baits: frozenset[str],
    kmer_size: int,
    fragments: Iterable[tuple[SeqRecord, ...]],
) -> CandidateReadRecount:
    metagenome_ids: list[str] = []
    fragment_ids: list[str] = []
    mates: list[str] = []
    read_lengths: list[int] = []
    bait_counts: list[int] = []
    candidate_sequence_ids: list[str] = []
    fasta_records: list[SeqRecord] = []
    seen_fragment_ids: set[str] = set()
    seen_signatures: set[tuple[str, ...]] = set()
    deacon_returned_read_count = 0
    duplicate_fragment_count = 0
    zero_bait_read_count = 0

    for records in fragments:
        if len(records) not in (1, 2):
            raise CandidateReadRecountError.invalid_fragment_layout()
        current_fragment_ids = tuple(fragment_id(record.id or "") for record in records)
        if len(set(current_fragment_ids)) != 1:
            raise CandidateReadRecountError.mismatched_fragment_identities()
        current_fragment_id = current_fragment_ids[0]
        if current_fragment_id in seen_fragment_ids:
            raise CandidateReadRecountError.duplicate_fragment_identity(current_fragment_id)
        seen_fragment_ids.add(current_fragment_id)
        deacon_returned_read_count += len(records)

        canonical_sequences = tuple(canonical(str(record.seq)) for record in records)
        if canonical_sequences in seen_signatures:
            duplicate_fragment_count += 1
            continue
        seen_signatures.add(canonical_sequences)

        bait_matches = tuple(
            matched_baits(sequence, baits, kmer_size) for sequence in canonical_sequences
        )
        for mate_index, (sequence, read_baits) in enumerate(
            zip(canonical_sequences, bait_matches, strict=True),
            start=1,
        ):
            bait_count = len(read_baits)
            if bait_count == 0:
                zero_bait_read_count += 1
                continue
            candidate_sequence_id = f"candidate_read_{len(fasta_records) + 1:06d}"
            metagenome_ids.append(metagenome_id)
            fragment_ids.append(current_fragment_id)
            mates.append("" if len(records) == 1 else str(mate_index))
            read_lengths.append(len(sequence))
            bait_counts.append(bait_count)
            candidate_sequence_ids.append(candidate_sequence_id)
            fasta_records.append(SeqRecord(Seq(sequence), id=candidate_sequence_id, description=""))

    counts = pl.DataFrame(
        (
            metagenome_ids,
            fragment_ids,
            mates,
            read_lengths,
            bait_counts,
            candidate_sequence_ids,
        ),
        schema=COUNT_SCHEMA,
        orient="col",
    ).lazy()
    status = RecountStatus(
        metagenome_id=metagenome_id,
        deacon_returned_read_count=deacon_returned_read_count,
        duplicate_fragment_count=duplicate_fragment_count,
        zero_bait_read_count=zero_bait_read_count,
        candidate_read_count=len(fasta_records),
    )
    return CandidateReadRecount(counts, tuple(fasta_records), status)


def construct_status_relation(status: RecountStatus) -> pl.LazyFrame:
    return pl.DataFrame(
        (
            (
                "metagenome_id",
                "deacon_returned_read_count",
                "duplicate_fragment_count",
                "zero_bait_read_count",
                "candidate_read_count",
            ),
            (
                status.metagenome_id,
                str(status.deacon_returned_read_count),
                str(status.duplicate_fragment_count),
                str(status.zero_bait_read_count),
                str(status.candidate_read_count),
            ),
        ),
        schema={"metric": pl.String, "value": pl.String},
        orient="col",
    ).lazy()


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    baits = load_bait_set(args.baits, args.kmer_size)
    recount = construct_candidate_read_recount(
        metagenome_id=args.metagenome_id,
        baits=baits,
        kmer_size=args.kmer_size,
        fragments=iter_read_fragments(args.reads),
    )
    recount.counts.sink_csv(args.counts_out, separator="\t", maintain_order=True)
    SeqIO.write(recount.fasta_records, args.fasta_out, "fasta")
    construct_status_relation(recount.status).sink_csv(args.status_out, separator="\t")


if __name__ == "__main__":
    main()
