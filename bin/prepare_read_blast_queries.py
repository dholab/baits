#!/usr/bin/env python3
"""Combine per-metagenome bait counts and prepare unique BLAST queries."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
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
    "query_id",
)
QUERY_PATTERN = re.compile(r"sequence_[0-9a-f]{64}")


class BlastQueryPreparationError(ValueError):
    """Raised when streamed read counts and unique queries disagree."""

    @classmethod
    def missing_counts(cls) -> BlastQueryPreparationError:
        return cls("Read bait count directory contains no tables")

    @classmethod
    def malformed_counts(cls, path: Path) -> BlastQueryPreparationError:
        return cls(f"Read bait counts are malformed: {path}")

    @classmethod
    def duplicate_metagenome(cls, metagenome_id: str) -> BlastQueryPreparationError:
        return cls(f"Read bait counts repeat metagenome_id {metagenome_id}")

    @classmethod
    def malformed_query(cls, identifier: str) -> BlastQueryPreparationError:
        return cls(f"Deduplicated calibration query is malformed: {identifier}")

    @classmethod
    def query_disagreement(cls) -> BlastQueryPreparationError:
        return cls("Read bait counts and deduplicated BLAST queries disagree")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--unique-fasta", type=Path, required=True)
    parser.add_argument("--counts-out", type=Path, required=True)
    parser.add_argument("--queries-out", type=Path, required=True)
    return parser.parse_args(argv)


def discover_count_tables(root: Path) -> tuple[Path, ...]:
    tables = tuple(sorted(root.glob("*.read_bait_counts.tsv")))
    if not tables:
        raise BlastQueryPreparationError.missing_counts()
    return tables


def valid_count_row(row: dict[str, str | None], ordinal: int) -> bool:
    return (
        bool(row["metagenome_id"])
        and row["read_id"] == f"read_{ordinal:012d}"
        and (row["read_length"] or "").isascii()
        and (row["read_length"] or "").isdecimal()
        and int(row["read_length"] or "0") > 0
        and (row["bait_count"] or "").isascii()
        and (row["bait_count"] or "").isdecimal()
        and int(row["bait_count"] or "0") > 0
        and QUERY_PATTERN.fullmatch(row["query_id"] or "") is not None
    )


def combine_read_bait_counts(
    paths: Sequence[Path],
    output: Path,
) -> frozenset[str]:
    metagenome_ids: set[str] = set()
    query_ids: set[str] = set()
    with output.open("w", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=COUNT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for path in paths:
            table_metagenome, table_queries = copy_count_table(path, writer)
            query_ids.update(table_queries)
            if table_metagenome is not None:
                if table_metagenome in metagenome_ids:
                    raise BlastQueryPreparationError.duplicate_metagenome(table_metagenome)
                metagenome_ids.add(table_metagenome)
    return frozenset(query_ids)


def copy_count_table(
    path: Path,
    writer: csv.DictWriter[str],
) -> tuple[str | None, frozenset[str]]:
    table_metagenome: str | None = None
    query_ids: set[str] = set()
    try:
        with path.open(newline="") as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != COUNT_FIELDS:
                raise BlastQueryPreparationError.malformed_counts(path)
            for ordinal, row in enumerate(reader, start=1):
                if not valid_count_row(row, ordinal):
                    raise BlastQueryPreparationError.malformed_counts(path)
                metagenome_id = row["metagenome_id"] or ""
                if table_metagenome is None:
                    table_metagenome = metagenome_id
                elif table_metagenome != metagenome_id:
                    raise BlastQueryPreparationError.malformed_counts(path)
                query_ids.add(row["query_id"] or "")
                writer.writerow(row)
    except (OSError, csv.Error, KeyError, ValueError) as error:
        if isinstance(error, BlastQueryPreparationError):
            raise
        raise BlastQueryPreparationError.malformed_counts(path) from error
    return table_metagenome, frozenset(query_ids)


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


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    counted_query_ids = combine_read_bait_counts(
        discover_count_tables(args.counts_dir),
        args.counts_out,
    )
    query_ids = write_unique_queries(args.unique_fasta, args.queries_out)
    if counted_query_ids != query_ids:
        raise BlastQueryPreparationError.query_disagreement()


if __name__ == "__main__":
    main()
