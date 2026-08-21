#!/usr/bin/env python3
"""Gather calibration-read BLAST shards into one evidence table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

BLAST_FIELDS = (
    "qseqid", "qlen", "saccver", "staxids", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "stitle",
)
PARAMETER_FIELDS = ("parameter", "value")
PARAMETER_COLUMN_COUNT = len(PARAMETER_FIELDS)


class BlastShardGatherError(ValueError):
    """Raised when sharded BLAST evidence cannot be gathered safely."""

    @classmethod
    def unexpected_hits_header(cls, path: Path) -> BlastShardGatherError:
        return cls(f"Read BLAST shard has an unexpected header: {path}")

    @classmethod
    def unexpected_parameters_header(cls, path: Path) -> BlastShardGatherError:
        return cls(f"Read BLAST search parameters have an unexpected header: {path}")

    @classmethod
    def malformed_parameters(cls, path: Path) -> BlastShardGatherError:
        return cls(f"Read BLAST search parameters are malformed: {path}")

    @classmethod
    def different_parameters(cls) -> BlastShardGatherError:
        return cls("Read BLAST shards reported different search parameters")

    @classmethod
    def different_shard_counts(cls) -> BlastShardGatherError:
        return cls("Read BLAST hit and search-parameter shard counts differ")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast-hits", type=Path, nargs="+", required=True)
    parser.add_argument("--search-parameters", type=Path, nargs="+", required=True)
    parser.add_argument("--hits-output", type=Path, required=True)
    parser.add_argument("--parameters-output", type=Path, required=True)
    return parser.parse_args(argv)


def gather_hits(paths: Sequence[Path], output: Path) -> None:
    expected_header = "\t".join(BLAST_FIELDS)
    with output.open("w", newline="") as destination:
        destination.write(f"{expected_header}\n")
        for path in sorted(paths, key=lambda candidate: candidate.name):
            with path.open(newline="") as source:
                if source.readline().rstrip("\r\n") != expected_header:
                    raise BlastShardGatherError.unexpected_hits_header(path)
                for line in source:
                    destination.write(line if line.endswith("\n") else f"{line}\n")


def read_search_parameters(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        if tuple(next(reader, ())) != PARAMETER_FIELDS:
            raise BlastShardGatherError.unexpected_parameters_header(path)
        rows: list[tuple[str, str]] = []
        for row in reader:
            if len(row) != PARAMETER_COLUMN_COUNT:
                raise BlastShardGatherError.malformed_parameters(path)
            rows.append((row[0], row[1]))
    return rows


def gather_search_parameters(paths: Sequence[Path], output: Path) -> None:
    ordered_paths = sorted(paths, key=lambda candidate: candidate.name)
    expected_rows = read_search_parameters(ordered_paths[0])
    for path in ordered_paths[1:]:
        if read_search_parameters(path) != expected_rows:
            raise BlastShardGatherError.different_parameters()

    with output.open("w", newline="") as destination:
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(PARAMETER_FIELDS)
        writer.writerows(expected_rows)
        writer.writerow(("query_file_count", len(ordered_paths)))


def main(argv: Sequence[str] | None = None) -> None:
    arguments = parse_args(argv)
    if len(arguments.blast_hits) != len(arguments.search_parameters):
        raise BlastShardGatherError.different_shard_counts()
    gather_hits(arguments.blast_hits, arguments.hits_output)
    gather_search_parameters(arguments.search_parameters, arguments.parameters_output)


if __name__ == "__main__":
    main()
