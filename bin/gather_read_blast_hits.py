#!/usr/bin/env python3
"""Gather calibration-read BLAST shards into one match table."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

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


class BlastShardGatherError(ValueError):
    """Raised when sharded calibration BLAST matches cannot be gathered."""

    @classmethod
    def missing_shards(cls) -> BlastShardGatherError:
        return cls("Calibration BLAST shard directory contains no match tables")

    @classmethod
    def unexpected_header(cls, path: Path) -> BlastShardGatherError:
        return cls(f"Calibration BLAST shard has an unexpected header: {path}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast-hits-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def discover_shards(root: Path) -> tuple[Path, ...]:
    shards = tuple(sorted(root.glob("*.read_blast_hits.tsv")))
    if not shards:
        raise BlastShardGatherError.missing_shards()
    return shards


def gather_hits(paths: Sequence[Path], output: Path) -> None:
    expected_header = "\t".join(BLAST_FIELDS)
    with output.open("w", newline="") as destination:
        destination.write(f"{expected_header}\n")
        for path in paths:
            with path.open(newline="") as source:
                if source.readline().rstrip("\r\n") != expected_header:
                    raise BlastShardGatherError.unexpected_header(path)
                for line in source:
                    destination.write(line if line.endswith("\n") else f"{line}\n")


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    gather_hits(discover_shards(args.blast_hits_dir), args.output)


if __name__ == "__main__":
    main()
