#!/usr/bin/env python3
"""Resolve a flat optimization read set directory into manifest evidence."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence

MANIFEST_COLUMNS = ("id", "design_id", "metagenome_id", "read_1", "read_2")
PAIRED_READ_COUNT = 2
FASTQ_PATTERN = re.compile(
    r"^(?P<metagenome_id>[A-Za-z0-9][A-Za-z0-9._-]*?)(?P<mate>_R[12])?\.(?:fastq|fq)(?:\.gz)?$",
)


class OptimizationReadSetError(ValueError):
    """Raised when optimization read set input cannot be resolved safely."""

    @classmethod
    def nested_directory(cls, path: Path) -> OptimizationReadSetError:
        return cls(f"Optimization read set contains a nested directory: {path.name}")

    @classmethod
    def no_accepted_fastq_files(cls) -> OptimizationReadSetError:
        return cls("Optimization read set contains no accepted FASTQ files")

    @classmethod
    def duplicate_files(cls, metagenome_id: str, role: str) -> OptimizationReadSetError:
        return cls(f"Optimization read set has duplicate files for {metagenome_id}: {role}")

    @classmethod
    def mixed_layout(cls, metagenome_id: str) -> OptimizationReadSetError:
        return cls(f"Optimization read set mixes single-end and paired files for {metagenome_id}")

    @classmethod
    def missing_mate(cls, metagenome_id: str, mate: str) -> OptimizationReadSetError:
        return cls(f"Optimization read set is missing {mate} for {metagenome_id}")

    @classmethod
    def invalid_files(cls, metagenome_id: str) -> OptimizationReadSetError:
        return cls(f"Optimization read set has invalid files for {metagenome_id}")

    @classmethod
    def not_a_directory(cls, path: Path) -> OptimizationReadSetError:
        return cls(f"Optimization read set is not a directory: {path}")

    @classmethod
    def invalid_read_count(cls) -> OptimizationReadSetError:
        return cls("Optimization read set must contain one or two reads")


@dataclass(frozen=True)
class ResolvedReadSet:
    read_set_id: str
    design_id: str
    metagenome_id: str
    reads: tuple[Path, ...]

    def __post_init__(self) -> None:
        if len(self.reads) not in (1, PAIRED_READ_COUNT):
            raise OptimizationReadSetError.invalid_read_count()


@dataclass(frozen=True)
class OptimizationReadFile:
    metagenome_id: str
    role: str
    path: Path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve one flat optimization read set directory.")
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    return parser.parse_args(argv)


def construct_optimization_read_sets(
    design_id: str,
    children: Sequence[Path],
) -> tuple[ResolvedReadSet, ...]:
    nested = next((path for path in children if path.is_dir()), None)
    if nested is not None:
        raise OptimizationReadSetError.nested_directory(nested)

    read_files = tuple(
        OptimizationReadFile(
            metagenome_id=match.group("metagenome_id"),
            role="single" if match.group("mate") is None else match.group("mate")[1:],
            path=path,
        )
        for path in children
        if path.is_file() and (match := FASTQ_PATTERN.fullmatch(path.name)) is not None
    )
    if not read_files:
        raise OptimizationReadSetError.no_accepted_fastq_files()

    ordered_read_files = sorted(
        read_files,
        key=lambda read_file: (read_file.metagenome_id, read_file.role),
    )
    return tuple(
        construct_read_set(design_id, metagenome_id, tuple(group))
        for metagenome_id, group in groupby(
            ordered_read_files,
            key=lambda read_file: read_file.metagenome_id,
        )
    )


def construct_read_set(
    design_id: str,
    metagenome_id: str,
    read_files: tuple[OptimizationReadFile, ...],
) -> ResolvedReadSet:
    roles = tuple(read_file.role for read_file in read_files)
    if len(set(roles)) != len(roles):
        duplicate_role = next(role for role in roles if roles.count(role) > 1)
        raise OptimizationReadSetError.duplicate_files(metagenome_id, duplicate_role)
    if "single" in roles and len(roles) != 1:
        raise OptimizationReadSetError.mixed_layout(metagenome_id)
    if roles == ("R1",):
        raise OptimizationReadSetError.missing_mate(metagenome_id, "R2")
    if roles == ("R2",):
        raise OptimizationReadSetError.missing_mate(metagenome_id, "R1")
    if roles == ("single",):
        reads = (read_files[0].path,)
    elif roles == ("R1", "R2"):
        reads = (read_files[0].path, read_files[1].path)
    else:
        raise OptimizationReadSetError.invalid_files(metagenome_id)
    return ResolvedReadSet(f"{design_id}__{metagenome_id}", design_id, metagenome_id, reads)


def construct_optimization_read_manifest(
    read_sets: Sequence[ResolvedReadSet],
) -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "id": [read_set.read_set_id for read_set in read_sets],
            "design_id": [read_set.design_id for read_set in read_sets],
            "metagenome_id": [read_set.metagenome_id for read_set in read_sets],
            "read_1": [read_set.reads[0].name for read_set in read_sets],
            "read_2": [
                read_set.reads[1].name
                if len(read_set.reads) == PAIRED_READ_COUNT
                else ""
                for read_set in read_sets
            ],
        },
        schema=dict.fromkeys(MANIFEST_COLUMNS, pl.String),
    ).lazy().select(MANIFEST_COLUMNS)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if not args.directory.is_dir():
        raise OptimizationReadSetError.not_a_directory(args.directory)
    read_sets = construct_optimization_read_sets(args.design_id, tuple(args.directory.iterdir()))
    construct_optimization_read_manifest(read_sets).sink_csv(args.manifest_out, separator="\t")


if __name__ == "__main__":
    main()
