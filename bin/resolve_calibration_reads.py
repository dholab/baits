#!/usr/bin/env python3
"""Resolve a flat calibration-read directory into manifest evidence."""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence

MANIFEST_COLUMNS = ("id", "design_id", "metagenome_id", "read")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
SAFE_SOURCE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


class CalibrationReadError(ValueError):
    """Raised when calibration-read input cannot be resolved safely."""

    @classmethod
    def nested_directory(cls, path: Path) -> CalibrationReadError:
        return cls(f"Calibration reads contain a nested directory: {path.name}")

    @classmethod
    def no_accepted_fastq_files(cls) -> CalibrationReadError:
        return cls("Calibration reads contain no accepted FASTQ files")

    @classmethod
    def not_a_directory(cls, path: Path) -> CalibrationReadError:
        return cls(f"Calibration reads are not a directory: {path}")

@dataclass(frozen=True)
class ResolvedReadSource:
    read_source_id: str
    design_id: str
    metagenome_id: str
    read: Path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve one flat calibration-read directory.")
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    return parser.parse_args(argv)


def construct_calibration_read_sources(
    design_id: str,
    children: Sequence[Path],
) -> tuple[ResolvedReadSource, ...]:
    nested = next((path for path in children if path.is_dir()), None)
    if nested is not None:
        raise CalibrationReadError.nested_directory(nested)

    read_files = tuple(
        path
        for path in children
        if path.is_file() and any(path.name.endswith(suffix) for suffix in FASTQ_SUFFIXES)
    )
    if not read_files:
        raise CalibrationReadError.no_accepted_fastq_files()

    base_source_ids = tuple(source_id(path) for path in read_files)
    base_counts = Counter(base_source_ids)
    provisional_ids = tuple(
        base_source_id if base_counts[base_source_id] == 1 else path.name
        for base_source_id, path in zip(base_source_ids, read_files, strict=True)
    )
    duplicate_ids = frozenset(
        source for source, count in Counter(provisional_ids).items() if count > 1
    )
    source_ids = tuple(
        safe_source_id(path, provisional_id, duplicate_ids)
        for provisional_id, path in zip(provisional_ids, read_files, strict=True)
    )
    if len(set(source_ids)) != len(source_ids):
        message = "Calibration reads could not be assigned unique source IDs"
        raise CalibrationReadError(message)
    return tuple(
        ResolvedReadSource(f"{design_id}__{source}", design_id, source, path)
        for source, path in sorted(zip(source_ids, read_files, strict=True))
    )


def source_id(path: Path) -> str:
    for suffix in FASTQ_SUFFIXES:
        if path.name.endswith(suffix):
            return path.name.removesuffix(suffix)
    message = f"Not an accepted FASTQ path: {path}"
    raise AssertionError(message)


def safe_source_id(path: Path, provisional_id: str, duplicate_ids: frozenset[str]) -> str:
    if provisional_id not in duplicate_ids and SAFE_SOURCE_ID.fullmatch(provisional_id):
        return provisional_id
    digest = hashlib.sha256(path.name.encode()).hexdigest()
    return f"read_{digest}"


def construct_calibration_read_manifest(
    read_sources: Sequence[ResolvedReadSource],
) -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "id": [source.read_source_id for source in read_sources],
            "design_id": [source.design_id for source in read_sources],
            "metagenome_id": [source.metagenome_id for source in read_sources],
            "read": [source.read.name for source in read_sources],
        },
        schema=dict.fromkeys(MANIFEST_COLUMNS, pl.String),
    ).lazy().select(MANIFEST_COLUMNS)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if not args.directory.is_dir():
        raise CalibrationReadError.not_a_directory(args.directory)
    read_sources = construct_calibration_read_sources(
        args.design_id,
        tuple(args.directory.iterdir()),
    )
    construct_calibration_read_manifest(read_sources).sink_csv(
        args.manifest_out,
        separator="\t",
    )


if __name__ == "__main__":
    main()
