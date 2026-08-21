#!/usr/bin/env python3
"""Resolve calibration-read filenames into manifest evidence."""

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
    def no_accepted_fastq_files(cls) -> CalibrationReadError:
        return cls("Calibration reads contain no accepted FASTQ files")


@dataclass(frozen=True)
class ResolvedReadSource:
    read_source_id: str
    design_id: str
    metagenome_id: str
    read_name: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve calibration-read filenames.")
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--names", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    return parser.parse_args(argv)


def construct_calibration_read_sources(
    design_id: str,
    names: Sequence[str],
) -> tuple[ResolvedReadSource, ...]:
    read_names = tuple(
        name for name in names if any(name.endswith(suffix) for suffix in FASTQ_SUFFIXES)
    )
    if not read_names:
        raise CalibrationReadError.no_accepted_fastq_files()
    base_source_ids = tuple(source_id(name) for name in read_names)
    base_counts = Counter(base_source_ids)
    provisional_ids = tuple(
        base_source_id if base_counts[base_source_id] == 1 else name
        for base_source_id, name in zip(base_source_ids, read_names, strict=True)
    )
    duplicate_ids = frozenset(
        source for source, count in Counter(provisional_ids).items() if count > 1
    )
    source_ids = tuple(
        safe_source_id(name, provisional_id, duplicate_ids)
        for provisional_id, name in zip(provisional_ids, read_names, strict=True)
    )
    if len(set(source_ids)) != len(source_ids):
        message = "Calibration reads could not be assigned unique source IDs"
        raise CalibrationReadError(message)
    return tuple(
        ResolvedReadSource(f"{design_id}__{source}", design_id, source, name)
        for source, name in sorted(zip(source_ids, read_names, strict=True))
    )


def source_id(name: str) -> str:
    for suffix in FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name.removesuffix(suffix)
    message = f"Not an accepted FASTQ filename: {name}"
    raise AssertionError(message)


def safe_source_id(name: str, provisional_id: str, duplicate_ids: frozenset[str]) -> str:
    if provisional_id not in duplicate_ids and SAFE_SOURCE_ID.fullmatch(provisional_id):
        return provisional_id
    digest = hashlib.sha256(name.encode()).hexdigest()
    return f"read_{digest}"


def construct_calibration_read_manifest(
    read_sources: Sequence[ResolvedReadSource],
) -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "id": [source.read_source_id for source in read_sources],
            "design_id": [source.design_id for source in read_sources],
            "metagenome_id": [source.metagenome_id for source in read_sources],
            "read": [source.read_name for source in read_sources],
        },
        schema=dict.fromkeys(MANIFEST_COLUMNS, pl.String),
    ).lazy().select(MANIFEST_COLUMNS)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    read_sources = construct_calibration_read_sources(
        args.design_id,
        args.names.read_text().splitlines(),
    )
    construct_calibration_read_manifest(read_sources).sink_csv(
        args.manifest_out,
        separator="\t",
    )


if __name__ == "__main__":
    main()
