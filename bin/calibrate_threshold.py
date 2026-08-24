#!/usr/bin/env python3
"""Derive the recommended threshold from scored calibration reads."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

READ_FIELDS = (
    "metagenome_id",
    "read_id",
    "read_length",
    "bait_count",
    "query_id",
    "best_target_bit_score",
    "best_other_bit_score",
)
TIE_TOLERANCE = Decimal("0.1")
MAX_COUNT = 2**63 - 1


class ThresholdCalibrationError(ValueError):
    """Raised when scored calibration reads cannot support threshold derivation."""


@dataclass(frozen=True)
class CalibrationRead:
    metagenome_id: str
    read_id: str
    read_length: int
    bait_count: int
    query_id: str
    best_target_bit_score: Decimal | None
    best_other_bit_score: Decimal | None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-id", required=True)
    parser.add_argument("--calibration-reads", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def parse_positive_integer(value: str | None, field: str) -> int:
    if value is None or not value.isascii() or not value.isdecimal():
        message = f"Calibration reads contain an invalid {field}"
        raise ThresholdCalibrationError(message)
    parsed = int(value)
    if parsed <= 0 or parsed > MAX_COUNT:
        message = f"Calibration reads contain an invalid {field}"
        raise ThresholdCalibrationError(message)
    return parsed


def parse_score(value: str | None, field: str) -> Decimal | None:
    if value in {None, ""}:
        return None
    try:
        score = Decimal(value)
    except InvalidOperation as error:
        message = f"Calibration reads contain an invalid {field}"
        raise ThresholdCalibrationError(message) from error
    if not score.is_finite():
        message = f"Calibration reads contain an invalid {field}"
        raise ThresholdCalibrationError(message)
    return score


def read_calibration_reads(path: Path) -> tuple[CalibrationRead, ...]:
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != READ_FIELDS:
                message = "Calibration reads have an unexpected schema"
                raise ThresholdCalibrationError(message)
            reads = tuple(
                CalibrationRead(
                    metagenome_id=row["metagenome_id"] or "",
                    read_id=row["read_id"] or "",
                    read_length=parse_positive_integer(row["read_length"], "read_length"),
                    bait_count=parse_positive_integer(row["bait_count"], "bait_count"),
                    query_id=row["query_id"] or "",
                    best_target_bit_score=parse_score(
                        row["best_target_bit_score"],
                        "best_target_bit_score",
                    ),
                    best_other_bit_score=parse_score(
                        row["best_other_bit_score"],
                        "best_other_bit_score",
                    ),
                )
                for row in reader
            )
    except (OSError, csv.Error, KeyError) as error:
        message = f"Could not parse calibration reads: {path}"
        raise ThresholdCalibrationError(message) from error
    if not reads:
        message = "Threshold derivation requires at least one calibration read"
        raise ThresholdCalibrationError(message)
    invalid_identity = next(
        (
            read
            for read in reads
            if not read.metagenome_id or not read.read_id or not read.query_id
        ),
        None,
    )
    if invalid_identity is not None:
        message = "Calibration reads contain a blank identity"
        raise ThresholdCalibrationError(message)
    identities = {(read.metagenome_id, read.read_id) for read in reads}
    if len(identities) != len(reads):
        message = "Calibration read identity is duplicated"
        raise ThresholdCalibrationError(message)
    return reads


def category_counts(read: CalibrationRead) -> tuple[int, int, int, int]:
    target = read.best_target_bit_score
    other = read.best_other_bit_score
    if target is None and other is None:
        return 0, 0, 0, 1
    if target is None:
        return 0, 1, 0, 0
    if other is None:
        return 1, 0, 0, 0
    if abs(target - other) <= TIE_TOLERANCE:
        return 0, 0, 1, 0
    return (1, 0, 0, 0) if target > other else (0, 1, 0, 0)


def recommended_threshold(reads: Sequence[CalibrationRead]) -> int | None:
    counts_by_bait_count: dict[int, list[int]] = {}
    for read in reads:
        totals = counts_by_bait_count.setdefault(read.bait_count, [0, 0, 0, 0])
        for index, count in enumerate(category_counts(read)):
            totals[index] += count
    running = [0, 0, 0, 0]
    supported: int | None = None
    maximum = max(counts_by_bait_count)
    for threshold in range(maximum, 0, -1):
        exact = counts_by_bait_count.get(threshold, [0, 0, 0, 0])
        running = [total + count for total, count in zip(running, exact, strict=True)]
        target, other, tied, no_hit = running
        if target > 0 and other == 0 and tied == 0 and no_hit == 0:
            supported = threshold
    return supported


def write_recommended_threshold(path: Path, design_id: str, threshold: int | None) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("design_id", "recommended_threshold"))
        writer.writerow((design_id, "" if threshold is None else threshold))


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    reads = read_calibration_reads(args.calibration_reads)
    write_recommended_threshold(
        args.output,
        args.design_id,
        recommended_threshold(reads),
    )


if __name__ == "__main__":
    main()
