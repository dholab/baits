#!/usr/bin/env python3
"""Publish one recommended-threshold row per design requesting calibration."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

THRESHOLD_FIELDS = ("design_id", "recommended_threshold")


class CalibrationThresholdReportError(ValueError):
    """Raised when private threshold rows cannot form the public relation."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-ids", type=Path, required=True)
    parser.add_argument("--thresholds-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def read_design_ids(path: Path) -> tuple[str, ...]:
    design_ids = tuple(path.read_text().splitlines())
    if any(not value or value != value.strip() or "\t" in value for value in design_ids):
        message = "Calibration design IDs are malformed"
        raise CalibrationThresholdReportError(message)
    if len(set(design_ids)) != len(design_ids):
        message = "Calibration design IDs are duplicated"
        raise CalibrationThresholdReportError(message)
    return tuple(sorted(design_ids))


def discover_thresholds(root: Path) -> tuple[Path, ...]:
    if not root.exists():
        return ()
    paths = tuple(
        directory / "recommended_threshold.tsv"
        for directory in sorted(path for path in root.iterdir() if path.is_dir())
    )
    missing = next((path for path in paths if not path.is_file()), None)
    if missing is not None:
        message = f"Calibration threshold input is missing: {missing}"
        raise CalibrationThresholdReportError(message)
    return paths


def read_threshold(path: Path) -> tuple[str, str]:
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != THRESHOLD_FIELDS:
                message = f"Calibration threshold schema is malformed: {path}"
                raise CalibrationThresholdReportError(message)
            rows = tuple(reader)
    except (OSError, csv.Error) as error:
        message = f"Could not read calibration threshold: {path}"
        raise CalibrationThresholdReportError(message) from error
    if len(rows) != 1:
        message = f"Calibration threshold must contain one row: {path}"
        raise CalibrationThresholdReportError(message)
    design_id = rows[0]["design_id"] or ""
    threshold = rows[0]["recommended_threshold"] or ""
    if not design_id or (threshold and (not threshold.isascii() or not threshold.isdecimal() or int(threshold) <= 0)):
        message = f"Calibration threshold row is malformed: {path}"
        raise CalibrationThresholdReportError(message)
    return design_id, threshold


def construct_thresholds(
    design_ids: Sequence[str],
    paths: Sequence[Path],
) -> tuple[tuple[str, str], ...]:
    requested = set(design_ids)
    thresholds: dict[str, str] = {}
    for path in paths:
        design_id, threshold = read_threshold(path)
        if design_id not in requested:
            message = f"Threshold belongs to an unrequested design: {design_id}"
            raise CalibrationThresholdReportError(message)
        if design_id in thresholds:
            message = f"Threshold design is duplicated: {design_id}"
            raise CalibrationThresholdReportError(message)
        thresholds[design_id] = threshold
    return tuple((design_id, thresholds.get(design_id, "")) for design_id in design_ids)


def write_thresholds(path: Path, rows: Sequence[tuple[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(THRESHOLD_FIELDS)
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    design_ids = read_design_ids(args.design_ids)
    write_thresholds(
        args.output,
        construct_thresholds(design_ids, discover_thresholds(args.thresholds_dir)),
    )


if __name__ == "__main__":
    main()
