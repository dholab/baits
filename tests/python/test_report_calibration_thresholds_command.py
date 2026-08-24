from pathlib import Path

import pytest
from report_calibration_thresholds import CalibrationThresholdReportError, main


def arguments(root: Path, output: Path) -> list[str]:
    return [
        "--design-ids",
        str(root / "design_ids.txt"),
        "--thresholds-dir",
        str(root / "thresholds"),
        "--output",
        str(output),
    ]


def write_threshold(root: Path, slot: str, design_id: str, threshold: str) -> None:
    directory = root / "thresholds" / slot
    directory.mkdir(parents=True)
    (directory / "recommended_threshold.tsv").write_text(
        f"design_id\trecommended_threshold\n{design_id}\t{threshold}\n",
    )


def test_reports_requested_designs_with_blank_missing_or_unsupported_thresholds(
    tmp_path: Path,
) -> None:
    (tmp_path / "design_ids.txt").write_text("design-c\ndesign-a\ndesign-b\n")
    write_threshold(tmp_path, "design0001", "design-a", "2")
    write_threshold(tmp_path, "design0002", "design-c", "")
    output = tmp_path / "calibration_thresholds.tsv"

    main(arguments(tmp_path, output))

    assert output.read_text() == (
        "design_id\trecommended_threshold\n"
        "design-a\t2\n"
        "design-b\t\n"
        "design-c\t\n"
    )


def test_rejects_a_threshold_for_an_unrequested_design(tmp_path: Path) -> None:
    (tmp_path / "design_ids.txt").write_text("design-a\n")
    write_threshold(tmp_path, "design0001", "other", "2")

    with pytest.raises(CalibrationThresholdReportError, match="unrequested"):
        main(arguments(tmp_path, tmp_path / "calibration_thresholds.tsv"))
