import csv
import gzip
import hashlib
import io
from pathlib import Path

import pytest
import report_calibration_results as reporting

PRIVATE_READ_HEADER = "\t".join(reporting.PRIVATE_READ_SCHEMA) + "\n"
PRIVATE_MATCH_HEADER = "\t".join(reporting.PRIVATE_MATCH_SCHEMA) + "\n"
QUERY_A = f"sequence_{hashlib.sha256(b'AAAAA').hexdigest()}"
QUERY_B = f"sequence_{hashlib.sha256(b'CCCCC').hexdigest()}"
QUERY_C = f"sequence_{hashlib.sha256(b'GGGGG').hexdigest()}"


def write_design_tables(
    root: Path,
    slot: str,
    read_rows: list[str],
    match_rows: list[str],
) -> None:
    read_dir = root / "reads" / slot
    match_dir = root / "matches" / slot
    read_dir.mkdir(parents=True)
    match_dir.mkdir(parents=True)
    (read_dir / "calibration_reads.tsv").write_text(
        PRIVATE_READ_HEADER + "".join(read_rows),
    )
    (match_dir / "read_blast_hits.tsv").write_text(
        PRIVATE_MATCH_HEADER + "".join(match_rows),
    )


def report_args(root: Path, output_dir: Path) -> list[str]:
    return [
        "--design-ids",
        str(root / "design_ids.txt"),
        "--calibration-reads-dir",
        str(root / "reads"),
        "--calibration-matches-dir",
        str(root / "matches"),
        "--reads-out",
        str(output_dir / "calibration_reads.tsv"),
        "--matches-out",
        str(output_dir / "calibration_matches.tsv.gz"),
    ]


def read_gzip(path: Path) -> str:
    with gzip.open(path, "rt") as handle:
        return handle.read()


def match_row(query_id: str, accession: str, bitscore: str) -> str:
    return (
        f"{query_id}\t100\t{accession}\t88456\t100\t100\t0\t0\t1\t100\t1\t100\t"
        f"1e-20\t{bitscore}\t100\ttitle\n"
    )


def test_reports_deterministic_calibration_rows_for_multiple_designs(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design-b\ndesign-a\n")
    write_design_tables(
        inputs,
        "design0001",
        [f"sample-b\tread_000000000001\t100\t2\t{QUERY_B}\t90\t80\n"],
        [match_row(QUERY_B, "B.1", "90")],
    )
    write_design_tables(
        inputs,
        "design0002",
        [
            f"sample-a\tread_000000000002\t100\t1\t{QUERY_A}\t100\t\n",
            f"sample-a\tread_000000000001\t100\t1\t{QUERY_A}\t100\t\n",
            f"sample-a\tread_000000000003\t100\t1\t{QUERY_C}\t\t\n",
        ],
        [match_row(QUERY_A, "A.1", "100")],
    )
    output = tmp_path / "results"
    output.mkdir()

    reporting.main(report_args(inputs, output))

    assert (output / "calibration_reads.tsv").read_text() == (
        "design_id\tmetagenome_id\tread_id\tread_length\tbait_count\tquery_id\tbest_target_bit_score\tbest_other_bit_score\n"
        f"design-a\tsample-a\tread_000000000001\t100\t1\t{QUERY_A}\t100\t\n"
        f"design-a\tsample-a\tread_000000000002\t100\t1\t{QUERY_A}\t100\t\n"
        f"design-a\tsample-a\tread_000000000003\t100\t1\t{QUERY_C}\t\t\n"
        f"design-b\tsample-b\tread_000000000001\t100\t2\t{QUERY_B}\t90\t80\n"
    )
    matches = read_gzip(output / "calibration_matches.tsv.gz")
    assert matches.splitlines()[0].split("\t") == reporting.PUBLIC_MATCH_COLUMNS
    rows = list(csv.DictReader(io.StringIO(matches), delimiter="\t"))
    assert [
        (row["design_id"], row["query_id"], row["subject_accession"], row["bit_score"])
        for row in rows
    ] == [
        ("design-a", QUERY_A, "A.1", "100.0"),
        ("design-b", QUERY_B, "B.1", "90.0"),
    ]


def test_reports_headers_when_there_are_no_completed_calibration_designs(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    (inputs / "reads").mkdir(parents=True)
    (inputs / "matches").mkdir()
    (inputs / "design_ids.txt").write_text("")
    output = tmp_path / "results"
    output.mkdir()

    reporting.main(report_args(inputs, output))

    assert (output / "calibration_reads.tsv").read_text() == (
        "\t".join(reporting.PUBLIC_READ_COLUMNS) + "\n"
    )
    assert read_gzip(output / "calibration_matches.tsv.gz") == (
        "\t".join(reporting.PUBLIC_MATCH_COLUMNS) + "\n"
    )


def test_rejects_a_match_with_an_unknown_query(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design\n")
    write_design_tables(
        inputs,
        "design0001",
        [f"sample\tread_000000000001\t100\t1\t{QUERY_A}\t100\t\n"],
        [match_row("unknown", "A.1", "100")],
    )
    output = tmp_path / "results"
    output.mkdir()

    with pytest.raises(reporting.CalibrationReportError, match="unknown query"):
        reporting.main(report_args(inputs, output))
