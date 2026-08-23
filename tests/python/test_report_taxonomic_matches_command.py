import gzip
from pathlib import Path

import pytest
import report_taxonomic_matches as reporting

PRIVATE_CANDIDATE_HEADER = (
    "candidate_kmer_id\tkmer\tsource_copy_count\tbackground_occurrences\toutcome\n"
)
PRIVATE_MATCH_HEADER = "\t".join(reporting.PRIVATE_MATCH_SCHEMA) + "\n"


def write_design_tables(
    root: Path,
    slot: str,
    candidate_rows: list[str],
    match_rows: list[str],
) -> None:
    manifest_dir = root / "manifests" / slot
    match_dir = root / "matches" / slot
    manifest_dir.mkdir(parents=True)
    match_dir.mkdir(parents=True)
    (manifest_dir / "candidate_kmers.tsv").write_text(
        PRIVATE_CANDIDATE_HEADER + "".join(candidate_rows),
    )
    (match_dir / "taxonomic_blast_hits.tsv").write_text(
        PRIVATE_MATCH_HEADER + "".join(match_rows),
    )


def report_args(root: Path, output: Path) -> list[str]:
    return [
        "--design-ids",
        str(root / "design_ids.txt"),
        "--candidate-manifests-dir",
        str(root / "manifests"),
        "--taxonomic-matches-dir",
        str(root / "matches"),
        "--matches-out",
        str(output),
    ]


def read_report(path: Path) -> str:
    with gzip.open(path, "rt") as handle:
        return handle.read()


def test_reports_deterministic_taxonomic_matches_for_multiple_designs(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design-b\ndesign-a\n")
    write_design_tables(
        inputs,
        "design0001",
        ["candidate-2\tCCCCC\t1\t0\toff_target_match\n"],
        ["candidate-2\tB.1\t9\tB species\t100\t5\t0\t0\t1\t5\t20\t24\tb title\n"],
    )
    write_design_tables(
        inputs,
        "design0002",
        ["candidate-1\tAAAAA\t1\t0\toff_target_match\n"],
        [
            "candidate-1\tZ.1\t8\tZ species\t100\t5\t0\t0\t1\t5\t30\t34\tz title\n",
            "candidate-1\tA.1\t7\tA species\t100\t5\t0\t0\t1\t5\t10\t14\ta title\n",
        ],
    )
    output = tmp_path / "taxonomic_matches.tsv.gz"

    reporting.main(report_args(inputs, output))

    assert read_report(output) == (
        "design_id\tkmer\tsubject_accession\ttaxids\tscientific_names\tpercent_identity\talignment_length\tmismatches\tgap_opens\tquery_start\tquery_end\tsubject_start\tsubject_end\tsubject_title\n"
        "design-a\tAAAAA\tA.1\t7\tA species\t100.0\t5\t0\t0\t1\t5\t10\t14\ta title\n"
        "design-a\tAAAAA\tZ.1\t8\tZ species\t100.0\t5\t0\t0\t1\t5\t30\t34\tz title\n"
        "design-b\tCCCCC\tB.1\t9\tB species\t100.0\t5\t0\t0\t1\t5\t20\t24\tb title\n"
    )


def test_reports_only_the_header_when_no_design_reaches_taxonomic_screening(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    (inputs / "manifests").mkdir(parents=True)
    (inputs / "matches").mkdir()
    (inputs / "design_ids.txt").write_text("")
    output = tmp_path / "taxonomic_matches.tsv.gz"

    reporting.main(report_args(inputs, output))

    assert read_report(output) == "\t".join(reporting.PUBLIC_MATCH_COLUMNS) + "\n"


def test_rejects_a_match_with_an_unknown_candidate_identifier(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design\n")
    write_design_tables(
        inputs,
        "design0001",
        ["candidate-1\tAAAAA\t1\t0\tselected\n"],
        ["unknown\tA.1\t7\tA species\t100\t5\t0\t0\t1\t5\t10\t14\ta title\n"],
    )

    with pytest.raises(reporting.TaxonomicMatchReportError, match="unknown candidate"):
        reporting.main(report_args(inputs, tmp_path / "taxonomic_matches.tsv.gz"))


def test_rejects_an_off_target_outcome_without_a_raw_match(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design\n")
    write_design_tables(
        inputs,
        "design0001",
        ["candidate-1\tAAAAA\t1\t0\toff_target_match\n"],
        [],
    )

    with pytest.raises(reporting.TaxonomicMatchReportError, match="do not support off-target"):
        reporting.main(report_args(inputs, tmp_path / "taxonomic_matches.tsv.gz"))


def test_rejects_an_unexpected_private_table_header(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design\n")
    write_design_tables(
        inputs,
        "design0001",
        ["candidate-1\tAAAAA\t1\t0\tselected\n"],
        [],
    )
    match_table = inputs / "matches" / "design0001" / "taxonomic_blast_hits.tsv"
    match_table.write_text(match_table.read_text().replace("qseqid", "wrong", 1))

    with pytest.raises(reporting.TaxonomicMatchReportError, match="unexpected schema"):
        reporting.main(report_args(inputs, tmp_path / "taxonomic_matches.tsv.gz"))
