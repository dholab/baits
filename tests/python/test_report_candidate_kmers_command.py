from pathlib import Path

import report_candidate_kmers as reporting

PRIVATE_CANDIDATE_HEADER = (
    "candidate_kmer_id\tkmer\tsource_copy_count\tbackground_occurrences\toutcome\n"
)
PRIVATE_OCCURRENCE_HEADER = "candidate_kmer_id\tsource_sequence_id\tstart\tquery_group\n"


def write_design_tables(
    root: Path,
    slot: str,
    candidate_rows: list[str],
    occurrence_rows: list[str],
) -> None:
    manifest_dir = root / "manifests" / slot
    occurrence_dir = root / "occurrences" / slot
    manifest_dir.mkdir(parents=True)
    occurrence_dir.mkdir(parents=True)
    (manifest_dir / "candidate_kmers.tsv").write_text(
        PRIVATE_CANDIDATE_HEADER + "".join(candidate_rows),
    )
    (occurrence_dir / "candidate_kmer_occurrences.tsv").write_text(
        PRIVATE_OCCURRENCE_HEADER + "".join(occurrence_rows),
    )


def report_args(root: Path, output_dir: Path) -> list[str]:
    return [
        "--design-ids",
        str(root / "design_ids.txt"),
        "--candidate-manifests-dir",
        str(root / "manifests"),
        "--candidate-occurrences-dir",
        str(root / "occurrences"),
        "--candidates-out",
        str(output_dir / "candidate_kmers.tsv"),
        "--occurrences-out",
        str(output_dir / "candidate_kmer_occurrences.tsv"),
    ]


def test_reports_deterministic_public_candidate_tables_for_multiple_designs(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design-b\ndesign-a\n")
    write_design_tables(
        inputs,
        "design0001",
        [
            "candidate-2\tCCCCC\t1\t0\tselected\n",
            "candidate-1\tAAAAA\t1\t3\tbackground_match\n",
        ],
        [
            "candidate-2\tsource-z\t2\tgroup-z\n",
            "candidate-1\tsource-y\t0\t\n",
        ],
    )
    write_design_tables(
        inputs,
        "design0002",
        [
            "candidate-2\tTTTTT\t1\t0\tlow_complexity\n",
            "candidate-1\tGGGGG\t2\t0\tselected\n",
        ],
        [
            "candidate-1\tsource-b\t10\tgroup-b\n",
            "candidate-2\tsource-a\t1\tgroup-a\n",
            "candidate-1\tsource-a\t2\tgroup-a\n",
        ],
    )
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    reporting.main(report_args(inputs, output_dir))

    assert (output_dir / "candidate_kmers.tsv").read_text() == (
        "design_id\tkmer\tbackground_occurrences\toutcome\n"
        "design-a\tGGGGG\t0\tselected\n"
        "design-a\tTTTTT\t0\tlow_complexity\n"
        "design-b\tAAAAA\t3\tbackground_match\n"
        "design-b\tCCCCC\t0\tselected\n"
    )
    assert (output_dir / "candidate_kmer_occurrences.tsv").read_text() == (
        "design_id\tkmer\tsource_sequence_id\tstart\tquery_group\n"
        "design-a\tGGGGG\tsource-a\t2\tgroup-a\n"
        "design-a\tGGGGG\tsource-b\t10\tgroup-b\n"
        "design-a\tTTTTT\tsource-a\t1\tgroup-a\n"
        "design-b\tAAAAA\tsource-y\t0\t\n"
        "design-b\tCCCCC\tsource-z\t2\tgroup-z\n"
    )


def test_reports_headers_when_a_completed_design_has_no_candidates(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("empty-design\n")
    write_design_tables(inputs, "design0001", [], [])
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    reporting.main(report_args(inputs, output_dir))

    assert (output_dir / "candidate_kmers.tsv").read_text() == (
        "design_id\tkmer\tbackground_occurrences\toutcome\n"
    )
    assert (output_dir / "candidate_kmer_occurrences.tsv").read_text() == (
        "design_id\tkmer\tsource_sequence_id\tstart\tquery_group\n"
    )
