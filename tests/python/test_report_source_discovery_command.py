import gzip
from pathlib import Path

import report_source_discovery as reporting

PRIVATE_LOCUS_HEADER = (
    "candidate_locus_id\tsource_sequence_id\tquery_id\tquery_group\tassembly_sequence_id\t"
    "assembly_start\tassembly_end\tstrand\tpercent_identity\tquery_coverage\n"
)


def write_discovery_tables(
    root: Path,
    slot: str,
    locus_rows: list[str],
    match_rows: list[str],
) -> None:
    locus_dir = root / "candidate_loci" / slot
    match_dir = root / "source_matches" / slot
    locus_dir.mkdir(parents=True)
    match_dir.mkdir(parents=True)
    (locus_dir / "candidate_loci.tsv").write_text(PRIVATE_LOCUS_HEADER + "".join(locus_rows))
    (match_dir / "query_blast_hits.txt").write_text("".join(match_rows))


def report_args(root: Path, output_dir: Path) -> list[str]:
    return [
        "--design-ids",
        str(root / "design_ids.txt"),
        "--candidate-loci-dir",
        str(root / "candidate_loci"),
        "--source-matches-dir",
        str(root / "source_matches"),
        "--candidate-loci-out",
        str(output_dir / "candidate_loci.tsv"),
        "--source-matches-out",
        str(output_dir / "source_matches.tsv.gz"),
    ]


def test_reports_deterministic_source_discovery_tables_for_multiple_designs(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "design_ids.txt").write_text("design-b\ndesign-a\n")
    write_discovery_tables(
        inputs,
        "design0001",
        ["locus-b\tsource-b\tquery-b\tB\tcontig-b\t4\t44\t+\t100.0\t100.0\n"],
        ["query-b\t41\t1\t41\tcontig-b\t5\t45\t41\t100.0\n"],
    )
    write_discovery_tables(
        inputs,
        "design0002",
        ["locus-a\tsource-a\tquery-a\tA\tcontig-a\t10\t50\t-\t99.5\t98.0\n"],
        ["query-a\t42\t2\t41\tcontig-a\t51\t12\t40\t99.5\n"],
    )
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    reporting.main(report_args(inputs, output_dir))

    assert (output_dir / "candidate_loci.tsv").read_text() == (
        "design_id\tcandidate_locus_id\tsource_sequence_id\tquery_id\tquery_group\t"
        "assembly_sequence_id\tassembly_start\tassembly_end\tstrand\tpercent_identity\t"
        "query_coverage\n"
        "design-a\tlocus-a\tsource-a\tquery-a\tA\tcontig-a\t10\t50\t-\t99.5\t98.0\n"
        "design-b\tlocus-b\tsource-b\tquery-b\tB\tcontig-b\t4\t44\t+\t100.0\t100.0\n"
    )
    with gzip.open(output_dir / "source_matches.tsv.gz", "rt") as handle:
        assert handle.read() == (
            "design_id\tquery_id\tquery_length\tquery_start\tquery_end\tassembly_sequence_id\t"
            "assembly_start\tassembly_end\talignment_length\tpercent_identity\n"
            "design-a\tquery-a\t42\t2\t41\tcontig-a\t51\t12\t40\t99.5\n"
            "design-b\tquery-b\t41\t1\t41\tcontig-b\t5\t45\t41\t100.0\n"
        )


def test_reports_headers_when_no_design_uses_source_discovery(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    (inputs / "candidate_loci").mkdir(parents=True)
    (inputs / "source_matches").mkdir()
    (inputs / "design_ids.txt").write_text("")
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    reporting.main(report_args(inputs, output_dir))

    assert (output_dir / "candidate_loci.tsv").read_text().count("\n") == 1
    with gzip.open(output_dir / "source_matches.tsv.gz", "rt") as handle:
        assert handle.read().count("\n") == 1
