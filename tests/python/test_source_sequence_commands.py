import subprocess
from pathlib import Path

import polars as pl
import pytest
from Bio import SeqIO
from extract_source_sequences import (
    all_query_groups_have_sources,
    construct_candidate_loci,
    construct_provisional_source_sequences,
)
from prepare_queries import construct_prepared_queries

PROJECT = Path(__file__).resolve().parents[2]


def run_command(*args: str) -> None:
    subprocess.run(args, cwd=PROJECT, check=True)


def run_failing_command(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=PROJECT, text=True, capture_output=True, check=False)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    return [(record.id, str(record.seq)) for record in SeqIO.parse(path, "fasta")]


def read_tsv(path: Path) -> list[dict[str, str]]:
    return pl.read_csv(path, separator="\t", infer_schema=False).to_dicts()


def write_common_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    queries = tmp_path / "queries.fasta"
    queries.write_text(">q1\nAACCGGTTAACC\n>q2\nAACCGGTT\n>q3\nCCCCAAAA\n")
    rules = tmp_path / "rules.tsv"
    rules.write_text(
        "query_id\tquery_group\tquery_start\tquery_end\tmin_identity\tmin_query_coverage\n"
        "q1\talpha\t2\t11\t98\t90\n"
        "q2\tbeta\t1\t0\t98\t90\n",
    )
    prepared = tmp_path / "prepared.fasta"
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nACCGGTTAACNNNAACCGGTT\n>contig2\nAACCGGTT\n>contig3\nGTTAACCGGT\n")
    return queries, rules, prepared, assembly


def test_construct_candidate_loci_extracts_and_orders_assembly_intervals() -> None:
    query_rules = pl.LazyFrame(
        {
            "query_id": ["q1"],
            "query_group": ["alpha"],
            "min_identity": [98.0],
            "min_query_coverage": [90.0],
        },
    )
    prepared_queries = pl.LazyFrame(
        {"query_id": ["q1"], "prepared_query_length": [4]},
    )
    target_assembly = pl.LazyFrame(
        {
            "assembly_sequence_id": ["forward", "reverse"],
            "assembly_sequence": ["ACCG", "NRYT"],
            "assembly_sequence_length": [4, 4],
        },
    )
    blast_hits = pl.LazyFrame(
        {
            "query_id": ["q1", "q1", "q1", "q1"],
            "query_length": [4, 4, 4, 4],
            "query_start": [1, 1, 1, 1],
            "query_end": [4, 4, 3, 4],
            "assembly_sequence_id": ["reverse", "forward", "forward", "forward"],
            "assembly_start_raw": [4, 1, 1, 1],
            "assembly_end_raw": [1, 4, 3, 4],
            "alignment_length": [4, 4, 3, 4],
            "percent_identity": [99.0, 98.0, 100.0, 97.9],
        },
    )

    candidate_loci = construct_candidate_loci(
        blast_hits,
        query_rules,
        prepared_queries,
        target_assembly,
    ).collect()

    assert candidate_loci.select(
        "candidate_locus_id",
        "assembly_sequence_id",
        "strand",
        "sequence",
    ).to_dicts() == [
        {
            "candidate_locus_id": "candidate_locus_000001",
            "assembly_sequence_id": "forward",
            "strand": "+",
            "sequence": "ACCG",
        },
        {
            "candidate_locus_id": "candidate_locus_000002",
            "assembly_sequence_id": "reverse",
            "strand": "-",
            "sequence": "ARYN",
        },
    ]


def test_construct_candidate_loci_rejects_missing_blast_fields() -> None:
    query_rules = pl.LazyFrame(
        {
            "query_id": ["q1"],
            "query_group": ["alpha"],
            "min_identity": [98.0],
            "min_query_coverage": [90.0],
        },
    )
    prepared_queries = pl.LazyFrame(
        {"query_id": ["q1"], "prepared_query_length": [4]},
    )
    target_assembly = pl.LazyFrame(
        {
            "assembly_sequence_id": ["contig"],
            "assembly_sequence": ["ACCG"],
            "assembly_sequence_length": [4],
        },
    )
    blast_hits = pl.LazyFrame(
        {
            "query_id": ["q1"],
            "query_length": [4],
            "query_start": [1],
            "query_end": [4],
            "assembly_sequence_id": ["contig"],
            "assembly_start_raw": [1],
            "assembly_end_raw": [4],
            "alignment_length": [None],
            "percent_identity": [100.0],
        },
    )

    with pytest.raises(ValueError, match="missing BLAST field"):
        construct_candidate_loci(
            blast_hits,
            query_rules,
            prepared_queries,
            target_assembly,
        )


def test_construct_provisional_source_sequences_deduplicates_within_query_groups() -> None:
    candidate_loci = pl.LazyFrame(
        {
            "candidate_locus_id": [
                "candidate_locus_000001",
                "candidate_locus_000002",
                "candidate_locus_000003",
            ],
            "query_group": ["alpha", "alpha", "beta"],
            "sequence": ["ACCG", "ACCG", "ACCG"],
        },
    )

    provisional_source_sequences = construct_provisional_source_sequences(
        candidate_loci,
    ).collect()

    assert provisional_source_sequences.to_dicts() == [
        {
            "source_sequence_id": "source_sequence_000001",
            "query_group": "alpha",
            "sequence": "ACCG",
        },
        {
            "source_sequence_id": "source_sequence_000002",
            "query_group": "beta",
            "sequence": "ACCG",
        },
    ]


def test_all_query_groups_have_sources_detects_a_missing_group() -> None:
    query_rules = pl.LazyFrame(
        {"query_id": ["q1", "q2", "q3"], "query_group": ["alpha", "alpha", "beta"]},
    )
    provisional_source_sequences = pl.LazyFrame(
        {
            "source_sequence_id": ["source_sequence_000001"],
            "query_group": ["alpha"],
            "sequence": ["ACCG"],
        },
    )

    assert not all_query_groups_have_sources(query_rules, provisional_source_sequences)
    assert all_query_groups_have_sources(
        query_rules.filter(pl.col("query_group") == "alpha"),
        provisional_source_sequences,
    )


def test_construct_prepared_queries_applies_one_based_closed_intervals() -> None:
    query_rules = pl.LazyFrame(
        {
            "query_id": ["q1", "q2"],
            "query_start": [2, 1],
            "query_end": [5, 0],
        },
    )
    representative_queries = pl.LazyFrame(
        {
            "query_id": ["q1", "q2"],
            "representative_query_sequence": ["AACCGG", "TTAA"],
            "representative_query_length": [6, 4],
        },
    )

    prepared_queries = construct_prepared_queries(
        query_rules,
        representative_queries,
    ).collect()

    assert prepared_queries.to_dicts() == [
        {"query_id": "q1", "sequence": "ACCG"},
        {"query_id": "q2", "sequence": "TTAA"},
    ]


def test_construct_prepared_queries_rejects_a_start_after_the_query() -> None:
    query_rules = pl.LazyFrame(
        {"query_id": ["q1"], "query_start": [5], "query_end": [0]},
    )
    representative_queries = pl.LazyFrame(
        {
            "query_id": ["q1"],
            "representative_query_sequence": ["ACCG"],
            "representative_query_length": [4],
        },
    )

    with pytest.raises(ValueError, match="Query interval out of range for q1"):
        construct_prepared_queries(query_rules, representative_queries)


def test_missing_candidate_locus_omits_provisional_source_sequence_outputs(tmp_path: Path) -> None:
    queries, rules, prepared, assembly = write_common_inputs(tmp_path)
    run_command(
        "bin/prepare_queries.py",
        "--representative-queries",
        str(queries),
        "--query-rules",
        str(rules),
        "--output-fasta",
        str(prepared),
    )
    hits = tmp_path / "hits.tsv"
    hits.write_text("q1\t10\t1\t10\tcontig1\t1\t10\t10\t100\n")
    source = tmp_path / "source.fasta"
    query_groups = tmp_path / "source_sequence_query_groups.tsv"
    loci = tmp_path / "loci.tsv"
    run_command(
        "bin/extract_source_sequences.py",
        "--prepared-queries",
        str(prepared),
        "--target-assembly",
        str(assembly),
        "--query-rules",
        str(rules),
        "--blast-hits",
        str(hits),
        "--source-sequences-out",
        str(source),
        "--source-sequence-query-groups-out",
        str(query_groups),
        "--candidate-loci-out",
        str(loci),
    )

    assert source.read_text() == ""
    assert query_groups.read_text() == "source_sequence_id\tquery_group\n"
    assert loci.read_text().splitlines()[1].startswith("candidate_locus_000001\tsource_sequence_000001\tq1\talpha")


def test_extract_source_sequences_rejects_malformed_blast_hits(tmp_path: Path) -> None:
    queries, rules, prepared, assembly = write_common_inputs(tmp_path)
    run_command(
        "bin/prepare_queries.py",
        "--representative-queries", str(queries),
        "--query-rules", str(rules),
        "--output-fasta", str(prepared),
    )
    invalid_rows = [
        ("unknown query ID", "unknown\t10\t1\t10\tcontig1\t1\t10\t10\t100"),
        ("unknown assembly ID", "q1\t10\t1\t10\tunknown\t1\t10\t10\t100"),
        ("query coordinates", "q1\t10\t0\t10\tcontig1\t1\t10\t10\t100"),
        ("assembly coordinates", "q1\t10\t1\t10\tcontig1\t1\t99\t10\t100"),
        ("percent identity", "q1\t10\t1\t10\tcontig1\t1\t10\t10\t101"),
        ("alignment length", "q1\t10\t1\t10\tcontig1\t1\t10\t0\t100"),
        ("inconsistent alignment", "q1\t10\t1\t10\tcontig1\t1\t10\t9\t100"),
    ]
    for expected_error, row in invalid_rows:
        hits = tmp_path / "hits.tsv"
        hits.write_text(f"{row}\n")
        result = run_failing_command(
            "bin/extract_source_sequences.py",
            "--prepared-queries", str(prepared),
            "--target-assembly", str(assembly),
            "--query-rules", str(rules),
            "--blast-hits", str(hits),
            "--source-sequences-out", str(tmp_path / "source.fasta"),
            "--source-sequence-query-groups-out", str(tmp_path / "source_sequence_query_groups.tsv"),
            "--candidate-loci-out", str(tmp_path / "loci.tsv"),
        )
        assert result.returncode != 0
        assert expected_error in result.stderr


def test_extract_source_sequences_orders_tied_loci_independent_of_blast_row_order(tmp_path: Path) -> None:
    queries, rules, prepared, assembly = write_common_inputs(tmp_path)
    run_command(
        "bin/prepare_queries.py",
        "--representative-queries", str(queries),
        "--query-rules", str(rules),
        "--output-fasta", str(prepared),
    )
    rows = [
        "q1\t10\t1\t10\tcontig1\t1\t10\t10\t98",
        "q1\t10\t1\t10\tcontig1\t1\t10\t10\t99",
    ]
    outputs: list[str] = []
    for index, ordered_rows in enumerate((rows, list(reversed(rows)))):
        hits = tmp_path / f"hits_{index}.tsv"
        hits.write_text("\n".join(ordered_rows) + "\n")
        loci = tmp_path / f"loci_{index}.tsv"
        run_command(
            "bin/extract_source_sequences.py",
            "--prepared-queries", str(prepared),
            "--target-assembly", str(assembly),
            "--query-rules", str(rules),
            "--blast-hits", str(hits),
            "--source-sequences-out", str(tmp_path / f"source_{index}.fasta"),
            "--source-sequence-query-groups-out", str(tmp_path / f"source_sequence_query_groups_{index}.tsv"),
            "--candidate-loci-out", str(loci),
        )
        outputs.append(loci.read_text())
    assert outputs[0] == outputs[1]


def test_curated_normalization_writes_blank_query_groups(tmp_path: Path) -> None:
    curated = tmp_path / "curated.fasta"
    curated.write_text(">sourceA description\nacgun\n>sourceB\nTTAA\n")
    source = tmp_path / "source.fasta"
    query_groups = tmp_path / "source_sequence_query_groups.tsv"

    run_command(
        "bin/normalize_curated_source_sequences.py",
        "--source-sequences",
        str(curated),
        "--source-sequences-out",
        str(source),
        "--source-sequence-query-groups-out",
        str(query_groups),
    )

    assert read_fasta(source) == [("sourceA", "ACGTN"), ("sourceB", "TTAA")]
    assert query_groups.read_text().splitlines() == ["source_sequence_id\tquery_group", "sourceA\t", "sourceB\t"]


def test_duplicate_representative_query_ids_are_malformed_input(tmp_path: Path) -> None:
    queries = tmp_path / "queries.fasta"
    queries.write_text(">q1\nAAAA\n>q1\nTTTT\n")
    rules = tmp_path / "rules.tsv"
    rules.write_text("query_id\tquery_group\tquery_start\tquery_end\tmin_identity\tmin_query_coverage\nq1\talpha\t1\t0\t98\t90\n")

    result = run_failing_command(
        "bin/prepare_queries.py",
        "--representative-queries",
        str(queries),
        "--query-rules",
        str(rules),
        "--output-fasta",
        str(tmp_path / "prepared.fasta"),
    )

    assert result.returncode != 0
    assert "Duplicate representative query FASTA record ID: q1" in result.stderr
