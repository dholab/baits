import hashlib
import json
import subprocess
from base64 import b64encode
from pathlib import Path

import polars as pl
import pytest
from Bio import SeqIO
from extract_source_sequences import (
    construct_candidate_loci,
    construct_provisional_source_sequences,
    construct_query_group_status,
)
from prepare_queries import construct_prepared_queries
from write_provenance import (
    ProvenanceInputError,
    ProvenanceRequest,
    construct_provenance,
    main,
)

PROJECT = Path(__file__).resolve().parents[2]


def run_command(*args: str) -> None:
    subprocess.run(args, cwd=PROJECT, check=True)


def run_failing_command(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=PROJECT, text=True, capture_output=True, check=False)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    return [(record.id, str(record.seq)) for record in SeqIO.parse(path, "fasta")]


def base64_text(value: str) -> str:
    return b64encode(value.encode()).decode()


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


def test_construct_query_group_status_reports_missing_provisional_source_sequences() -> None:
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

    status = construct_query_group_status(
        query_rules,
        provisional_source_sequences,
    ).collect()

    assert status.to_dicts() == [
        {
            "query_group": "alpha",
            "status": "PASS",
            "provisional_source_sequence_count": 1,
        },
        {
            "query_group": "beta",
            "status": "NO_CANDIDATE_LOCUS",
            "provisional_source_sequence_count": 0,
        },
    ]


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
    status = tmp_path / "status.tsv"

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
        "--discovery-status-out",
        str(status),
    )

    assert not source.exists()
    assert not query_groups.exists()
    assert "beta\tNO_CANDIDATE_LOCUS\t0" in status.read_text().splitlines()
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
            "--discovery-status-out", str(tmp_path / "status.tsv"),
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
            "--discovery-status-out", str(tmp_path / f"status_{index}.tsv"),
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
    assert "Duplicate Representative Query FASTA record ID: q1" in result.stderr


def construct_source_sequence_provenance(
    tmp_path: Path,
) -> tuple[tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame], tuple[Path, Path, Path, Path]]:
    background = tmp_path / "background.fasta"
    background.write_text(">b1\nTTTT\n")
    queries = tmp_path / "representative_queries.fasta"
    queries.write_text(">q1\nACGT\n")
    rules = tmp_path / "query_rules.tsv"
    rules.write_text("query_id\tquery_group\tquery_start\tquery_end\tmin_identity\tmin_query_coverage\nq1\talpha\t1\t0\t98\t90\n")
    assembly = tmp_path / "target_assembly.fasta"
    assembly.write_text(">contig\nACGT\n")
    provenance = construct_provenance(ProvenanceRequest(
        '[{"input_role":"design","input_id":"design1","attribute":"source_sequence_origin","value":"query_guided_discovery"},{"input_role":"target_taxon","input_id":"design1","attribute":"target_taxid","value":"88456"}]',
        ["interference_background", "representative_queries", "query_rules", "target_assembly"],
        ["interference_background", "representative_queries", "query_rules", "target_assembly"],
        ["file", "file", "file", "file"],
        [background, queries, rules, assembly],
        '[{"parameter":"target_taxid","value":"88456"}]',
        '[{"component":"biopython","version":"1.87"},{"component":"blast","version":"2.16.0"},{"component":"biopython","version":"1.87"}]',
    ))
    return provenance, (background, queries, rules, assembly)


def test_provenance_constructs_ordered_source_sequence_rows(tmp_path: Path) -> None:
    (inputs, parameters, software), (_background, queries, _rules, _assembly) = construct_source_sequence_provenance(tmp_path)

    assert inputs.to_dicts() == sorted(
        inputs.to_dicts(),
        key=lambda row: (row["input_role"], row["input_id"], row["attribute"]),
    )
    assert parameters.to_dicts() == [{"parameter": "target_taxid", "value": "88456"}]
    assert {row["component"] for row in software.to_dicts()} == {"biopython", "blast", "python"}
    assert {
        (row["input_role"], row["input_id"], row["attribute"], row["value"])
        for row in inputs.to_dicts()
    } >= {
        ("design", "design1", "source_sequence_origin", "query_guided_discovery"),
        ("representative_queries", "representative_queries", "sha256", hashlib.sha256(queries.read_bytes()).hexdigest()),
        ("target_assembly", "target_assembly", "filename", "target_assembly.fasta"),
    }
    assert all(not Path(row["value"]).is_absolute() for row in inputs.to_dicts())


def test_provenance_constructs_directory_basename(tmp_path: Path) -> None:
    directory = tmp_path / "reference"
    directory.mkdir()
    inputs, _parameters, _software = construct_provenance(ProvenanceRequest(
        "[]", ["reference"], ["reference"], ["directory"], [directory], "[]", "[]",
    ))
    assert inputs.to_dicts() == [{"input_role": "reference", "input_id": "reference", "attribute": "directory", "value": "reference"}]


@pytest.mark.parametrize(
    ("input_facts", "parameters", "versions"),
    [
        ("{}", "[]", "[]"),
        ('[{"input_role":"role","input_id":"id","attribute":"a","value":"v","extra":"x"}]', "[]", "[]"),
        ('[{"input_role":"role","input_id":"id","attribute":"a","value":null}]', "[]", "[]"),
        ('[{"input_role":" ","input_id":"id","attribute":"a","value":"v"}]', "[]", "[]"),
        ("[]", '[{"parameter":"","value":"v"}]', "[]"),
        ("[]", "[]", '[{"component":"tool","version":" "}]'),
    ],
)
def test_provenance_rejects_malformed_or_blank_json_rows(
    tmp_path: Path, input_facts: str, parameters: str, versions: str,
) -> None:
    path = tmp_path / "input.txt"
    path.write_text("input\n")
    with pytest.raises(ProvenanceInputError):
        construct_provenance(ProvenanceRequest(input_facts, ["role"], ["id"], ["file"], [path], parameters, versions))


def test_provenance_rejects_bad_file_descriptors(tmp_path: Path) -> None:
    path = tmp_path / "input.txt"
    path.write_text("input\n")
    directory = tmp_path / "directory"
    directory.mkdir()
    with pytest.raises(ProvenanceInputError, match="equal lengths"):
        construct_provenance(ProvenanceRequest("[]", ["role"], ["id", "other"], ["file"], [path], "[]", "[]"))
    with pytest.raises(ProvenanceInputError, match="directory"):
        construct_provenance(ProvenanceRequest("[]", ["role"], ["id"], ["directory"], [path], "[]", "[]"))
    with pytest.raises(ProvenanceInputError, match="regular file"):
        construct_provenance(ProvenanceRequest("[]", ["role"], ["id"], ["file"], [directory], "[]", "[]"))
    with pytest.raises(ProvenanceInputError, match="nonblank"):
        construct_provenance(ProvenanceRequest("[]", [" "], ["id"], ["file"], [path], "[]", "[]"))
    with pytest.raises(ProvenanceInputError, match="basename"):
        construct_provenance(ProvenanceRequest("[]", ["role"], ["id"], ["directory"], [Path("/")], "[]", "[]"))


@pytest.mark.parametrize(
    ("input_facts", "parameters", "versions", "message"),
    [
        ('[{"input_role":"role","input_id":"id","attribute":"a","value":"v"},{"input_role":"role","input_id":"id","attribute":"a","value":"other"}]', "[]", "[]", "input facts"),
        ("[]", '[{"parameter":"p","value":"one"},{"parameter":"p","value":"two"}]', "[]", "parameters"),
        ("[]", "[]", '[{"component":"tool","version":"one"},{"component":"tool","version":"two"}]', "software versions"),
    ],
)
def test_provenance_rejects_conflicting_keys(
    tmp_path: Path, input_facts: str, parameters: str, versions: str, message: str,
) -> None:
    path = tmp_path / "input.txt"
    path.write_text("input\n")
    with pytest.raises(ProvenanceInputError, match=message):
        construct_provenance(ProvenanceRequest(input_facts, ["role"], ["id"], ["file"], [path], parameters, versions))


def test_provenance_deduplicates_identical_software_facts(tmp_path: Path) -> None:
    path = tmp_path / "input.txt"
    path.write_text("input\n")
    _inputs, _parameters, software = construct_provenance(ProvenanceRequest(
        "[]", ["role"], ["id"], ["file"], [path], "[]",
        '[{"component":"tool","version":"1"},{"component":"tool","version":"1"}]',
    ))
    assert [row for row in software.to_dicts() if row["component"] == "tool"] == [
        {"component": "tool", "version": "1"},
    ]


def test_provenance_does_not_write_partial_outputs_after_failure(tmp_path: Path) -> None:
    path = tmp_path / "input.txt"
    path.write_text("input\n")
    outputs = [tmp_path / name for name in ("inputs.tsv", "parameters.tsv", "software_versions.tsv")]
    with pytest.raises(ProvenanceInputError):
        main([
            "--input-facts-base64", base64_text('[{"input_role":"role","input_id":"id","attribute":"a","value":"one"},{"input_role":"role","input_id":"id","attribute":"a","value":"two"}]'),
            "--input-file-roles-base64", base64_text('["role"]'), "--input-file-ids-base64", base64_text('["id"]'),
            "--input-file-kinds-base64", base64_text('["file"]'), "--input-files", str(path),
            "--parameters-base64", base64_text("[]"), "--software-versions-base64", base64_text("[]"),
            "--inputs-out", str(outputs[0]), "--parameters-out", str(outputs[1]), "--software-versions-out", str(outputs[2]),
        ])
    assert not any(output.exists() for output in outputs)


@pytest.mark.parametrize("value", ["two\tcolumns", "two\nrows", "carriage\rreturn"])
def test_provenance_rejects_tsv_control_characters(tmp_path: Path, value: str) -> None:
    path = tmp_path / "input.fasta"
    path.write_text(">input\nACGT\n")
    facts = json.dumps([
        {"input_role": "design", "input_id": "design", "attribute": "origin", "value": value},
    ])

    with pytest.raises(ProvenanceInputError, match="single-line TSV"):
        construct_provenance(ProvenanceRequest(facts, ["role"], ["id"], ["file"], [path], "[]", "[]"))
