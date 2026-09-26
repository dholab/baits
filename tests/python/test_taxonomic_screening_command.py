from pathlib import Path

import apply_taxonomic_screening as screening
import polars as pl
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

KMER_REJECTED = "AAAAA"
KMERS = ("AACGT", "ACGTA", "AGTCC", "CCGTA", "CGTAC")


def candidate_frame() -> pl.LazyFrame:
    return pl.LazyFrame(
        {
            "_candidate_order": range(1, len(KMERS) + 1),
            "candidate_kmer_id": [
                f"candidate_kmer_{number + 1:06d}"
                for number in range(1, len(KMERS) + 1)
            ],
            "kmer": KMERS,
        },
        schema={
            "_candidate_order": pl.UInt32,
            "candidate_kmer_id": pl.String,
            "kmer": pl.String,
        },
    )


def manifest_frame() -> pl.LazyFrame:
    rows = [
        {
            "candidate_kmer_id": "candidate_kmer_000001",
            "kmer": KMER_REJECTED,
            "source_copy_count": 1,
            "background_occurrences": 0,
            "outcome": "low_complexity",
            "_manifest_order": 1,
        },
        *[
            {
                "candidate_kmer_id": f"candidate_kmer_{number + 1:06d}",
                "kmer": kmer,
                "source_copy_count": 1,
                "background_occurrences": 0,
                "outcome": "selected",
                "_manifest_order": number + 1,
            }
            for number, kmer in enumerate(KMERS, start=1)
        ],
    ]
    return pl.LazyFrame(rows)


def hit_frame(rows: list[tuple[str, bool]]) -> pl.LazyFrame:
    return pl.LazyFrame(
        {
            "qseqid": [candidate_kmer_id for candidate_kmer_id, _ in rows],
            "_on_target": [on_target for _, on_target in rows],
        },
        schema={"qseqid": pl.String, "_on_target": pl.Boolean},
    )


def write_inputs(tmp_path: Path, *, all_rejected: bool = False) -> tuple[Path, Path, Path]:
    candidates = tmp_path / "selected_candidates.fasta"
    candidates.write_text(
        "".join(
            f">candidate_kmer_{number + 1:06d}\n{kmer}\n"
            for number, kmer in enumerate(KMERS, start=1)
        ),
    )
    manifest = tmp_path / "candidate_kmers.tsv"
    manifest_frame().collect().sort("_manifest_order").drop("_manifest_order").write_csv(
        manifest,
        separator="\t",
        quote_style="never",
    )
    hits = tmp_path / "blast_hits.tsv"
    blast_rows = (
        [
            ("candidate_kmer_000002", "123"),
            ("candidate_kmer_000003", "N/A"),
            ("candidate_kmer_000004", "88456;123"),
            ("candidate_kmer_000005", "999"),
            ("candidate_kmer_000006", "N/A"),
        ]
        if all_rejected
        else [
            ("candidate_kmer_000002", "88456"),
            ("candidate_kmer_000002", "88456"),
            ("candidate_kmer_000004", "123"),
            ("candidate_kmer_000005", "88456"),
            ("candidate_kmer_000005", "88456;123"),
            ("candidate_kmer_000006", "N/A"),
        ]
    )
    hits.write_text(
        "\t".join(screening.BLAST_FIELDS)
        + "\n"
        + "".join(
            f"{candidate_kmer_id}\tACCESSION.1\t{taxids}\tname\t100\t5\t0\t0\t1\t5\t10\t14\ttitle\n"
            for candidate_kmer_id, taxids in blast_rows
        ),
    )
    return candidates, manifest, hits


def command_args(tmp_path: Path, candidates: Path, manifest: Path, hits: Path) -> list[str]:
    return [
        "--target-taxid",
        "88456",
        "--kmer-size",
        "5",
        "--candidates",
        str(candidates),
        "--blast-hits",
        str(hits),
        "--manifest-in",
        str(manifest),
        "--manifest-out",
        str(tmp_path / "screened_candidate_kmers.tsv"),
        "--selected-candidates-out",
        str(tmp_path / "screened_candidates.fasta"),
    ]


def test_taxonomic_screening_preserves_prior_outcomes_and_marks_off_target_matches() -> None:
    hits = hit_frame(
        [
            ("candidate_kmer_000005", False),
            ("candidate_kmer_000002", True),
            ("candidate_kmer_000004", False),
            ("candidate_kmer_000005", True),
            ("candidate_kmer_000002", True),
            ("candidate_kmer_000006", False),
        ],
    )

    result = screening.construct_screening_result(candidate_frame(), manifest_frame(), hits)
    manifest, survivors = pl.collect_all([result.manifest, result.survivors])

    assert survivors.select("candidate_kmer_id", "kmer").rows() == [
        ("candidate_kmer_000002", "AACGT"),
        ("candidate_kmer_000003", "ACGTA"),
    ]
    assert manifest.select("kmer", "outcome").rows() == [
        (KMER_REJECTED, "low_complexity"),
        ("AACGT", "selected"),
        ("ACGTA", "selected"),
        ("AGTCC", "off_target_match"),
        ("CCGTA", "off_target_match"),
        ("CGTAC", "off_target_match"),
    ]


def test_no_database_hits_allow_every_candidate_to_pass() -> None:
    result = screening.construct_screening_result(candidate_frame(), manifest_frame(), hit_frame([]))
    assert result.manifest.collect().filter(pl.col("outcome") == "selected").height == len(KMERS)


def test_manifest_and_candidate_fasta_relations_must_be_bijective() -> None:
    missing_candidate = candidate_frame().filter(
        pl.col("candidate_kmer_id") != "candidate_kmer_000006",
    )
    with pytest.raises(screening.TaxonomicScreeningError, match="absent from the candidate FASTA"):
        screening.construct_screening_result(missing_candidate, manifest_frame(), hit_frame([]))

    mismatched_sequence = candidate_frame().with_columns(
        pl.when(pl.col("candidate_kmer_id") == "candidate_kmer_000002")
        .then(pl.lit("AAAAC"))
        .otherwise(pl.col("kmer"))
        .alias("kmer"),
    )
    with pytest.raises(screening.TaxonomicScreeningError, match="absent from the candidate FASTA"):
        screening.construct_screening_result(mismatched_sequence, manifest_frame(), hit_frame([]))


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("pident", "99", "non-exact full-length"),
        ("length", "4", "non-exact full-length"),
        ("mismatch", "1", "non-exact full-length"),
        ("gapopen", "1", "non-exact full-length"),
        ("qstart", "2", "non-exact full-length"),
        ("qend", "4", "non-exact full-length"),
        ("sstart", "0", "non-exact full-length"),
        ("send", "10", "non-exact full-length"),
        ("length", "", "invalid integer field"),
        ("length", "999999999999999999999999", "invalid integer field"),
        ("staxids", "88456;", "malformed staxids"),
    ],
)
def test_blast_ingestion_rejects_invalid_exact_matches(
    tmp_path: Path,
    field: str,
    value: str,
    message: str,
) -> None:
    row = dict.fromkeys(screening.BLAST_FIELDS, "value")
    row.update(
        {
            "qseqid": "candidate_kmer_000002",
            "staxids": "88456",
            "pident": "100",
            "length": "5",
            "mismatch": "0",
            "gapopen": "0",
            "qstart": "1",
            "qend": "5",
            "sstart": "10",
            "send": "14",
            field: value,
        },
    )
    path = tmp_path / "hits.tsv"
    pl.DataFrame([row], schema=dict.fromkeys(screening.BLAST_FIELDS, pl.String)).write_csv(
        path,
        separator="\t",
        quote_style="never",
    )
    with pytest.raises(screening.TaxonomicScreeningError, match=message):
        screening.scan_blast_hits(path, 5, 88456)


@pytest.mark.parametrize("alias_kind", ["same_path", "symlink", "hardlink"])
@pytest.mark.parametrize("output_flag", ["--manifest-out", "--selected-candidates-out"])
@pytest.mark.parametrize(
    "input_index",
    [0, 1, 2, 3],
    ids=["candidates", "manifest", "hits", "screening_scope"],
)
def test_command_rejects_output_aliasing_an_input_before_writing(
    tmp_path: Path,
    alias_kind: str,
    output_flag: str,
    input_index: int,
) -> None:
    command_inputs = write_inputs(tmp_path)
    scope = tmp_path / "taxonomic_screening_taxids.tsv"
    scope.write_text("taxid\n88456\n")
    protected_inputs = (*command_inputs, scope)
    originals = [path.read_bytes() for path in protected_inputs]
    args = command_args(tmp_path, *command_inputs)
    args.extend(["--taxonomic-screening-taxids", str(scope)])
    destination = protected_inputs[input_index]
    if alias_kind != "same_path":
        destination = tmp_path / "aliased_output"
        if alias_kind == "symlink":
            destination.symlink_to(protected_inputs[input_index])
        else:
            destination.hardlink_to(protected_inputs[input_index])
    args[args.index(output_flag) + 1] = str(destination)

    with pytest.raises(screening.TaxonomicScreeningError, match="alias"):
        screening.main(args)

    assert [path.read_bytes() for path in protected_inputs] == originals
    assert not (tmp_path / "screened_candidate_kmers.tsv").exists()
    assert not (tmp_path / "screened_candidates.fasta").exists()


@pytest.mark.parametrize("alias_kind", ["same_path", "symlink", "hardlink"])
def test_command_rejects_aliased_outputs_before_writing(tmp_path: Path, alias_kind: str) -> None:
    inputs = write_inputs(tmp_path)
    args = command_args(tmp_path, *inputs)
    manifest_output = tmp_path / "screened_candidate_kmers.tsv"
    fasta_output = tmp_path / "screened_candidates.fasta"
    if alias_kind == "same_path":
        args[args.index("--selected-candidates-out") + 1] = str(manifest_output)
    elif alias_kind == "symlink":
        fasta_output.symlink_to(manifest_output)
    else:
        manifest_output.write_text("existing output\n")
        fasta_output.hardlink_to(manifest_output)

    with pytest.raises(screening.TaxonomicScreeningError, match="alias"):
        screening.main(args)

    if alias_kind == "hardlink":
        assert manifest_output.read_text() == "existing output\n"
        assert fasta_output.read_text() == "existing output\n"
    else:
        assert not manifest_output.exists()
        assert not fasta_output.exists()


def test_command_writes_candidate_outcomes_and_selected_baits(tmp_path: Path) -> None:
    candidates, manifest, hits = write_inputs(tmp_path)
    screening.main(command_args(tmp_path, candidates, manifest, hits))

    assert (tmp_path / "screened_candidates.fasta").read_text() == (
        ">candidate_kmer_000002\nAACGT\n>candidate_kmer_000003\nACGTA\n"
    )
    assert pl.read_csv(tmp_path / "screened_candidate_kmers.tsv", separator="\t").get_column(
        "outcome",
    ).to_list() == [
        "low_complexity",
        "selected",
        "selected",
        "off_target_match",
        "off_target_match",
        "off_target_match",
    ]

    terminal = tmp_path / "terminal"
    terminal.mkdir()
    candidates, manifest, hits = write_inputs(terminal, all_rejected=True)
    screening.main(command_args(terminal, candidates, manifest, hits))
    assert (terminal / "screened_candidates.fasta").read_text() == ""
    terminal_manifest = pl.read_csv(terminal / "screened_candidate_kmers.tsv", separator="\t")
    assert terminal_manifest.filter(pl.col("outcome") == "off_target_match").height == len(
        KMERS,
    )


def test_command_treats_supplied_screening_scope_as_authoritative(tmp_path: Path) -> None:
    candidates, manifest, hits = write_inputs(tmp_path)
    scope = tmp_path / "taxonomic_screening_taxids.tsv"
    scope.write_text("taxid\n44417\n")
    hits.write_text(hits.read_text().replace("88456", "44417"))

    args = command_args(tmp_path, candidates, manifest, hits)
    args.extend(["--taxonomic-screening-taxids", str(scope)])
    screening.main(args)

    assert pl.read_csv(tmp_path / "screened_candidate_kmers.tsv", separator="\t").select(
        "kmer",
        "outcome",
    ).rows() == [
        (KMER_REJECTED, "low_complexity"),
        ("AACGT", "selected"),
        ("ACGTA", "selected"),
        ("AGTCC", "off_target_match"),
        ("CCGTA", "off_target_match"),
        ("CGTAC", "off_target_match"),
    ]


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("target_taxid\n88456\n", "exactly one header"),
        ("taxid\n", "canonical positive taxids"),
        ("taxid\n0\n", "canonical positive taxids"),
        ("taxid\n88456\n88456\n", "unique taxids"),
    ],
)
def test_taxonomic_screening_scope_rejects_invalid_rows(
    tmp_path: Path,
    content: str,
    message: str,
) -> None:
    scope = tmp_path / "taxonomic_screening_taxids.tsv"
    scope.write_text(content)

    with pytest.raises(screening.TaxonomicScreeningError, match=message):
        screening.load_taxonomic_screening_taxids(scope, 88456)


def test_blast_ingestion_accepts_reverse_strand_exact_subject_span(tmp_path: Path) -> None:
    row = dict(
        zip(
            screening.BLAST_FIELDS,
            (
                "candidate_kmer_000002",
                "ACCESSION.1",
                "88456",
                "name",
                "100",
                "5",
                "0",
                "0",
                "1",
                "5",
                "14",
                "10",
                "title",
            ),
            strict=True,
        ),
    )
    path = tmp_path / "reverse.tsv"
    pl.DataFrame([row], schema=dict.fromkeys(screening.BLAST_FIELDS, pl.String)).write_csv(
        path,
        separator="\t",
        quote_style="never",
    )
    assert screening.scan_blast_hits(path, 5, 88456).collect().height == 1


def test_manifest_ingestion_rejects_overflowing_source_count(tmp_path: Path) -> None:
    path = tmp_path / "manifest.tsv"
    manifest_frame().collect().sort("_manifest_order").drop("_manifest_order").with_columns(
        pl.when(pl.col("candidate_kmer_id") == "candidate_kmer_000001")
        .then(pl.lit("999999999999999999999999"))
        .otherwise(pl.col("source_copy_count").cast(pl.String))
        .alias("source_copy_count"),
    ).write_csv(path, separator="\t", quote_style="never")
    with pytest.raises(screening.TaxonomicScreeningError, match="invalid source_copy_count"):
        screening.scan_manifest(path, 5)


def test_fasta_ingestion_requires_canonical_candidate_ids_and_sequences() -> None:
    with pytest.raises(screening.TaxonomicScreeningError, match="malformed"):
        screening.construct_candidates(
            [
                SeqRecord(
                    Seq(""),
                    id="candidate_kmer_000001",
                    description="candidate_kmer_000001",
                ),
            ],
            5,
        )
    with pytest.raises(screening.TaxonomicScreeningError, match="invalid candidate identifier"):
        screening.construct_candidates(
            [SeqRecord(Seq("AACGT"), id="candidate_000002", description="candidate_000002")],
            5,
        )
    with pytest.raises(screening.TaxonomicScreeningError, match="noncanonical"):
        screening.construct_candidates(
            [
                SeqRecord(
                    Seq("TACGG"),
                    id="candidate_kmer_000001",
                    description="candidate_kmer_000001",
                ),
            ],
            5,
        )
