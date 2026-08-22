from pathlib import Path

import apply_taxonomic_screening as screening
import polars as pl
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

KMER_REJECTED = "AAAAA"
KMERS = ("AACGT", "ACGTA", "AGTCC", "CCGTA", "CGTAC")


def bait_frame() -> pl.LazyFrame:
    return pl.LazyFrame(
        {
            "_bait_order": range(1, len(KMERS) + 1),
            "bait_id": [f"bait_{number:06d}" for number in range(1, len(KMERS) + 1)],
            "kmer": KMERS,
        },
        schema={"_bait_order": pl.UInt32, "bait_id": pl.String, "kmer": pl.String},
    )


def manifest_frame() -> pl.LazyFrame:
    rows = [
        {
            "candidate_kmer_id": "candidate_kmer_000001",
            "bait_id": "",
            "kmer": KMER_REJECTED,
            "source_copy_count": 1,
            "background_occurrences": 0,
            "outcome": "low_complexity",
            "_manifest_order": 1,
        },
        *[
            {
                "candidate_kmer_id": f"candidate_kmer_{number + 1:06d}",
                "bait_id": f"bait_{number:06d}",
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
            "qseqid": [bait_id for bait_id, _ in rows],
            "_on_target": [on_target for _, on_target in rows],
        },
        schema={"qseqid": pl.String, "_on_target": pl.Boolean},
    )


def write_inputs(tmp_path: Path, *, all_rejected: bool = False) -> tuple[Path, Path, Path]:
    baits = tmp_path / "baits.fasta"
    baits.write_text(
        "".join(
            f">bait_{number:06d}\n{kmer}\n"
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
            ("bait_000001", "123"),
            ("bait_000002", "N/A"),
            ("bait_000003", "88456;123"),
            ("bait_000004", "999"),
            ("bait_000005", "N/A"),
        ]
        if all_rejected
        else [
            ("bait_000001", "88456"),
            ("bait_000001", "88456"),
            ("bait_000003", "123"),
            ("bait_000004", "88456"),
            ("bait_000004", "88456;123"),
            ("bait_000005", "N/A"),
        ]
    )
    hits.write_text(
        "\t".join(screening.BLAST_FIELDS)
        + "\n"
        + "".join(
            f"{bait_id}\tACCESSION.1\t{taxids}\tname\t100\t5\t0\t0\t1\t5\t10\t14\ttitle\n"
            for bait_id, taxids in blast_rows
        ),
    )
    return baits, manifest, hits


def command_args(tmp_path: Path, baits: Path, manifest: Path, hits: Path) -> list[str]:
    return [
        "--target-taxid",
        "88456",
        "--kmer-size",
        "5",
        "--baits",
        str(baits),
        "--blast-hits",
        str(hits),
        "--manifest-in",
        str(manifest),
        "--manifest-out",
        str(tmp_path / "screened_candidate_kmers.tsv"),
        "--baits-out",
        str(tmp_path / "selected_baits.fasta"),
        "--decisions-out",
        str(tmp_path / "screening_decisions.tsv"),
    ]


def test_taxonomic_screening_preserves_prior_outcomes_and_marks_off_target_matches() -> None:
    hits = hit_frame(
        [
            ("bait_000004", False),
            ("bait_000001", True),
            ("bait_000003", False),
            ("bait_000004", True),
            ("bait_000001", True),
            ("bait_000005", False),
        ],
    )

    result = screening.construct_screening_result(bait_frame(), manifest_frame(), hits)
    manifest, survivors = pl.collect_all([result.manifest, result.survivors])

    assert survivors.select("bait_id", "kmer").rows() == [
        ("bait_000001", "AACGT"),
        ("bait_000002", "ACGTA"),
    ]
    assert manifest.select("kmer", "outcome").rows() == [
        (KMER_REJECTED, "low_complexity"),
        ("AACGT", "selected"),
        ("ACGTA", "selected"),
        ("AGTCC", "off_target_match"),
        ("CCGTA", "off_target_match"),
        ("CGTAC", "off_target_match"),
    ]


def test_no_database_hits_allow_every_bait_to_pass() -> None:
    result = screening.construct_screening_result(bait_frame(), manifest_frame(), hit_frame([]))
    assert result.decisions.collect().select("outcome").rows() == [("selected",)] * len(KMERS)


def test_manifest_and_bait_relations_must_be_bijective() -> None:
    missing_bait = bait_frame().filter(pl.col("bait_id") != "bait_000005")
    with pytest.raises(screening.TaxonomicScreeningError, match="absent from the bait FASTA"):
        screening.construct_screening_result(missing_bait, manifest_frame(), hit_frame([]))

    mismatched_sequence = bait_frame().with_columns(
        pl.when(pl.col("bait_id") == "bait_000001")
        .then(pl.lit("AAAAC"))
        .otherwise(pl.col("kmer"))
        .alias("kmer"),
    )
    with pytest.raises(screening.TaxonomicScreeningError, match="absent from the bait FASTA"):
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
            "qseqid": "bait_000001",
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


def test_command_writes_candidate_outcomes_and_selected_baits(tmp_path: Path) -> None:
    baits, manifest, hits = write_inputs(tmp_path)
    screening.main(command_args(tmp_path, baits, manifest, hits))

    assert (tmp_path / "selected_baits.fasta").read_text() == (
        ">bait_000001\nAACGT\n>bait_000002\nACGTA\n"
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
    baits, manifest, hits = write_inputs(terminal, all_rejected=True)
    screening.main(command_args(terminal, baits, manifest, hits))
    assert not (terminal / "selected_baits.fasta").exists()
    terminal_manifest = pl.read_csv(terminal / "screened_candidate_kmers.tsv", separator="\t")
    assert terminal_manifest.filter(pl.col("outcome") == "off_target_match").get_column(
        "bait_id",
    ).fill_null("").to_list() == [""] * len(KMERS)


def test_blast_ingestion_accepts_reverse_strand_exact_subject_span(tmp_path: Path) -> None:
    row = dict(
        zip(
            screening.BLAST_FIELDS,
            (
                "bait_000001",
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


def test_fasta_ingestion_requires_canonical_sequential_baits() -> None:
    with pytest.raises(screening.TaxonomicScreeningError, match="malformed"):
        screening.construct_baits(
            [SeqRecord(Seq(""), id="bait_000001", description="bait_000001")],
            5,
        )
    with pytest.raises(screening.TaxonomicScreeningError, match="sequential"):
        screening.construct_baits(
            [SeqRecord(Seq("AACGT"), id="bait_000002", description="bait_000002")],
            5,
        )
    with pytest.raises(screening.TaxonomicScreeningError, match="noncanonical"):
        screening.construct_baits(
            [SeqRecord(Seq("TACGG"), id="bait_000001", description="bait_000001")],
            5,
        )
