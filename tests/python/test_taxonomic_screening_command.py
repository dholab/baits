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
            "status": "REJECT_LOW_COMPLEXITY",
            "rejection_reason": "low_complexity",
            "taxonomic_screening_status": "NOT_APPLICABLE",
            "on_target_hits": "",
            "off_target_hits": "",
            "_manifest_order": 1,
        },
        *[
            {
                "candidate_kmer_id": f"candidate_kmer_{number + 1:06d}",
                "bait_id": f"bait_{number:06d}",
                "kmer": kmer,
                "source_copy_count": 1,
                "background_occurrences": 0,
                "status": "PASS",
                "rejection_reason": "none",
                "taxonomic_screening_status": "NOT_RUN",
                "on_target_hits": "",
                "off_target_hits": "",
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


def write_inputs(tmp_path: Path, *, all_rejected: bool = False) -> tuple[Path, Path, Path, Path]:
    baits = tmp_path / "locally_filtered_baits.fasta"
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
    status = tmp_path / "bait_set_status.tsv"
    status.write_text(
        "metric\tvalue\n"
        "design_id\tdesign\n"
        "source_sequence_origin\tcurated_input\n"
        "candidate_kmer_count\t6\n"
        "locally_filtered_bait_count\t5\n"
        "taxonomic_screening_status\tNOT_RUN\n"
        "taxonomically_screened_bait_count\t\n"
        "deepest_bait_set\tlocally_filtered\n"
        "deacon_index_source\t\n",
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
    return baits, manifest, status, hits


def command_args(tmp_path: Path, baits: Path, manifest: Path, status: Path, hits: Path) -> list[str]:
    return [
        "--design-id",
        "design",
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
        "--bait-set-status-in",
        str(status),
        "--manifest-out",
        str(tmp_path / "screened_candidate_kmers.tsv"),
        "--baits-out",
        str(tmp_path / "taxonomically_screened_baits.fasta"),
        "--decisions-out",
        str(tmp_path / "screening_decisions.tsv"),
        "--screening-status-out",
        str(tmp_path / "screening_status.tsv"),
        "--bait-set-status-out",
        str(tmp_path / "screened_bait_set_status.tsv"),
        "--terminal-bait-set-status-out",
        str(tmp_path / "terminal_bait_set_status.tsv"),
    ]


def test_screening_transformation_preserves_bait_order_and_counts_blast_rows() -> None:
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
    manifest, decisions, survivors = pl.collect_all(
        [result.manifest, result.decisions, result.survivors],
    )

    assert decisions.rows() == [
        ("bait_000001", "PASS", 2, 0),
        ("bait_000002", "PASS", 0, 0),
        ("bait_000003", "REJECT_OFF_TARGET_HIT", 0, 1),
        ("bait_000004", "REJECT_OFF_TARGET_HIT", 1, 1),
        ("bait_000005", "REJECT_OFF_TARGET_HIT", 0, 1),
    ]
    assert survivors.select("bait_id", "kmer").rows() == [
        ("bait_000001", "AACGT"),
        ("bait_000002", "ACGTA"),
    ]
    assert manifest.select("status", "taxonomic_screening_status").rows() == [
        ("REJECT_LOW_COMPLEXITY", "NOT_APPLICABLE"),
        ("PASS", "PASS"),
        ("PASS", "PASS"),
        ("REJECT_OFF_TARGET_HIT", "REJECT_OFF_TARGET_HIT"),
        ("REJECT_OFF_TARGET_HIT", "REJECT_OFF_TARGET_HIT"),
        ("REJECT_OFF_TARGET_HIT", "REJECT_OFF_TARGET_HIT"),
    ]


def test_no_database_hits_allow_every_bait_to_pass() -> None:
    result = screening.construct_screening_result(bait_frame(), manifest_frame(), hit_frame([]))
    assert result.decisions.collect().select(
        "taxonomic_screening_status",
        "on_target_hits",
        "off_target_hits",
    ).rows() == [("PASS", 0, 0)] * len(KMERS)


def test_manifest_and_bait_relations_must_be_bijective() -> None:
    missing_bait = bait_frame().filter(pl.col("bait_id") != "bait_000005")
    with pytest.raises(screening.TaxonomicScreeningError, match="absent from the Bait FASTA"):
        screening.construct_screening_result(missing_bait, manifest_frame(), hit_frame([]))

    mismatched_sequence = bait_frame().with_columns(
        pl.when(pl.col("bait_id") == "bait_000001")
        .then(pl.lit("AAAAC"))
        .otherwise(pl.col("kmer"))
        .alias("kmer"),
    )
    with pytest.raises(screening.TaxonomicScreeningError, match="absent from the Bait FASTA"):
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
def test_blast_ingestion_rejects_invalid_exact_hit_evidence(
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


def test_command_writes_canonical_outputs_and_terminal_evidence(tmp_path: Path) -> None:
    baits, manifest, status, hits = write_inputs(tmp_path)
    screening.main(command_args(tmp_path, baits, manifest, status, hits))

    assert pl.read_csv(tmp_path / "screening_decisions.tsv", separator="\t").rows() == [
        ("bait_000001", "PASS", 2, 0),
        ("bait_000002", "PASS", 0, 0),
        ("bait_000003", "REJECT_OFF_TARGET_HIT", 0, 1),
        ("bait_000004", "REJECT_OFF_TARGET_HIT", 1, 1),
        ("bait_000005", "REJECT_OFF_TARGET_HIT", 0, 1),
    ]
    assert (tmp_path / "taxonomically_screened_baits.fasta").read_text() == (
        ">bait_000001\nAACGT\n>bait_000002\nACGTA\n"
    )
    assert not (tmp_path / "terminal_bait_set_status.tsv").exists()

    terminal = tmp_path / "terminal"
    terminal.mkdir()
    baits, manifest, status, hits = write_inputs(terminal, all_rejected=True)
    screening.main(command_args(terminal, baits, manifest, status, hits))
    assert not (terminal / "taxonomically_screened_baits.fasta").exists()
    assert pl.read_csv(terminal / "screening_status.tsv", separator="\t").rows()[-2:] == [
        ("taxonomic_screening_status", "NO_BAITS"),
        ("taxonomically_screened_bait_count", "0"),
    ]
    assert (terminal / "terminal_bait_set_status.tsv").read_text() == (
        terminal / "screened_bait_set_status.tsv"
    ).read_text()

    baits, manifest, status, hits = write_inputs(terminal)
    screening.main(command_args(terminal, baits, manifest, status, hits))
    assert (terminal / "taxonomically_screened_baits.fasta").exists()
    assert not (terminal / "terminal_bait_set_status.tsv").exists()


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
