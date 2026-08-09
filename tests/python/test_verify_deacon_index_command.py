import csv
import hashlib
import subprocess
from pathlib import Path

import polars as pl
import pytest
import verify_deacon_index as verification

PROJECT = Path(__file__).resolve().parents[2]
MANIFEST_HEADER = (
    "candidate_kmer_id\tbait_id\tkmer\tsource_copy_count\tbackground_occurrences\t"
    "status\trejection_reason\ttaxonomic_screening_status\ton_target_hits\toff_target_hits\n"
)


def test_verification_accepts_reordered_roundtrip_sequences_with_incidental_ids() -> None:
    baits = pl.LazyFrame(
        {
            "record_id": ["bait_000001", "bait_000002"],
            "sequence": ["ACGTA", "CGTAC"],
        },
    )
    roundtrip = pl.LazyFrame(
        {
            "record_id": ["returned_2", "returned_1"],
            "sequence": ["CGTAC", "ACGTA"],
        },
    )
    background = pl.LazyFrame(
        schema=verification.FASTA_SCHEMA,
    )
    manifest = pl.LazyFrame(
        {
            "candidate_kmer_id": ["candidate_kmer_000001", "candidate_kmer_000002"],
            "bait_id": ["bait_000001", "bait_000002"],
            "kmer": ["ACGTA", "CGTAC"],
            "status": ["PASS", "PASS"],
            "rejection_reason": ["none", "none"],
            "taxonomic_screening_status": ["NOT_RUN", "NOT_RUN"],
        },
    )
    status = verification.BaitSetStatus(
        design_id="design",
        source_sequence_origin=verification.SourceSequenceOrigin.CURATED_INPUT,
        candidate_kmer_count=2,
        locally_filtered_bait_count=2,
        taxonomic_screening_status=verification.TaxonomicScreeningStatus.NOT_RUN,
        taxonomically_screened_bait_count=None,
        deepest_bait_set=verification.BaitSetSource.LOCALLY_FILTERED,
        deacon_index_source=None,
    )

    result = verification.construct_verification_result(
        verification.VerificationRelations(
            baits=baits,
            roundtrip=roundtrip,
            background=background,
            manifest=manifest,
        ),
        status,
        kmer_size=5,
        deacon_window=1,
    )

    assert result.bait_set_source is verification.BaitSetSource.LOCALLY_FILTERED
    assert result.bait_count == 2
    assert result.roundtrip_count == 2


def read_key_values(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [(row["metric"], row["value"]) for row in rows]


def test_bait_set_status_rejects_an_overflowing_count(tmp_path: Path) -> None:
    status = tmp_path / "bait_set_status.tsv"
    status.write_text(
        "metric\tvalue\n"
        "design_id\tdesign\n"
        "source_sequence_origin\tcurated_input\n"
        "candidate_kmer_count\t9223372036854775808\n"
        "locally_filtered_bait_count\t2\n"
        "taxonomic_screening_status\tNOT_RUN\n"
        "taxonomically_screened_bait_count\t\n"
        "deepest_bait_set\tlocally_filtered\n"
        "deacon_index_source\t\n",
    )

    with pytest.raises(verification.DeaconVerificationError, match="candidate_kmer_count"):
        verification.read_bait_set_status(status, "design")


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("source_copy_count", "", "source_copy_count"),
        ("source_copy_count", "9223372036854775808", "source_copy_count"),
        ("background_occurrences", "-1", "background_occurrences"),
        ("on_target_hits", "9223372036854775808", "candidate state"),
    ],
)
def test_manifest_ingestion_rejects_invalid_numeric_evidence(
    tmp_path: Path,
    field: str,
    value: str,
    message: str,
) -> None:
    row = {
        "candidate_kmer_id": "candidate_kmer_000001",
        "bait_id": "bait_000001",
        "kmer": "ACGTA",
        "source_copy_count": "1",
        "background_occurrences": "0",
        "status": "PASS",
        "rejection_reason": "none",
        "taxonomic_screening_status": "PASS",
        "on_target_hits": "1",
        "off_target_hits": "0",
    }
    row[field] = value
    manifest = tmp_path / "candidate_kmers.tsv"
    pl.DataFrame([row], schema=dict.fromkeys(verification.MANIFEST_FIELDS, pl.String)).write_csv(
        manifest,
        separator="\t",
        quote_style="never",
    )

    with pytest.raises(verification.DeaconVerificationError, match=message):
        verification.scan_manifest(manifest, 5)


def test_manifest_ingestion_rejects_a_noncanonical_kmer(tmp_path: Path) -> None:
    manifest = tmp_path / "candidate_kmers.tsv"
    manifest.write_text(
        MANIFEST_HEADER
        + "candidate_kmer_000001\tbait_000001\tTACGT\t1\t0\tPASS\tnone\tNOT_RUN\t\t\n",
    )

    with pytest.raises(verification.DeaconVerificationError, match="noncanonical"):
        verification.scan_manifest(manifest, 5)


def test_roundtrip_fasta_rejects_a_duplicate_sequence(tmp_path: Path) -> None:
    roundtrip = tmp_path / "bait_roundtrip.fasta"
    roundtrip.write_text(">first\nACGTA\n>second\nACGTA\n")

    with pytest.raises(verification.DeaconVerificationError, match="sequences must be unique"):
        verification.read_fasta(roundtrip, "bait round-trip FASTA", 5)


def test_verifier_publishes_compact_evidence_and_finalizes_the_bait_set(tmp_path: Path) -> None:
    baits = tmp_path / "locally_filtered_baits.fasta"
    baits.write_text(">bait_000001\nACGTA\n>bait_000002\nCGTAC\n")
    roundtrip = tmp_path / "bait_roundtrip.fasta"
    roundtrip.write_text(">returned_2\nCGTAC\n>returned_1\nACGTA\n")
    background_roundtrip = tmp_path / "background_roundtrip.fasta"
    background_roundtrip.write_text("")
    manifest = tmp_path / "candidate_kmers.tsv"
    manifest.write_text(
        MANIFEST_HEADER
        + "candidate_kmer_000001\t\tAAAAA\t1\t0\tREJECT_LOW_COMPLEXITY\tlow_complexity\tNOT_APPLICABLE\t\t\n"
        + "candidate_kmer_000002\tbait_000001\tACGTA\t1\t0\tPASS\tnone\tNOT_RUN\t\t\n"
        + "candidate_kmer_000003\tbait_000002\tCGTAC\t1\t0\tPASS\tnone\tNOT_RUN\t\t\n",
    )
    status = tmp_path / "bait_set_status.tsv"
    status.write_text(
        "metric\tvalue\n"
        "design_id\tdesign\n"
        "source_sequence_origin\tcurated_input\n"
        "candidate_kmer_count\t3\n"
        "locally_filtered_bait_count\t2\n"
        "taxonomic_screening_status\tNOT_RUN\n"
        "taxonomically_screened_bait_count\t\n"
        "deepest_bait_set\tlocally_filtered\n"
        "deacon_index_source\t\n",
    )
    summary = tmp_path / "verification_summary.tsv"
    report = tmp_path / "verification_report.md"
    final_status = tmp_path / "final_bait_set_status.tsv"

    subprocess.run(
        [
            str(PROJECT / "bin/verify_deacon_index.py"),
            "--design-id", "design",
            "--kmer-size", "5",
            "--deacon-window", "1",
            "--baits", str(baits),
            "--bait-roundtrip", str(roundtrip),
            "--background-roundtrip", str(background_roundtrip),
            "--manifest", str(manifest),
            "--bait-set-status-in", str(status),
            "--bait-set-status-out", str(final_status),
            "--summary-out", str(summary),
            "--report-out", str(report),
        ],
        cwd=PROJECT,
        check=True,
    )

    assert read_key_values(summary) == [
        ("design_id", "design"),
        ("kmer_size", "5"),
        ("deacon_window", "1"),
        ("deacon_index_entropy_threshold", "0"),
        ("deacon_filter_absolute_threshold", "1"),
        ("deacon_filter_relative_threshold", "0"),
        ("bait_roundtrip_record_count", "2"),
        ("bait_roundtrip_sha256", hashlib.sha256(roundtrip.read_bytes()).hexdigest()),
        ("bait_sequence_sets_equal", "true"),
        ("interference_background_roundtrip_record_count", "0"),
        ("interference_background_roundtrip_sha256", hashlib.sha256(b"").hexdigest()),
        (
            "conclusion",
            "The Deacon index reproduces the bait set and retains no interference background records.",
        ),
    ]
    assert read_key_values(final_status) == [
        ("design_id", "design"),
        ("source_sequence_origin", "curated_input"),
        ("candidate_kmer_count", "3"),
        ("locally_filtered_bait_count", "2"),
        ("taxonomic_screening_status", "NOT_RUN"),
        ("taxonomically_screened_bait_count", ""),
        ("deepest_bait_set", "locally_filtered"),
        ("deacon_index_source", "locally_filtered"),
    ]
    assert "- Source: locally filtered bait set" in report.read_text()


def test_verifier_finalizes_a_taxonomically_screened_bait_set(tmp_path: Path) -> None:
    baits = tmp_path / "taxonomically_screened_baits.fasta"
    baits.write_text(">bait_000001\nACGTA\n")
    roundtrip = tmp_path / "bait_roundtrip.fasta"
    roundtrip.write_text(">returned\nACGTA\n")
    background_roundtrip = tmp_path / "background_roundtrip.fasta"
    background_roundtrip.write_text("")
    manifest = tmp_path / "candidate_kmers.tsv"
    manifest.write_text(
        MANIFEST_HEADER
        + "candidate_kmer_000001\tbait_000001\tACGTA\t1\t0\tPASS\tnone\tPASS\t1\t0\n"
        + "candidate_kmer_000002\tbait_000002\tCGTAC\t1\t0\tREJECT_OFF_TARGET_HIT\toff_target_exact_match\tREJECT_OFF_TARGET_HIT\t0\t1\n",
    )
    status = tmp_path / "bait_set_status.tsv"
    status.write_text(
        "metric\tvalue\n"
        "design_id\tdesign\n"
        "source_sequence_origin\tcurated_input\n"
        "candidate_kmer_count\t2\n"
        "locally_filtered_bait_count\t2\n"
        "taxonomic_screening_status\tPASS\n"
        "taxonomically_screened_bait_count\t1\n"
        "deepest_bait_set\ttaxonomically_screened\n"
        "deacon_index_source\t\n",
    )
    final_status = tmp_path / "final_bait_set_status.tsv"

    verification.main(
        [
            "--design-id",
            "design",
            "--kmer-size",
            "5",
            "--deacon-window",
            "1",
            "--baits",
            str(baits),
            "--bait-roundtrip",
            str(roundtrip),
            "--background-roundtrip",
            str(background_roundtrip),
            "--manifest",
            str(manifest),
            "--bait-set-status-in",
            str(status),
            "--bait-set-status-out",
            str(final_status),
            "--summary-out",
            str(tmp_path / "summary.tsv"),
            "--report-out",
            str(tmp_path / "report.md"),
        ],
    )

    assert read_key_values(final_status) == [
        ("design_id", "design"),
        ("source_sequence_origin", "curated_input"),
        ("candidate_kmer_count", "2"),
        ("locally_filtered_bait_count", "2"),
        ("taxonomic_screening_status", "PASS"),
        ("taxonomically_screened_bait_count", "1"),
        ("deepest_bait_set", "taxonomically_screened"),
        ("deacon_index_source", "taxonomically_screened"),
    ]


def test_verifier_rejects_retained_interference_background_records(tmp_path: Path) -> None:
    baits = tmp_path / "baits.fasta"
    baits.write_text(">bait_000001\nACGTA\n")
    roundtrip = tmp_path / "bait_roundtrip.fasta"
    roundtrip.write_text(">bait_000001\nACGTA\n")
    background_roundtrip = tmp_path / "background_roundtrip.fasta"
    background_roundtrip.write_text(">retained_background_record\nACGTACGT\n")
    manifest = tmp_path / "candidate_kmers.tsv"
    manifest.write_text(
        MANIFEST_HEADER
        + "candidate_kmer_000001\tbait_000001\tACGTA\t1\t0\tPASS\tnone\tNOT_RUN\t\t\n",
    )
    status = tmp_path / "bait_set_status.tsv"
    status.write_text(
        "metric\tvalue\n"
        "design_id\tdesign\n"
        "source_sequence_origin\tcurated_input\n"
        "candidate_kmer_count\t1\n"
        "locally_filtered_bait_count\t1\n"
        "taxonomic_screening_status\tNOT_RUN\n"
        "taxonomically_screened_bait_count\t\n"
        "deepest_bait_set\tlocally_filtered\n"
        "deacon_index_source\t\n",
    )

    result = subprocess.run(
        [
            str(PROJECT / "bin/verify_deacon_index.py"),
            "--design-id", "design",
            "--kmer-size", "5",
            "--deacon-window", "1",
            "--baits", str(baits),
            "--bait-roundtrip", str(roundtrip),
            "--background-roundtrip", str(background_roundtrip),
            "--manifest", str(manifest),
            "--bait-set-status-in", str(status),
            "--bait-set-status-out", str(tmp_path / "final.tsv"),
            "--summary-out", str(tmp_path / "summary.tsv"),
            "--report-out", str(tmp_path / "report.md"),
        ],
        cwd=PROJECT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert "retained interference background records" in result.stderr


def test_verifier_rejects_a_bait_roundtrip_sequence_set_mismatch(tmp_path: Path) -> None:
    roundtrip = tmp_path / "bait_roundtrip.fasta"
    roundtrip.write_text(">bait_000001\nACGTTGCATGTCAGTACGATCGTAGCTAGCA\n")
    background_roundtrip = tmp_path / "background_roundtrip.fasta"
    background_roundtrip.write_text("")

    result = subprocess.run(
        [
            str(PROJECT / "bin/verify_deacon_index.py"),
            "--design-id", "design",
            "--kmer-size", "31",
            "--deacon-window", "1",
            "--baits", str(PROJECT / "tests/data/taxonomic_locally_filtered_baits.fasta"),
            "--bait-roundtrip", str(roundtrip),
            "--background-roundtrip", str(background_roundtrip),
            "--manifest", str(PROJECT / "tests/data/taxonomic_candidate_kmers.tsv"),
            "--bait-set-status-in", str(PROJECT / "tests/data/taxonomic_bait_set_status.tsv"),
            "--bait-set-status-out", str(tmp_path / "final.tsv"),
            "--summary-out", str(tmp_path / "summary.tsv"),
            "--report-out", str(tmp_path / "report.md"),
        ],
        cwd=PROJECT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert "Bait sequence set differs" in result.stderr
    assert not (tmp_path / "summary.tsv").exists()
