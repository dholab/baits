import subprocess
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
MANIFEST_HEADER = (
    "candidate_kmer_id\tbait_id\tkmer\tsource_copy_count\tbackground_occurrences\toutcome\n"
)


def write_inputs(
    tmp_path: Path,
    *,
    selected_kmer: str = "ACGTA",
    background: str = "",
    off_target_bait_id: str = "",
) -> list[str]:
    baits = tmp_path / "baits.fasta"
    baits.write_text(">private-bait\nACGTA\n")
    bait_roundtrip = tmp_path / "bait-roundtrip.fasta"
    bait_roundtrip.write_text(">deacon-output\nACGTA\n")
    background_roundtrip = tmp_path / "background-roundtrip.fasta"
    background_roundtrip.write_text(background)
    manifest = tmp_path / "candidate-kmers.tsv"
    manifest.write_text(
        MANIFEST_HEADER
        + f"candidate-a\tprivate-bait\t{selected_kmer}\t1\t0\tselected\n"
        + "candidate-b\t\tAAAAA\t2\t0\tlow_complexity\n"
        + "candidate-c\t\tAACGT\t1\t3\tbackground_match\n"
        + f"candidate-d\t{off_target_bait_id}\tACCGT\t1\t0\toff_target_match\n",
    )
    return [
        "--kmer-size",
        "5",
        "--baits",
        str(baits),
        "--bait-roundtrip",
        str(bait_roundtrip),
        "--background-roundtrip",
        str(background_roundtrip),
        "--manifest",
        str(manifest),
    ]


def run_verifier(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(PROJECT / "bin/verify_deacon_index.py"), *args],
        cwd=PROJECT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_verifier_accepts_a_valid_candidate_manifest_and_round_trips(tmp_path: Path) -> None:
    result = run_verifier(write_inputs(tmp_path))

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""


def test_verifier_rejects_a_bait_and_manifest_mismatch(tmp_path: Path) -> None:
    result = run_verifier(write_inputs(tmp_path, selected_kmer="AACGA"))

    assert result.returncode != 0
    assert "Selected candidate k-mers do not equal the bait FASTA" in result.stderr


def test_verifier_rejects_a_nonempty_background_roundtrip(tmp_path: Path) -> None:
    result = run_verifier(write_inputs(tmp_path, background=">retained\nACGTACGT\n"))

    assert result.returncode != 0
    assert "Background round-trip FASTA must be empty" in result.stderr


def test_verifier_rejects_a_bait_roundtrip_sequence_mismatch(tmp_path: Path) -> None:
    args = write_inputs(tmp_path)
    (tmp_path / "bait-roundtrip.fasta").write_text(">deacon-output\nAACGT\n")

    result = run_verifier(args)

    assert result.returncode != 0
    assert "Bait sequence set differs after the Deacon round trip" in result.stderr


def test_verifier_rejects_an_unknown_candidate_outcome(tmp_path: Path) -> None:
    args = write_inputs(tmp_path)
    manifest = tmp_path / "candidate-kmers.tsv"
    manifest.write_text(manifest.read_text().replace("\tselected\n", "\tlegacy_selected\n"))

    result = run_verifier(args)

    assert result.returncode != 0
    assert "unknown outcome" in result.stderr


def test_verifier_rejects_a_bait_identifier_on_an_off_target_candidate(tmp_path: Path) -> None:
    result = run_verifier(write_inputs(tmp_path, off_target_bait_id="removed-bait"))

    assert result.returncode != 0
    assert "incoherent candidate outcome" in result.stderr
