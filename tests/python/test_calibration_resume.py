import csv
import subprocess
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]


def run_command(command: list[str], cwd: Path) -> None:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"


def task_statuses(trace: Path, process_name: str) -> dict[str, str]:
    with trace.open(newline="") as lines:
        rows = csv.DictReader(lines, delimiter="\t")
        return {
            row["name"].rsplit(" (", maxsplit=1)[1].removesuffix(")"): row["status"]
            for row in rows
            if process_name in row["name"]
        }


def process_status_counts(trace: Path, process_name: str) -> Counter[str]:
    return Counter(task_statuses(trace, process_name).values())


def write_fastq(path: Path, read_id: str, sequence: str) -> None:
    path.write_text(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def test_adding_calibration_fastq_reuses_existing_per_stream_tasks(tmp_path: Path) -> None:
    sequence = "ACGTTGCATGTCAGTACGATCGTAG" * 4
    source = tmp_path / "source.fasta"
    source.write_text(f">target\n{sequence}\n")
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    write_fastq(reads / "sample_b.fastq", "sample_b", sequence)
    write_fastq(reads / "sample_c.fastq", "sample_c", sequence)
    taxonomic_reference_db = tmp_path / "taxonomic_reference_db"
    taxonomic_reference_db.mkdir()
    run_command(
        [
            "makeblastdb",
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-taxid",
            "88456",
            "-in",
            str(source),
            "-out",
            str(taxonomic_reference_db / "taxonomic"),
        ],
        tmp_path,
    )

    work = tmp_path / "work"
    common_command = [
        "nextflow",
        "run",
        str(PROJECT_ROOT / "main.nf"),
        "-profile",
        "test",
        "-work-dir",
        str(work),
        "--id",
        "resume_probe",
        "--source_sequences",
        str(source),
        "--target_taxid",
        "88456",
        "--background",
        str(PROJECT_ROOT / "tests/data/interference_background.fasta"),
        "--calibration_reads",
        str(reads),
        "--taxon_ref_db",
        str(taxonomic_reference_db),
    ]
    run_command([*common_command, "-with-trace", str(tmp_path / "trace1.txt")], tmp_path)

    write_fastq(reads / "sample_a.fastq", "sample_a", sequence)
    run_command(
        [
            *common_command,
            "-resume",
            "-with-trace",
            str(tmp_path / "trace2.txt"),
        ],
        tmp_path,
    )

    trace = tmp_path / "trace2.txt"
    assert process_status_counts(trace, "RESOLVE_CALIBRATION_READS") == Counter(
        {"COMPLETED": 1},
    )
    expected_per_stream_statuses = {
        "resume_probe__sample_a": "COMPLETED",
        "resume_probe__sample_b": "CACHED",
        "resume_probe__sample_c": "CACHED",
    }
    assert task_statuses(trace, "DEACON_RETRIEVE_CALIBRATION_READS") == (
        expected_per_stream_statuses
    )
    assert task_statuses(trace, "COUNT_READ_BAITS") == expected_per_stream_statuses
