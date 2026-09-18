import shutil
import subprocess
from pathlib import Path

import polars as pl
from Bio import SeqIO

PROJECT = Path(__file__).resolve().parents[2]


def test_screening_preserves_upstream_artifacts_on_resume(tmp_path: Path) -> None:
    makeblastdb = shutil.which("makeblastdb")
    nextflow = shutil.which("nextflow")
    assert makeblastdb, "Run this integration test in the Pixi dev environment"
    assert nextflow, "Run this integration test in the Pixi dev environment"
    source = tmp_path / "source.fasta"
    source.write_text(">source\n" + "ACGTTGCATGTCAGTACGATCGTAG" * 4 + "\n")
    database = tmp_path / "database"
    database.mkdir()
    subprocess.run(
        [
            makeblastdb, "-dbtype", "nucl", "-parse_seqids", "-taxid", "123",
            "-in", str(source), "-out", str(database / "reference"),
        ],
        check=True, capture_output=True, text=True,
    )
    config = tmp_path / "trace.config"
    config.write_text("trace.fields = 'process,status,workdir'\n")
    command = [
        nextflow, "run", str(PROJECT / "main.nf"),
        "-profile", "test", "-ansi-log", "false", "-c", str(config),
        "--id", "screening_terminal", "--source_sequences", str(source),
        "--target_taxid", "88456", "--background", str(PROJECT / "tests/data/background.fasta"),
        "--taxon_ref_db", str(database),
    ]
    originals = None
    upstream_dir = None
    for execution in ("first", "resumed"):
        trace = tmp_path / f"{execution}.tsv"
        result = subprocess.run(
            [*command, "-with-trace", str(trace), *(["-resume"] if execution == "resumed" else [])],
            cwd=tmp_path, capture_output=True, text=True, check=False, timeout=180,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        tasks = pl.read_csv(trace, separator="\t")
        upstream = tasks.filter(pl.col("process").str.ends_with(":APPLY_COMPLEXITY_FILTER"))
        assert upstream.height == 1
        task_dir = Path(upstream.item(0, "workdir"))
        fasta = task_dir / "selected_candidates.fasta"
        manifest = task_dir / "candidate_kmers.tsv"
        selected = pl.read_csv(manifest, separator="\t").filter(pl.col("outcome") == "selected")
        records = {record.id: str(record.seq) for record in SeqIO.parse(fasta, "fasta")}
        assert records
        assert records == dict(selected.select("candidate_kmer_id", "kmer").iter_rows())
        artifacts = (fasta.read_bytes(), manifest.read_bytes())
        if execution == "first":
            assert upstream.item(0, "status") == "COMPLETED"
            upstream_dir = task_dir
            originals = artifacts
        else:
            assert upstream.item(0, "status") == "CACHED"
            assert task_dir == upstream_dir
            assert artifacts == originals
        assert (tmp_path / "results/artifacts/screening_terminal/baits.fasta").read_bytes() == b""
