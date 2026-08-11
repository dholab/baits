import gzip
from pathlib import Path

import pytest
from resolve_optimization_reads import (
    OptimizationReadSetError,
    ResolvedReadSet,
    construct_optimization_read_manifest,
    construct_optimization_read_sets,
    main,
)


def test_construct_optimization_read_sets_orders_valid_layouts(tmp_path: Path) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    alpha_read_1 = reads / "alpha_R1.fastq"
    alpha_read_1.write_text("@fragment/1\nACGTA\n+\nIIIII\n")
    alpha_read_2 = reads / "alpha_R2.fastq.gz"
    alpha_read_2.write_bytes(gzip.compress(b"@fragment/2\nTGCAT\n+\nIIIII\n", mtime=0))
    beta_read = reads / "beta.fq"
    beta_read.write_text("@single\nGATCC\n+\nIIIII\n")
    (reads / "README.txt").write_text("ignored\n")

    read_sets = construct_optimization_read_sets("design", tuple(reads.iterdir()))

    assert read_sets == (
        ResolvedReadSet("design__alpha", "design", "alpha", (alpha_read_1, alpha_read_2)),
        ResolvedReadSet("design__beta", "design", "beta", (beta_read,)),
    )
    assert construct_optimization_read_manifest(read_sets).collect().to_dicts() == [
        {
            "id": "design__alpha",
            "design_id": "design",
            "metagenome_id": "alpha",
            "read_1": "alpha_R1.fastq",
            "read_2": "alpha_R2.fastq.gz",
        },
        {
            "id": "design__beta",
            "design_id": "design",
            "metagenome_id": "beta",
            "read_1": "beta.fq",
            "read_2": "",
        },
    ]


def test_construct_optimization_read_sets_rejects_mixed_layouts(tmp_path: Path) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    (reads / "sample.fastq").touch()
    (reads / "sample_R1.fastq").touch()
    (reads / "sample_R2.fastq").touch()

    with pytest.raises(OptimizationReadSetError, match="mixes single-end and paired"):
        construct_optimization_read_sets("design", tuple(reads.iterdir()))


def test_construct_optimization_read_sets_rejects_a_missing_mate(tmp_path: Path) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    (reads / "sample_R1.fastq").touch()

    with pytest.raises(OptimizationReadSetError, match="missing R2"):
        construct_optimization_read_sets("design", tuple(reads.iterdir()))


def test_construct_optimization_read_sets_rejects_duplicate_roles(tmp_path: Path) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    (reads / "sample.fastq").touch()
    (reads / "sample.fastq.gz").touch()

    with pytest.raises(OptimizationReadSetError, match="duplicate files"):
        construct_optimization_read_sets("design", tuple(reads.iterdir()))


def test_construct_optimization_read_sets_rejects_nested_directories(tmp_path: Path) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    (reads / "nested").mkdir()

    with pytest.raises(OptimizationReadSetError, match="nested directory"):
        construct_optimization_read_sets("design", tuple(reads.iterdir()))


def test_construct_optimization_read_sets_rejects_no_accepted_fastq_files(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    (reads / "README.txt").touch()

    with pytest.raises(OptimizationReadSetError, match="no accepted FASTQ files"):
        construct_optimization_read_sets("design", tuple(reads.iterdir()))


def test_construct_optimization_read_sets_ignores_unrelated_regular_files(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    read = reads / "sample.fastq"
    read.touch()
    (reads / "README.txt").touch()

    assert construct_optimization_read_sets("design", tuple(reads.iterdir())) == (
        ResolvedReadSet("design__sample", "design", "sample", (read,)),
    )


def test_construct_optimization_read_sets_preserves_dots_and_dashes_in_ids(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    read = reads / "sample.2026-08_R1.fq"
    read.touch()
    mate = reads / "sample.2026-08_R2.fq"
    mate.touch()

    assert construct_optimization_read_sets("design", tuple(reads.iterdir())) == (
        ResolvedReadSet(
            "design__sample.2026-08",
            "design",
            "sample.2026-08",
            (read, mate),
        ),
    )


def test_main_writes_exact_optimization_read_headers(tmp_path: Path) -> None:
    reads = tmp_path / "optimization_reads"
    reads.mkdir()
    (reads / "sample.fastq").touch()
    manifest = tmp_path / "optimization_reads.tsv"

    main(
        (
            "--design-id",
            "design",
            "--directory",
            str(reads),
            "--manifest-out",
            str(manifest),
        ),
    )

    assert manifest.read_text() == 'id\tdesign_id\tmetagenome_id\tread_1\tread_2\ndesign__sample\tdesign\tsample\tsample.fastq\t""\n'


def test_main_rejects_a_regular_file_as_an_optimization_read_directory(
    tmp_path: Path,
) -> None:
    read_file = tmp_path / "sample.fastq"
    read_file.touch()

    with pytest.raises(OptimizationReadSetError, match="is not a directory"):
        main(
            (
                "--design-id",
                "design",
                "--directory",
                str(read_file),
                "--manifest-out",
                str(tmp_path / "optimization_reads.tsv"),
            ),
        )
