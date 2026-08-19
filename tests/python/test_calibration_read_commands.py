import gzip
import re
from pathlib import Path

import pytest
from resolve_calibration_reads import (
    CalibrationReadError,
    ResolvedReadSource,
    construct_calibration_read_manifest,
    construct_calibration_read_sources,
    main,
)


def test_construct_calibration_read_sources_treats_every_fastq_as_an_independent_source(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    alpha_read_1 = reads / "alpha_R1.fastq"
    alpha_read_1.write_text("@fragment/1\nACGTA\n+\nIIIII\n")
    alpha_read_2 = reads / "alpha_R2.fastq.gz"
    alpha_read_2.write_bytes(gzip.compress(b"@fragment/2\nTGCAT\n+\nIIIII\n", mtime=0))
    beta_read = reads / "beta.fq"
    beta_read.write_text("@single\nGATCC\n+\nIIIII\n")
    (reads / "README.txt").write_text("ignored\n")

    read_sources = construct_calibration_read_sources("design", tuple(reads.iterdir()))

    assert read_sources == (
        ResolvedReadSource("design__alpha_R1", "design", "alpha_R1", alpha_read_1),
        ResolvedReadSource("design__alpha_R2", "design", "alpha_R2", alpha_read_2),
        ResolvedReadSource("design__beta", "design", "beta", beta_read),
    )
    assert construct_calibration_read_manifest(read_sources).collect().to_dicts() == [
        {
            "id": "design__alpha_R1",
            "design_id": "design",
            "metagenome_id": "alpha_R1",
            "read": "alpha_R1.fastq",
        },
        {
            "id": "design__alpha_R2",
            "design_id": "design",
            "metagenome_id": "alpha_R2",
            "read": "alpha_R2.fastq.gz",
        },
        {
            "id": "design__beta",
            "design_id": "design",
            "metagenome_id": "beta",
            "read": "beta.fq",
        },
    ]


def test_construct_calibration_read_sources_accepts_mate_like_names_without_pairing(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    (reads / "sample.fastq").touch()
    (reads / "sample_R1.fastq").touch()
    (reads / "sample_R2.fastq").touch()

    assert construct_calibration_read_sources("design", tuple(reads.iterdir())) == (
        ResolvedReadSource(
            "design__sample",
            "design",
            "sample",
            reads / "sample.fastq",
        ),
        ResolvedReadSource(
            "design__sample_R1",
            "design",
            "sample_R1",
            reads / "sample_R1.fastq",
        ),
        ResolvedReadSource(
            "design__sample_R2",
            "design",
            "sample_R2",
            reads / "sample_R2.fastq",
        ),
    )


def test_construct_calibration_read_sources_disambiguates_equal_stems(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    plain = reads / "sample.fastq"
    compressed = reads / "sample.fastq.gz"
    plain.touch()
    compressed.touch()

    assert construct_calibration_read_sources("design", tuple(reads.iterdir())) == (
        ResolvedReadSource("design__sample.fastq", "design", "sample.fastq", plain),
        ResolvedReadSource(
            "design__sample.fastq.gz",
            "design",
            "sample.fastq.gz",
            compressed,
        ),
    )


def test_construct_calibration_read_sources_keeps_final_keys_unique(tmp_path: Path) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    for name in ("sample.fastq", "sample.fastq.gz", "sample.fastq.fastq"):
        (reads / name).touch()

    sources = construct_calibration_read_sources("design", tuple(reads.iterdir()))

    source_ids = [source.read_source_id for source in sources]
    assert len(set(source_ids)) == 3
    assert all(re.fullmatch(r"design__[A-Za-z0-9][A-Za-z0-9._-]*", value) for value in source_ids)


def test_construct_calibration_read_sources_makes_unsafe_filenames_safe(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration reads"
    reads.mkdir()
    filenames = ("-leading.fastq", "sample one.fastq", "sample;touch PWN.fastq")
    for name in filenames:
        (reads / name).touch()

    sources = construct_calibration_read_sources("design", tuple(reads.iterdir()))

    assert {source.read.name for source in sources} == set(filenames)
    assert all(
        re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", source.metagenome_id)
        for source in sources
    )


def test_construct_calibration_read_sources_rejects_nested_directories(tmp_path: Path) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    (reads / "nested").mkdir()

    with pytest.raises(CalibrationReadError, match="nested directory"):
        construct_calibration_read_sources("design", tuple(reads.iterdir()))


def test_construct_calibration_read_sources_rejects_no_accepted_fastq_files(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    (reads / "README.txt").touch()

    with pytest.raises(CalibrationReadError, match="no accepted FASTQ files"):
        construct_calibration_read_sources("design", tuple(reads.iterdir()))


def test_construct_calibration_read_sources_ignores_unrelated_regular_files(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    read = reads / "sample.fastq"
    read.touch()
    (reads / "README.txt").touch()

    assert construct_calibration_read_sources("design", tuple(reads.iterdir())) == (
        ResolvedReadSource("design__sample", "design", "sample", read),
    )


def test_construct_calibration_read_sources_preserves_dots_and_dashes_in_ids(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    read = reads / "sample.2026-08_R1.fq"
    read.touch()
    assert construct_calibration_read_sources("design", tuple(reads.iterdir())) == (
        ResolvedReadSource(
            "design__sample.2026-08_R1",
            "design",
            "sample.2026-08_R1",
            read,
        ),
    )


def test_main_writes_exact_calibration_read_headers(tmp_path: Path) -> None:
    reads = tmp_path / "calibration_reads"
    reads.mkdir()
    (reads / "sample.fastq").touch()
    manifest = tmp_path / "calibration_reads.tsv"

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

    assert manifest.read_text() == "id\tdesign_id\tmetagenome_id\tread\ndesign__sample\tdesign\tsample\tsample.fastq\n"


def test_main_rejects_a_regular_file_as_a_calibration_read_directory(
    tmp_path: Path,
) -> None:
    read_file = tmp_path / "sample.fastq"
    read_file.touch()

    with pytest.raises(CalibrationReadError, match="not a directory"):
        main(
            (
                "--design-id",
                "design",
                "--directory",
                str(read_file),
                "--manifest-out",
                str(tmp_path / "calibration_reads.tsv"),
            ),
        )
