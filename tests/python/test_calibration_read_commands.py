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
) -> None:
    read_sources = construct_calibration_read_sources(
        "design",
        ("alpha_R1.fastq", "alpha_R2.fastq.gz", "beta.fq", "README.txt"),
    )

    assert read_sources == (
        ResolvedReadSource("design__alpha_R1", "design", "alpha_R1", "alpha_R1.fastq"),
        ResolvedReadSource("design__alpha_R2", "design", "alpha_R2", "alpha_R2.fastq.gz"),
        ResolvedReadSource("design__beta", "design", "beta", "beta.fq"),
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
) -> None:
    assert construct_calibration_read_sources(
        "design",
        ("sample.fastq", "sample_R1.fastq", "sample_R2.fastq"),
    ) == (
        ResolvedReadSource(
            "design__sample",
            "design",
            "sample",
            "sample.fastq",
        ),
        ResolvedReadSource(
            "design__sample_R1",
            "design",
            "sample_R1",
            "sample_R1.fastq",
        ),
        ResolvedReadSource(
            "design__sample_R2",
            "design",
            "sample_R2",
            "sample_R2.fastq",
        ),
    )


def test_construct_calibration_read_sources_disambiguates_equal_stems(
) -> None:
    assert construct_calibration_read_sources(
        "design",
        ("sample.fastq", "sample.fastq.gz"),
    ) == (
        ResolvedReadSource("design__sample.fastq", "design", "sample.fastq", "sample.fastq"),
        ResolvedReadSource(
            "design__sample.fastq.gz",
            "design",
            "sample.fastq.gz",
            "sample.fastq.gz",
        ),
    )


def test_construct_calibration_read_sources_keeps_final_keys_unique() -> None:
    sources = construct_calibration_read_sources(
        "design",
        ("sample.fastq", "sample.fastq.gz", "sample.fastq.fastq"),
    )

    source_ids = [source.read_source_id for source in sources]
    assert len(set(source_ids)) == 3
    assert all(re.fullmatch(r"design__[A-Za-z0-9][A-Za-z0-9._-]*", value) for value in source_ids)


def test_construct_calibration_read_sources_makes_unsafe_filenames_safe(
) -> None:
    filenames = ("-leading.fastq", "sample one.fastq", "sample;touch PWN.fastq")
    sources = construct_calibration_read_sources("design", filenames)

    assert {source.read_name for source in sources} == set(filenames)
    assert all(
        re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", source.metagenome_id)
        for source in sources
    )


def test_construct_calibration_read_sources_rejects_no_accepted_fastq_files(
) -> None:
    with pytest.raises(CalibrationReadError, match="no accepted FASTQ files"):
        construct_calibration_read_sources("design", ("README.txt",))


def test_construct_calibration_read_sources_ignores_unrelated_regular_files(
) -> None:
    assert construct_calibration_read_sources(
        "design",
        ("sample.fastq", "README.txt"),
    ) == (
        ResolvedReadSource("design__sample", "design", "sample", "sample.fastq"),
    )


def test_construct_calibration_read_sources_preserves_dots_and_dashes_in_ids(
) -> None:
    assert construct_calibration_read_sources("design", ("sample.2026-08_R1.fq",)) == (
        ResolvedReadSource(
            "design__sample.2026-08_R1",
            "design",
            "sample.2026-08_R1",
            "sample.2026-08_R1.fq",
        ),
    )


def test_main_writes_exact_calibration_read_headers(tmp_path: Path) -> None:
    names = tmp_path / "calibration_read_names.txt"
    names.write_text("sample.fastq\n")
    manifest = tmp_path / "calibration_reads.tsv"

    main(
        (
            "--design-id",
            "design",
            "--names",
            str(names),
            "--manifest-out",
            str(manifest),
        ),
    )

    assert manifest.read_text() == "id\tdesign_id\tmetagenome_id\tread\ndesign__sample\tdesign\tsample\tsample.fastq\n"


def test_main_resolves_calibration_read_names_without_accessing_fastq_files(
    tmp_path: Path,
) -> None:
    names = tmp_path / "calibration_read_names.txt"
    names.write_text("sample_R2.fastq.gz\nsample_R1.fastq\n")
    manifest = tmp_path / "calibration_reads.tsv"

    main(
        (
            "--design-id",
            "design",
            "--names",
            str(names),
            "--manifest-out",
            str(manifest),
        ),
    )

    assert manifest.read_text() == (
        "id\tdesign_id\tmetagenome_id\tread\n"
        "design__sample_R1\tdesign\tsample_R1\tsample_R1.fastq\n"
        "design__sample_R2\tdesign\tsample_R2\tsample_R2.fastq.gz\n"
    )
