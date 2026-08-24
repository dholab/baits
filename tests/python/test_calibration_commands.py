import hashlib
from decimal import Decimal
from pathlib import Path

import pytest
from calibrate_threshold import (
    CalibrationRead,
    ThresholdCalibrationError,
    category_counts,
    read_calibration_reads,
    recommended_threshold,
)
from calibrate_threshold import main as calibrate_threshold_main
from prepare_read_blast_queries import BlastQueryPreparationError
from prepare_read_blast_queries import main as prepare_read_blast_queries_main
from score_calibration_reads import (
    BLAST_FIELDS,
    CalibrationReadScoreError,
    construct_read_bait_counts,
    construct_read_hits,
    construct_scored_reads,
    load_calibration_target_taxids,
)
from score_calibration_reads import main as score_calibration_reads_main

QUERY_A = f"sequence_{hashlib.sha256(b'AAAAA').hexdigest()}"
QUERY_B = f"sequence_{hashlib.sha256(b'CCCCC').hexdigest()}"
READ_HEADER = "metagenome_id\tread_id\tread_length\tbait_count\tquery_id\n"
SCORED_HEADER = (
    "metagenome_id\tread_id\tread_length\tbait_count\tquery_id\t"
    "best_target_bit_score\tbest_other_bit_score\n"
)
BLAST_HEADER = "\t".join(BLAST_FIELDS) + "\n"


def write_calibration_target_taxids(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "calibration_target_taxids.tsv"
    path.write_text("taxid\n" + rows)
    return path


def test_load_calibration_target_taxids_uses_a_supplied_scope(tmp_path: Path) -> None:
    scope = write_calibration_target_taxids(tmp_path, "44417\n88456\n")
    assert load_calibration_target_taxids(scope, "88456") == frozenset(
        {"44417", "88456"},
    )
    assert load_calibration_target_taxids(None, "88456") == frozenset({"88456"})


@pytest.mark.parametrize(
    ("content", "error"),
    [
        ("target_taxid\n88456\n", "exactly one header"),
        ("taxid\n88456\n88456\n", "unique"),
        ("taxid\n088456\n", "canonical positive taxids"),
        ("taxid\n44417\n", "must include target_taxid 88456"),
    ],
)
def test_load_calibration_target_taxids_rejects_invalid_scopes(
    tmp_path: Path,
    content: str,
    error: str,
) -> None:
    scope = tmp_path / "calibration_target_taxids.tsv"
    scope.write_text(content)
    with pytest.raises(CalibrationReadScoreError, match=error):
        load_calibration_target_taxids(scope, "88456")


def write_count_table(root: Path, name: str, rows: str) -> Path:
    path = root / f"{name}.read_bait_counts.tsv"
    path.write_text(READ_HEADER + rows)
    return path


def preparation_arguments(tmp_path: Path, unique_fasta: Path) -> list[str]:
    return [
        "--counts-dir",
        str(tmp_path / "counts"),
        "--unique-fasta",
        str(unique_fasta),
        "--counts-out",
        str(tmp_path / "read_bait_counts.tsv"),
        "--queries-out",
        str(tmp_path / "read_queries.fasta"),
    ]


def test_prepare_read_blast_queries_preserves_observations_and_unique_queries(
    tmp_path: Path,
) -> None:
    counts = tmp_path / "counts"
    counts.mkdir()
    write_count_table(
        counts,
        "design__sample_a",
        f"sample_a\tread_000000000001\t5\t1\t{QUERY_A}\n"
        f"sample_a\tread_000000000002\t5\t1\t{QUERY_A}\n",
    )
    write_count_table(
        counts,
        "design__sample_b",
        f"sample_b\tread_000000000001\t5\t1\t{QUERY_B}\n",
    )
    unique = tmp_path / "unique.fasta"
    unique.write_text(f">{QUERY_A}\nAAAAA\n>{QUERY_B}\nCCCCC\n")

    prepare_read_blast_queries_main(preparation_arguments(tmp_path, unique))

    assert (tmp_path / "read_bait_counts.tsv").read_text().splitlines() == [
        READ_HEADER.rstrip(),
        f"sample_a\tread_000000000001\t5\t1\t{QUERY_A}",
        f"sample_a\tread_000000000002\t5\t1\t{QUERY_A}",
        f"sample_b\tread_000000000001\t5\t1\t{QUERY_B}",
    ]
    assert (tmp_path / "read_queries.fasta").read_text() == (
        f">{QUERY_A}\nAAAAA\n>{QUERY_B}\nCCCCC\n"
    )


def test_prepare_read_blast_queries_emits_only_empty_primary_outputs(
    tmp_path: Path,
) -> None:
    counts = tmp_path / "counts"
    counts.mkdir()
    write_count_table(counts, "design__sample", "")
    unique = tmp_path / "unique.fasta"
    unique.write_text("")

    prepare_read_blast_queries_main(preparation_arguments(tmp_path, unique))

    assert (tmp_path / "read_bait_counts.tsv").read_text() == READ_HEADER
    assert (tmp_path / "read_queries.fasta").read_text() == ""
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "counts",
        "read_bait_counts.tsv",
        "read_queries.fasta",
        "unique.fasta",
    ]


def test_prepare_read_blast_queries_rejects_query_disagreement(tmp_path: Path) -> None:
    counts = tmp_path / "counts"
    counts.mkdir()
    write_count_table(
        counts,
        "design__sample",
        f"sample\tread_000000000001\t5\t1\t{QUERY_A}\n",
    )
    unique = tmp_path / "unique.fasta"
    unique.write_text(f">{QUERY_B}\nCCCCC\n")

    with pytest.raises(BlastQueryPreparationError, match="queries disagree"):
        prepare_read_blast_queries_main(preparation_arguments(tmp_path, unique))


def write_read_bait_counts(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "read_bait_counts.tsv"
    path.write_text(READ_HEADER + rows)
    return path


def write_blast_matches(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "read_blast_hits.tsv"
    path.write_text(BLAST_HEADER + rows)
    return path


def blast_row(query: str, taxids: str, bitscore: str) -> str:
    return (
        f"{query}\t100\taccession\t{taxids}\t100\t100\t0\t0\t1\t100\t1\t100\t0\t"
        f"{bitscore}\t100\ttitle\n"
    )


def scored_rows(
    tmp_path: Path,
    read_rows: str,
    match_rows: str,
    *,
    target_scope: frozenset[str] | None = None,
) -> list[dict[str, object]]:
    reads = construct_read_bait_counts(write_read_bait_counts(tmp_path, read_rows))
    hits = construct_read_hits(
        write_blast_matches(tmp_path, match_rows),
        query_ids=reads.select("query_id").unique(),
        target_taxid="88456",
        calibration_target_taxids=target_scope,
    )
    return construct_scored_reads(reads, hits).collect().to_dicts()


def test_scoring_projects_each_shared_query_back_to_every_read(tmp_path: Path) -> None:
    rows = scored_rows(
        tmp_path,
        f"sample\tread_000000000001\t100\t2\t{QUERY_A}\n"
        f"sample\tread_000000000002\t100\t2\t{QUERY_A}\n"
        f"sample\tread_000000000003\t100\t1\t{QUERY_B}\n",
        blast_row(QUERY_A, "88456", "100") + blast_row(QUERY_B, "123", "90"),
    )
    assert [
        (row["read_id"], row["best_target_bit_score"], row["best_other_bit_score"])
        for row in rows
    ] == [
        ("read_000000000001", "100", ""),
        ("read_000000000002", "100", ""),
        ("read_000000000003", "", "90"),
    ]


def test_scoring_preserves_best_scores_and_blank_no_hit_rows(tmp_path: Path) -> None:
    rows = scored_rows(
        tmp_path,
        f"sample\tread_000000000001\t100\t1\t{QUERY_A}\n"
        f"sample\tread_000000000002\t100\t1\t{QUERY_B}\n",
        blast_row(QUERY_A, "88456", "99")
        + blast_row(QUERY_A, "123", "80")
        + blast_row(QUERY_A, "88456;123", "95")
        + blast_row(QUERY_A, "88456", "100")
        + blast_row(QUERY_A, "123", "90"),
    )
    assert rows[0]["best_target_bit_score"] == "100"
    assert rows[0]["best_other_bit_score"] == "95"
    assert rows[1]["best_target_bit_score"] == ""
    assert rows[1]["best_other_bit_score"] == ""


def test_scoring_uses_calibration_target_scope_and_unknown_taxonomy(tmp_path: Path) -> None:
    rows = scored_rows(
        tmp_path,
        f"sample\tread_000000000001\t100\t1\t{QUERY_A}\n",
        blast_row(QUERY_A, "44417", "100") + blast_row(QUERY_A, "N/A", "90"),
        target_scope=frozenset({"44417", "88456"}),
    )
    assert rows[0]["best_target_bit_score"] == "100"
    assert rows[0]["best_other_bit_score"] == "90"


@pytest.mark.parametrize("taxids", ["88456;", "0", "-1", "088456"])
def test_scoring_rejects_malformed_taxids(tmp_path: Path, taxids: str) -> None:
    reads = construct_read_bait_counts(
        write_read_bait_counts(
            tmp_path,
            f"sample\tread_000000000001\t100\t1\t{QUERY_A}\n",
        ),
    )
    with pytest.raises(CalibrationReadScoreError, match="malformed staxids"):
        construct_read_hits(
            write_blast_matches(tmp_path, blast_row(QUERY_A, taxids, "100")),
            query_ids=reads.select("query_id").unique(),
            target_taxid="88456",
        )


def test_score_calibration_reads_main_writes_no_classification_column(
    tmp_path: Path,
) -> None:
    reads = write_read_bait_counts(
        tmp_path,
        f"sample\tread_000000000001\t100\t1\t{QUERY_A}\n",
    )
    matches = write_blast_matches(tmp_path, "")
    output = tmp_path / "calibration_reads.tsv"

    score_calibration_reads_main(
        [
            "--read-bait-counts",
            str(reads),
            "--blast-matches",
            str(matches),
            "--target-taxid",
            "88456",
            "--output",
            str(output),
        ],
    )

    assert output.read_text() == (
        SCORED_HEADER
        + f"sample\tread_000000000001\t100\t1\t{QUERY_A}\t\t\n"
    )


def calibration_read(
    read_id: str,
    bait_count: int,
    target: str = "",
    other: str = "",
) -> CalibrationRead:
    return CalibrationRead(
        "sample",
        read_id,
        100,
        bait_count,
        f"query_{read_id}",
        None if target == "" else Decimal(target),
        None if other == "" else Decimal(other),
    )


def test_threshold_categories_are_derived_only_from_scores() -> None:
    assert category_counts(calibration_read("target", 1, target="100")) == (1, 0, 0, 0)
    assert category_counts(calibration_read("other", 1, other="100")) == (0, 1, 0, 0)
    assert category_counts(calibration_read("tie", 1, target="100.1", other="100")) == (0, 0, 1, 0)
    assert category_counts(calibration_read("target_above_tie", 1, target="100.1001", other="100")) == (1, 0, 0, 0)
    assert category_counts(calibration_read("no_hit", 1)) == (0, 0, 0, 1)


def test_recommended_threshold_uses_the_smallest_supported_count() -> None:
    reads = (
        calibration_read("bad", 1, other="100"),
        calibration_read("target_low", 2, target="100"),
        calibration_read("target_high", 3, target="100"),
    )
    assert recommended_threshold(reads) == 2
    assert recommended_threshold(tuple(reversed(reads))) == 2
    assert recommended_threshold((calibration_read("bad", 3, other="100"),)) is None


def test_calibrate_threshold_main_writes_only_the_recommended_value(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "calibration_reads.tsv"
    reads.write_text(
        SCORED_HEADER
        + f"sample\tread_000000000001\t100\t1\t{QUERY_A}\t100\t\n",
    )
    output = tmp_path / "recommended_threshold.tsv"

    calibrate_threshold_main(
        [
            "--design-id",
            "design",
            "--calibration-reads",
            str(reads),
            "--output",
            str(output),
        ],
    )

    assert read_calibration_reads(reads)[0].query_id == QUERY_A
    assert output.read_text() == "design_id\trecommended_threshold\ndesign\t1\n"


def test_threshold_rejects_duplicate_read_identity(tmp_path: Path) -> None:
    reads = tmp_path / "calibration_reads.tsv"
    row = f"sample\tread_000000000001\t100\t1\t{QUERY_A}\t100\t\n"
    reads.write_text(SCORED_HEADER + row + row)
    with pytest.raises(ThresholdCalibrationError, match="duplicated"):
        read_calibration_reads(reads)
