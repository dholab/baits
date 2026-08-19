from pathlib import Path

import pytest
from calibrate_threshold import (
    CalibrationStatus,
    DeaconThresholdCalibrationError,
    PreparationSummary,
    ThresholdCalibration,
    construct_classified_reads,
    construct_preparation_summary,
    construct_threshold_calibration,
)
from calibrate_threshold import main as calibrate_threshold_main
from classify_read_hits import (
    CandidateReadClassificationError,
    construct_candidate_read_classifications,
    construct_candidate_read_counts,
    construct_read_hits,
    load_calibration_target_taxids,
)
from classify_read_hits import main as classify_read_hits_main
from prepare_read_blast_queries import BlastQueryPreparationError
from prepare_read_blast_queries import main as prepare_read_blast_queries_main

CANDIDATE_HEADER = (
    "metagenome_id\tread_id\tread_length\tbait_count\trepresentative_id\n"
)
CLASSIFIED_HEADER = (
    "metagenome_id\tread_id\tread_length\tbait_count\trepresentative_id\t"
    "classification\tbest_target_bit_score\tbest_non_target_bit_score\n"
)
BLAST_HEADER = (
    "qseqid\tqlen\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\t"
    "sstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n"
)


def write_calibration_target_taxids(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "calibration_target_taxids.tsv"
    path.write_text("taxid\n" + rows)
    return path


def test_load_calibration_target_taxids_uses_a_supplied_scope(tmp_path: Path) -> None:
    scope = write_calibration_target_taxids(tmp_path, "44417\n88456\n")

    assert load_calibration_target_taxids(scope, target_taxid="88456") == frozenset(
        {"44417", "88456"},
    )


def test_load_calibration_target_taxids_defaults_to_the_bait_set_target() -> None:
    assert load_calibration_target_taxids(None, target_taxid="88456") == frozenset(
        {"88456"},
    )


@pytest.mark.parametrize(
    ("content", "error"),
    [
        ("target_taxid\n88456\n", "exactly one header"),
        ("taxid\n88456\n88456\n", "unique"),
        ("taxid\n088456\n", "canonical positive taxid"),
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

    with pytest.raises(CandidateReadClassificationError, match=error):
        load_calibration_target_taxids(scope, target_taxid="88456")


def write_preparation_inputs(
    tmp_path: Path,
    *,
    rows: str,
    candidate_count: int,
    unique_fasta: str,
) -> tuple[Path, Path, Path]:
    counts = tmp_path / "sample.counted_reads.tsv"
    counts.write_text(CANDIDATE_HEADER + rows)
    status = tmp_path / "sample.candidate_read_status.tsv"
    status.write_text(
        "metric\tvalue\n"
        "metagenome_id\tsample\n"
        f"deacon_returned_read_count\t{candidate_count}\n"
        f"candidate_read_count\t{candidate_count}\n",
    )
    fasta = tmp_path / "unique.fasta"
    fasta.write_text(unique_fasta)
    return counts, status, fasta


def preparation_arguments(
    tmp_path: Path,
    counts: Path,
    status: Path,
    fasta: Path,
) -> list[str]:
    return [
        "--design-id",
        "design",
        "--counts",
        str(counts),
        "--statuses",
        str(status),
        "--unique-fasta",
        str(fasta),
        "--candidate-counts-out",
        str(tmp_path / "candidate_read_counts.tsv"),
        "--query-out",
        str(tmp_path / "read_queries.fasta"),
        "--summary-out",
        str(tmp_path / "preparation_summary.tsv"),
        "--terminal-blast-hits-out",
        str(tmp_path / "read_blast_hits.tsv"),
        "--terminal-search-parameters-out",
        str(tmp_path / "read_blast_search_parameters.tsv"),
        "--terminal-classified-reads-out",
        str(tmp_path / "classified_reads.tsv"),
        "--terminal-read-counts-out",
        str(tmp_path / "threshold_read_counts.tsv"),
        "--terminal-summary-out",
        str(tmp_path / "threshold_summary.tsv"),
    ]


def test_prepare_read_blast_queries_preserves_rows_and_uses_unique_queries(
    tmp_path: Path,
) -> None:
    representative = (
        "sequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219"
    )
    counts, status, fasta = write_preparation_inputs(
        tmp_path,
        rows=(
            f"sample\tread_000000000001\t5\t1\t{representative}\n"
            f"sample\tread_000000000002\t5\t1\t{representative}\n"
        ),
        candidate_count=2,
        unique_fasta=f">{representative}\nACGTA\n",
    )

    prepare_read_blast_queries_main(preparation_arguments(tmp_path, counts, status, fasta))

    assert (tmp_path / "candidate_read_counts.tsv").read_text().splitlines() == [
        CANDIDATE_HEADER.rstrip(),
        f"sample\tread_000000000001\t5\t1\t{representative}",
        f"sample\tread_000000000002\t5\t1\t{representative}",
    ]
    assert (tmp_path / "read_queries.fasta").read_text() == f">{representative}\nACGTA\n"
    summary = (tmp_path / "preparation_summary.tsv").read_text()
    assert "candidate_read_count\t2" in summary
    assert "duplicate_sequence_count\t1" in summary
    assert "read_blast_query_count\t1" in summary
    assert not (tmp_path / "threshold_summary.tsv").exists()


def test_prepare_read_blast_queries_emits_complete_terminal_outputs(
    tmp_path: Path,
) -> None:
    counts, status, fasta = write_preparation_inputs(
        tmp_path,
        rows="",
        candidate_count=0,
        unique_fasta="",
    )

    prepare_read_blast_queries_main(preparation_arguments(tmp_path, counts, status, fasta))

    assert not (tmp_path / "read_queries.fasta").exists()
    assert (tmp_path / "read_blast_hits.tsv").read_text() == BLAST_HEADER
    assert (tmp_path / "read_blast_search_parameters.tsv").read_text() == "parameter\tvalue\n"
    assert (tmp_path / "classified_reads.tsv").read_text() == CLASSIFIED_HEADER
    assert (tmp_path / "threshold_read_counts.tsv").read_text().endswith("1\t0\t0\t0\t0\n")
    assert "calibration_status\tNO_CANDIDATE_READS" in (
        tmp_path / "threshold_summary.tsv"
    ).read_text()


def test_prepare_read_blast_queries_rejects_a_query_identifier_disagreement(
    tmp_path: Path,
) -> None:
    counts, status, fasta = write_preparation_inputs(
        tmp_path,
        rows=(
            "sample\tread_000000000001\t5\t1\t"
            "sequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219\n"
        ),
        candidate_count=1,
        unique_fasta=(
            ">sequence_17b80bc751e1f35c75d6ada07267a5fd9981b7b3655cffcaf94a636e93eb27ba\n"
            "CCCCC\n"
        ),
    )

    with pytest.raises(BlastQueryPreparationError, match="queries disagree"):
        prepare_read_blast_queries_main(preparation_arguments(tmp_path, counts, status, fasta))


def test_prepare_read_blast_queries_matches_counts_and_statuses_by_source(
    tmp_path: Path,
) -> None:
    representative = (
        "sequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219"
    )
    count_paths = []
    status_paths = []
    for source in ("sample", "sample.cat"):
        counts = tmp_path / f"design__{source}.counted_reads.tsv"
        counts.write_text(
            CANDIDATE_HEADER
            + f"{source}\tread_000000000001\t5\t1\t{representative}\n",
        )
        status = tmp_path / f"design__{source}.candidate_read_status.tsv"
        status.write_text(
            "metric\tvalue\n"
            f"metagenome_id\t{source}\n"
            "deacon_returned_read_count\t1\n"
            "candidate_read_count\t1\n",
        )
        count_paths.append(counts)
        status_paths.append(status)
    unique = tmp_path / "unique.fasta"
    unique.write_text(f">{representative}\nACGTA\n")
    arguments = preparation_arguments(tmp_path, count_paths[0], status_paths[0], unique)
    counts_index = arguments.index("--counts") + 1
    arguments[counts_index:counts_index + 1] = [str(path) for path in count_paths]
    statuses_index = arguments.index("--statuses") + 1
    arguments[statuses_index:statuses_index + 1] = [
        str(path) for path in reversed(status_paths)
    ]

    prepare_read_blast_queries_main(arguments)

    assert (tmp_path / "candidate_read_counts.tsv").read_text().count("read_000000000001") == 2


def write_classification_candidates(tmp_path: Path, rows: str) -> Path:
    candidates = tmp_path / "candidate_read_counts.tsv"
    candidates.write_text(CANDIDATE_HEADER + rows)
    return candidates


def write_classification_hits(tmp_path: Path, rows: str) -> Path:
    hits = tmp_path / "read_blast_hits.tsv"
    hits.write_text(BLAST_HEADER + rows)
    return hits


def blast_row(query: str, taxids: str, bitscore: str, *, qlen: str = "100") -> str:
    return (
        f"{query}\t{qlen}\taccession\t{taxids}\t100\t100\t0\t0\t1\t100\t1\t100\t0\t"
        f"{bitscore}\t100\ttitle\n"
    )


def classification_rows(
    tmp_path: Path,
    candidate_rows: str,
    hit_rows: str,
    *,
    calibration_target_taxids: frozenset[str] | None = None,
) -> list[dict[str, object]]:
    candidates = construct_candidate_read_counts(
        write_classification_candidates(tmp_path, candidate_rows),
    )
    hits = construct_read_hits(
        write_classification_hits(tmp_path, hit_rows),
        candidate_representatives=candidates.select("representative_id").unique(),
        target_taxid="88456",
        calibration_target_taxids=calibration_target_taxids,
    )
    return construct_candidate_read_classifications(candidates, hits).collect().to_dicts()


def test_classification_projects_one_query_result_to_every_candidate_read(
    tmp_path: Path,
) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tread_000000000001\t100\t2\tquery_a\n"
        "sample\tread_000000000002\t100\t2\tquery_a\n"
        "sample\tread_000000000003\t100\t1\tquery_b\n",
        blast_row("query_a", "88456", "100") + blast_row("query_b", "123", "90"),
    )

    assert [(row["read_id"], row["classification"]) for row in rows] == [
        ("read_000000000001", "TARGET"),
        ("read_000000000002", "TARGET"),
        ("read_000000000003", "NON_TARGET"),
    ]


def test_classification_preserves_target_tie_non_target_and_no_hit_order(
    tmp_path: Path,
) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tread_000000000001\t100\t4\tq1\n"
        "sample\tread_000000000002\t100\t3\tq2\n"
        "sample\tread_000000000003\t100\t2\tq3\n"
        "sample\tread_000000000004\t100\t1\tq4\n",
        blast_row("q1", "88456", "100")
        + blast_row("q2", "123", "100")
        + blast_row("q3", "88456;123", "100"),
    )

    assert [row["classification"] for row in rows] == [
        "TARGET",
        "NON_TARGET",
        "TIED",
        "NO_HIT",
    ]


def test_classification_uses_the_supplied_target_scope(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tread_000000000001\t100\t1\tq1\n",
        blast_row("q1", "44417", "100"),
        calibration_target_taxids=frozenset({"44417", "88456"}),
    )

    assert rows[0]["classification"] == "TARGET"


def test_classification_keeps_the_26th_co_best_subject_decisive(
    tmp_path: Path,
) -> None:
    target_hits = "".join(
        blast_row("q1", "88456", "100").replace("accession", f"target_{index:02d}")
        for index in range(1, 26)
    )
    rows = classification_rows(
        tmp_path,
        "sample\tread_000000000001\t100\t1\tq1\n",
        target_hits
        + blast_row("q1", "123", "100").replace("accession", "non_target_26"),
    )

    assert rows[0]["classification"] == "TIED"
    assert rows[0]["best_target_bit_score"] == "100"
    assert rows[0]["best_non_target_bit_score"] == "100"


def test_classification_uses_the_decimal_tie_boundary(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tread_000000000001\t100\t1\tq1\n"
        "sample\tread_000000000002\t100\t1\tq2\n",
        blast_row("q1", "88456", "100.1")
        + blast_row("q1", "123", "100")
        + blast_row("q2", "88456", "100.1001")
        + blast_row("q2", "123", "100"),
    )

    assert [row["classification"] for row in rows] == ["TIED", "TARGET"]


def test_classification_treats_unknown_taxonomy_as_non_target(
    tmp_path: Path,
) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tread_000000000001\t100\t1\tq1\n",
        blast_row("q1", "N/A", "100"),
    )

    assert rows[0]["classification"] == "NON_TARGET"


def test_candidate_counts_reject_duplicate_read_identity_but_share_queries(
    tmp_path: Path,
) -> None:
    repeated_query = construct_candidate_read_counts(
        write_classification_candidates(
            tmp_path,
            "sample\tread_000000000001\t100\t1\tq1\n"
            "sample\tread_000000000002\t100\t1\tq1\n",
        ),
    )
    assert repeated_query.collect().height == 2

    with pytest.raises(CandidateReadClassificationError, match="identity is duplicated"):
        construct_candidate_read_counts(
            write_classification_candidates(
                tmp_path,
                "sample\tread_000000000001\t100\t1\tq1\n"
                "sample\tread_000000000001\t100\t1\tq2\n",
            ),
        )


@pytest.mark.parametrize("taxids", ["88456;", "0", "-1", "088456"])
def test_read_hits_reject_malformed_taxids(tmp_path: Path, taxids: str) -> None:
    candidates = construct_candidate_read_counts(
        write_classification_candidates(
            tmp_path,
            "sample\tread_000000000001\t100\t1\tq1\n",
        ),
    )

    with pytest.raises(CandidateReadClassificationError, match="malformed staxids"):
        construct_read_hits(
            write_classification_hits(tmp_path, blast_row("q1", taxids, "100")),
            candidate_representatives=candidates.select("representative_id"),
            target_taxid="88456",
        )


def test_classify_read_hits_main_writes_exact_header_and_blank_scores(
    tmp_path: Path,
) -> None:
    candidates = write_classification_candidates(
        tmp_path,
        "sample\tread_000000000001\t100\t1\tq1\n",
    )
    hits = write_classification_hits(tmp_path, "")
    output = tmp_path / "classified_reads.tsv"

    classify_read_hits_main(
        [
            "--candidate-read-counts",
            str(candidates),
            "--blast-hits",
            str(hits),
            "--target-taxid",
            "88456",
            "--output",
            str(output),
        ],
    )

    assert output.read_text() == (
        CLASSIFIED_HEADER
        + "sample\tread_000000000001\t100\t1\tq1\tNO_HIT\t\t\n"
    )


def write_classified_reads(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "classified_reads.tsv"
    path.write_text(CLASSIFIED_HEADER + rows)
    return path


def write_preparation_summary(tmp_path: Path, candidate_count: int) -> Path:
    path = tmp_path / "preparation_summary.tsv"
    path.write_text(
        "metric\tvalue\n"
        "design_id\tdesign\n"
        f"deacon_returned_read_count\t{candidate_count}\n"
        f"candidate_read_count\t{candidate_count}\n"
        "duplicate_sequence_count\t1\n"
        "read_blast_query_count\t4\n",
    )
    return path


def calibration(tmp_path: Path, rows: str) -> ThresholdCalibration:
    classified = construct_classified_reads(write_classified_reads(tmp_path, rows))
    candidate_count = len(rows.strip().splitlines())
    preparation = construct_preparation_summary(
        write_preparation_summary(tmp_path, candidate_count),
        "design",
    )
    return construct_threshold_calibration(classified, preparation)


def classified_row(
    read_id: str,
    bait_count: int,
    classification: str,
) -> str:
    target = "100" if classification in {"TARGET", "TIED"} else ""
    non_target = "100" if classification in {"NON_TARGET", "TIED"} else ""
    return (
        f"sample\t{read_id}\t100\t{bait_count}\tquery_{read_id}\t{classification}\t"
        f"{target}\t{non_target}\n"
    )


def test_threshold_calibration_builds_the_complete_cumulative_curve(
    tmp_path: Path,
) -> None:
    rows = (
        classified_row("bad_1", 1, "NON_TARGET")
        + classified_row("bad_2", 1, "NON_TARGET")
        + classified_row("target_1", 2, "TARGET")
        + classified_row("target_2", 2, "TARGET")
        + classified_row("target_3", 2, "TARGET")
        + classified_row("tie", 2, "TIED")
        + classified_row("no_hit", 2, "NO_HIT")
        + classified_row("target_high", 3, "TARGET")
    )

    result = calibration(tmp_path, rows)

    assert result.read_counts.collect().to_dicts() == [
        {
            "threshold": 1,
            "target_read_count": 4,
            "non_target_read_count": 2,
            "tied_read_count": 1,
            "no_hit_read_count": 1,
        },
        {
            "threshold": 2,
            "target_read_count": 4,
            "non_target_read_count": 0,
            "tied_read_count": 1,
            "no_hit_read_count": 1,
        },
        {
            "threshold": 3,
            "target_read_count": 1,
            "non_target_read_count": 0,
            "tied_read_count": 0,
            "no_hit_read_count": 0,
        },
        {
            "threshold": 4,
            "target_read_count": 0,
            "non_target_read_count": 0,
            "tied_read_count": 0,
            "no_hit_read_count": 0,
        },
    ]
    assert result.conclusion.status is CalibrationStatus.RECOMMENDED_THRESHOLD
    assert result.conclusion.recommended_threshold == 3


def test_threshold_calibration_reports_no_supported_threshold(
    tmp_path: Path,
) -> None:
    result = calibration(
        tmp_path,
        classified_row("target", 2, "TARGET")
        + classified_row("bad", 3, "NON_TARGET"),
    )

    assert result.conclusion.status is CalibrationStatus.NO_SUPPORTED_THRESHOLD
    assert result.conclusion.recommended_threshold is None


def test_threshold_calibration_is_invariant_to_input_order(tmp_path: Path) -> None:
    rows = [
        classified_row("target", 3, "TARGET"),
        classified_row("bad", 2, "NON_TARGET"),
        classified_row("tie", 1, "TIED"),
    ]
    forward = calibration(tmp_path, "".join(rows))
    reverse = calibration(tmp_path, "".join(reversed(rows)))

    assert forward.read_counts.collect().to_dicts() == reverse.read_counts.collect().to_dicts()
    assert forward.conclusion == reverse.conclusion


def test_threshold_calibration_rejects_candidate_count_disagreement(
    tmp_path: Path,
) -> None:
    classified = construct_classified_reads(
        write_classified_reads(tmp_path, classified_row("target", 1, "TARGET")),
    )
    preparation = PreparationSummary("design", 2, 2, 1, 1)

    with pytest.raises(DeaconThresholdCalibrationError, match="count disagrees"):
        construct_threshold_calibration(classified, preparation)


def test_calibrate_threshold_main_writes_exact_outputs(tmp_path: Path) -> None:
    classified = write_classified_reads(
        tmp_path,
        classified_row("target", 1, "TARGET"),
    )
    preparation = write_preparation_summary(tmp_path, 1)
    counts = tmp_path / "threshold_read_counts.tsv"
    summary = tmp_path / "threshold_summary.tsv"

    calibrate_threshold_main(
        [
            "--design-id",
            "design",
            "--classified-reads",
            str(classified),
            "--preparation-summary",
            str(preparation),
            "--read-counts-out",
            str(counts),
            "--summary-out",
            str(summary),
        ],
    )

    assert counts.read_text().splitlines() == [
        "threshold\ttarget_read_count\tnon_target_read_count\ttied_read_count\tno_hit_read_count",
        "1\t1\t0\t0\t0",
        "2\t0\t0\t0\t0",
    ]
    assert "recommended_threshold\t1" in summary.read_text()
