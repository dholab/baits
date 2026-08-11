from pathlib import Path

import polars as pl
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from calibrate_deacon_threshold import (
    CalibrationStatus,
    DeaconThresholdCalibrationError,
    ThresholdCalibration,
    construct_classified_reads,
    construct_preparation_summary,
    construct_threshold_calibration,
)
from calibrate_deacon_threshold import main as calibrate_deacon_threshold_main
from classify_read_hits import (
    CandidateReadClassificationError,
    construct_candidate_read_classifications,
    construct_candidate_read_counts,
    construct_whole_read_hits,
)
from classify_read_hits import main as classify_read_hits_main
from count_read_baits import (
    CandidateReadRecountError,
    construct_candidate_read_recount,
    iter_read_fragments,
    load_bait_set,
)
from count_read_baits import main as count_read_baits_main
from prepare_read_blast_queries import (
    BlastQueryPreparationError,
    construct_blast_query_preparation,
    construct_metagenome_candidate_evidence,
)
from prepare_read_blast_queries import main as prepare_read_blast_queries_main


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n" for name, sequence in records))


def read_record(identifier: str, sequence: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id=identifier)


def test_construct_candidate_read_recount_deduplicates_before_excluding_zero_bait_mates() -> None:
    recount = construct_candidate_read_recount(
        metagenome_id="sample",
        baits=frozenset({"ACGTA"}),
        kmer_size=5,
        fragments=(
            (read_record("f1/1", "ACGTA"), read_record("f1/2", "TTTTT")),
            (read_record("f2/1", "ACGTA"), read_record("f2/2", "TTTTT")),
            (read_record("f3/1", "TTTTT"), read_record("f3/2", "ACGTA")),
            (read_record("f4/1", "TTTTT"), read_record("f4/2", "CCCCC")),
        ),
    )

    assert recount.counts.collect().to_dicts() == [
        {
            "metagenome_id": "sample",
            "fragment_id": "f1",
            "mate": "1",
            "read_length": 5,
            "distinct_bait_count": 1,
            "fragment_distinct_bait_count": 1,
            "candidate_sequence_id": "candidate_read_000001",
        },
        {
            "metagenome_id": "sample",
            "fragment_id": "f3",
            "mate": "2",
            "read_length": 5,
            "distinct_bait_count": 1,
            "fragment_distinct_bait_count": 1,
            "candidate_sequence_id": "candidate_read_000002",
        },
    ]
    assert [(record.id, str(record.seq)) for record in recount.fasta_records] == [
        ("candidate_read_000001", "ACGTA"),
        ("candidate_read_000002", "ACGTA"),
    ]
    assert recount.status.deacon_returned_read_count == 8
    assert recount.status.duplicate_fragment_count == 1
    assert recount.status.zero_bait_read_count == 4
    assert recount.status.candidate_read_count == 2


def test_construct_candidate_read_recount_records_the_paired_fragment_bait_union() -> None:
    recount = construct_candidate_read_recount(
        metagenome_id="sample",
        baits=frozenset({"ACGTA", "CCCCC"}),
        kmer_size=5,
        fragments=((read_record("fragment/1", "ACGTA"), read_record("fragment/2", "CCCCC")),),
    )

    assert recount.counts.collect().select(
        "distinct_bait_count",
        "fragment_distinct_bait_count",
    ).rows() == [(1, 2), (1, 2)]


def test_construct_candidate_read_recount_deduplicates_shared_baits_across_mates() -> None:
    recount = construct_candidate_read_recount(
        metagenome_id="sample",
        baits=frozenset({"ACGTA", "CCCCC"}),
        kmer_size=5,
        fragments=((read_record("fragment/1", "ACGTACCCCC"), read_record("fragment/2", "ACGTA")),),
    )

    assert recount.counts.collect().select(
        "distinct_bait_count",
        "fragment_distinct_bait_count",
    ).rows() == [(2, 2), (1, 2)]


def test_construct_candidate_read_recount_counts_repeated_bait_once() -> None:
    recount = construct_candidate_read_recount(
        metagenome_id="sample",
        baits=frozenset({"ACGTA"}),
        kmer_size=5,
        fragments=((read_record("fragment", "ACGTAACGTA"),),),
    )

    assert recount.counts.collect().select("distinct_bait_count").item() == 1


def test_construct_candidate_read_recount_matches_reverse_complements() -> None:
    recount = construct_candidate_read_recount(
        metagenome_id="sample",
        baits=frozenset({"ACGTA"}),
        kmer_size=5,
        fragments=((read_record("fragment", "TACGT"),),),
    )

    assert recount.counts.collect().select("distinct_bait_count").item() == 1
    assert [(record.id, str(record.seq)) for record in recount.fasta_records] == [
        ("candidate_read_000001", "ACGTA"),
    ]


def test_iter_read_fragments_rejects_unequal_paired_file_record_counts(tmp_path: Path) -> None:
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    write_fastq(read_1, [("f1/1", "ACGTA"), ("f2/1", "CCCCC")])
    write_fastq(read_2, [("f1/2", "ACGTA")])

    with pytest.raises(CandidateReadRecountError, match="different record counts"):
        tuple(iter_read_fragments((read_1, read_2)))


def test_construct_candidate_read_recount_rejects_mismatched_paired_fragment_identities() -> None:
    with pytest.raises(CandidateReadRecountError, match="different fragment identifiers"):
        construct_candidate_read_recount(
            metagenome_id="sample",
            baits=frozenset({"ACGTA"}),
            kmer_size=5,
            fragments=((read_record("left/1", "ACGTA"), read_record("right/2", "ACGTA")),),
        )


def test_construct_candidate_read_recount_rejects_repeated_fragment_identities() -> None:
    with pytest.raises(CandidateReadRecountError, match="fragment identifier is duplicated"):
        construct_candidate_read_recount(
            metagenome_id="sample",
            baits=frozenset({"ACGTA"}),
            kmer_size=5,
            fragments=(
                (read_record("fragment/1", "ACGTA"),),
                (read_record("fragment/1", "CCCCC"),),
            ),
        )


def test_construct_candidate_read_recount_excludes_reads_shorter_than_k() -> None:
    recount = construct_candidate_read_recount(
        metagenome_id="sample",
        baits=frozenset({"ACGTA"}),
        kmer_size=5,
        fragments=((read_record("fragment", "ACGT"),),),
    )

    assert recount.counts.collect().is_empty()
    assert recount.status.zero_bait_read_count == 1


@pytest.mark.parametrize(
    ("bait_fasta", "kmer_size", "error"),
    [
        ("", 5, "must not be empty"),
        (">one\nACGTA\n>two\nTACGT\n", 5, "must be unique"),
        (">one\nACGTA\n>two\nACGT\n", 5, "one k-mer size"),
        (">one\nACGTA\n", 4, "disagrees with kmer_size"),
    ],
)
def test_load_bait_set_rejects_malformed_bait_sets(
    tmp_path: Path,
    bait_fasta: str,
    kmer_size: int,
    error: str,
) -> None:
    baits = tmp_path / "baits.fasta"
    baits.write_text(bait_fasta)

    with pytest.raises(CandidateReadRecountError, match=error):
        load_bait_set(baits, kmer_size)


def test_count_read_baits_main_writes_exact_headers_and_readable_fasta(tmp_path: Path) -> None:
    baits = tmp_path / "baits.fasta"
    baits.write_text(">bait\nACGTA\n")
    reads = tmp_path / "sample.fastq"
    write_fastq(reads, [("fragment", "ACGTA")])
    counts = tmp_path / "counts.tsv"
    fasta = tmp_path / "candidate_reads.fasta"
    status = tmp_path / "status.tsv"

    count_read_baits_main([
        "--metagenome-id", "sample",
        "--baits", str(baits),
        "--kmer-size", "5",
        "--reads", str(reads),
        "--counts-out", str(counts),
        "--fasta-out", str(fasta),
        "--status-out", str(status),
    ])

    assert counts.read_text().splitlines()[0] == (
        "metagenome_id\tfragment_id\tmate\tread_length\tdistinct_bait_count\t"
        "fragment_distinct_bait_count\tcandidate_sequence_id"
    )
    assert status.read_text().splitlines()[0] == "metric\tvalue"
    assert [(record.id, str(record.seq)) for record in SeqIO.parse(fasta, "fasta")] == [
        ("candidate_read_000001", "ACGTA"),
    ]


def write_candidate_evidence(
    directory: Path,
    metagenome_id: str,
    rows: list[str],
    fasta: str,
    *,
    candidate_read_count: int | None = None,
) -> tuple[Path, Path, Path]:
    count_header = (
        "metagenome_id\tfragment_id\tmate\tread_length\tdistinct_bait_count\t"
        "fragment_distinct_bait_count\tcandidate_sequence_id\n"
    )
    counts = directory / f"{metagenome_id}.counts.tsv"
    counts.write_text(count_header + "".join(rows))
    fasta_path = directory / f"{metagenome_id}.fasta"
    fasta_path.write_text(fasta)
    status = directory / f"{metagenome_id}.status.tsv"
    status.write_text(
        "metric\tvalue\n"
        f"metagenome_id\t{metagenome_id}\n"
        "deacon_returned_read_count\t4\n"
        "duplicate_fragment_count\t1\n"
        "zero_bait_read_count\t1\n"
        f"candidate_read_count\t{len(rows) if candidate_read_count is None else candidate_read_count}\n",
    )
    return counts, fasta_path, status


def test_construct_blast_query_preparation_deduplicates_reverse_complements(tmp_path: Path) -> None:
    alpha = construct_metagenome_candidate_evidence(*write_candidate_evidence(tmp_path, "alpha", ["alpha\tf2\t1\t5\t2\t2\tread\n"], ">read\nACGTA\n"))
    beta = construct_metagenome_candidate_evidence(*write_candidate_evidence(tmp_path, "beta", ["beta\tf1\t\t5\t1\t1\tread\n"], ">read\nTACGT\n"))
    preparation = construct_blast_query_preparation(design_id="design", evidence=(beta, alpha))

    assert preparation.candidate_read_counts.collect().to_dicts() == [
        {"metagenome_id": "alpha", "fragment_id": "f2", "mate": "1", "read_length": 5, "distinct_bait_count": 2, "fragment_distinct_bait_count": 2, "representative_id": "representative_000001"},
        {"metagenome_id": "beta", "fragment_id": "f1", "mate": "", "read_length": 5, "distinct_bait_count": 1, "fragment_distinct_bait_count": 1, "representative_id": "representative_000001"},
    ]
    assert [(record.id, str(record.seq)) for record in preparation.query_records] == [("representative_000001", "ACGTA")]
    assert preparation.summary.candidate_read_count == 2
    assert preparation.summary.whole_read_blast_query_count == 1


def test_construct_blast_query_preparation_joins_repeated_sequence_ids_by_metagenome(tmp_path: Path) -> None:
    alpha = construct_metagenome_candidate_evidence(*write_candidate_evidence(
        tmp_path,
        "alpha",
        ["alpha\tf\t\t5\t1\t1\tshared\n"],
        ">shared\nAAAAA\n",
    ))
    beta = construct_metagenome_candidate_evidence(*write_candidate_evidence(
        tmp_path,
        "beta",
        ["beta\tf\t\t5\t1\t1\tshared\n"],
        ">shared\nCCCCC\n",
    ))

    preparation = construct_blast_query_preparation(design_id="design", evidence=(alpha, beta))

    assert preparation.candidate_read_counts.collect().select(
        "metagenome_id",
        "representative_id",
    ).to_dicts() == [
        {"metagenome_id": "alpha", "representative_id": "representative_000001"},
        {"metagenome_id": "beta", "representative_id": "representative_000002"},
    ]


def test_construct_blast_query_preparation_preserves_empty_metagenome_status(tmp_path: Path) -> None:
    empty = construct_metagenome_candidate_evidence(*write_candidate_evidence(tmp_path, "empty", [], ""))
    nonempty = construct_metagenome_candidate_evidence(*write_candidate_evidence(tmp_path, "full", ["full\tf\t\t5\t1\t1\tread\n"], ">read\nAAAAA\n"))
    preparation = construct_blast_query_preparation(design_id="design", evidence=(empty, nonempty))

    assert preparation.summary.deacon_returned_read_count == 8
    assert preparation.candidate_read_counts.collect().get_column("metagenome_id").to_list() == ["full"]


def test_construct_blast_query_preparation_emits_complete_terminal_evidence(tmp_path: Path) -> None:
    evidence = construct_metagenome_candidate_evidence(*write_candidate_evidence(tmp_path, "sample", [], ""))
    preparation = construct_blast_query_preparation(design_id="design", evidence=(evidence,))

    assert preparation.query_records == ()
    assert preparation.terminal_evidence is not None
    assert preparation.terminal_evidence.threshold_read_counts.collect().to_dicts() == [{"threshold": 1, "target_read_count": 0, "non_target_read_count": 0, "tied_read_count": 0, "no_hit_read_count": 0}]
    assert preparation.terminal_evidence.summary.collect().filter(pl.col("metric") == "conclusion").item(0, "value") == "The optimization read set contains no candidate reads."


@pytest.mark.parametrize("fasta", [">other\nACGTA\n", ">read\nACGTA\n>extra\nTTTTT\n", ">read\nACGTA\n>read\nACGTA\n"])
def test_construct_metagenome_candidate_evidence_rejects_fasta_identity_disagreements(tmp_path: Path, fasta: str) -> None:
    paths = write_candidate_evidence(tmp_path, "sample", ["sample\tf\t\t5\t1\t1\tread\n"], fasta)
    with pytest.raises(BlastQueryPreparationError, match="FASTA"):
        construct_metagenome_candidate_evidence(*paths)


def test_construct_metagenome_candidate_evidence_rejects_count_status_disagreement(tmp_path: Path) -> None:
    paths = write_candidate_evidence(tmp_path, "sample", ["other\tf\t\t5\t1\t1\tread\n"], ">read\nACGTA\n")
    with pytest.raises(BlastQueryPreparationError, match="different metagenomes"):
        construct_metagenome_candidate_evidence(*paths)


def test_construct_blast_query_preparation_rejects_duplicate_status_metagenomes(tmp_path: Path) -> None:
    paths = write_candidate_evidence(tmp_path, "sample", [], "")
    evidence = construct_metagenome_candidate_evidence(*paths)
    with pytest.raises(BlastQueryPreparationError, match="duplicate"):
        construct_blast_query_preparation(design_id="design", evidence=(evidence, evidence))


def test_construct_metagenome_candidate_evidence_rejects_duplicate_candidate_read_identity(tmp_path: Path) -> None:
    paths = write_candidate_evidence(
        tmp_path,
        "sample",
        [
            "sample\tfragment\t1\t5\t1\t1\tfirst\n",
            "sample\tfragment\t1\t5\t1\t1\tsecond\n",
        ],
        ">first\nACGTA\n>second\nTTTTT\n",
    )

    with pytest.raises(BlastQueryPreparationError, match="Candidate-read identity is duplicated"):
        construct_metagenome_candidate_evidence(*paths)


@pytest.mark.parametrize(
    "row",
    [
        "sample\tf\t\tbad\t1\t1\tread\n",
        "sample\tf\t\t-1\t1\t1\tread\n",
        f"sample\tf\t\t{2**63}\t1\t1\tread\n",
        "sample\t\t\t5\t1\t1\tread\n",
        "sample\tf\t\t5\t1\t1\tread\nsample\tf\t\t5\t1\t1\tread2\n",
    ],
)
def test_construct_metagenome_candidate_evidence_rejects_invalid_counts(tmp_path: Path, row: str) -> None:
    paths = write_candidate_evidence(tmp_path, "sample", [row], ">read\nACGTA\n")
    with pytest.raises(BlastQueryPreparationError):
        construct_metagenome_candidate_evidence(*paths)


@pytest.mark.parametrize(
    "header",
    [
        "renamed\tfragment_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\tcandidate_sequence_id\n",
        "fragment_id\tmetagenome_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\tcandidate_sequence_id\n",
    ],
)
def test_construct_metagenome_candidate_evidence_rejects_invalid_count_headers(
    tmp_path: Path,
    header: str,
) -> None:
    counts, fasta, status = write_candidate_evidence(
        tmp_path,
        "sample",
        ["sample\tf\t\t5\t1\t1\tread\n"],
        ">read\nACGTA\n",
    )
    counts.write_text(header + "sample\tf\t\t5\t1\t1\tread\n")

    with pytest.raises(BlastQueryPreparationError, match="Candidate-read counts are malformed"):
        construct_metagenome_candidate_evidence(counts, fasta, status)


def test_construct_metagenome_candidate_evidence_rejects_renamed_status_header(tmp_path: Path) -> None:
    counts, fasta, status = write_candidate_evidence(
        tmp_path,
        "sample",
        ["sample\tf\t\t5\t1\t1\tread\n"],
        ">read\nACGTA\n",
    )
    status.write_text(status.read_text().replace("metric\tvalue", "renamed\tvalue", 1))

    with pytest.raises(BlastQueryPreparationError, match="Candidate-read status is malformed"):
        construct_metagenome_candidate_evidence(counts, fasta, status)


def test_construct_metagenome_candidate_evidence_rejects_status_count_and_sequence_length(tmp_path: Path) -> None:
    paths = write_candidate_evidence(tmp_path, "sample", ["sample\tf\t\t5\t1\t1\tread\n"], ">read\nACGT\n", candidate_read_count=2)
    with pytest.raises(BlastQueryPreparationError, match="length"):
        construct_metagenome_candidate_evidence(*paths)


def test_construct_metagenome_candidate_evidence_rejects_status_count_disagreement(tmp_path: Path) -> None:
    paths = write_candidate_evidence(tmp_path, "sample", ["sample\tf\t\t5\t1\t1\tread\n"], ">read\nACGTA\n", candidate_read_count=2)
    with pytest.raises(BlastQueryPreparationError, match="status count"):
        construct_metagenome_candidate_evidence(*paths)


def test_prepare_read_blast_queries_main_writes_headers_and_branch_files(tmp_path: Path) -> None:
    counts, fasta, status = write_candidate_evidence(tmp_path, "sample", [], "")
    query = tmp_path / "queries.fasta"
    terminal_summary = tmp_path / "terminal.tsv"

    prepare_read_blast_queries_main([
            "--design-id", "design",
            "--counts", str(counts), "--fastas", str(fasta), "--statuses", str(status),
            "--candidate-counts-out", str(tmp_path / "published.tsv"), "--query-out", str(query),
            "--summary-out", str(tmp_path / "summary.tsv"), "--terminal-blast-hits-out", str(tmp_path / "hits.tsv"),
            "--terminal-classified-reads-out", str(tmp_path / "classified.tsv"), "--terminal-read-counts-out", str(tmp_path / "counts.tsv"),
            "--terminal-curve-out", str(tmp_path / "curve.tsv"), "--terminal-summary-out", str(terminal_summary),
    ])
    assert (tmp_path / "published.tsv").read_text().splitlines()[0].endswith("representative_id")
    assert terminal_summary.exists()
    assert not query.exists()

    nonempty_counts, nonempty_fasta, nonempty_status = write_candidate_evidence(
        tmp_path,
        "nonempty",
        ["nonempty\tf\t\t5\t1\t1\tread\n"],
        ">read\nACGTA\n",
    )
    prepare_read_blast_queries_main([
        "--design-id", "design",
        "--counts", str(nonempty_counts), "--fastas", str(nonempty_fasta), "--statuses", str(nonempty_status),
        "--candidate-counts-out", str(tmp_path / "nonempty-published.tsv"), "--query-out", str(query),
        "--summary-out", str(tmp_path / "nonempty-summary.tsv"), "--terminal-summary-out", str(terminal_summary),
    ])
    assert query.exists()
    assert not terminal_summary.exists()


def write_classification_candidates(tmp_path: Path, rows: str) -> Path:
    candidates = tmp_path / "candidate_read_counts.tsv"
    candidates.write_text(
        "metagenome_id\tfragment_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\trepresentative_id\n"
        + rows,
    )
    return candidates


def write_classification_hits(tmp_path: Path, rows: str) -> Path:
    hits = tmp_path / "whole_read_blast_hits.tsv"
    hits.write_text(
        "qseqid\tqlen\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\t"
        "sstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n"
        + rows,
    )
    return hits


def classification_rows(tmp_path: Path, candidate_rows: str, hit_rows: str) -> list[dict[str, object]]:
    candidates = construct_candidate_read_counts(write_classification_candidates(tmp_path, candidate_rows))
    hits = construct_whole_read_hits(
        write_classification_hits(tmp_path, hit_rows),
        candidate_representatives=candidates.select("representative_id").unique(),
        target_taxid="88456",
    )
    return construct_candidate_read_classifications(candidates, hits).collect().to_dicts()


def blast_row(query: str, taxids: str, bitscore: str, *, qlen: str = "100") -> str:
    return f"{query}\t{qlen}\taccession\t{taxids}\t100\t100\t0\t0\t1\t100\t1\t100\t0\t{bitscore}\t100\ttitle\n"


def test_construct_candidate_read_classifications_preserves_candidate_order(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tf1\t1\t100\t4\t4\tr1\n"
        "sample\tf2\t1\t100\t3\t3\tr2\n"
        "sample\tf3\t1\t100\t2\t2\tr3\n"
        "sample\tf4\t1\t100\t1\t1\tr4\n",
        blast_row("r1", "88456", "100")
        + blast_row("r2", "123", "100")
        + blast_row("r3", "88456;123", "100"),
    )

    assert [(row["fragment_id"], row["classification"]) for row in rows] == [
        ("f1", "TARGET"),
        ("f2", "NON_TARGET"),
        ("f3", "TIED"),
        ("f4", "NO_HIT"),
    ]


def test_construct_candidate_read_classifications_uses_decimal_tie_boundary(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tf1\t\t100\t1\t1\tr1\n"
        "sample\tf2\t\t100\t1\t1\tr2\n",
        blast_row("r1", "88456", "100.1")
        + blast_row("r1", "123", "100")
        + blast_row("r2", "88456", "100.1001")
        + blast_row("r2", "123", "100"),
    )

    assert [row["classification"] for row in rows] == ["TIED", "TARGET"]


def test_construct_candidate_read_classifications_preserves_source_decimal_precision(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tf\t\t100\t1\t1\tr1\n",
        blast_row("r1", "88456", "100.100000000001")
        + blast_row("r1", "123", "100.000000000000"),
    )

    assert rows[0]["classification"] == "TARGET"
    assert rows[0]["best_target_bit_score"] == "100.100000000001"
    assert rows[0]["best_non_target_bit_score"] == "100"


def test_construct_whole_read_hits_uses_independent_best_hsp_scores(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tf\t\t100\t1\t1\tr1\n",
        blast_row("r1", "88456", "90")
        + blast_row("r1", "88456", "100")
        + blast_row("r1", "123", "80")
        + blast_row("r1", "123", "99"),
    )

    assert rows[0]["best_target_bit_score"] == "100"
    assert rows[0]["best_non_target_bit_score"] == "99"


def test_construct_whole_read_hits_treats_na_as_non_target_evidence(tmp_path: Path) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tf1\t\t100\t1\t1\tr1\n"
        "sample\tf2\t\t100\t1\t1\tr2\n",
        blast_row("r1", "N/A", "100"),
    )

    assert [(row["classification"], row["best_non_target_bit_score"]) for row in rows] == [
        ("NON_TARGET", "100"),
        ("NO_HIT", ""),
    ]


@pytest.mark.parametrize("taxids", ["88456,123", "88456; 123"])
def test_construct_whole_read_hits_accepts_taxid_delimiters(tmp_path: Path, taxids: str) -> None:
    rows = classification_rows(
        tmp_path,
        "sample\tf\t\t100\t1\t1\tr1\n",
        blast_row("r1", taxids, "100"),
    )

    assert rows[0]["classification"] == "TIED"


@pytest.mark.parametrize("taxids", ["88456;", "0", "-1", "088456"])
def test_construct_whole_read_hits_rejects_malformed_taxids(tmp_path: Path, taxids: str) -> None:
    candidates = construct_candidate_read_counts(
        write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n"),
    )

    with pytest.raises(CandidateReadClassificationError, match="malformed staxids"):
        construct_whole_read_hits(
            write_classification_hits(tmp_path, blast_row("r1", taxids, "100")),
            candidate_representatives=candidates.select("representative_id"),
            target_taxid="88456",
        )


def test_construct_whole_read_hits_rejects_unknown_query(tmp_path: Path) -> None:
    candidates = construct_candidate_read_counts(
        write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n"),
    )

    with pytest.raises(CandidateReadClassificationError, match="unknown query"):
        construct_whole_read_hits(
            write_classification_hits(tmp_path, blast_row("unknown", "88456", "100")),
            candidate_representatives=candidates.select("representative_id"),
            target_taxid="88456",
        )


@pytest.mark.parametrize("count", ["bad", "", str(2**63)])
def test_construct_candidate_read_counts_rejects_malformed_integer_evidence(tmp_path: Path, count: str) -> None:
    with pytest.raises(CandidateReadClassificationError, match="malformed"):
        construct_candidate_read_counts(
            write_classification_candidates(tmp_path, f"sample\tf\t\t{count}\t1\t1\tr1\n"),
        )


def test_construct_candidate_read_counts_rejects_duplicate_identity_but_accepts_representatives(tmp_path: Path) -> None:
    repeated = construct_candidate_read_counts(
        write_classification_candidates(
            tmp_path,
            "sample\tf1\t\t100\t1\t1\tr1\n"
            "sample\tf2\t\t100\t1\t1\tr1\n",
        ),
    )
    assert repeated.collect().height == 2

    with pytest.raises(CandidateReadClassificationError, match="identity is duplicated"):
        construct_candidate_read_counts(
            write_classification_candidates(
                tmp_path,
                "sample\tf\t\t100\t1\t1\tr1\n"
                "sample\tf\t\t100\t1\t1\tr2\n",
            ),
        )


@pytest.mark.parametrize(
    "header",
    [
        "renamed\tfragment_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\trepresentative_id\n",
        "fragment_id\tmetagenome_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\trepresentative_id\n",
    ],
)
def test_construct_candidate_read_counts_rejects_invalid_headers(tmp_path: Path, header: str) -> None:
    candidates = write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n")
    candidates.write_text(header + "sample\tf\t\t100\t1\t1\tr1\n")

    with pytest.raises(CandidateReadClassificationError, match="Candidate-read counts are malformed"):
        construct_candidate_read_counts(candidates)


@pytest.mark.parametrize(
    "header",
    [
        "renamed\tqlen\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n",
        "qlen\tqseqid\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n",
    ],
)
def test_construct_whole_read_hits_rejects_invalid_headers(tmp_path: Path, header: str) -> None:
    candidates = construct_candidate_read_counts(
        write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n"),
    )
    hits = write_classification_hits(tmp_path, blast_row("r1", "88456", "100"))
    hits.write_text(header + blast_row("r1", "88456", "100"))

    with pytest.raises(CandidateReadClassificationError, match="Whole-read BLAST hits are malformed"):
        construct_whole_read_hits(
            hits,
            candidate_representatives=candidates.select("representative_id"),
            target_taxid="88456",
        )


@pytest.mark.parametrize(
    ("field", "error"),
    [("bitscore", "invalid bitscore"), ("qlen", "malformed")],
)
def test_construct_whole_read_hits_rejects_malformed_numeric_evidence(
    tmp_path: Path,
    field: str,
    error: str,
) -> None:
    values = {"bitscore": "100", "qlen": "100"}
    values[field] = "bad"
    candidates = construct_candidate_read_counts(
        write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n"),
    )

    with pytest.raises(CandidateReadClassificationError, match=error):
        construct_whole_read_hits(
            write_classification_hits(tmp_path, blast_row("r1", "88456", values["bitscore"], qlen=values["qlen"])),
            candidate_representatives=candidates.select("representative_id"),
            target_taxid="88456",
        )


@pytest.mark.parametrize("bitscore", ["NaN", "inf"])
def test_construct_whole_read_hits_rejects_nonfinite_bitscore(tmp_path: Path, bitscore: str) -> None:
    candidates = construct_candidate_read_counts(
        write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n"),
    )

    with pytest.raises(CandidateReadClassificationError, match="invalid bitscore"):
        construct_whole_read_hits(
            write_classification_hits(tmp_path, blast_row("r1", "88456", bitscore)),
            candidate_representatives=candidates.select("representative_id"),
            target_taxid="88456",
        )


@pytest.mark.parametrize("target_taxid", ["0", "-1", "bad", "088456"])
def test_construct_whole_read_hits_rejects_invalid_target_taxid(tmp_path: Path, target_taxid: str) -> None:
    with pytest.raises(CandidateReadClassificationError, match="target_taxid"):
        construct_whole_read_hits(
            write_classification_hits(tmp_path, ""),
            candidate_representatives=pl.LazyFrame({"representative_id": []}, schema={"representative_id": pl.String}),
            target_taxid=target_taxid,
        )


def test_classify_read_hits_main_writes_exact_header_and_blank_scores(tmp_path: Path) -> None:
    candidates = write_classification_candidates(tmp_path, "sample\tf\t\t100\t1\t1\tr1\n")
    hits = write_classification_hits(tmp_path, "")
    output = tmp_path / "classified_reads.tsv"

    classify_read_hits_main([
        "--candidate-read-counts", str(candidates),
        "--blast-hits", str(hits),
        "--target-taxid", "88456",
        "--output", str(output),
    ])

    assert output.read_text() == (
        "metagenome_id\tfragment_id\tmate\tread_length\tdistinct_bait_count\t"
        "fragment_distinct_bait_count\trepresentative_id\tclassification\t"
        "best_target_bit_score\tbest_non_target_bit_score\n"
        "sample\tf\t\t100\t1\t1\tr1\tNO_HIT\t\t\n"
    )


CLASSIFIED_HEADER = (
    "metagenome_id\tfragment_id\tmate\tread_length\tdistinct_bait_count\t"
    "fragment_distinct_bait_count\trepresentative_id\tclassification\t"
    "best_target_bit_score\tbest_non_target_bit_score\n"
)


def write_classified_reads(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "classified_reads.tsv"
    path.write_text(CLASSIFIED_HEADER + rows)
    return path


def write_preparation_summary(
    tmp_path: Path,
    *,
    rows: str,
) -> Path:
    path = tmp_path / "preparation_summary.tsv"
    path.write_text("metric\tvalue\n" + rows)
    return path


def preparation_rows(candidate_count: int) -> str:
    return (
        "design_id\tdesign\n"
        "deacon_returned_read_count\t9\n"
        "duplicate_fragment_count\t2\n"
        "zero_bait_read_count\t3\n"
        f"candidate_read_count\t{candidate_count}\n"
        "whole_read_blast_query_count\t4\n"
    )


def calibration(tmp_path: Path, classified_rows: str) -> ThresholdCalibration:
    classified = construct_classified_reads(write_classified_reads(tmp_path, classified_rows))
    row_count = len(classified_rows.strip().splitlines())
    preparation = construct_preparation_summary(
        write_preparation_summary(tmp_path, rows=preparation_rows(row_count)),
        "design",
    )
    return construct_threshold_calibration(classified, preparation)


def test_construct_threshold_calibration_builds_the_worked_multi_metagenome_relation(tmp_path: Path) -> None:
    result = calibration(
        tmp_path,
        "alpha\tf1\t1\t100\t4\t4\tr1\tTARGET\t100\t\n"
        "alpha\tf2\t1\t100\t2\t2\tr2\tNON_TARGET\t\t100\n"
        "beta\tf3\t1\t100\t3\t3\tr3\tTIED\t100\t100\n"
        "beta\tf4\t1\t100\t5\t5\tr4\tNO_HIT\t\t\n",
    )

    assert result.read_counts.collect().to_dicts() == [
        {"threshold": 1, "target_read_count": 1, "non_target_read_count": 1, "tied_read_count": 1, "no_hit_read_count": 1},
        {"threshold": 2, "target_read_count": 1, "non_target_read_count": 1, "tied_read_count": 1, "no_hit_read_count": 1},
        {"threshold": 3, "target_read_count": 1, "non_target_read_count": 0, "tied_read_count": 1, "no_hit_read_count": 1},
        {"threshold": 4, "target_read_count": 1, "non_target_read_count": 0, "tied_read_count": 0, "no_hit_read_count": 1},
        {"threshold": 5, "target_read_count": 0, "non_target_read_count": 0, "tied_read_count": 0, "no_hit_read_count": 1},
        {"threshold": 6, "target_read_count": 0, "non_target_read_count": 0, "tied_read_count": 0, "no_hit_read_count": 0},
    ]
    assert result.curve.collect().to_dicts() == [
        {"threshold": 1, "retained_metagenome_count": 2, "retained_fragment_count": 4},
        {"threshold": 2, "retained_metagenome_count": 2, "retained_fragment_count": 4},
        {"threshold": 3, "retained_metagenome_count": 2, "retained_fragment_count": 3},
        {"threshold": 4, "retained_metagenome_count": 2, "retained_fragment_count": 2},
        {"threshold": 5, "retained_metagenome_count": 1, "retained_fragment_count": 1},
        {"threshold": 6, "retained_metagenome_count": 0, "retained_fragment_count": 0},
    ]
    assert result.conclusion.status is CalibrationStatus.RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD
    assert result.conclusion.recommended_threshold == 4
    assert result.summary.collect().to_dicts()[-8:] == [
        {"metric": "target_classified_read_count", "value": "1"},
        {"metric": "non_target_classified_read_count", "value": "1"},
        {"metric": "tied_read_count", "value": "1"},
        {"metric": "no_hit_read_count", "value": "1"},
        {"metric": "calibration_status", "value": "RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD"},
        {"metric": "recommended_deacon_absolute_threshold", "value": "4"},
        {"metric": "specificity_floor", "value": ""},
        {"metric": "conclusion", "value": "The recommended Deacon absolute threshold is 4 for this optimization read set."},
    ]


def test_construct_threshold_calibration_recommends_one_when_one_is_sufficient(tmp_path: Path) -> None:
    result = calibration(tmp_path, "sample\ttarget\t\t100\t1\t1\tr1\tTARGET\t100\t\n")

    assert result.conclusion.status is CalibrationStatus.RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD
    assert result.conclusion.recommended_threshold == 1


def test_construct_threshold_calibration_uses_shared_fragment_wide_bait_counts(tmp_path: Path) -> None:
    result = calibration(
        tmp_path,
        "sample\tfragment\t1\t100\t1\t2\tr1\tTARGET\t100\t\n"
        "sample\tfragment\t2\t100\t2\t2\tr2\tNON_TARGET\t\t100\n",
    )

    assert result.read_counts.collect().to_dicts() == [
        {"threshold": 1, "target_read_count": 1, "non_target_read_count": 1, "tied_read_count": 0, "no_hit_read_count": 0},
        {"threshold": 2, "target_read_count": 1, "non_target_read_count": 1, "tied_read_count": 0, "no_hit_read_count": 0},
        {"threshold": 3, "target_read_count": 0, "non_target_read_count": 0, "tied_read_count": 0, "no_hit_read_count": 0},
    ]


def test_construct_threshold_calibration_reports_the_first_specificity_floor(tmp_path: Path) -> None:
    result = calibration(
        tmp_path,
        "sample\ttarget\t\t100\t2\t2\tr1\tTARGET\t100\t\n"
        "sample\tnon-target\t\t100\t3\t3\tr2\tNON_TARGET\t\t100\n",
    )

    assert result.conclusion.status is CalibrationStatus.SPECIFICITY_FLOOR
    assert result.conclusion.specificity_floor == 4
    assert result.conclusion.conclusion == "The specificity floor is 4; no target-classified read remains."


def test_construct_threshold_calibration_reports_all_no_hit_reads(tmp_path: Path) -> None:
    result = calibration(tmp_path, "sample\tf\t\t100\t2\t2\tr\tNO_HIT\t\t\n")

    assert result.read_counts.collect().get_column("threshold").to_list() == [1, 2, 3]
    assert result.conclusion.status is CalibrationStatus.NO_CLASSIFIED_READS
    assert result.conclusion.recommended_threshold is None
    assert result.conclusion.specificity_floor is None
    assert result.conclusion.conclusion == "Every candidate read is a no-hit read; no threshold is supported."


def test_construct_threshold_calibration_allows_no_hit_reads_at_a_recommendation(tmp_path: Path) -> None:
    result = calibration(
        tmp_path,
        "sample\ttarget\t\t100\t3\t3\tr1\tTARGET\t100\t\n"
        "sample\tnon-target\t\t100\t2\t2\tr2\tNON_TARGET\t\t100\n"
        "sample\tno-hit\t\t100\t3\t3\tr3\tNO_HIT\t\t\n",
    )

    assert result.conclusion.status is CalibrationStatus.RECOMMENDED_DEACON_ABSOLUTE_THRESHOLD
    assert result.conclusion.recommended_threshold == 3
    assert result.read_counts.collect().filter(pl.col("threshold") == 3).item(0, "no_hit_read_count") == 1


@pytest.mark.parametrize(
    "rows",
    [
        (
            "sample\tf\t1\t100\t1\t2\tr1\tTARGET\t\t\n"
            "sample\tf\t2\t100\t1\t3\tr2\tTARGET\t\t\n"
        ),
        (
            "sample\tf\t\t100\t1\t1\tr1\tTARGET\t\t\n"
            "sample\tf\t\t100\t1\t1\tr2\tTARGET\t\t\n"
        ),
        "sample\tf\t\t100\t1\t1\tr\tUNKNOWN\t\t\n",
    ],
)
def test_construct_classified_reads_rejects_fragment_consistency_identity_and_classification(
    tmp_path: Path,
    rows: str,
) -> None:
    with pytest.raises(DeaconThresholdCalibrationError):
        construct_classified_reads(write_classified_reads(tmp_path, rows))


@pytest.mark.parametrize(
    "row",
    [
        "\tf\t\t100\t1\t1\tr\tTARGET\t\t\n",
        "sample\tf\t\t\t1\t1\tr\tTARGET\t\t\n",
        "sample\tf\t\tbad\t1\t1\tr\tTARGET\t\t\n",
        "sample\tf\t\t100\tbad\t1\tr\tTARGET\t\t\n",
        "sample\tf\t\t100\t0\t1\tr\tTARGET\t\t\n",
        "sample\tf\t\t100\t1\t0\tr\tTARGET\t\t\n",
        f"sample\tf\t\t{2**63}\t1\t1\tr\tTARGET\t\t\n",
        f"sample\tf\t\t100\t{2**63}\t{2**63}\tr\tTARGET\t\t\n",
    ],
)
def test_construct_classified_reads_rejects_invalid_read_and_bait_evidence(tmp_path: Path, row: str) -> None:
    with pytest.raises(DeaconThresholdCalibrationError):
        construct_classified_reads(write_classified_reads(tmp_path, row))


@pytest.mark.parametrize(
    "header",
    [
        "renamed\tfragment_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\trepresentative_id\tclassification\tbest_target_bit_score\tbest_non_target_bit_score\n",
        "fragment_id\tmetagenome_id\tmate\tread_length\tdistinct_bait_count\tfragment_distinct_bait_count\trepresentative_id\tclassification\tbest_target_bit_score\tbest_non_target_bit_score\n",
    ],
)
def test_construct_classified_reads_rejects_invalid_headers(tmp_path: Path, header: str) -> None:
    classified = write_classified_reads(tmp_path, "sample\tf\t\t100\t1\t1\tr\tTARGET\t100\t\n")
    classified.write_text(header + "sample\tf\t\t100\t1\t1\tr\tTARGET\t100\t\n")

    with pytest.raises(DeaconThresholdCalibrationError, match="Classified-read evidence is malformed"):
        construct_classified_reads(classified)


def test_construct_preparation_summary_rejects_renamed_header(tmp_path: Path) -> None:
    preparation = write_preparation_summary(tmp_path, rows=preparation_rows(1))
    preparation.write_text(preparation.read_text().replace("metric\tvalue", "renamed\tvalue", 1))

    with pytest.raises(DeaconThresholdCalibrationError, match="preparation summary is malformed"):
        construct_preparation_summary(preparation, "design")


@pytest.mark.parametrize(
    "rows",
    [
        "design_id\tother\n" + preparation_rows(1).removeprefix("design_id\tdesign\n"),
        preparation_rows(1).replace("candidate_read_count\t1\n", "candidate_read_count\t1\ncandidate_read_count\t1\n"),
        preparation_rows(1).replace("whole_read_blast_query_count\t4\n", ""),
        preparation_rows(1).replace("zero_bait_read_count\t3\n", "zero_bait_read_count\tbad\n"),
        preparation_rows(1).replace("duplicate_fragment_count\t2\n", "duplicate_fragment_count\t-1\n"),
        preparation_rows(1).replace("whole_read_blast_query_count\t4\n", f"whole_read_blast_query_count\t{2**63}\n"),
    ],
)
def test_construct_preparation_summary_rejects_invalid_metrics(tmp_path: Path, rows: str) -> None:
    with pytest.raises(DeaconThresholdCalibrationError):
        construct_preparation_summary(write_preparation_summary(tmp_path, rows=rows), "design")


def test_construct_threshold_calibration_rejects_preparation_candidate_count_disagreement(tmp_path: Path) -> None:
    classified = construct_classified_reads(write_classified_reads(tmp_path, "sample\tf\t\t100\t1\t1\tr\tTARGET\t\t\n"))
    preparation = construct_preparation_summary(
        write_preparation_summary(tmp_path, rows=preparation_rows(2)),
        "design",
    )

    with pytest.raises(DeaconThresholdCalibrationError):
        construct_threshold_calibration(classified, preparation)


def test_construct_threshold_calibration_is_invariant_to_classified_read_order(tmp_path: Path) -> None:
    rows = (
        "beta\tf2\t\t100\t3\t3\tr2\tNON_TARGET\t\t100\n"
        "alpha\tf1\t\t100\t2\t2\tr1\tTARGET\t100\t\n"
        "beta\tf3\t\t100\t4\t4\tr3\tNO_HIT\t\t\n"
    )
    reversed_rows = "".join(reversed(rows.splitlines(keepends=True)))
    first = calibration(tmp_path, rows)
    other_directory = tmp_path / "permuted"
    other_directory.mkdir()
    second = calibration(other_directory, reversed_rows)

    assert first.read_counts.collect().to_dicts() == second.read_counts.collect().to_dicts()
    assert first.curve.collect().to_dicts() == second.curve.collect().to_dicts()
    assert first.summary.collect().to_dicts() == second.summary.collect().to_dicts()


def test_calibrate_deacon_threshold_main_writes_exact_output_contract(tmp_path: Path) -> None:
    classified = write_classified_reads(tmp_path, "sample\tf\t\t100\t1\t1\tr\tNO_HIT\t\t\n")
    preparation = write_preparation_summary(tmp_path, rows=preparation_rows(1))
    read_counts = tmp_path / "read_counts.tsv"
    curve = tmp_path / "curve.tsv"
    summary = tmp_path / "summary.tsv"

    calibrate_deacon_threshold_main([
        "--design-id", "design",
        "--classified-reads", str(classified),
        "--preparation-summary", str(preparation),
        "--read-counts-out", str(read_counts),
        "--curve-out", str(curve),
        "--summary-out", str(summary),
    ])

    assert read_counts.read_text().splitlines()[0] == "threshold\ttarget_read_count\tnon_target_read_count\ttied_read_count\tno_hit_read_count"
    assert curve.read_text().splitlines()[0] == "threshold\tretained_metagenome_count\tretained_fragment_count"
    assert summary.read_text().splitlines() == [
        "metric\tvalue",
        "design_id\tdesign",
        "deacon_returned_read_count\t9",
        "duplicate_fragment_count\t2",
        "zero_bait_read_count\t3",
        "candidate_read_count\t1",
        "whole_read_blast_query_count\t4",
        "target_classified_read_count\t0",
        "non_target_classified_read_count\t0",
        "tied_read_count\t0",
        "no_hit_read_count\t1",
        "calibration_status\tNO_CLASSIFIED_READS",
        "recommended_deacon_absolute_threshold\t",
        "specificity_floor\t",
        "conclusion\tEvery candidate read is a no-hit read; no threshold is supported.",
    ]
