from pathlib import Path

from Bio import SeqIO
from count_read_baits import canonical
from count_read_baits import main as count_read_baits_main


def test_candidate_read_stream_ignores_headers_and_preserves_observations(
    tmp_path: Path,
) -> None:
    baits = tmp_path / "baits.fasta"
    baits.write_text(">bait_1\nACGTA\n>bait_2\nCCCCC\n")
    reads = tmp_path / "candidate_reads.fasta"
    reads.write_text(
        ">repeated\nACGTA\n"
        ">repeated\nTACGT\n"
        ">another\nacgta\n"
        ">another\nCCCCC\n",
    )
    counts = tmp_path / "counts.tsv"
    candidates = tmp_path / "candidates.fasta"
    status = tmp_path / "status.tsv"

    count_read_baits_main(
        [
            "--metagenome-id",
            "sample",
            "--baits",
            str(baits),
            "--kmer-size",
            "5",
            "--read",
            str(reads),
            "--counts-out",
            str(counts),
            "--fasta-out",
            str(candidates),
            "--status-out",
            str(status),
        ],
    )

    assert counts.read_text().splitlines() == [
        "metagenome_id\tread_id\tread_length\tbait_count\trepresentative_id",
        "sample\tread_000000000001\t5\t1\tsequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219",
        "sample\tread_000000000002\t5\t1\tsequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219",
        "sample\tread_000000000003\t5\t1\tsequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219",
        "sample\tread_000000000004\t5\t1\tsequence_17b80bc751e1f35c75d6ada07267a5fd9981b7b3655cffcaf94a636e93eb27ba",
    ]
    assert [(record.id, str(record.seq)) for record in SeqIO.parse(candidates, "fasta")] == [
        (
            "sequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219",
            "ACGTA",
        ),
        (
            "sequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219",
            "ACGTA",
        ),
        (
            "sequence_8be635b3b461bfc5d03b08aa3ee554d044abedfa1ee7dc70ee5fbf7f665ae219",
            "ACGTA",
        ),
        (
            "sequence_17b80bc751e1f35c75d6ada07267a5fd9981b7b3655cffcaf94a636e93eb27ba",
            "CCCCC",
        ),
    ]
    assert status.read_text().splitlines() == [
        "metric\tvalue",
        "metagenome_id\tsample",
        "deacon_returned_read_count\t4",
        "candidate_read_count\t4",
    ]


def test_canonical_supports_the_full_iupac_dna_alphabet() -> None:
    assert canonical("ARY") == canonical("RYT") == "ARY"
