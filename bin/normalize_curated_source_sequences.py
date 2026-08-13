#!/usr/bin/env python3
"""Normalize curated source sequences and write their query groups."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from collections.abc import Sequence


class CuratedSourceSequenceError(ValueError):
    """Raised when curated source sequences do not satisfy their contract."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sequences", type=Path, required=True)
    parser.add_argument("--source-sequences-out", type=Path, required=True)
    parser.add_argument("--source-sequence-query-groups-out", type=Path, required=True)
    return parser.parse_args(argv)


def normalize_curated_source_sequences(path: Path) -> tuple[SeqRecord, ...]:
    records = tuple(SeqIO.parse(path, "fasta"))
    if not records:
        message = "Curated source sequence FASTA must contain at least one record"
        raise CuratedSourceSequenceError(message)
    duplicate_ids = sorted(
        identifier
        for identifier, count in Counter(record.id for record in records).items()
        if count > 1
    )
    if duplicate_ids:
        message = f"Duplicate curated source sequence FASTA record ID: {duplicate_ids[0]}"
        raise CuratedSourceSequenceError(message)
    return tuple(
        SeqRecord(
            Seq(str(record.seq).upper().replace("U", "T")),
            id=record.id,
            description="",
        )
        for record in records
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    source_sequences = normalize_curated_source_sequences(args.source_sequences)
    SeqIO.write(source_sequences, args.source_sequences_out, "fasta")
    pl.DataFrame(
        {
            "source_sequence_id": [record.id for record in source_sequences],
            "query_group": ["" for _record in source_sequences],
        },
        schema={"source_sequence_id": pl.String, "query_group": pl.String},
    ).write_csv(
        args.source_sequence_query_groups_out,
        separator="\t",
        quote_style="never",
    )


if __name__ == "__main__":
    main()
