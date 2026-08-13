#!/usr/bin/env python3
"""Prepare representative queries for candidate locus discovery."""

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

QUERY_RULE_SCHEMA = {
    "query_id": pl.String,
    "query_group": pl.String,
    "query_start": pl.Int64,
    "query_end": pl.Int64,
    "min_identity": pl.Float64,
    "min_query_coverage": pl.Float64,
}


class QueryInputError(ValueError):
    """Raised when representative queries or Query Rules are invalid."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representative-queries", type=Path, required=True)
    parser.add_argument("--query-rules", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    return parser.parse_args(argv)


def scan_query_rules(path: Path) -> pl.LazyFrame:
    try:
        query_rules = pl.scan_csv(
            path,
            separator="\t",
            schema_overrides=QUERY_RULE_SCHEMA,
        )
        columns = query_rules.collect_schema().names()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Query Rules: {error}"
        raise QueryInputError(message) from error
    if columns != list(QUERY_RULE_SCHEMA):
        expected = ", ".join(QUERY_RULE_SCHEMA)
        message = f"Query Rules columns must be exactly: {expected}"
        raise QueryInputError(message)

    query_id = pl.col("query_id").fill_null("")
    invalid_rule = (
        pl.when((query_id == "") | (pl.len().over("query_id") > 1))
        .then(pl.concat_str(pl.lit("missing or duplicated query_id: "), query_id))
        .when(
            pl.col("query_start").is_null()
            | pl.col("query_end").is_null()
            | (pl.col("query_start") < 1)
            | (pl.col("query_end") < 0)
            | ((pl.col("query_end") != 0) & (pl.col("query_end") < pl.col("query_start"))),
        )
        .then(pl.concat_str(pl.lit("Query interval out of range for "), query_id))
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias("_invalid_rule")
    )
    parsed = query_rules.with_columns(invalid_rule)
    try:
        validation = parsed.select(
            pl.len().alias("row_count"),
            pl.col("_invalid_rule").drop_nulls().first().alias("first_error"),
        )
        validation = validation.collect()
    except pl.exceptions.PolarsError as error:
        message = f"Could not parse Query Rules: {error}"
        raise QueryInputError(message) from error
    if validation.item(0, "row_count") == 0:
        message = "Query Rules must contain at least one row"
        raise QueryInputError(message)
    if first_error := validation.item(0, "first_error"):
        raise QueryInputError(first_error)
    return parsed.drop("_invalid_rule")


def construct_prepared_queries(
    query_rules: pl.LazyFrame,
    representative_queries: pl.LazyFrame,
) -> pl.LazyFrame:
    joined = query_rules.with_row_index("_query_rule_order").join(
        representative_queries,
        on="query_id",
        how="left",
    ).with_columns(
        pl.when(pl.col("representative_query_sequence").is_null())
        .then(
            pl.concat_str(
                pl.lit("Missing representative query FASTA record: "),
                pl.col("query_id"),
            ),
        )
        .when(
            (pl.col("query_start") > pl.col("representative_query_length"))
            | (
                (pl.col("query_end") != 0)
                & (pl.col("query_end") > pl.col("representative_query_length"))
            ),
        )
        .then(
            pl.concat_str(
                pl.lit("Query interval out of range for "),
                pl.col("query_id"),
            ),
        )
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias("_invalid_query"),
    )
    first_error = (
        joined.filter(pl.col("_invalid_query").is_not_null())
        .select("_invalid_query")
        .limit(1)
        .collect()
    )
    if first_error.height:
        raise QueryInputError(first_error.item(0, "_invalid_query"))
    return (
        joined.filter(pl.col("_invalid_query").is_null())
        .with_columns(
            pl.when(pl.col("query_end") == 0)
            .then(pl.col("representative_query_length"))
            .otherwise(pl.col("query_end"))
            .alias("_resolved_query_end"),
        )
        .with_columns(
            pl.col("representative_query_sequence")
            .str.slice(
                pl.col("query_start") - 1,
                pl.col("_resolved_query_end") - pl.col("query_start") + 1,
            )
            .alias("sequence"),
        )
        .sort("_query_rule_order")
        .select("query_id", "sequence")
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    fasta_records = tuple(SeqIO.parse(args.representative_queries, "fasta"))
    duplicate_id = next(
        (
            identifier
            for identifier, count in Counter(record.id for record in fasta_records).items()
            if count > 1
        ),
        None,
    )
    if duplicate_id:
        message = f"Duplicate representative query FASTA record ID: {duplicate_id}"
        raise QueryInputError(message)
    records = SeqIO.to_dict(fasta_records)
    sequences = [str(record.seq).upper().replace("U", "T") for record in records.values()]
    representative_queries = pl.LazyFrame(
        {
            "query_id": list(records),
            "representative_query_sequence": sequences,
            "representative_query_length": list(map(len, sequences)),
        },
        schema={
            "query_id": pl.String,
            "representative_query_sequence": pl.String,
            "representative_query_length": pl.Int64,
        },
    )
    prepared_queries = construct_prepared_queries(
        scan_query_rules(args.query_rules),
        representative_queries,
    ).collect()
    SeqIO.write(
        (
            SeqRecord(Seq(sequence), id=query_id, description="")
            for query_id, sequence in prepared_queries.iter_rows()
        ),
        args.output_fasta,
        "fasta",
    )


if __name__ == "__main__":
    main()
