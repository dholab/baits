#!/usr/bin/env python3
"""Construct Provisional Source Sequences from Candidate Loci."""

from __future__ import annotations

import argparse
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
BLAST_SCHEMA = {
    "query_id": pl.String,
    "query_length": pl.Int64,
    "query_start": pl.Int64,
    "query_end": pl.Int64,
    "assembly_sequence_id": pl.String,
    "assembly_start_raw": pl.Int64,
    "assembly_end_raw": pl.Int64,
    "alignment_length": pl.Int64,
    "percent_identity": pl.Float64,
}
CANDIDATE_LOCUS_COLUMNS = [
    "candidate_locus_id",
    "source_sequence_id",
    "query_id",
    "query_group",
    "assembly_sequence_id",
    "assembly_start",
    "assembly_end",
    "strand",
    "percent_identity",
    "query_coverage",
]


class QueryRulesError(ValueError):
    """Raised when Query Rules do not satisfy their tabular contract."""


class SequenceInputError(ValueError):
    """Raised when a sequence input does not satisfy its contract."""


class BlastHitInputError(ValueError):
    """Raised when BLAST output cannot construct Candidate Loci."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-queries", type=Path, required=True)
    parser.add_argument("--target-assembly", type=Path, required=True)
    parser.add_argument("--query-rules", type=Path, required=True)
    parser.add_argument("--blast-hits", type=Path, required=True)
    parser.add_argument("--source-sequences-out", type=Path, required=True)
    parser.add_argument("--source-sequence-query-groups-out", type=Path, required=True)
    parser.add_argument("--candidate-loci-out", type=Path, required=True)
    parser.add_argument("--discovery-status-out", type=Path, required=True)
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
        raise QueryRulesError(message) from error
    if columns != list(QUERY_RULE_SCHEMA):
        expected = ", ".join(QUERY_RULE_SCHEMA)
        message = f"Query Rules columns must be exactly: {expected}"
        raise QueryRulesError(message)

    query_id = pl.col("query_id").fill_null("")
    invalid_rule = (
        pl.when((query_id == "") | (pl.len().over("query_id") > 1))
        .then(pl.concat_str(pl.lit("missing or duplicated query_id: "), query_id))
        .when(pl.col("query_group").fill_null("") == "")
        .then(pl.concat_str(pl.lit("missing Query Group for "), query_id))
        .when(
            pl.col("min_identity").is_null()
            | pl.col("min_query_coverage").is_null()
            | ~pl.col("min_identity").is_between(0, 100, closed="both")
            | ~pl.col("min_query_coverage").is_between(0, 100, closed="both"),
        )
        .then(pl.concat_str(pl.lit("query thresholds out of range for "), query_id))
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
        raise QueryRulesError(message) from error
    if validation.item(0, "row_count") == 0:
        message = "Query Rules must contain at least one row"
        raise QueryRulesError(message)
    if first_error := validation.item(0, "first_error"):
        raise QueryRulesError(first_error)
    return parsed.drop("_invalid_rule")


def scan_blast_hits(path: Path) -> pl.LazyFrame:
    if path.stat().st_size == 0:
        return pl.LazyFrame(schema=BLAST_SCHEMA)
    return pl.scan_csv(
        path,
        separator="\t",
        quote_char=None,
        has_header=False,
        comment_prefix="#",
        schema=BLAST_SCHEMA,
    )


def blast_hit_error() -> pl.Expr:
    query_id = pl.col("query_id")
    assembly_id = pl.col("assembly_sequence_id")
    invalid_query_coordinates = (
        (pl.col("query_start") < 1)
        | (pl.col("query_start") > pl.col("query_end"))
        | (pl.col("query_end") > pl.col("query_length"))
    )
    invalid_assembly_coordinates = (
        (pl.col("assembly_start_raw") < 1)
        | (pl.col("assembly_end_raw") < 1)
        | (
            pl.max_horizontal("assembly_start_raw", "assembly_end_raw")
            > pl.col("assembly_sequence_length")
        )
    )
    inconsistent_alignment_length = (
        (pl.col("alignment_length") < pl.col("query_end") - pl.col("query_start") + 1)
        | (
            pl.col("alignment_length")
            < (pl.col("assembly_end_raw") - pl.col("assembly_start_raw")).abs() + 1
        )
    )
    return pl.coalesce(
        pl.when(pl.any_horizontal([pl.col(column).is_null() for column in BLAST_SCHEMA]))
        .then(pl.lit("missing BLAST field")),
        pl.when(pl.col("query_group").is_null() | pl.col("prepared_query_length").is_null())
        .then(pl.concat_str(pl.lit("unknown query ID: "), query_id)),
        pl.when(pl.col("assembly_sequence_length").is_null())
        .then(pl.concat_str(pl.lit("unknown assembly ID: "), assembly_id)),
        pl.when(pl.col("query_length") <= 0)
        .then(pl.concat_str(pl.lit("invalid query length for "), query_id)),
        pl.when(pl.col("query_length") != pl.col("prepared_query_length"))
        .then(pl.concat_str(pl.lit("inconsistent query length for "), query_id)),
        pl.when(invalid_query_coordinates)
        .then(pl.concat_str(pl.lit("invalid query coordinates for "), query_id)),
        pl.when(invalid_assembly_coordinates)
        .then(pl.concat_str(pl.lit("invalid assembly coordinates for "), assembly_id)),
        pl.when(pl.col("alignment_length") <= 0)
        .then(pl.concat_str(pl.lit("invalid alignment length for "), query_id)),
        pl.when(
            ~pl.col("percent_identity").is_finite()
            | ~pl.col("percent_identity").is_between(0, 100),
        )
        .then(pl.concat_str(pl.lit("invalid percent identity for "), query_id)),
        pl.when(inconsistent_alignment_length)
        .then(pl.concat_str(pl.lit("inconsistent alignment length for "), query_id)),
    )


def construct_candidate_loci(
    blast_hits: pl.LazyFrame,
    query_rules: pl.LazyFrame,
    prepared_queries: pl.LazyFrame,
    target_assembly: pl.LazyFrame,
) -> pl.LazyFrame:
    missing_query = (
        query_rules.join(prepared_queries, on="query_id", how="anti")
        .select("query_id")
        .limit(1)
        .collect()
    )
    if missing_query.height:
        query_id = missing_query.item(0, "query_id")
        message = f"Missing Representative Query FASTA record: {query_id}"
        raise SequenceInputError(message)

    contextualized = (
        blast_hits.with_row_index("_row")
        .join(
            query_rules.select(
                "query_id",
                "query_group",
                "min_identity",
                "min_query_coverage",
            ),
            on="query_id",
            how="left",
        )
        .join(prepared_queries, on="query_id", how="left")
        .join(target_assembly, on="assembly_sequence_id", how="left")
        .with_columns(blast_hit_error().alias("_blast_hit_error"))
    )
    first_error = (
        contextualized.filter(pl.col("_blast_hit_error").is_not_null())
        .sort("_row")
        .select("_blast_hit_error")
        .limit(1)
        .collect()
    )
    if first_error.height:
        raise BlastHitInputError(first_error.item(0, "_blast_hit_error"))

    normalized_interval = (
        pl.col("assembly_sequence")
        .str.slice(
            pl.col("assembly_start") - 1,
            pl.col("assembly_end") - pl.col("assembly_start") + 1,
        )
        .str.to_uppercase()
        .str.replace_all("U", "T")
    )
    return (
        contextualized.filter(pl.col("_blast_hit_error").is_null())
        .with_columns(
            pl.min_horizontal("assembly_start_raw", "assembly_end_raw").alias("assembly_start"),
            pl.max_horizontal("assembly_start_raw", "assembly_end_raw").alias("assembly_end"),
            pl.when(pl.col("assembly_start_raw") <= pl.col("assembly_end_raw"))
            .then(pl.lit("+"))
            .otherwise(pl.lit("-"))
            .alias("strand"),
            ((pl.col("query_end") - pl.col("query_start") + 1) / pl.col("query_length") * 100).alias(
                "query_coverage",
            ),
        )
        .filter(
            (pl.col("percent_identity") >= pl.col("min_identity"))
            & (pl.col("query_coverage") >= pl.col("min_query_coverage")),
        )
        .with_columns(normalized_interval.alias("_forward_sequence"))
        .with_columns(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("_forward_sequence"))
            .otherwise(
                pl.col("_forward_sequence").map_elements(
                    lambda sequence: str(Seq(sequence).reverse_complement()),
                    return_dtype=pl.String,
                ),
            )
            .alias("sequence"),
        )
        .sort(
            "query_group",
            "sequence",
            "query_id",
            "assembly_sequence_id",
            "assembly_start",
            "assembly_end",
            "strand",
            "percent_identity",
            "query_coverage",
        )
        .with_row_index("_candidate_locus_number", offset=1)
        .with_columns(
            pl.concat_str(
                pl.lit("candidate_locus_"),
                pl.col("_candidate_locus_number").cast(pl.String).str.pad_start(6, "0"),
            ).alias("candidate_locus_id"),
        )
        .select(
            "candidate_locus_id",
            "query_id",
            "query_group",
            "assembly_sequence_id",
            "assembly_start",
            "assembly_end",
            "strand",
            "percent_identity",
            "query_coverage",
            "sequence",
        )
    )


def construct_provisional_source_sequences(
    candidate_loci: pl.LazyFrame,
) -> pl.LazyFrame:
    return (
        candidate_loci.select("query_group", "sequence")
        .unique()
        .sort("query_group", "sequence")
        .with_row_index("_source_sequence_number", offset=1)
        .with_columns(
            pl.concat_str(
                pl.lit("source_sequence_"),
                pl.col("_source_sequence_number").cast(pl.String).str.pad_start(6, "0"),
            ).alias("source_sequence_id"),
        )
        .select("source_sequence_id", "query_group", "sequence")
    )


def construct_query_group_status(
    query_rules: pl.LazyFrame,
    provisional_source_sequences: pl.LazyFrame,
) -> pl.LazyFrame:
    return (
        query_rules.select("query_group")
        .unique()
        .join(
            provisional_source_sequences.group_by("query_group").len(),
            on="query_group",
            how="left",
        )
        .with_columns(
            pl.col("len").fill_null(0).alias("provisional_source_sequence_count"),
        )
        .with_columns(
            pl.when(pl.col("provisional_source_sequence_count") > 0)
            .then(pl.lit("PASS"))
            .otherwise(pl.lit("NO_CANDIDATE_LOCUS"))
            .alias("status"),
        )
        .select("query_group", "status", "provisional_source_sequence_count")
        .sort("query_group")
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    query_rules = scan_query_rules(args.query_rules).collect().lazy()
    try:
        prepared_records = SeqIO.to_dict(SeqIO.parse(args.prepared_queries, "fasta"))
        assembly_records = SeqIO.to_dict(SeqIO.parse(args.target_assembly, "fasta"))
    except ValueError as error:
        message = str(error)
        raise SequenceInputError(message) from error
    prepared_queries = pl.LazyFrame(
        {
            "query_id": list(prepared_records),
            "prepared_query_length": [len(record.seq) for record in prepared_records.values()],
        },
        schema={"query_id": pl.String, "prepared_query_length": pl.Int64},
    )
    assembly_sequences = [
        str(record.seq).upper().replace("U", "T")
        for record in assembly_records.values()
    ]
    target_assembly = pl.LazyFrame(
        {
            "assembly_sequence_id": list(assembly_records),
            "assembly_sequence": assembly_sequences,
            "assembly_sequence_length": list(map(len, assembly_sequences)),
        },
        schema={
            "assembly_sequence_id": pl.String,
            "assembly_sequence": pl.String,
            "assembly_sequence_length": pl.Int64,
        },
    )
    candidate_loci = (
        construct_candidate_loci(
            scan_blast_hits(args.blast_hits),
            query_rules,
            prepared_queries,
            target_assembly,
        )
        .collect()
        .lazy()
    )
    provisional_source_sequences = (
        construct_provisional_source_sequences(candidate_loci).collect().lazy()
    )
    status = construct_query_group_status(
        query_rules,
        provisional_source_sequences,
    ).collect()
    candidate_loci_output = (
        candidate_loci.join(
            provisional_source_sequences,
            on=["query_group", "sequence"],
            how="left",
        )
        .with_columns(
            pl.col("percent_identity").cast(pl.String).str.replace(r"\.0$", ""),
            pl.col("query_coverage").cast(pl.String).str.replace(r"\.0$", ""),
        )
        .select(CANDIDATE_LOCUS_COLUMNS)
    )

    candidate_loci_output.collect().write_csv(
        args.candidate_loci_out,
        separator="\t",
        quote_style="never",
    )
    status.write_csv(
        args.discovery_status_out,
        separator="\t",
        quote_style="never",
    )
    if status.get_column("provisional_source_sequence_count").gt(0).all():
        source_sequences = provisional_source_sequences.collect()
        SeqIO.write(
            (
                SeqRecord(Seq(sequence), id=source_sequence_id, description="")
                for source_sequence_id, _query_group, sequence in source_sequences.iter_rows()
            ),
            args.source_sequences_out,
            "fasta",
        )
        source_sequences.select("source_sequence_id", "query_group").write_csv(
            args.source_sequence_query_groups_out,
            separator="\t",
            quote_style="never",
        )


if __name__ == "__main__":
    main()
