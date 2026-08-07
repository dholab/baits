# dholab/baits

`dholab/baits` will build auditable Bait Sets of Clean K-mers for biological targets. It will remove Candidate K-mers that have exact matches in configured Interference Backgrounds. When a Taxonomic Reference Database is supplied, it will then apply Taxonomic Exact-Match Screening.

This repository implements Source Sequence acquisition, Candidate K-mer filtering, and optional Taxonomic Exact-Match Screening. It enumerates canonical Candidate K-mers, cancels exact Interference Background matches, and removes low-complexity k-mers with Deacon to produce a Locally Filtered Bait Set when survivors remain. When supplied a Taxonomic Reference Database, it screens those Baits with full-length exact BLAST matches under the Target Taxon. Final Deacon index verification and threshold calibration remain for later commits.

For one curated execution, supply a Curated Source Sequence FASTA, the Target Taxon, and one Interference Background FASTA:

```bash
nextflow run dholab/baits \
    --source_sequences sources.fasta \
    --target_taxid 88456 \
    --interference_background interference.fasta
```

Add `--taxonomic_reference_db taxonomic_reference_db` to run Taxonomic Exact-Match Screening.

For one or more executions, use `--input` with a CSV samplesheet. See `assets/schema_input.json` for its columns. Query-guided rows refer to a separate query-rules TSV defined by `assets/schema_query_rules.json`.

Successful curated and query-guided Source Sequence acquisition publishes `results/<id>/01_source_sequences/source_sequences.fasta`. Query-guided discovery also publishes `candidate_loci.tsv`, `query_blast_hits.tsv`, and `discovery_status.tsv`. The internal `source_sequence_query_groups.tsv` table is generated for downstream stages but is not published.

Whenever Source Sequence acquisition continues to filtering, `results/<id>/02_candidate_kmers/` contains the complete `candidate_kmers.tsv` manifest, `candidate_kmer_occurrences.tsv`, and `filtering_status.tsv`. These tables remain published for biological terminal results. A nonempty local result additionally publishes `results/<id>/04_bait_sets/locally_filtered_baits.fasta`; internal Meryl and Deacon artifacts and the Bait Set status draft are not published.

Taxonomic Exact-Match Screening is optional and uses the supplied Taxonomic Reference Database directory. BLAST searches use `-task blastn -word_size <kmer_size> -ungapped -perc_identity 100 -qcov_hsp_perc 100 -dust no`; an exact match outside the Target Taxon or without a usable taxonomic assignment rejects that Bait. Screening publishes its hits, decisions, status, nonempty Taxonomically Screened Bait Set, and the database's unparsed `blastdbcmd -info` self-report.

If query-guided discovery finds no accepted Candidate Locus for one or more configured Query Groups, the run succeeds as a terminal scientific result for that design. It publishes BLAST hits, Candidate Loci, and Query Group status, but no Provisional Source Sequence FASTA.

`assets/cyclospora_rrna_query_rules.tsv` contains the real query configuration from the Cyclospora rRNA study. It is not a test fixture and does not include the reference sequence files required to run the study.

Run the pipeline with Nextflow:

```bash
nextflow run dholab/baits --help
nextflow run dholab/baits --version
```

## Development

Pixi provides the locked development and test environments. It is not required to run the pipeline.

Install the development environments:

```bash
pixi install --locked
```

Inspect the local checkout:

```bash
pixi run nextflow run . --help
pixi run nextflow run . --version
```

Run the project tests:

```bash
pixi run --environment dev nf-test test
pixi run --environment dev test-python
```

The project uses nf-core-compatible components, nf-schema for parameter validation, and nf-test for workflow tests.
