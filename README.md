# dholab/baits

`dholab/baits` will build auditable Bait Sets of Clean K-mers for biological targets. It will remove Candidate K-mers that have exact matches in configured Interference Backgrounds. When a Taxonomic Reference Database is supplied, it will then apply Taxonomic Exact-Match Screening.

This repository implements Source Sequence acquisition, Candidate K-mer filtering, optional Taxonomic Exact-Match Screening, verified Deacon Index construction, and evidence-qualified Deacon threshold calibration. It enumerates canonical Candidate K-mers, cancels exact Interference Background matches, and removes low-complexity k-mers with Deacon to produce a Locally Filtered Bait Set when survivors remain. When supplied a Taxonomic Reference Database, it screens those Baits with full-length exact BLAST matches under the Target Taxon.

For one curated execution, supply a Curated Source Sequence FASTA, the Target Taxon, and one Interference Background FASTA:

```bash
nextflow run dholab/baits \
    --source_sequences sources.fasta \
    --target_taxid 88456 \
    --interference_background interference.fasta
```

Add `--taxonomic_reference_db taxonomic_reference_db` to run Taxonomic Exact-Match Screening.

To calibrate the Deacon absolute threshold, also supply `--optimization_read_set optimization_reads`. This must be a flat directory containing single-end `<metagenome>.fastq[.gz]` files or complete `<metagenome>_R1.fastq[.gz]` and `<metagenome>_R2.fastq[.gz]` pairs. Calibration runs only for a Taxonomically Screened Bait Set backed by its verified Deacon Index and the Taxonomic Reference Database.

For one or more executions, use `--input` with a CSV samplesheet. See `assets/schema_input.json` for its columns. Query-guided rows refer to a separate query-rules TSV defined by `assets/schema_query_rules.json`.

Successful curated and query-guided Source Sequence acquisition publishes `results/<id>/01_source_sequences/source_sequences.fasta`. Query-guided discovery also publishes `candidate_loci.tsv`, `query_blast_hits.tsv`, and `discovery_status.tsv`. The internal `source_sequence_query_groups.tsv` table is generated for downstream stages but is not published.

Whenever Source Sequence acquisition continues to filtering, `results/<id>/02_candidate_kmers/` contains the complete `candidate_kmers.tsv` manifest, `candidate_kmer_occurrences.tsv`, and `filtering_status.tsv`. These tables remain published for biological terminal results. A nonempty local result additionally publishes `results/<id>/04_bait_sets/locally_filtered_baits.fasta`.

Taxonomic Exact-Match Screening is optional and uses the supplied Taxonomic Reference Database directory. BLAST searches use `-task blastn -word_size <kmer_size> -ungapped -perc_identity 100 -qcov_hsp_perc 100 -dust no`; an exact match outside the Target Taxon or without a usable taxonomic assignment rejects that Bait. Screening publishes its hits, decisions, status, nonempty Taxonomically Screened Bait Set, and the database's unparsed `blastdbcmd -info` self-report.

The pipeline builds one Deacon Index from the deepest justified nonempty Bait Set. It verifies that the index reproduces the Bait sequence set at `-a 1 -r 0` and retains no Interference Background records. The final Bait Set status is published in `04_bait_sets/`; the index, machine-readable verification summary, and readable report are published in `05_deacon_index/`.

Calibration publishes per-read and Read Fragment-wide distinct-Bait counts, whole-read BLAST evidence, classifications, threshold counts and retention curves, and `threshold_summary.tsv` in `06_calibration/`. For paired-end data, the fragment-wide count is the union across mates and is the value compared with Deacon `-a`. The summary reports a Recommended Deacon Absolute Threshold when supported, or the evidence-limited status `SPECIFICITY_FLOOR`, `NO_CLASSIFIED_READS`, or `NO_CANDIDATE_READS`. Optimization Read Set file identities and SHA-256 digests are included in `07_provenance/inputs.tsv`; fixed Deacon, BLAST, and tie-tolerance settings are added to `parameters.tsv` when calibration runs.

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
