# dholab/baits

`dholab/baits` will build auditable Bait Sets of Clean K-mers for biological targets. It will remove Candidate K-mers that have exact matches in configured Interference Backgrounds. It will then apply Global Exact-Match Validation. The name refers to the CLEAN algorithm from radio astronomy.

This repository validates Target Sequence inputs but does not yet contain an analysis workflow.

For one curated execution, supply a Curated Target Sequence FASTA, the Target Taxon, and one Interference Background FASTA:

```bash
nextflow run dholab/baits \
    --target_sequences targets.fasta \
    --target_taxid 88456 \
    --interference_background interference.fasta
```

For one or more executions, use `--input` with a CSV samplesheet. See `assets/schema_input.json` for its columns. Query-guided rows refer to a separate query-rules TSV defined by `assets/schema_query_rules.json`.

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
```

The project uses nf-core-compatible components, nf-schema for parameter validation, and nf-test for workflow tests.
