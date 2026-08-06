# dholab/baits

`dholab/baits` will build auditable Bait Sets of Clean K-mers for biological targets. It will remove Candidate K-mers that have exact matches in configured Interference Backgrounds. It will then apply Global Exact-Match Validation. The name refers to the CLEAN algorithm from radio astronomy.

This repository is an initial structural scaffold. It does not yet contain an analysis workflow or define biological input files.

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

Run the scaffold tests:

```bash
pixi run --environment dev nf-test test
```

The project uses nf-core-compatible components, nf-schema for parameter validation, and nf-test for workflow tests.
