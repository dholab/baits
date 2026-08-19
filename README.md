# dholab/baits

`dholab/baits` is an nf-core-style bioinformatics pipeline. It takes biological sequences from a particular taxon and boils them down to the set of k-mers--"baits"--that most confidently identify that taxon. Baits are like the *in silico* equivalent of PCR primers or hybrid capture probes. Use them with [Deacon](https://github.com/bede/deacon) to enrich target reads from shotgun reads, remove a known contaminant, or separate out closely related strains.

## Overview

`dholab/baits` is a generalization of the idea developed in our study [*Putative Cyclospora cayetanensis detection in wastewater metagenomic datasets*](https://github.com/dholab/cyclospora-in-wastewater-metagenomics). The corrected analysis has 1,670 baits before taxonomic exact-match screening; the final bait set remains pending a new `core_nt` run. With `dholab/baits` and the right choice of reference databases, users can identify baits for any taxon.

`dholab/baits` offers three major operations for each target taxon. First, all k-mer subsequences from the provided source sequences are checked against a "background", which can be any reference in FASTA format. For the aforementioned *Cyclospora* study, the ribosomal rRNA databases SILVA and Rfam together with NCBI GenBank *Cyclospora* sequences comprised the background. Any k-mers from the source sequences that are present in the background are rejected. Give `dholab/baits` source sequences and a background, and it will find the source k-mers that distinguish it from the background.

K-mers that are absent from the background are filtered for low complexity. The survivors can then go through the second operation: an exact, taxonomically informed screening against NCBI Core NT. Core NT can be swapped for another taxonomic BLAST database, and this step can also be skipped. If a database is provided, k-mers matching sequences from any taxon other than the target are rejected.

Third, `dholab/baits` can use provided sequencing reads to estimate how many bait matches are required before off-target hits go to zero. For many sequencing reads, a single bait match may still allow off-target hits. In the *Cyclospora* study, we found that requiring 20 bait matches reduced off-target hits to zero in the tested reads. `dholab/baits` can help you calibrate your own threshold for your own baits and your own data.

In summary, `dholab/baits` allows you to bring your own source sequences, your own background reference, your own BLAST-compatible taxonomic database, and your own calibration reads. In return, it gives you baits for your target, the evidence behind each screening decision, and the number of bait matches you should require to call a detection "good".

## Quick Start

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/), then launch the pipeline from the directory where you want its results and work files:

```bash
nextflow run dholab/baits \
    -r main \
    -profile docker \
    --source_sequences source_sequences.fasta \
    --target_taxid 88456 \
    --background background.fasta
```

Nextflow fetches and caches the pipeline under `$NXF_HOME/assets` (normally `~/.nextflow/assets`), so users do not need to clone this repository. Here, `88456` is the NCBI taxonomy ID for *C. cayetanensis*; substitute the ID for your target taxon. The bait FASTA is written under `results/<id>/04_bait_sets/`, and a verified Deacon index and its report are written under `results/<id>/05_deacon_index/`.

The example follows the current `main` branch. For reproducible analyses, replace `main` with a release tag or immutable commit revision.

## Pipeline Configuration

A basic run needs source sequences, an NCBI taxonomy ID, and a background FASTA. These are the parameters most users will control:

| Parameter | Default | Purpose |
|---|---|---|
| `--input` | — | Design CSV for one or more designs; replaces all direct-run inputs below. |
| `--source_sequences` | — | Direct-run source-sequence FASTA. |
| `--target_taxid` | — | Direct-run NCBI taxonomy ID for the target taxon. |
| `--background` | — | Direct-run FASTA containing sequences that the baits must not match. |
| `--id` | FASTA basename | Direct-run name used for the design and its results directory. |
| `--taxon_ref_db` | not run | Directory containing an optional taxonomic BLAST database. |
| `--calibration_reads` | not run | Direct-run flat FASTQ directory used to calibrate a threshold. Requires `--taxon_ref_db`. |
| `--calibration_target_taxids` | `target_taxid` only | Direct-run `taxid` TSV declaring target-compatible calibration assignments. Requires `--calibration_reads` and must include `target_taxid`. |
| `--max_blast_targets` | `100` | Maximum subject sequences retained per calibration-read BLAST query. |
| `--kmer_size` | `31` | Length of each candidate k-mer. |
| `--deacon_window` | `1` | Deacon minimizer window length. |
| `--entropy_threshold` | `0.6` | Minimum scaled entropy used for low-complexity filtering. |

`--input` selects design-CSV mode and cannot be combined with `--id`, `--source_sequences`, `--target_taxid`, `--background`, `--calibration_reads`, or `--calibration_target_taxids`. Put those values in each applicable CSV row instead. Run-level parameters such as `--taxon_ref_db`, `--kmer_size`, `--deacon_window`, and `--entropy_threshold` still apply to every design in the run.

Every process has a portable resource baseline. These requests are conservative starting points, not claims that one allocation is optimal for every dataset or executor:

| Process label | CPUs | Memory | Time |
|---|---:|---:|---:|
| `process_low` | 1 | 2 GB | 2 hours |
| `process_medium` | 2 | 4 GB | 8 hours |
| `process_high` | 2 | 4 GB | 24 hours |

BLAST searches, Meryl counting, and Deacon read filtering consume their allocated CPUs. `BLAST_MAKEBLASTDB` remains single-core because `makeblastdb` does not provide a corresponding thread-count option. A site configuration can override any baseline without modifying the pipeline:

```groovy
process {
    withLabel: process_medium {
        cpus = 4
        memory = '8 GB'
        time = '12h'
    }
}
```

Use Nextflow trace and report output from representative runs to tune requests for the actual source, background, taxonomic reference database, calibration reads, and executor. Keep queue, accounting, scratch, transfer, and container-cache policy in the site configuration.

For a complete direct run with taxonomic screening and threshold calibration:

```bash
nextflow run dholab/baits \
    -r main \
    -profile docker \
    --id my_design \
    --source_sequences source_sequences.fasta \
    --target_taxid 88456 \
    --background background.fasta \
    --taxon_ref_db /path/to/blast/database \
    --calibration_reads reads/ \
    --calibration_target_taxids calibration_target_taxids.tsv
```

The scope file has exactly one column. Taxids must be canonical positive decimal identifiers and unique:

```tsv
taxid
44417
88456
```

Use `--input designs.csv` instead of the direct source, taxon, and background parameters for multiple designs or query-guided source discovery. Copy and edit [`assets/example_designs.csv`](assets/example_designs.csv), then run:

```bash
nextflow run dholab/baits \
    -r main \
    -profile docker \
    --input designs.csv \
    --taxon_ref_db /path/to/blast/database
```

Run `nextflow run dholab/baits -r main --help` for the full parameter list. The design and query-rule schemas live under [`assets/`](assets/). See [Curated sequences and targeted designs](docs/designs.md) for a practical guide to the two input systems.

## Containers and Dependencies

Remote runs require Nextflow and a supported container runtime on the host. Each process declares a pinned container image, which Nextflow and the selected runtime fetch automatically. No Pixi environment or source checkout is required.

The commands above select the bundled `docker` profile, which enables Docker while retaining Nextflow's local executor. On a cluster with [Apptainer](https://apptainer.org/docs/admin/latest/installation.html), select `-profile apptainer` instead. Executor, queue, and storage settings can be supplied separately through the site's Nextflow configuration.

No runtime profile is selected automatically. Without one, tasks run in the host environment using Nextflow's default local executor; this is primarily useful for contributors working inside the Pixi environment or for users who manage every dependency themselves.

The main software components are:

| Component | Role |
|---|---|
| [Nextflow](https://www.nextflow.io/) | Host-side workflow runner that fetches the pipeline and orchestrates local or cluster execution. |
| [Meryl](https://github.com/marbl/meryl) | Counts source k-mers and finds exact occurrences in the background. |
| [Deacon](https://github.com/bede/deacon) | Filters low-complexity k-mers and builds the final read-filtering index. |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | Locates source sequences, screens baits taxonomically, and classifies calibration reads. |
| [Biopython](https://biopython.org/) and [Polars](https://pola.rs/) | Parse sequence data and construct the pipeline's evidence tables. |

Reference databases and sequencing reads remain user-supplied. In particular, the pipeline does not download or manage a taxonomic BLAST database such as `core_nt`.

## Development

Pixi is intended for contributors working from a source checkout. It provides locked environments for development on macOS ARM and Linux. Install them with:

```bash
pixi install --locked
```

[`Containerfile`](Containerfile) packages the complete locked development environment and repository source for CI or container-based development:

```bash
docker build --file Containerfile --tag nrminor/baits:dev .
docker run --rm --interactive --tty nrminor/baits:dev
```

Releases publish this development monoimage as `nrminor/baits:<version>`. It is separate from the process-specific images selected by the pipeline's `docker` and `apptainer` profiles.

The repository follows nf-core's modular Nextflow conventions: `main.nf` is the entry point, `workflows/` and `subworkflows/` compose the analysis, and `modules/` contains individual processes. The pipeline reuses nf-core modules where their contracts fit and keeps project-specific processes under `modules/local/`. [nf-schema](https://nextflow-io.github.io/nf-schema/latest/) validates command-line parameters and design CSVs against the schemas under `assets/`.

Inspect the local checkout with:

```bash
pixi run nextflow run . --help
pixi run nextflow run . --version
```

Pure Python tests cover scientific transformations and input contracts. nf-test exercises modules, subworkflows, and complete pipeline paths; Ruff and ty provide static checks.

```bash
pixi run --environment dev test-python
pixi run --environment dev lint-python
pixi run --environment dev typecheck-python
pixi run --environment dev nf-test test
```

Pass a test file to nf-test when working on one pipeline path:

```bash
pixi run --environment dev nf-test test tests/filter_candidate_kmers.nf.test
```
