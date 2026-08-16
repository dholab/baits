# Curated sequences and targeted designs

`dholab/baits` has two data input systems, a simpler “curated source” system and a more sophisticated “designs” system. The curated source system is demonstrated in the project README; you give `dholab/baits` the source sequences to search, and it proceeds directly to k-mer-izing them and screening them against the background.

Sometimes, though, these sequences are non-trivial to curate. `dholab/baits` provides its designs input system to help. Using a designs CSV, users can provide raw representative sequences, locate them in a known target assembly, and apply trimming rules and labels in a “query rules” file. Users are free to mix these more sophisticated designs with the simpler curated input system depending on the needs of the target taxon as well as the limitations of the available assemblies.

## The two input systems

Each run obtains its source sequences through one of these systems:

```text
Curated sequence input

curated_source_sequences
          │
          ▼
accepted source sequences
          │
          ▼
candidate k-mers


Design input

representative_queries ── query_rules
          │
          ▼
prepared reference sequences
          │
          ▼
megablast against target_assembly
          │
          ▼
accepted assembly intervals
          │
          ▼
provisional source sequences
          │
          ▼
candidate k-mers
```

Both input systems then use the same filtering path:

```text
candidate k-mers
      │
      ├── cancel exact matches in background
      ├── remove low-complexity candidates
      └── optionally screen against a taxonomic reference database
```

(Like everywhere else in the pipeline, the candidate k-mer length is controlled by `--kmer_size` and defaults to 31.)

The designs CSV schema enforces one source-sequence path per row. A curated-sequence row must leave `representative_queries`, `query_rules`, and `target_assembly` empty. A design row must leave `curated_source_sequences` empty. See [`assets/schema_input.json`](../assets/schema_input.json) for the machine-readable contract.

| Input system | Required fields in addition to `id`, `target_taxid`, and `background` |
|---|---|
| Curated sequences | `curated_source_sequences` |
| Designs | `representative_queries`, `query_rules`, and `target_assembly` |

`calibration_reads` is optional with either input system. When calibration uses a target-compatible scope broader than `target_taxid`, `calibration_target_taxids` names a one-column TSV with the exact header `taxid`. Its positive unique taxids must include `target_taxid`. The scope changes calibration-read classification only; bait design and taxonomic exact-match screening continue to use `target_taxid`. `background` is required because background cancellation is part of every run.

## What each input contributes

| Input | Purpose | Supplies candidate k-mer sequence? |
|---|---|---:|
| `curated_source_sequences` | Provides source sequences that a curator has already accepted | Yes, directly |
| `representative_queries` | Provides representative sequences for the intended sequence classes | No |
| `query_rules` | Crops each representative sequence and sets its locus-acceptance rules | No |
| `target_assembly` | Provides the assembly intervals found by the representative sequences | Yes, after discovery |
| `background` | Declares the sequence collection whose exact k-mers must be cancelled | No; it is subtractive |
| `target_taxid` | Identifies the target taxon for later taxonomic interpretation | No |
| `calibration_target_taxids` | Optionally declares taxids treated as target-compatible only during calibration | No |

These inputs describe roles, not necessarily disjoint collections of organisms. The important distinction is whether an input generates candidate sequence or rejects it.

### Curated source sequences

`curated_source_sequences` is a FASTA file of sequences that you accept as direct sources of candidate k-mers. `dholab/baits` normalizes these sequences and enumerates their canonical k-mers. It does not perform a design step or search an assembly first.

Use curated sequence input when you can answer this question directly:

> Which exact sequences should generate candidate k-mers?

For example:

```csv
id,target_taxid,curated_source_sequences,representative_queries,query_rules,target_assembly,background,calibration_reads,calibration_target_taxids
curated_rrna,88456,data/source_sequences.fasta,,,,data/background.fasta,,
```

Supplying a curated source sequence is a scientific assertion. A record is curated because a curator accepted it, not because `dholab/baits` proved that it belongs to the target taxon.

### Representative sequences in a design

In the designs system, `representative_queries` is a FASTA file of known sequences that should map to the intended sequence classes in a target assembly. This can be an imperfect assembly contig or a consensus sequence that can be used to locate the closest homolog in an established (and ideally cleaner) reference assembly.

`dholab/baits` does not directly call candidate k-mers from the representative sequence. It uses that sequence to find an assembly interval, then enumerates k-mers from the interval. This keeps the resulting bait set tied to the target assembly rather than to an external accession used only to locate it. Want to call baits directly from your representative sequences instead? Use them as curated source sequences instead of through the designs input system.

Every `query_id` in the `query_rules` TSV must match one FASTA record identifier exactly. The FASTA can contain additional records, but records without design rules are not used.

### Target assembly

`target_assembly` is the positive search space for a design. `dholab/baits` builds a nucleotide BLAST database from this FASTA and searches it with the prepared representative sequences.

Each accepted BLAST interval is extracted from the assembly. Reverse-strand intervals are reverse-complemented so that their orientation agrees with the representative sequence. These extracted sequences become provisional source sequences.

The target assembly is not automatically a source of every possible k-mer. The design selects only assembly intervals that satisfy its declared rules. This is more targeted than enumerating every k-mer from an entire genome or transcriptome.

`dholab/baits` trusts the assembly identity supplied by the user. `target_taxid` does not verify that the assembly belongs to that taxon.

### Background

`background` is a FASTA collection of sequences whose exact canonical k-mers must not survive. It answers a different question:

> Which known sequence collection must the candidate k-mers distinguish themselves from?

The background does not locate source sequences. It cancels a candidate k-mer whenever that k-mer occurs exactly in the background.

A background can contain conserved references, related organisms, host sequence, expected environmental sequence, or target-derived sequence that you deliberately consider uninformative. The role is operational rather than purely taxonomic. However, if the background contains an entire source locus exactly, it can cancel every candidate from that locus.

Do not use the target assembly itself as the background unless cancelling all target-assembly k-mers is truly the intended rule.

## How design rules work

A design uses a `query_rules` TSV with these columns:

| Column | Meaning |
|---|---|
| `query_id` | FASTA record identifier for one representative sequence |
| `query_group` | User-supplied label for representative sequences that serve the same intended sequence class |
| `query_start` | First representative-sequence base to use, with 1-based inclusive coordinates |
| `query_end` | Last representative-sequence base to use, with 1-based inclusive coordinates; `0` means the end of the record |
| `min_identity` | Minimum BLAST percent identity for an accepted interval |
| `min_query_coverage` | Minimum percentage of the prepared reference-sequence span covered by an accepted interval |

Rules and representative sequences must also satisfy these constraints:

- Each `query_id` is unique in the rules and in the `representative_queries` FASTA.
- Each rule has a matching `representative_queries` FASTA record.
- `query_start` is at least 1 and does not exceed the FASTA record length.
- `query_end` is `0` or is at least `query_start` and no greater than the record length.
- Identity and coverage thresholds are numbers from 0 through 100.

The selected interval is applied before BLAST:

```text
representative FASTA record
          │
          ├── take query_start through query_end
          │   (`query_end = 0` means through the record end)
          ▼
prepared reference sequence
```

An assembly hit passes when both conditions are true:

```text
percent identity ≥ min_identity
reference-sequence span coverage ≥ min_query_coverage
```

Coverage is calculated from the reference-coordinate span reported by BLAST, divided by the prepared reference-sequence length. The accepted source sequence is the aligned assembly interval. `dholab/baits` does not extend that interval to an inferred gene boundary.

### Required sequence classes

Each `query_group` represents one required sequence class. Several `query_id` values can belong to the same group when alternative representative sequences can locate that class.

`dholab/baits` accepts all assembly intervals that pass their design rule. It does not select only the best hit or require exactly one locus. Overlapping or nested hits are not merged only because their coordinates overlap. Passing intervals are deduplicated when they produce the same extracted sequence in the same required group. The same sequence can remain as a separate source record in two different groups because the groups express different biological requirements. Later candidate-k-mer enumeration deduplicates identical canonical k-mer sequences.

Source-sequence construction continues only when every required group produces at least one provisional source sequence:

```text
every required group has an accepted sequence  → continue to candidate-k-mer generation
any required group has no accepted sequence    → terminal discovery result
```

The terminal result is not a workflow error. The run succeeds and publishes the BLAST hits, candidate-locus table, and discovery status for review. It does not publish a source-sequence FASTA or start candidate-k-mer filtering for that incomplete design. Automation can identify this result in `discovery_status.tsv` by the `NO_CANDIDATE_LOCUS` status.

## Cyclospora rRNA example

The Cyclospora analysis sought source sequences for four mature rRNA classes in a Cyclospora assembly:

```text
18S
28S
5.8S
5S
```

The bundled [`cyclospora_rrna_query_rules.tsv`](../assets/cyclospora_rrna_query_rules.tsv) contains one representative sequence for 18S, 28S, and 5.8S, plus six alternatives for 5S. All accepted hits require at least 98% identity and 90% coverage of the prepared reference sequence.

| Reference sequence | Group | Prepared interval | Minimum identity | Minimum coverage |
|---|---|---:|---:|---:|
| `AF111183.1` | `18S` | 1–1795 | 98% | 90% |
| `MPGL01000046.1` | `28S` | 11–3500 | 98% | 90% |
| `XR_003297357.1` | `5.8S` | Full record | 98% | 90% |
| Six `XR_...` records | `5S` | Full record | 98% | 90% |

The six 5S reference sequences are alternatives for one required group. They do not create six separate completeness requirements. At least one accepted 5S source sequence is required, just as at least one accepted sequence is required for each of 18S, 28S, and 5.8S.

```text
known rRNA accessions
        │
        ▼
representative sequences and design rules
        │
        ▼
Cyclospora target assembly
        │
        ▼
assembly-derived 18S, 28S, 5.8S, and 5S intervals
        │
        ▼
provisional source sequences
        │
        ▼
5,561 candidate 31-mers
        │
        ▼
background cancellation against the declared background
built from SILVA, Rfam, and GenBank-derived inputs
```

The representative accessions located the intended rRNA classes. The target assembly supplied the sequence that generated candidate k-mers. The background then removed exact 31-mers already present in the declared reference collection. These are three separate scientific roles.

This separation is what makes the design both targeted and auditable:

- Required groups state which sequence classes must be present.
- Design rules state what counts as an acceptable assembly locus.
- Candidate-locus evidence records where each provisional source came from.
- The background states which declared sequence collection each candidate must survive.

## Inspect design results before using the bait set

The designs system publishes these files under `results/<id>/01_source_sequences/`:

| File | Review question |
|---|---|
| `query_blast_hits.tsv` | What did each prepared representative sequence align to? |
| `candidate_loci.tsv` | Which assembly intervals passed, and with what identity and coverage? |
| `discovery_status.tsv` | Did every required group produce an accepted source sequence? |
| `source_sequences.fasta` | Which assembly-derived sequences generated candidate k-mers? |

Do not review only the final bait FASTA. First confirm that the candidate loci and provisional source sequences represent the intended biology. A high-scoring BLAST hit is evidence for a locus; it is not automatic curation.

## Practical guidance

Use curated sequence input when:

- You already know the exact source sequences.
- You do not need an assembly search to establish them.
- You are willing to accept responsibility for their curation.

Use the designs system when:

- You have a trusted target assembly.
- You know which sequence classes you want.
- Representative sequences can locate those classes.
- You want the assembly, rather than the external reference accession, to supply the final source sequence.

Before trusting a target assembly, record its accession or source, release, checksum, and taxonomic assignment. Review known contamination or assembly-quality warnings. An accepted assembly match can show that a region resembles a representative sequence, but it cannot prove that the assembly record has correct provenance.

When writing design rules:

1. Give each required sequence class a clear group.
2. Add alternative representative sequences to one group when they serve the same biological requirement.
3. Crop representative sequences to the region that should drive locus discovery.
4. Choose identity and coverage thresholds that reject partial or remote homologues without excluding expected target variation.
5. Inspect every passing candidate locus rather than assuming that BLAST selected one correct answer.

When constructing the background:

1. State the background-cancellation question before collecting FASTA files.
2. Record each source and release used to build the combined background.
3. Expect conserved loci such as rRNA to lose many candidate k-mers.
4. Check whether target-derived records in the background express an intentional cancellation rule.
5. Do not describe the background as exhaustive evidence about all non-target organisms.

## Common mistakes

**Using representative sequences as if they were source sequences.** In the designs system, the assembly intervals are the provisional source sequences. The representative sequences only locate them.

**Enumerating the whole target assembly.** This changes a targeted locus design into a genome-wide design and can create a much larger, less interpretable candidate set.

**Using the target assembly as the background.** This usually cancels the target k-mers that the design just generated.

**Treating `query_end = 0` as base zero.** It means “through the end of this FASTA record.” Other coordinates are 1-based and inclusive.

**Assuming one hit per representative sequence.** Every passing interval is retained. Repeats and paralogues can therefore produce several provisional source sequences.

**Assuming `target_taxid` verifies the design.** It does not validate the target assembly or representative sequences. It is used later for taxonomic interpretation.

**Treating the background as universal specificity evidence.** A bait survives only the evidence sources and rules enabled in its run. It can still match sequences absent from that background.

Start from [`assets/example_designs.csv`](../assets/example_designs.csv), then replace every placeholder path with a file that expresses the intended scientific role. The `query_rules` schema is available at [`assets/schema_query_rules.json`](../assets/schema_query_rules.json).
