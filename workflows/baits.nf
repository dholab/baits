include { NORMALIZE_CURATED_SOURCE_SEQUENCES } from '../modules/local/normalize_curated_source_sequences/main'
include { WRITE_PROVENANCE } from '../modules/local/write_provenance/main'
include { DISCOVER_SOURCE_SEQUENCES } from '../subworkflows/local/discover_source_sequences/main'
include { FILTER_CANDIDATE_KMERS } from '../subworkflows/local/filter_candidate_kmers/main'

workflow BAITS_MAIN {
    take:
    ch_design_context
    ch_curated_inputs
    ch_discovery_inputs
    ch_curated_provenance_source_inputs
    ch_discovery_provenance_source_inputs
    ch_kmer_size
    ch_deacon_window
    ch_entropy_threshold

    main:

    // Source Sequence acquisition
    NORMALIZE_CURATED_SOURCE_SEQUENCES(ch_curated_inputs)
    DISCOVER_SOURCE_SEQUENCES(ch_discovery_inputs)

    // Converge curated and discovered Source Sequences
    ch_source_sequences = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.source_sequences
        .mix(DISCOVER_SOURCE_SEQUENCES.out.source_sequences)

    // Filter Candidate K-mers
    ch_filtering_inputs = ch_source_sequences
        .join(ch_design_context)
        .map { meta, source_sequences, source_sequence_query_groups, target_taxid, interference_background ->
            tuple(meta, source_sequences, source_sequence_query_groups, interference_background)
        }
    FILTER_CANDIDATE_KMERS(ch_filtering_inputs, ch_kmer_size, ch_deacon_window, ch_entropy_threshold)

    // Construct immutable source-input and design provenance facts
    ch_curated_provenance_facts = ch_design_context
        .join(ch_curated_provenance_source_inputs)
        .join(NORMALIZE_CURATED_SOURCE_SEQUENCES.out.reported_biopython)
        .map { meta, target_taxid, interference_background, source_input_roles, source_input_files, component, version ->
            tuple(
                meta,
                [
                    [input_role: 'design', input_id: meta.id, attribute: 'source_sequence_origin', value: meta.source_sequence_origin],
                    [input_role: 'target_taxon', input_id: meta.id, attribute: 'target_taxid', value: target_taxid],
                ],
                ['interference_background'] + source_input_roles,
                ['interference_background'] + source_input_roles,
                ['file'] * (source_input_files.size() + 1),
                [interference_background] + source_input_files,
                [[parameter: 'target_taxid', value: target_taxid]],
                [[component: component, version: version]],
            )
        }

    ch_discovery_versions = DISCOVER_SOURCE_SEQUENCES.out.reported_biopython
        .join(DISCOVER_SOURCE_SEQUENCES.out.reported_blast)
    ch_discovery_provenance_facts = ch_design_context
        .join(ch_discovery_provenance_source_inputs)
        .join(ch_discovery_versions)
        .map { meta, target_taxid, interference_background, source_input_roles, source_input_files, biopython_component, biopython_version, blast_component, blast_version ->
            tuple(
                meta,
                [
                    [input_role: 'design', input_id: meta.id, attribute: 'source_sequence_origin', value: meta.source_sequence_origin],
                    [input_role: 'target_taxon', input_id: meta.id, attribute: 'target_taxid', value: target_taxid],
                ],
                ['interference_background'] + source_input_roles,
                ['interference_background'] + source_input_roles,
                ['file'] * (source_input_files.size() + 1),
                [interference_background] + source_input_files,
                [[parameter: 'target_taxid', value: target_taxid]],
                [
                    [component: biopython_component, version: biopython_version],
                    [component: blast_component, version: blast_version],
                ],
            )
        }

    ch_discovery_terminal_provenance_facts = ch_discovery_provenance_facts
        .join(DISCOVER_SOURCE_SEQUENCES.out.discovery_terminal)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, discovery_terminal ->
            tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions)
        }
    ch_filtering_provenance_facts = ch_curated_provenance_facts
        .mix(
            ch_discovery_provenance_facts
                .join(DISCOVER_SOURCE_SEQUENCES.out.source_sequences)
                .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, source_sequences, source_sequence_query_groups ->
                    tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions)
                },
        )

    // Add facts from the filtering path that actually ran
    ch_without_deacon_provenance_facts = ch_filtering_provenance_facts
        .join(FILTER_CANDIDATE_KMERS.out.terminal_without_deacon)
        .join(FILTER_CANDIDATE_KMERS.out.reported_meryl)
        .join(FILTER_CANDIDATE_KMERS.out.reported_biopython)
        .join(FILTER_CANDIDATE_KMERS.out.reported_polars)
        .combine(ch_kmer_size)
        .combine(ch_deacon_window)
        .combine(ch_entropy_threshold)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, filtering_status, meryl_component, meryl_version, biopython_component, biopython_version, polars_component, polars_version, kmer_size, deacon_window, entropy_threshold ->
            tuple(
                meta,
                input_facts,
                input_file_roles,
                input_file_ids,
                input_file_kinds,
                input_files,
                parameters + [
                    [parameter: 'kmer_size', value: kmer_size.toString()],
                    [parameter: 'deacon_window', value: deacon_window.toString()],
                    [parameter: 'entropy_threshold', value: entropy_threshold.toString()],
                ],
                software_versions + [
                    [component: meryl_component, version: meryl_version],
                    [component: biopython_component, version: biopython_version],
                    [component: polars_component, version: polars_version],
                ],
            )
        }
    ch_with_deacon_provenance_facts = ch_filtering_provenance_facts
        .join(FILTER_CANDIDATE_KMERS.out.after_deacon)
        .join(FILTER_CANDIDATE_KMERS.out.reported_meryl)
        .join(FILTER_CANDIDATE_KMERS.out.reported_deacon)
        .join(FILTER_CANDIDATE_KMERS.out.reported_biopython)
        .join(FILTER_CANDIDATE_KMERS.out.reported_polars)
        .combine(ch_kmer_size)
        .combine(ch_deacon_window)
        .combine(ch_entropy_threshold)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, manifest, filtering_status, meryl_component, meryl_version, deacon_component, deacon_version, biopython_component, biopython_version, polars_component, polars_version, kmer_size, deacon_window, entropy_threshold ->
            tuple(
                meta,
                input_facts,
                input_file_roles,
                input_file_ids,
                input_file_kinds,
                input_files,
                parameters + [
                    [parameter: 'kmer_size', value: kmer_size.toString()],
                    [parameter: 'deacon_window', value: deacon_window.toString()],
                    [parameter: 'entropy_threshold', value: entropy_threshold.toString()],
                ],
                software_versions + [
                    [component: meryl_component, version: meryl_version],
                    [component: deacon_component, version: deacon_version],
                    [component: biopython_component, version: biopython_version],
                    [component: polars_component, version: polars_version],
                ],
            )
        }

    WRITE_PROVENANCE(
        ch_discovery_terminal_provenance_facts
            .mix(ch_without_deacon_provenance_facts)
            .mix(ch_with_deacon_provenance_facts),
    )

    emit:
    source_sequences = ch_source_sequences
    candidate_loci = DISCOVER_SOURCE_SEQUENCES.out.candidate_loci
    query_blast_hits = DISCOVER_SOURCE_SEQUENCES.out.blast_hits
    discovery_status = DISCOVER_SOURCE_SEQUENCES.out.discovery_status
    provenance = WRITE_PROVENANCE.out.provenance
    candidate_kmers = FILTER_CANDIDATE_KMERS.out.manifest
    candidate_kmer_occurrences = FILTER_CANDIDATE_KMERS.out.occurrences
    filtering_status = FILTER_CANDIDATE_KMERS.out.filtering_status
    locally_filtered_baits = FILTER_CANDIDATE_KMERS.out.baits
    versions = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.versions_biopython
        .mix(DISCOVER_SOURCE_SEQUENCES.out.versions)
        .mix(WRITE_PROVENANCE.out.versions_python)
        .mix(FILTER_CANDIDATE_KMERS.out.versions)
}
