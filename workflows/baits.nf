include { NORMALIZE_CURATED_SOURCE_SEQUENCES } from '../modules/local/normalize_curated_source_sequences/main'
include { WRITE_PROVENANCE } from '../modules/local/write_provenance/main'
include { DISCOVER_SOURCE_SEQUENCES } from '../subworkflows/local/discover_source_sequences/main'
include { FILTER_CANDIDATE_KMERS } from '../subworkflows/local/filter_candidate_kmers/main'
include { TAXONOMIC_EXACT_MATCH_SCREENING } from '../subworkflows/local/taxonomic_exact_match_screening/main'
include { BUILD_VERIFY_DEACON_INDEX } from '../subworkflows/local/build_verify_deacon_index/main'

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
    ch_taxonomic_reference_db
    ch_taxonomic_screening_not_run

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

    // Screen the Locally Filtered Bait Set when a Taxonomic Reference Database is supplied
    ch_screening_inputs = FILTER_CANDIDATE_KMERS.out.baits
        .join(ch_design_context)
        .map { meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid, interference_background ->
            tuple(meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid)
        }
    TAXONOMIC_EXACT_MATCH_SCREENING(ch_screening_inputs, ch_taxonomic_reference_db, ch_kmer_size)

    // Build and verify one Deacon Index from the deepest justified Bait Set
    ch_locally_filtered_index_inputs = FILTER_CANDIDATE_KMERS.out.baits
        .join(ch_design_context)
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid, interference_background, not_run ->
            tuple(meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, interference_background)
        }
    ch_taxonomically_screened_index_inputs = TAXONOMIC_EXACT_MATCH_SCREENING.out.baits
        .join(ch_design_context)
        .map { meta, taxonomically_screened_baits, candidate_kmer_manifest, bait_set_status, target_taxid, interference_background ->
            tuple(meta, taxonomically_screened_baits, candidate_kmer_manifest, bait_set_status, interference_background)
        }
    // These keyed markers restore the Bait Set branch after index verification.
    ch_locally_filtered_index_keys = ch_locally_filtered_index_inputs.map { meta, baits, candidate_kmer_manifest, bait_set_status, interference_background -> tuple(meta, true) }
    ch_taxonomically_screened_index_keys = ch_taxonomically_screened_index_inputs.map { meta, baits, candidate_kmer_manifest, bait_set_status, interference_background -> tuple(meta, true) }
    ch_index_inputs = ch_locally_filtered_index_inputs.mix(ch_taxonomically_screened_index_inputs)
    BUILD_VERIFY_DEACON_INDEX(ch_index_inputs, ch_kmer_size, ch_deacon_window)

    ch_locally_filtered_indexes = BUILD_VERIFY_DEACON_INDEX.out.index
        .join(ch_locally_filtered_index_keys)
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index, selected -> tuple(meta, deacon_index) }
    ch_taxonomically_screened_indexes = BUILD_VERIFY_DEACON_INDEX.out.index
        .join(ch_taxonomically_screened_index_keys)
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index, selected -> tuple(meta, deacon_index) }

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

    // Finalize designs according to whether Taxonomic Exact-Match Screening ran
    ch_filtering_terminal_provenance_facts = ch_without_deacon_provenance_facts
        .mix(
            ch_with_deacon_provenance_facts
                .join(FILTER_CANDIDATE_KMERS.out.terminal_manifest)
                .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, manifest ->
                    tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions)
                },
        )
    ch_filtering_continuation_provenance_facts = ch_with_deacon_provenance_facts
        .join(FILTER_CANDIDATE_KMERS.out.baits)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, baits, manifest, bait_set_status ->
            tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions)
        }

    ch_screening_not_requested_provenance_facts = ch_filtering_terminal_provenance_facts
        .mix(ch_filtering_continuation_provenance_facts)
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, not_run ->
            tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions)
        }
    ch_filtering_terminal_with_reference_facts = ch_filtering_terminal_provenance_facts
        .combine(ch_taxonomic_reference_db)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, taxonomic_reference_db ->
            tuple(
                meta,
                input_facts,
                input_file_roles + ['taxonomic_reference_database'],
                input_file_ids + ['taxonomic_reference_database'],
                input_file_kinds + ['directory'],
                input_files + [taxonomic_reference_db],
                parameters,
                software_versions,
            )
        }

    ch_screened_provenance_facts = ch_filtering_continuation_provenance_facts
        .join(TAXONOMIC_EXACT_MATCH_SCREENING.out.screening_status)
        .join(TAXONOMIC_EXACT_MATCH_SCREENING.out.reported_blast)
        .join(TAXONOMIC_EXACT_MATCH_SCREENING.out.reported_biopython)
        .join(TAXONOMIC_EXACT_MATCH_SCREENING.out.reported_polars)
        .combine(ch_taxonomic_reference_db)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, screening_status, blast_component, blast_version, biopython_component, biopython_version, polars_component, polars_version, taxonomic_reference_db ->
            tuple(
                meta,
                input_facts,
                input_file_roles + ['taxonomic_reference_database'],
                input_file_ids + ['taxonomic_reference_database'],
                input_file_kinds + ['directory'],
                input_files + [taxonomic_reference_db],
                parameters,
                software_versions + [
                    [component: blast_component, version: blast_version],
                    [component: biopython_component, version: biopython_version],
                    [component: polars_component, version: polars_version],
                ],
            )
        }

    WRITE_PROVENANCE(
        ch_discovery_terminal_provenance_facts
            .mix(ch_screening_not_requested_provenance_facts)
            .mix(ch_filtering_terminal_with_reference_facts)
            .mix(ch_screened_provenance_facts),
    )

    emit:
    source_sequences = ch_source_sequences
    candidate_loci = DISCOVER_SOURCE_SEQUENCES.out.candidate_loci
    query_blast_hits = DISCOVER_SOURCE_SEQUENCES.out.blast_hits
    discovery_status = DISCOVER_SOURCE_SEQUENCES.out.discovery_status
    provenance = WRITE_PROVENANCE.out.provenance
    candidate_kmers = FILTER_CANDIDATE_KMERS.out.manifest
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, manifest, not_run -> tuple(meta, manifest) }
        .mix(FILTER_CANDIDATE_KMERS.out.terminal_manifest.combine(ch_taxonomic_reference_db).map { meta, manifest, taxonomic_reference_db -> tuple(meta, manifest) })
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.manifest)
    candidate_kmer_occurrences = FILTER_CANDIDATE_KMERS.out.occurrences
    filtering_status = FILTER_CANDIDATE_KMERS.out.filtering_status
    locally_filtered_baits = FILTER_CANDIDATE_KMERS.out.baits
    taxonomic_blast_hits = TAXONOMIC_EXACT_MATCH_SCREENING.out.hits
    screening_decisions = TAXONOMIC_EXACT_MATCH_SCREENING.out.decisions
    screening_status = TAXONOMIC_EXACT_MATCH_SCREENING.out.screening_status
    taxonomically_screened_baits = TAXONOMIC_EXACT_MATCH_SCREENING.out.baits
    bait_set_status = BUILD_VERIFY_DEACON_INDEX.out.bait_set_status
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.terminal_bait_set_status)
    locally_filtered_deacon_index = ch_locally_filtered_indexes
    taxonomically_screened_deacon_index = ch_taxonomically_screened_indexes
    index_verification_summary = BUILD_VERIFY_DEACON_INDEX.out.summary
    index_verification_report = BUILD_VERIFY_DEACON_INDEX.out.report
    taxonomic_reference_database = TAXONOMIC_EXACT_MATCH_SCREENING.out.reference_database_report
    versions = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.versions_biopython
        .mix(DISCOVER_SOURCE_SEQUENCES.out.versions)
        .mix(WRITE_PROVENANCE.out.versions_python)
        .mix(FILTER_CANDIDATE_KMERS.out.versions)
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.versions)
        .mix(BUILD_VERIFY_DEACON_INDEX.out.versions)
}
