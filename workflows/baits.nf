include { NORMALIZE_CURATED_SOURCE_SEQUENCES } from '../modules/local/normalize_curated_source_sequences/main'
include { WRITE_PROVENANCE } from '../modules/local/write_provenance/main'
include { DISCOVER_SOURCE_SEQUENCES } from '../subworkflows/local/discover_source_sequences/main'
include { FILTER_CANDIDATE_KMERS } from '../subworkflows/local/filter_candidate_kmers/main'
include { TAXONOMIC_EXACT_MATCH_SCREENING } from '../subworkflows/local/taxonomic_exact_match_screening/main'
include { BUILD_VERIFY_DEACON_INDEX } from '../subworkflows/local/build_verify_deacon_index/main'
include { CALIBRATE_THRESHOLD } from '../subworkflows/local/calibrate_threshold/main'

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
    ch_calibration_reads
    ch_calibration_target_scopes
    ch_calibration_provenance_facts
    ch_without_calibration_keys

    main:

    // Source sequence acquisition
    NORMALIZE_CURATED_SOURCE_SEQUENCES(ch_curated_inputs)
    DISCOVER_SOURCE_SEQUENCES(ch_discovery_inputs)

    // Converge curated and discovered source sequences
    ch_source_sequences = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.source_sequences
        .mix(DISCOVER_SOURCE_SEQUENCES.out.source_sequences)

    // Filter candidate k-mers
    ch_filtering_inputs = ch_source_sequences
        .join(ch_design_context)
        .map { meta, source_sequences, source_sequence_query_groups, target_taxid, interference_background ->
            tuple(meta, source_sequences, source_sequence_query_groups, interference_background)
        }
    FILTER_CANDIDATE_KMERS(ch_filtering_inputs, ch_kmer_size, ch_deacon_window, ch_entropy_threshold)

    // Screen the locally filtered bait set when a taxonomic reference database is supplied
    ch_screening_inputs = FILTER_CANDIDATE_KMERS.out.baits
        .join(ch_design_context)
        .map { meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid, interference_background ->
            tuple(meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid)
        }
    TAXONOMIC_EXACT_MATCH_SCREENING(ch_screening_inputs, ch_taxonomic_reference_db, ch_kmer_size)

    // Build and verify one Deacon index from the deepest justified nonempty bait set
    ch_locally_filtered_index_inputs = FILTER_CANDIDATE_KMERS.out.baits
        .join(ch_design_context)
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid, interference_background, not_run ->
            tuple(meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, interference_background, 'locally_filtered')
        }
    ch_taxonomically_screened_index_inputs = TAXONOMIC_EXACT_MATCH_SCREENING.out.baits
        .join(ch_design_context)
        .map { meta, taxonomically_screened_baits, candidate_kmer_manifest, bait_set_status, target_taxid, interference_background ->
            tuple(meta, taxonomically_screened_baits, candidate_kmer_manifest, bait_set_status, interference_background, 'taxonomically_screened')
        }
    ch_index_inputs = ch_locally_filtered_index_inputs.mix(ch_taxonomically_screened_index_inputs)
    BUILD_VERIFY_DEACON_INDEX(ch_index_inputs, ch_kmer_size, ch_deacon_window)

    ch_locally_filtered_indexes = BUILD_VERIFY_DEACON_INDEX.out.index
        .filter { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index, index_source -> index_source == 'locally_filtered' }
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index, index_source -> tuple(meta, deacon_index) }
    ch_taxonomically_screened_verified_indexes = BUILD_VERIFY_DEACON_INDEX.out.index
        .filter { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index, index_source -> index_source == 'taxonomically_screened' }
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index, index_source ->
            tuple(meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index)
        }
    ch_taxonomically_screened_indexes = ch_taxonomically_screened_verified_indexes
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index -> tuple(meta, deacon_index) }

    // Calibrate only taxonomically screened, verified Deacon Indexes
    ch_calibration_inputs = ch_calibration_reads
        .combine(ch_taxonomically_screened_verified_indexes, by: 0)
        .combine(ch_design_context, by: 0)
        .combine(ch_calibration_target_scopes, by: 0)
        .map { meta, read_meta, reads, taxonomically_screened_baits, candidate_kmer_manifest, bait_set_status, deacon_index, target_taxid, interference_background, calibration_target_taxids ->
            tuple(meta, read_meta, reads, taxonomically_screened_baits, bait_set_status, deacon_index, target_taxid, calibration_target_taxids)
        }
    CALIBRATE_THRESHOLD(ch_calibration_inputs, ch_taxonomic_reference_db, ch_kmer_size)

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

    // Finalize designs according to whether taxonomic exact-match screening ran
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

    ch_provenance_facts = ch_discovery_terminal_provenance_facts
        .mix(ch_screening_not_requested_provenance_facts)
        .mix(ch_filtering_terminal_with_reference_facts)
        .mix(ch_screened_provenance_facts)

    ch_without_calibration_provenance_facts = ch_provenance_facts
        .join(ch_without_calibration_keys)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, without_calibration ->
            tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions)
        }
    ch_with_calibration_provenance_facts = ch_provenance_facts
        .join(ch_calibration_provenance_facts)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, calibration_input_facts, calibration_file_roles, calibration_file_ids, calibration_file_kinds, calibration_files, calibration_versions ->
            tuple(
                meta,
                input_facts + calibration_input_facts,
                input_file_roles + calibration_file_roles,
                input_file_ids + calibration_file_ids,
                input_file_kinds + calibration_file_kinds,
                input_files + calibration_files,
                parameters,
                software_versions + calibration_versions,
            )
        }

    ch_calibrated_design_keys = CALIBRATE_THRESHOLD.out.summary
        .map { meta, summary -> tuple(meta, true) }

    ch_final_calibration_provenance_facts = ch_with_calibration_provenance_facts
        .join(ch_calibrated_design_keys, remainder: true)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, parameters, software_versions, calibration_ran ->
            def calibrationParameters = calibration_ran ? [
                [parameter: 'deacon_filter_absolute_threshold', value: '1'],
                [parameter: 'deacon_filter_relative_threshold', value: '0'],
                [parameter: 'read_blast_dust', value: 'no'],
                [parameter: 'read_blast_evalue', value: '1e-10'],
                [parameter: 'max_blast_targets', value: params.max_blast_targets.toString()],
                [parameter: 'read_blast_task', value: 'blastn'],
                [parameter: 'read_classification_tie_tolerance', value: '0.1'],
            ] : []
            tuple(
                meta,
                input_facts,
                input_file_roles,
                input_file_ids,
                input_file_kinds,
                input_files,
                parameters + calibrationParameters,
                software_versions,
            )
        }

    WRITE_PROVENANCE(
        ch_without_calibration_provenance_facts
            .mix(ch_final_calibration_provenance_facts),
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
    candidate_read_counts = CALIBRATE_THRESHOLD.out.candidate_read_counts
    read_blast_hits = CALIBRATE_THRESHOLD.out.read_blast_hits
    read_blast_search_parameters = CALIBRATE_THRESHOLD.out.read_blast_search_parameters
    classified_reads = CALIBRATE_THRESHOLD.out.classified_reads
    threshold_read_counts = CALIBRATE_THRESHOLD.out.threshold_read_counts
    threshold_summary = CALIBRATE_THRESHOLD.out.summary
    versions = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.versions_biopython
        .mix(DISCOVER_SOURCE_SEQUENCES.out.versions)
        .mix(WRITE_PROVENANCE.out.versions_python)
        .mix(FILTER_CANDIDATE_KMERS.out.versions)
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.versions)
        .mix(BUILD_VERIFY_DEACON_INDEX.out.versions)
        .mix(CALIBRATE_THRESHOLD.out.versions)
}
