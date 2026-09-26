include { NORMALIZE_CURATED_SOURCE_SEQUENCES } from '../modules/local/normalize_curated_source_sequences/main'
include { DISCOVER_SOURCE_SEQUENCES } from '../subworkflows/local/discover_source_sequences/main'
include { REPORT_SOURCE_DISCOVERY_RESULTS } from '../subworkflows/local/report_source_discovery/main'
include { FILTER_CANDIDATE_KMERS } from '../subworkflows/local/filter_candidate_kmers/main'
include { TAXONOMIC_EXACT_MATCH_SCREENING } from '../subworkflows/local/taxonomic_exact_match_screening/main'
include { BUILD_VERIFY_DEACON_INDEX as BUILD_VERIFY_DEACON_INDEX_WITHOUT_TAXONOMY } from '../subworkflows/local/build_verify_deacon_index/main'
include { BUILD_VERIFY_DEACON_INDEX as BUILD_VERIFY_DEACON_INDEX_WITH_TAXONOMY } from '../subworkflows/local/build_verify_deacon_index/main'
include { CALIBRATE_THRESHOLD } from '../subworkflows/local/calibrate_threshold/main'

workflow BAITS_MAIN {
    take:
    ch_design_context
    ch_curated_inputs
    ch_discovery_inputs
    ch_kmer_size
    ch_deacon_window
    ch_entropy_threshold
    ch_taxonomic_reference_db
    ch_taxonomic_screening_not_run
    ch_taxonomic_screening_scopes
    ch_calibration_reads
    ch_calibration_target_scopes

    main:

    // Source sequence acquisition
    NORMALIZE_CURATED_SOURCE_SEQUENCES(ch_curated_inputs)
    DISCOVER_SOURCE_SEQUENCES(ch_discovery_inputs)
    REPORT_SOURCE_DISCOVERY_RESULTS(
        DISCOVER_SOURCE_SEQUENCES.out.candidate_loci,
        DISCOVER_SOURCE_SEQUENCES.out.blast_hits,
    )

    // Converge curated and discovered source sequences
    ch_source_sequences = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.source_sequences
        .mix(DISCOVER_SOURCE_SEQUENCES.out.source_sequences)

    // Filter candidate k-mers
    ch_filtering_inputs = ch_source_sequences
        .join(ch_design_context)
        .map { meta, source_sequences, source_sequence_query_groups, target_taxid, background ->
            tuple(meta, source_sequences, source_sequence_query_groups, background)
        }
    FILTER_CANDIDATE_KMERS(ch_filtering_inputs, ch_kmer_size, ch_deacon_window, ch_entropy_threshold)

    // Screen selected candidates when a taxonomic reference database is supplied
    ch_screening_inputs = FILTER_CANDIDATE_KMERS.out.candidates
        .join(ch_design_context)
        .combine(ch_taxonomic_screening_scopes, by: 0)
        .map { meta, candidates, candidate_kmer_manifest, target_taxid, background, taxonomic_screening_taxids ->
            tuple(meta, candidates, candidate_kmer_manifest, target_taxid, taxonomic_screening_taxids)
        }
    TAXONOMIC_EXACT_MATCH_SCREENING(ch_screening_inputs, ch_taxonomic_reference_db, ch_kmer_size)

    // Build and verify one Deacon index from the selected candidates
    ch_index_inputs_without_taxonomy = FILTER_CANDIDATE_KMERS.out.candidates
        .join(ch_design_context)
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, candidates, candidate_kmer_manifest, target_taxid, background, not_run ->
            tuple(meta, candidates, candidate_kmer_manifest, background)
        }
    ch_index_inputs_with_taxonomy = TAXONOMIC_EXACT_MATCH_SCREENING.out.candidates
        .join(ch_design_context)
        .map { meta, candidates, candidate_kmer_manifest, target_taxid, background ->
            tuple(meta, candidates, candidate_kmer_manifest, background)
        }
    BUILD_VERIFY_DEACON_INDEX_WITHOUT_TAXONOMY(ch_index_inputs_without_taxonomy, ch_kmer_size, ch_deacon_window)
    BUILD_VERIFY_DEACON_INDEX_WITH_TAXONOMY(ch_index_inputs_with_taxonomy, ch_kmer_size, ch_deacon_window)

    ch_indexes_without_taxonomy = BUILD_VERIFY_DEACON_INDEX_WITHOUT_TAXONOMY.out.index
        .map { meta, baits, candidate_kmer_manifest, deacon_index -> tuple(meta, deacon_index) }
    ch_verified_indexes_with_taxonomy = BUILD_VERIFY_DEACON_INDEX_WITH_TAXONOMY.out.index
    ch_indexes_with_taxonomy = ch_verified_indexes_with_taxonomy
        .map { meta, baits, candidate_kmer_manifest, deacon_index -> tuple(meta, deacon_index) }

    // Calibrate only taxonomically screened, verified Deacon Indexes
    ch_calibration_inputs = ch_calibration_reads
        .combine(ch_verified_indexes_with_taxonomy, by: 0)
        .combine(ch_design_context, by: 0)
        .combine(ch_calibration_target_scopes, by: 0)
        .map { meta, read_meta, reads, baits, candidate_kmer_manifest, deacon_index, target_taxid, background, calibration_target_taxids ->
            tuple(meta, read_meta, reads, baits, deacon_index, target_taxid, calibration_target_taxids)
        }
    CALIBRATE_THRESHOLD(ch_calibration_inputs, ch_taxonomic_reference_db, ch_kmer_size)

    emit:
    source_sequences = ch_source_sequences
    candidate_loci = REPORT_SOURCE_DISCOVERY_RESULTS.out.candidate_loci
    source_matches = REPORT_SOURCE_DISCOVERY_RESULTS.out.source_matches
    candidate_kmers = FILTER_CANDIDATE_KMERS.out.manifest
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, manifest, not_run -> tuple(meta, manifest) }
        .mix(FILTER_CANDIDATE_KMERS.out.manifest_without_candidates.combine(ch_taxonomic_reference_db).map { meta, manifest, taxonomic_reference_db -> tuple(meta, manifest) })
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.manifest)
    candidate_kmer_occurrences = FILTER_CANDIDATE_KMERS.out.occurrences
    taxonomic_blast_hits = TAXONOMIC_EXACT_MATCH_SCREENING.out.hits
    taxonomic_candidate_kmers = TAXONOMIC_EXACT_MATCH_SCREENING.out.manifest
    baits = FILTER_CANDIDATE_KMERS.out.candidates
        .combine(ch_taxonomic_screening_not_run)
        .map { meta, candidates, manifest, not_run -> tuple(meta, candidates) }
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.candidates.map { meta, candidates, manifest -> tuple(meta, candidates) })
    deacon_index = ch_indexes_without_taxonomy.mix(ch_indexes_with_taxonomy)
    calibration_reads = CALIBRATE_THRESHOLD.out.calibration_reads
    calibration_matches = CALIBRATE_THRESHOLD.out.calibration_matches
    recommended_threshold = CALIBRATE_THRESHOLD.out.recommended_threshold
    versions = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.versions_biopython
        .mix(DISCOVER_SOURCE_SEQUENCES.out.versions)
        .mix(REPORT_SOURCE_DISCOVERY_RESULTS.out.versions)
        .mix(FILTER_CANDIDATE_KMERS.out.versions)
        .mix(TAXONOMIC_EXACT_MATCH_SCREENING.out.versions)
        .mix(BUILD_VERIFY_DEACON_INDEX_WITHOUT_TAXONOMY.out.versions)
        .mix(BUILD_VERIFY_DEACON_INDEX_WITH_TAXONOMY.out.versions)
        .mix(CALIBRATE_THRESHOLD.out.versions)
}
