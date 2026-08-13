include { MERYL_COUNT as MERYL_COUNT_SOURCE } from '../../../modules/nf-core/meryl/count/main'
include { MERYL_COUNT as MERYL_COUNT_BACKGROUND } from '../../../modules/nf-core/meryl/count/main'
include { MERYL_INTERSECT } from '../../../modules/local/meryl_intersect/main'
include { MERYL_PRINT as MERYL_PRINT_BACKGROUND } from '../../../modules/local/meryl_print/main'
include { MERYL_PRINT as MERYL_PRINT_SOURCE } from '../../../modules/local/meryl_print/main'
include { BUILD_CANDIDATE_TABLES } from '../../../modules/local/build_candidate_tables/main'
include { DEACON_INDEX_BUILD as DEACON_INDEX_ENTROPY } from '../../../modules/local/deacon_index_build/main'
include { DEACON_INDEX_DUMP } from '../../../modules/local/deacon_index_dump/main'
include { APPLY_COMPLEXITY_FILTER } from '../../../modules/local/apply_complexity_filter/main'

workflow FILTER_CANDIDATE_KMERS {
    take:
    ch_filtering_inputs
    ch_kmer_size
    ch_deacon_window
    ch_entropy_threshold

    main:

    // Count source and background k-mers
    ch_source_count_inputs = ch_filtering_inputs.map { meta, source_sequences, source_sequence_query_groups, interference_background ->
        tuple(meta, source_sequences)
    }
    ch_background_count_inputs = ch_filtering_inputs.map { meta, source_sequences, source_sequence_query_groups, interference_background ->
        tuple(meta, interference_background)
    }
    ch_source_context = ch_filtering_inputs.map { meta, source_sequences, source_sequence_query_groups, interference_background ->
        tuple(meta, source_sequences, source_sequence_query_groups)
    }
    MERYL_COUNT_SOURCE(ch_source_count_inputs, ch_kmer_size)
    MERYL_COUNT_BACKGROUND(ch_background_count_inputs, ch_kmer_size)

    // Isolate source k-mers absent from the background
    ch_count_pairs = MERYL_COUNT_SOURCE.out.meryl_db
        .join(MERYL_COUNT_BACKGROUND.out.meryl_db)
        .map { meta, source_db, background_db -> tuple(meta, source_db, background_db) }
    MERYL_INTERSECT(ch_count_pairs.map { meta, source_db, background_db -> tuple(meta, background_db, source_db) })
    MERYL_PRINT_BACKGROUND(MERYL_INTERSECT.out.db)
    MERYL_PRINT_SOURCE(MERYL_COUNT_SOURCE.out.meryl_db)

    // Build candidate evidence and apply complexity filtering
    ch_candidate_inputs = ch_source_context
        .join(MERYL_PRINT_BACKGROUND.out.txt)
        .join(MERYL_PRINT_SOURCE.out.txt)
        .combine(ch_kmer_size)
        .map { meta, source_sequences, source_sequence_query_groups, background_intersection_counts, meryl_source_counts, kmer_size ->
            tuple(meta, source_sequences, source_sequence_query_groups, meryl_source_counts, background_intersection_counts, kmer_size)
        }
    BUILD_CANDIDATE_TABLES(ch_candidate_inputs)
    ch_entropy_inputs = BUILD_CANDIDATE_TABLES.out.complexity_candidates
        .combine(ch_kmer_size)
        .combine(ch_deacon_window)
        .combine(ch_entropy_threshold)
        .map { meta, complexity_candidates, kmer_size, deacon_window, entropy_threshold ->
            tuple(meta, complexity_candidates, kmer_size, deacon_window, entropy_threshold)
        }
    DEACON_INDEX_ENTROPY(ch_entropy_inputs)
    DEACON_INDEX_DUMP(DEACON_INDEX_ENTROPY.out.index)
    ch_complexity_inputs = BUILD_CANDIDATE_TABLES.out.tables
        .join(DEACON_INDEX_DUMP.out.fasta)
        .map { meta, manifest, occurrences, passing_kmers -> tuple(meta, meta.source_sequence_origin, manifest, passing_kmers) }
    APPLY_COMPLEXITY_FILTER(ch_complexity_inputs)

    emit:
    baits = APPLY_COMPLEXITY_FILTER.out.baits
    manifest = BUILD_CANDIDATE_TABLES.out.terminal_manifest
        .mix(APPLY_COMPLEXITY_FILTER.out.evidence.map { meta, manifest, filtering_status -> tuple(meta, manifest) })
    terminal_manifest = BUILD_CANDIDATE_TABLES.out.terminal_manifest
        .mix(APPLY_COMPLEXITY_FILTER.out.terminal_manifest)
    occurrences = BUILD_CANDIDATE_TABLES.out.tables.map { meta, manifest, occurrences -> tuple(meta, occurrences) }
    filtering_status = BUILD_CANDIDATE_TABLES.out.terminal_status
        .mix(APPLY_COMPLEXITY_FILTER.out.evidence.map { meta, manifest, filtering_status -> tuple(meta, filtering_status) })
    reported_meryl = MERYL_INTERSECT.out.reported_meryl
    reported_deacon = DEACON_INDEX_ENTROPY.out.reported_deacon
    reported_biopython = BUILD_CANDIDATE_TABLES.out.reported_biopython
    reported_polars = BUILD_CANDIDATE_TABLES.out.reported_polars
    terminal_without_deacon = BUILD_CANDIDATE_TABLES.out.terminal_status
    after_deacon = APPLY_COMPLEXITY_FILTER.out.evidence
    versions = MERYL_COUNT_SOURCE.out.versions_meryl
        .mix(MERYL_COUNT_BACKGROUND.out.versions_meryl)
        .mix(MERYL_INTERSECT.out.versions_meryl)
        .mix(MERYL_PRINT_BACKGROUND.out.versions_meryl)
        .mix(MERYL_PRINT_SOURCE.out.versions_meryl)
        .mix(DEACON_INDEX_ENTROPY.out.versions_deacon)
        .mix(DEACON_INDEX_DUMP.out.versions_deacon)
        .mix(BUILD_CANDIDATE_TABLES.out.versions_biopython)
        .mix(BUILD_CANDIDATE_TABLES.out.versions_polars)
        .mix(APPLY_COMPLEXITY_FILTER.out.versions_biopython)
        .mix(APPLY_COMPLEXITY_FILTER.out.versions_polars)
}
