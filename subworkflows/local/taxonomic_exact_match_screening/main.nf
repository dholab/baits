include { BLASTN_TAXONOMIC_SCREENING } from '../../../modules/local/blastn_taxonomic_screening/main'
include { APPLY_TAXONOMIC_SCREENING } from '../../../modules/local/apply_taxonomic_screening/main'

workflow TAXONOMIC_EXACT_MATCH_SCREENING {
    take:
    ch_screening_inputs
    ch_taxonomic_reference_db
    ch_kmer_size

    main:

    // Search locally filtered Baits against the shared Taxonomic Reference Database
    ch_blast_inputs = ch_screening_inputs
        .combine(ch_taxonomic_reference_db)
        .combine(ch_kmer_size)
        .map { meta, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid, taxonomic_reference_db, kmer_size ->
            tuple(meta, locally_filtered_baits, taxonomic_reference_db, kmer_size)
        }
    BLASTN_TAXONOMIC_SCREENING(ch_blast_inputs)

    // Apply the Allowed Taxonomic Scope
    ch_screening_interpreter_inputs = BLASTN_TAXONOMIC_SCREENING.out.hits
        .join(ch_screening_inputs)
        .combine(ch_kmer_size)
        .map { meta, blast_hits, locally_filtered_baits, candidate_kmer_manifest, bait_set_status, target_taxid, kmer_size ->
            tuple(meta, locally_filtered_baits, blast_hits, candidate_kmer_manifest, bait_set_status, target_taxid, kmer_size)
        }
    APPLY_TAXONOMIC_SCREENING(ch_screening_interpreter_inputs)

    emit:
    hits = BLASTN_TAXONOMIC_SCREENING.out.hits
    reference_database_report = BLASTN_TAXONOMIC_SCREENING.out.reference_database_report
    decisions = APPLY_TAXONOMIC_SCREENING.out.decisions
    manifest = APPLY_TAXONOMIC_SCREENING.out.manifest
    screening_status = APPLY_TAXONOMIC_SCREENING.out.screening_status
    bait_set_status = APPLY_TAXONOMIC_SCREENING.out.bait_set_status
    terminal_bait_set_status = APPLY_TAXONOMIC_SCREENING.out.terminal_bait_set_status
    baits = APPLY_TAXONOMIC_SCREENING.out.baits
    reported_blast = BLASTN_TAXONOMIC_SCREENING.out.reported_blast
    reported_biopython = APPLY_TAXONOMIC_SCREENING.out.reported_biopython
    reported_polars = APPLY_TAXONOMIC_SCREENING.out.reported_polars
    versions = BLASTN_TAXONOMIC_SCREENING.out.versions_blast
        .mix(APPLY_TAXONOMIC_SCREENING.out.versions_biopython)
        .mix(APPLY_TAXONOMIC_SCREENING.out.versions_polars)
}
