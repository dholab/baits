include { BLASTN_TAXONOMIC_SCREENING } from '../../../modules/local/blastn_taxonomic_screening/main'
include { APPLY_TAXONOMIC_SCREENING } from '../../../modules/local/apply_taxonomic_screening/main'

workflow TAXONOMIC_EXACT_MATCH_SCREENING {
    take:
    ch_screening_inputs
    ch_taxonomic_reference_db
    ch_kmer_size

    main:

    // Search selected candidates against the shared taxonomic reference database
    ch_blast_inputs = ch_screening_inputs
        .combine(ch_taxonomic_reference_db)
        .combine(ch_kmer_size)
        .map { meta, baits, candidate_kmer_manifest, target_taxid, taxonomic_reference_db, kmer_size ->
            tuple(meta, baits, taxonomic_reference_db, kmer_size)
        }
    BLASTN_TAXONOMIC_SCREENING(ch_blast_inputs)

    // Apply the allowed taxonomic scope
    ch_screening_interpreter_inputs = BLASTN_TAXONOMIC_SCREENING.out.hits
        .join(ch_screening_inputs)
        .combine(ch_kmer_size)
        .map { meta, blast_hits, baits, candidate_kmer_manifest, target_taxid, kmer_size ->
            tuple(meta, baits, blast_hits, candidate_kmer_manifest, target_taxid, kmer_size)
        }
    APPLY_TAXONOMIC_SCREENING(ch_screening_interpreter_inputs)

    ch_screened_baits = APPLY_TAXONOMIC_SCREENING.out.manifest
        .join(APPLY_TAXONOMIC_SCREENING.out.baits)
        .map { meta, manifest, baits -> tuple(meta, baits, manifest) }

    emit:
    hits = BLASTN_TAXONOMIC_SCREENING.out.hits
    reference_database_report = BLASTN_TAXONOMIC_SCREENING.out.reference_database_report
    decisions = APPLY_TAXONOMIC_SCREENING.out.decisions
    manifest = APPLY_TAXONOMIC_SCREENING.out.manifest
    completed = APPLY_TAXONOMIC_SCREENING.out.manifest
    baits = ch_screened_baits
    reported_blast = BLASTN_TAXONOMIC_SCREENING.out.reported_blast
    reported_biopython = APPLY_TAXONOMIC_SCREENING.out.reported_biopython
    reported_polars = APPLY_TAXONOMIC_SCREENING.out.reported_polars
    versions = BLASTN_TAXONOMIC_SCREENING.out.versions_blast
        .mix(APPLY_TAXONOMIC_SCREENING.out.versions_biopython)
        .mix(APPLY_TAXONOMIC_SCREENING.out.versions_polars)
}
