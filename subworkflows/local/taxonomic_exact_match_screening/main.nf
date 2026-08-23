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
        .map { meta, candidates, candidate_kmer_manifest, target_taxid, taxonomic_reference_db, kmer_size ->
            tuple(meta, candidates, taxonomic_reference_db, kmer_size)
        }
    BLASTN_TAXONOMIC_SCREENING(ch_blast_inputs)

    // Apply the allowed taxonomic scope
    ch_screening_interpreter_inputs = BLASTN_TAXONOMIC_SCREENING.out.hits
        .join(ch_screening_inputs)
        .combine(ch_kmer_size)
        .map { meta, blast_hits, candidates, candidate_kmer_manifest, target_taxid, kmer_size ->
            tuple(meta, candidates, blast_hits, candidate_kmer_manifest, target_taxid, kmer_size)
        }
    APPLY_TAXONOMIC_SCREENING(ch_screening_interpreter_inputs)

    ch_selected_candidates = APPLY_TAXONOMIC_SCREENING.out.manifest
        .join(APPLY_TAXONOMIC_SCREENING.out.candidates)
        .filter { meta, manifest, candidates -> candidates.size() > 0 }
        .map { meta, manifest, candidates -> tuple(meta, candidates, manifest) }

    emit:
    hits = BLASTN_TAXONOMIC_SCREENING.out.hits
    reference_database_report = BLASTN_TAXONOMIC_SCREENING.out.reference_database_report
    manifest = APPLY_TAXONOMIC_SCREENING.out.manifest
    completed = APPLY_TAXONOMIC_SCREENING.out.manifest
    candidates = ch_selected_candidates
    reported_blast = BLASTN_TAXONOMIC_SCREENING.out.reported_blast
    reported_biopython = APPLY_TAXONOMIC_SCREENING.out.reported_biopython
    reported_polars = APPLY_TAXONOMIC_SCREENING.out.reported_polars
    versions = BLASTN_TAXONOMIC_SCREENING.out.versions_blast
        .mix(APPLY_TAXONOMIC_SCREENING.out.versions_biopython)
        .mix(APPLY_TAXONOMIC_SCREENING.out.versions_polars)
}
