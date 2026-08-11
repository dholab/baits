include { DEACON_FILTER_READS as DEACON_FILTER_CANDIDATE_READS } from '../../../modules/local/deacon_filter_reads/main'
include { COUNT_READ_BAITS } from '../../../modules/local/count_read_baits/main'
include { PREPARE_READ_BLAST_QUERIES } from '../../../modules/local/prepare_read_blast_queries/main'
include { BLASTN_WHOLE_READ_CLASSIFICATION } from '../../../modules/local/blastn_whole_read_classification/main'
include { CLASSIFY_READ_HITS } from '../../../modules/local/classify_read_hits/main'
include { CALIBRATE_DEACON_THRESHOLD_EVIDENCE } from '../../../modules/local/calibrate_deacon_threshold_evidence/main'


workflow CALIBRATE_DEACON_THRESHOLD {
    take:
    ch_calibration_inputs
    ch_taxonomic_reference_db
    ch_kmer_size

    main:

    // Retain every read with at least one indexed Bait
    ch_deacon_inputs = ch_calibration_inputs.map { design_meta, read_meta, reads, baits, bait_set_status, deacon_index, target_taxid ->
        tuple(read_meta, deacon_index, reads, 1, 0)
    }
    ch_read_context = ch_calibration_inputs.map { design_meta, read_meta, reads, baits, bait_set_status, deacon_index, target_taxid ->
        tuple(read_meta, design_meta, baits, bait_set_status, target_taxid)
    }
    ch_design_context = ch_calibration_inputs
        .map { design_meta, read_meta, reads, baits, bait_set_status, deacon_index, target_taxid -> tuple(design_meta, target_taxid) }
        .unique()
    DEACON_FILTER_CANDIDATE_READS(ch_deacon_inputs)

    // Recount distinct Baits on each individual read
    ch_recount_inputs = DEACON_FILTER_CANDIDATE_READS.out.fastq_filtered
        .join(ch_read_context)
        .combine(ch_kmer_size)
        .map { read_meta, candidate_reads, design_meta, baits, bait_set_status, target_taxid, kmer_size ->
            tuple(design_meta, read_meta, baits, candidate_reads, kmer_size)
        }
    COUNT_READ_BAITS(ch_recount_inputs)

    // Aggregate metagenomes and prepare one whole-read BLAST Query FASTA per design
    ch_grouped_counts = COUNT_READ_BAITS.out.counts
        .map { design_meta, read_meta, counts -> tuple(design_meta, counts) }
        .groupTuple(by: 0)
    ch_grouped_fastas = COUNT_READ_BAITS.out.fasta
        .map { design_meta, read_meta, fasta -> tuple(design_meta, fasta) }
        .groupTuple(by: 0)
    ch_grouped_statuses = COUNT_READ_BAITS.out.status
        .map { design_meta, read_meta, status -> tuple(design_meta, status) }
        .groupTuple(by: 0)
    ch_preparation_inputs = ch_grouped_counts
        .join(ch_grouped_fastas)
        .join(ch_grouped_statuses)
    PREPARE_READ_BLAST_QUERIES(ch_preparation_inputs)

    // Classify each Candidate Read from its representative's best whole-read alignments
    ch_whole_read_blast_inputs = PREPARE_READ_BLAST_QUERIES.out.queries
        .combine(ch_taxonomic_reference_db)
        .map { design_meta, queries, taxonomic_reference_db ->
            tuple(design_meta, queries, taxonomic_reference_db)
        }
    BLASTN_WHOLE_READ_CLASSIFICATION(ch_whole_read_blast_inputs)
    ch_classification_inputs = PREPARE_READ_BLAST_QUERIES.out.candidate_read_counts
        .join(BLASTN_WHOLE_READ_CLASSIFICATION.out.hits)
        .join(ch_design_context)
        .map { design_meta, candidate_read_counts, blast_hits, target_taxid ->
            tuple(design_meta, candidate_read_counts, blast_hits, target_taxid)
        }
    CLASSIFY_READ_HITS(ch_classification_inputs)

    // Sweep every observed absolute threshold plus the all-reads-removed boundary
    ch_threshold_inputs = CLASSIFY_READ_HITS.out.classified_reads
        .join(PREPARE_READ_BLAST_QUERIES.out.preparation_summary)
    CALIBRATE_DEACON_THRESHOLD_EVIDENCE(ch_threshold_inputs)

    ch_candidate_reads = DEACON_FILTER_CANDIDATE_READS.out.fastq_filtered
        .join(ch_read_context)
        .map { read_meta, candidate_reads, design_meta, baits, bait_set_status, target_taxid ->
            tuple(design_meta, read_meta, candidate_reads)
        }
    ch_deacon_summaries = DEACON_FILTER_CANDIDATE_READS.out.log
        .join(ch_read_context)
        .map { read_meta, deacon_summary, design_meta, baits, bait_set_status, target_taxid ->
            tuple(design_meta, read_meta, deacon_summary)
        }

    emit:
    candidate_reads = ch_candidate_reads
    deacon_summaries = ch_deacon_summaries
    candidate_read_counts = PREPARE_READ_BLAST_QUERIES.out.candidate_read_counts
    whole_read_blast_hits = BLASTN_WHOLE_READ_CLASSIFICATION.out.hits
        .mix(PREPARE_READ_BLAST_QUERIES.out.terminal_blast_hits)
    classified_reads = CLASSIFY_READ_HITS.out.classified_reads
        .mix(PREPARE_READ_BLAST_QUERIES.out.terminal_classified_reads)
    threshold_read_counts = CALIBRATE_DEACON_THRESHOLD_EVIDENCE.out.read_counts
        .mix(PREPARE_READ_BLAST_QUERIES.out.terminal_read_counts)
    threshold_curve = CALIBRATE_DEACON_THRESHOLD_EVIDENCE.out.curve
        .mix(PREPARE_READ_BLAST_QUERIES.out.terminal_curve)
    summary = PREPARE_READ_BLAST_QUERIES.out.terminal_summary
        .mix(CALIBRATE_DEACON_THRESHOLD_EVIDENCE.out.summary)
    reported_deacon = DEACON_FILTER_CANDIDATE_READS.out.versions_deacon
    reported_blast = BLASTN_WHOLE_READ_CLASSIFICATION.out.reported_blast
    versions = DEACON_FILTER_CANDIDATE_READS.out.versions_deacon
        .mix(COUNT_READ_BAITS.out.versions_biopython)
        .mix(COUNT_READ_BAITS.out.versions_polars)
        .mix(PREPARE_READ_BLAST_QUERIES.out.versions_biopython)
        .mix(PREPARE_READ_BLAST_QUERIES.out.versions_polars)
        .mix(BLASTN_WHOLE_READ_CLASSIFICATION.out.versions_blast)
        .mix(CLASSIFY_READ_HITS.out.versions_python)
        .mix(CLASSIFY_READ_HITS.out.versions_polars)
        .mix(CALIBRATE_DEACON_THRESHOLD_EVIDENCE.out.versions_python)
        .mix(CALIBRATE_DEACON_THRESHOLD_EVIDENCE.out.versions_polars)
}
