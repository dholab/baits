include { DEACON_FILTER_READS as DEACON_RETRIEVE_CALIBRATION_READS } from '../../../modules/local/deacon_filter_reads/main'
include { COUNT_READ_BAITS } from '../../../modules/local/count_read_baits/main'
include { FIND_CONCATENATE } from '../../../modules/nf-core/find/concatenate/main'
include { SEQKIT_RMDUP } from '../../../modules/nf-core/seqkit/rmdup/main'
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT_READ_QUERIES } from '../../../modules/nf-core/seqkit/split2/main'
include { PREPARE_READ_BLAST_QUERIES } from '../../../modules/local/prepare_read_blast_queries/main'
include { BLASTN_CALIBRATION_READS } from '../../../modules/local/blastn_calibration_reads/main'
include { GATHER_READ_BLAST_HITS } from '../../../modules/local/gather_read_blast_hits/main'
include { SCORE_CALIBRATION_READS } from '../../../modules/local/score_calibration_reads/main'
include { DERIVE_RECOMMENDED_THRESHOLD } from '../../../modules/local/calibrate_threshold/main'


workflow CALIBRATE_THRESHOLD {
    take:
    ch_calibration_inputs
    ch_taxonomic_reference_db
    ch_kmer_size

    main:

    // Retain every read with at least one indexed bait
    ch_deacon_inputs = ch_calibration_inputs.map { design_meta, read_meta, reads, baits, deacon_index, target_taxid, calibration_target_taxids ->
        tuple(read_meta, deacon_index, reads, 1, 0)
    }
    ch_read_context = ch_calibration_inputs.map { design_meta, read_meta, reads, baits, deacon_index, target_taxid, calibration_target_taxids ->
        tuple(read_meta, design_meta, baits, target_taxid, calibration_target_taxids)
    }
    ch_design_context = ch_calibration_inputs
        .map { design_meta, read_meta, reads, baits, deacon_index, target_taxid, calibration_target_taxids ->
            tuple(design_meta, target_taxid, calibration_target_taxids)
        }
        .unique()
    DEACON_RETRIEVE_CALIBRATION_READS(ch_deacon_inputs)

    // Count baits on each individual read
    ch_recount_inputs = DEACON_RETRIEVE_CALIBRATION_READS.out.fasta_filtered
        .join(ch_read_context)
        .combine(ch_kmer_size)
        .map { read_meta, candidate_reads, design_meta, baits, target_taxid, calibration_target_taxids, kmer_size ->
            tuple(design_meta, read_meta, baits, candidate_reads, kmer_size)
        }
    COUNT_READ_BAITS(ch_recount_inputs)

    // Aggregate metagenomes and prepare one read BLAST query FASTA per design
    ch_grouped_counts = COUNT_READ_BAITS.out.counts
        .map { design_meta, read_meta, counts -> tuple(design_meta, counts) }
        .groupTuple(by: 0)
    ch_grouped_fastas = COUNT_READ_BAITS.out.fasta
        .map { design_meta, read_meta, fasta -> tuple(design_meta, fasta) }
        .groupTuple(by: 0)
    FIND_CONCATENATE(ch_grouped_fastas)
    SEQKIT_RMDUP(FIND_CONCATENATE.out.file_out)
    ch_seqkit_version = SEQKIT_RMDUP.out.versions_seqkit
        .map { process, component, version -> version }
        .first()
    ch_reported_seqkit = SEQKIT_RMDUP.out.fastx
        .combine(ch_seqkit_version)
        .map { meta, fasta, version -> tuple(meta, 'seqkit', version) }
    ch_preparation_inputs = ch_grouped_counts
        .join(SEQKIT_RMDUP.out.fastx)
    PREPARE_READ_BLAST_QUERIES(ch_preparation_inputs)

    // Split unique read queries across a bounded number of BLAST searches
    ch_split_inputs = PREPARE_READ_BLAST_QUERIES.out.queries
        .filter { design_meta, queries -> queries.size() > 0 }
        .map { design_meta, queries -> tuple([id: design_meta.id, single_end: true, design: design_meta], queries) }
    SEQKIT_SPLIT_READ_QUERIES(ch_split_inputs)
    ch_read_blast_inputs = SEQKIT_SPLIT_READ_QUERIES.out.reads
        .transpose()
        .map { split_meta, queries -> tuple(split_meta.design, queries) }
        .combine(ch_taxonomic_reference_db)
    BLASTN_CALIBRATION_READS(ch_read_blast_inputs)
    ch_grouped_blast_hits = BLASTN_CALIBRATION_READS.out.hits.groupTuple(by: 0)
    GATHER_READ_BLAST_HITS(ch_grouped_blast_hits)

    // Project each shared query's best target and other scores onto every read
    ch_scoring_inputs = PREPARE_READ_BLAST_QUERIES.out.counts
        .join(GATHER_READ_BLAST_HITS.out.hits)
        .join(ch_design_context)
        .map { design_meta, read_bait_counts, blast_hits, target_taxid, calibration_target_taxids ->
            tuple(design_meta, read_bait_counts, blast_hits, target_taxid, calibration_target_taxids)
        }
    SCORE_CALIBRATION_READS(ch_scoring_inputs)

    // Keep only the recommended value needed by the final per-design projection
    DERIVE_RECOMMENDED_THRESHOLD(SCORE_CALIBRATION_READS.out.reads)

    emit:
    prepared_reads = PREPARE_READ_BLAST_QUERIES.out.counts
    calibration_reads = SCORE_CALIBRATION_READS.out.reads
    calibration_matches = GATHER_READ_BLAST_HITS.out.hits
    recommended_threshold = DERIVE_RECOMMENDED_THRESHOLD.out.threshold
    reported_blast = BLASTN_CALIBRATION_READS.out.reported_blast.unique()
    reported_seqkit = ch_reported_seqkit
    versions = DEACON_RETRIEVE_CALIBRATION_READS.out.versions_deacon
        .mix(COUNT_READ_BAITS.out.versions_biopython)
        .mix(PREPARE_READ_BLAST_QUERIES.out.versions_biopython)
        .mix(FIND_CONCATENATE.out.versions_find)
        .mix(FIND_CONCATENATE.out.versions_pigz)
        .mix(FIND_CONCATENATE.out.versions_coreutils)
        .mix(SEQKIT_RMDUP.out.versions_seqkit)
        .mix(SEQKIT_SPLIT_READ_QUERIES.out.versions_seqkit)
        .mix(BLASTN_CALIBRATION_READS.out.versions_blast)
        .mix(GATHER_READ_BLAST_HITS.out.versions_python)
        .mix(SCORE_CALIBRATION_READS.out.versions_python)
        .mix(SCORE_CALIBRATION_READS.out.versions_polars)
        .mix(DERIVE_RECOMMENDED_THRESHOLD.out.versions_python)
}
