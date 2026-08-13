include { PREPARE_QUERIES } from '../../../modules/local/prepare_queries/main'
include { EXTRACT_SOURCE_SEQUENCES } from '../../../modules/local/extract_source_sequences/main'
include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLASTN_SOURCE_SEQUENCES } from '../../../modules/local/blastn_source_sequences/main'

workflow DISCOVER_SOURCE_SEQUENCES {
    take:
    ch_discovery_inputs

    main:

    // Prepare per-design inputs
    ch_query_inputs = ch_discovery_inputs.map { meta, representative_queries, query_rules, target_assembly ->
        tuple(meta, representative_queries, query_rules)
    }
    ch_assembly_inputs = ch_discovery_inputs.map { meta, representative_queries, query_rules, target_assembly ->
        tuple(meta, target_assembly)
    }
    ch_extraction_context = ch_discovery_inputs.map { meta, representative_queries, query_rules, target_assembly ->
        tuple(meta, query_rules, target_assembly)
    }

    PREPARE_QUERIES(ch_query_inputs)

    // Build and search the assembly database
    BLAST_MAKEBLASTDB(ch_assembly_inputs)

    ch_blast_inputs = PREPARE_QUERIES.out.fasta
        .join(BLAST_MAKEBLASTDB.out.db)
        .map { meta, prepared_queries, blast_database ->
            tuple(meta, prepared_queries, blast_database)
        }

    BLASTN_SOURCE_SEQUENCES(ch_blast_inputs)

    // Extract source sequences and evidence
    ch_extraction_inputs = PREPARE_QUERIES.out.fasta
        .join(ch_extraction_context)
        .join(BLASTN_SOURCE_SEQUENCES.out.txt)
        .map { meta, prepared_queries, query_rules, target_assembly, blast_hits ->
            tuple(meta, prepared_queries, target_assembly, query_rules, blast_hits)
        }

    EXTRACT_SOURCE_SEQUENCES(ch_extraction_inputs)

    emit:
    source_sequences = EXTRACT_SOURCE_SEQUENCES.out.source_sequences
    candidate_loci = EXTRACT_SOURCE_SEQUENCES.out.candidate_loci
    blast_hits = BLASTN_SOURCE_SEQUENCES.out.txt
    discovery_status = EXTRACT_SOURCE_SEQUENCES.out.discovery_status
    discovery_terminal = EXTRACT_SOURCE_SEQUENCES.out.discovery_terminal
    reported_biopython = EXTRACT_SOURCE_SEQUENCES.out.reported_biopython
    reported_blast = BLASTN_SOURCE_SEQUENCES.out.reported_blast
    versions = PREPARE_QUERIES.out.versions_biopython
        .mix(BLAST_MAKEBLASTDB.out.versions)
        .mix(BLASTN_SOURCE_SEQUENCES.out.versions_blast)
        .mix(EXTRACT_SOURCE_SEQUENCES.out.versions_biopython)
}
