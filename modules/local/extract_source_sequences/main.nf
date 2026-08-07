process EXTRACT_SOURCE_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_python:43cc0a959225cd71@sha256:abd18a716d4d3c5b158555116181c91dda4f1fa8edbcfb07cad291364ed436e4'

    input:
    tuple val(meta), path(prepared_queries), path(target_assembly), path(query_rules), path(blast_hits)

    output:
    tuple val(meta), path('source_sequences.fasta'), path('source_sequence_query_groups.tsv'), optional: true, emit: source_sequences
    tuple val(meta), path('candidate_loci.tsv')                                                    , emit: candidate_loci
    tuple val(meta), path('discovery_status.tsv')                                                  , emit: discovery_status
    tuple val(meta), path('discovery_terminal_status.tsv'), optional: true                         , emit: discovery_terminal
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val(meta), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), emit: reported_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_source_sequences.py \
        --prepared-queries ${prepared_queries} \
        --target-assembly ${target_assembly} \
        --query-rules ${query_rules} \
        --blast-hits ${blast_hits} \
        --source-sequences-out source_sequences.fasta \
        --source-sequence-query-groups-out source_sequence_query_groups.tsv \
        --candidate-loci-out candidate_loci.tsv \
        --discovery-status-out discovery_status.tsv \
        --terminal-status-out discovery_terminal_status.tsv
    """
}
