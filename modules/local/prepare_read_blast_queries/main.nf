process PREPARE_READ_BLAST_QUERIES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(counts, stageAs: 'counted_reads/*'), path(fastas, stageAs: 'candidate_fastas/*'), path(statuses, stageAs: 'candidate_statuses/*')

    output:
    tuple val(meta), path('candidate_read_counts.tsv'), emit: candidate_read_counts
    tuple val(meta), path('preparation_summary.tsv'), emit: preparation_summary
    tuple val(meta), path('read_queries.fasta'), optional: true, emit: queries
    tuple val(meta), path('read_blast_hits.tsv'), optional: true, emit: terminal_blast_hits
    tuple val(meta), path('classified_reads.tsv'), optional: true, emit: terminal_classified_reads
    tuple val(meta), path('threshold_read_counts.tsv'), optional: true, emit: terminal_read_counts
    tuple val(meta), path('threshold_summary.tsv'), optional: true, emit: terminal_summary
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    def sorted_counts = (counts instanceof List ? counts : [counts]).sort { path -> path.name }
    def sorted_fastas = (fastas instanceof List ? fastas : [fastas]).sort { path -> path.name }
    def sorted_statuses = (statuses instanceof List ? statuses : [statuses]).sort { path -> path.name }
    """
    prepare_read_blast_queries.py \
        --design-id ${meta.id} \
        --counts ${sorted_counts.join(' ')} \
        --fastas ${sorted_fastas.join(' ')} \
        --statuses ${sorted_statuses.join(' ')} \
        --candidate-counts-out candidate_read_counts.tsv \
        --query-out read_queries.fasta \
        --summary-out preparation_summary.tsv \
        --terminal-blast-hits-out read_blast_hits.tsv \
        --terminal-classified-reads-out classified_reads.tsv \
        --terminal-read-counts-out threshold_read_counts.tsv \
        --terminal-summary-out threshold_summary.tsv
    """
}
