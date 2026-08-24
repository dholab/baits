process PREPARE_READ_BLAST_QUERIES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(counts, stageAs: 'read_bait_counts/*'), path(unique_fasta)

    output:
    tuple val(meta), path('read_bait_counts.tsv'), emit: counts
    tuple val(meta), path('read_queries.fasta'), emit: queries
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_read_blast_queries.py \
        --counts-dir read_bait_counts \
        --unique-fasta ${unique_fasta} \
        --counts-out read_bait_counts.tsv \
        --queries-out read_queries.fasta
    """
}
