process PREPARE_QUERIES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_python:43cc0a959225cd71@sha256:abd18a716d4d3c5b158555116181c91dda4f1fa8edbcfb07cad291364ed436e4'

    input:
    tuple val(meta), path(representative_queries), path(query_rules)

    output:
    tuple val(meta), path('prepared_queries.fasta'), emit: fasta
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_queries.py \
        --representative-queries ${representative_queries} \
        --query-rules ${query_rules} \
        --output-fasta prepared_queries.fasta
    """
}
