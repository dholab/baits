process NORMALIZE_CURATED_SOURCE_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_python:43cc0a959225cd71@sha256:abd18a716d4d3c5b158555116181c91dda4f1fa8edbcfb07cad291364ed436e4'

    input:
    tuple val(meta), path(source_sequences)

    output:
    tuple val(meta), path('source_sequences.fasta'), path('source_sequence_query_groups.tsv'), emit: source_sequences
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val(meta), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), emit: reported_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    normalize_curated_source_sequences.py \
        --source-sequences ${source_sequences} \
        --source-sequences-out source_sequences.fasta \
        --source-sequence-query-groups-out source_sequence_query_groups.tsv
    """
}
