process APPLY_COMPLEXITY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), val(source_sequence_origin), path(manifest), path(passing_kmers)

    output:
    tuple val(meta), path('candidate_kmers.tsv'), path('filtering_status.tsv'), emit: evidence
    tuple val(meta), path('locally_filtered_baits.fasta'), path('candidate_kmers.tsv'), path('bait_set_status.tsv'), optional: true, emit: baits
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    apply_complexity_filter.py \
        --design-id ${meta.id} \
        --source-sequence-origin ${source_sequence_origin} \
        --manifest-in ${manifest} \
        --passing-kmers ${passing_kmers} \
        --manifest-out candidate_kmers.tsv \
        --baits-out locally_filtered_baits.fasta \
        --filtering-status-out filtering_status.tsv \
        --bait-set-status-out bait_set_status.tsv
    """
}
