process CALIBRATE_THRESHOLD_EVIDENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(classified_reads), path(preparation_summary)

    output:
    tuple val(meta), path('threshold_read_counts.tsv'), emit: read_counts
    tuple val(meta), path('threshold_summary.tsv'), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    calibrate_threshold.py \
        --design-id ${meta.id} \
        --classified-reads ${classified_reads} \
        --preparation-summary ${preparation_summary} \
        --read-counts-out threshold_read_counts.tsv \
        --summary-out threshold_summary.tsv
    """
}
