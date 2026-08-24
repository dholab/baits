process DERIVE_RECOMMENDED_THRESHOLD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(calibration_reads)

    output:
    tuple val(meta), path('recommended_threshold.tsv'), emit: threshold
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    calibrate_threshold.py \
        --design-id ${meta.id} \
        --calibration-reads ${calibration_reads} \
        --output recommended_threshold.tsv
    """
}
