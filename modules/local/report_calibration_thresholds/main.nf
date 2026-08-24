process REPORT_CALIBRATION_THRESHOLDS {
    tag "${recommended_thresholds.size()} completed thresholds"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple path(design_ids, stageAs: 'calibration_design_ids.txt'), path(recommended_thresholds, stageAs: 'recommended_thresholds/design????/*', arity: '0..*')

    output:
    path('calibration_thresholds.tsv'), emit: thresholds
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    report_calibration_thresholds.py \
        --design-ids ${design_ids} \
        --thresholds-dir recommended_thresholds \
        --output calibration_thresholds.tsv
    """

    stub:
    """
    printf 'design_id\trecommended_threshold\n' > calibration_thresholds.tsv
    """
}
