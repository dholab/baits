process RESOLVE_CALIBRATION_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(calibration_read_set)

    output:
    tuple val(meta), path('calibration_reads.tsv'), path(calibration_read_set), emit: manifest
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars
    tuple val(meta), val('python'), eval("python --version | sed 's/Python //'"), emit: reported_python
    tuple val(meta), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: reported_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    resolve_calibration_reads.py \
        --design-id ${meta.id} \
        --directory ${calibration_read_set} \
        --manifest-out calibration_reads.tsv
    """
}
