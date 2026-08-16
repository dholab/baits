process CLASSIFY_READ_HITS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(candidate_read_counts), path(blast_hits), val(target_taxid), path(calibration_target_taxids)

    output:
    tuple val(meta), path('classified_reads.tsv'), emit: classified_reads
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    def calibrationTargetTaxidsArg = calibration_target_taxids
        ? "--calibration-target-taxids ${calibration_target_taxids}"
        : ''
    """
    classify_read_hits.py \
        --candidate-read-counts ${candidate_read_counts} \
        --blast-hits ${blast_hits} \
        --target-taxid ${target_taxid} \
        ${calibrationTargetTaxidsArg} \
        --output classified_reads.tsv
    """
}
