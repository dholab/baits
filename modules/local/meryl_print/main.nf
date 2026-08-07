process MERYL_PRINT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/meryl:1.4.1--h4ac6f70_0@sha256:60ba02cde408b606fc1834ef3261c5abc33796d39bf9640dcee307c256501093'

    input:
    tuple val(meta), path(database)

    output:
    tuple val(meta), path('meryl_print.tsv'), emit: txt
    tuple val("${task.process}"), val('meryl'), eval("meryl --version |& sed 's/meryl //'"), topic: versions, emit: versions_meryl

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    meryl print ${database} > meryl_print.tsv
    """

    stub:
    """
    touch meryl_print.tsv
    """
}
