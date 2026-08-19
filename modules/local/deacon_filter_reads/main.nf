process DEACON_FILTER_READS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deacon:0.15.0--hdd79491_0@sha256:d90033841907baece8760175fd8c8de9609c09dbd801ca7afa4cfed0de97c30d'

    input:
    tuple val(meta), path(index, stageAs: 'deacon.idx'), path(read, stageAs: 'calibration_read'), val(abs_threshold), val(rel_threshold)

    output:
    tuple val(meta), path("${prefix}.fasta.gz"), emit: fasta_filtered
    tuple val(meta), path("${prefix}.json")  , emit: log
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.id}"
    """
    deacon \
        filter \
        --threads ${task.cpus} \
        --abs-threshold ${abs_threshold} \
        --rel-threshold ${rel_threshold} \
        --summary '${prefix}.json' \
        --fasta \
        --output '${prefix}.fasta.gz' \
        deacon.idx \
        calibration_read
    """

    stub:
    prefix = "${meta.id}"
    """
    gzip </dev/null > '${prefix}.fasta.gz'
    touch '${prefix}.json'
    """
}
