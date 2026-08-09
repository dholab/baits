process DEACON_FILTER_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deacon:0.15.0--hdd79491_0@sha256:d90033841907baece8760175fd8c8de9609c09dbd801ca7afa4cfed0de97c30d'

    input:
    tuple val(meta), path(index, stageAs: 'deacon_filter_input/deacon.idx'), path(fasta, stageAs: 'deacon_filter_input/fasta/*')

    output:
    tuple val(meta), path('filtered.fasta'), emit: fasta
    tuple val("${task.process}"), val('deacon'), eval("deacon --version | sed 's/^deacon //'"), topic: versions, emit: versions_deacon
    tuple val(meta), val('deacon'), eval("deacon --version | sed 's/^deacon //'"), emit: reported_deacon

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    deacon filter \
        -a 1 \
        -r 0 \
        ${index} \
        ${fasta} \
        -o filtered.fasta
    """

    stub:
    """
    touch filtered.fasta
    """
}
