process DEACON_INDEX_DUMP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deacon:0.15.0--hdd79491_0@sha256:d90033841907baece8760175fd8c8de9609c09dbd801ca7afa4cfed0de97c30d'

    input:
    tuple val(meta), path(index)

    output:
    tuple val(meta), path('passing_kmers.fasta'), emit: fasta
    tuple val("${task.process}"), val('deacon'), eval("deacon --version | sed 's/^deacon //'"), topic: versions, emit: versions_deacon

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    deacon index dump \
        -o passing_kmers.fasta \
        ${index}
    """

    stub:
    """
    touch passing_kmers.fasta
    """
}
