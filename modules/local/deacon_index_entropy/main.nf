process DEACON_INDEX_ENTROPY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deacon:0.15.0--hdd79491_0@sha256:d90033841907baece8760175fd8c8de9609c09dbd801ca7afa4cfed0de97c30d'

    input:
    tuple val(meta), path(fasta), val(kmer_size), val(deacon_window), val(entropy_threshold)

    output:
    tuple val(meta), path('entropy.idx'), emit: index
    tuple val("${task.process}"), val('deacon'), eval("deacon --version | sed 's/^deacon //'"), topic: versions, emit: versions_deacon
    tuple val(meta), val('deacon'), eval("deacon --version | sed 's/^deacon //'"), emit: reported_deacon

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    deacon index build -k ${kmer_size} -w ${deacon_window} -e ${entropy_threshold} -o entropy.idx ${fasta}
    """

    stub:
    """
    touch entropy.idx
    """
}
