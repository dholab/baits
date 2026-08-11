process DEACON_FILTER_READS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/deacon:0.15.0--hdd79491_0@sha256:d90033841907baece8760175fd8c8de9609c09dbd801ca7afa4cfed0de97c30d'

    input:
    tuple val(meta), path(index), path(reads), val(abs_threshold), val(rel_threshold)

    output:
    tuple val(meta), path("${prefix}*.fq.gz"), emit: fastq_filtered
    tuple val(meta), path("${prefix}.json")  , emit: log
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.id}"
    def read_type = reads.size() == 1 ? "> ${prefix}.fq" : "-o ${prefix}_1.fq -O ${prefix}_2.fq" // deacon's automatic compression does not work
    """
    deacon \
        filter \
        --threads ${task.cpus} \
        --abs-threshold ${abs_threshold} \
        --rel-threshold ${rel_threshold} \
        --summary ${prefix}.json \
        $index \
        $reads \
        ${read_type}

    gzip -f ${prefix}*.fq
    """

    stub:
    prefix = "${meta.id}"
    """
    gzip </dev/null > '${prefix}.fq.gz'
    touch ${prefix}.json
    """
}
