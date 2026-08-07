process BLASTN_SOURCE_SEQUENCES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/blast:2.16.0--h66d330f_5@sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c'

    input:
    tuple val(meta), path(fasta), path(db)

    output:
    tuple val(meta), path('query_blast_hits.txt'), emit: txt
    tuple val("${task.process}"), val('blast'), eval('blastn -version 2>&1 | sed "s/^.*blastn: //; s/ .*$//; s/+.*$//"'), topic: versions, emit: versions_blast
    tuple val(meta), val('blast'), eval('blastn -version 2>&1 | sed "s/^.*blastn: //; s/ .*$//; s/+.*$//"'), emit: reported_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    db_prefix=\$(for file in ${db}/*.nin ${db}/*.ndb; do [ -e "\$file" ] && printf '%s\n' "\${file%.*}" && break; done)

    blastn \
        -query ${fasta} \
        -db "\$db_prefix" \
        -task megablast \
        -outfmt '6 qseqid qlen qstart qend sseqid sstart send length pident' \
        -out query_blast_hits.txt
    """

    stub:
    """
    touch query_blast_hits.txt
    """
}
