process GATHER_READ_BLAST_HITS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(blast_hit_shards), path(search_parameter_shards)

    output:
    tuple val(meta), path('read_blast_hits.tsv'), emit: hits
    tuple val(meta), path('read_blast_search_parameters.tsv'), emit: search_parameters
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gather_read_blast_hits.py \
        --blast-hits ${blast_hit_shards.join(' ')} \
        --search-parameters ${search_parameter_shards.join(' ')} \
        --hits-output read_blast_hits.tsv \
        --parameters-output read_blast_search_parameters.tsv
    """

    stub:
    """
    printf 'qseqid\tqlen\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n' > read_blast_hits.tsv
    printf 'parameter\tvalue\nquery_file_count\t%s\n' '${blast_hit_shards.size()}' > read_blast_search_parameters.tsv
    """
}
