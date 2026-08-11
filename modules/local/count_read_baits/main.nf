process COUNT_READ_BAITS {
    tag "${read_meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(design_meta), val(read_meta), path(baits), path(reads, stageAs: 'candidate_read_input/*'), val(kmer_size)

    output:
    tuple val(design_meta), val(read_meta), path("${read_meta.id}.counted_reads.tsv"), emit: counts
    tuple val(design_meta), val(read_meta), path("${read_meta.id}.candidate_reads.fasta"), emit: fasta
    tuple val(design_meta), val(read_meta), path("${read_meta.id}.candidate_read_status.tsv"), emit: status
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    count_read_baits.py \
        --metagenome-id ${read_meta.metagenome_id} \
        --baits ${baits} \
        --kmer-size ${kmer_size} \
        --reads ${reads} \
        --counts-out ${read_meta.id}.counted_reads.tsv \
        --fasta-out ${read_meta.id}.candidate_reads.fasta \
        --status-out ${read_meta.id}.candidate_read_status.tsv
    """
}
