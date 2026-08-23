process APPLY_COMPLEXITY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(manifest, stageAs: 'filtering_input/candidate_kmers.tsv'), path(passing_kmers)

    output:
    tuple val(meta), path('candidate_kmers.tsv'), emit: manifest
    tuple val(meta), path('selected_candidates.fasta'), emit: candidates
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    apply_complexity_filter.py \
        --manifest-in ${manifest} \
        --passing-kmers ${passing_kmers} \
        --manifest-out candidate_kmers.tsv \
        --selected-candidates-out selected_candidates.fasta
    """
}
