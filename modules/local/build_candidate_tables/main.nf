process BUILD_CANDIDATE_TABLES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(source_sequences), path(source_sequence_query_groups), path(meryl_source_counts, stageAs: 'meryl_source_counts.tsv'), path(background_intersection_counts, stageAs: 'background_intersection_counts.tsv'), val(kmer_size)

    output:
    tuple val(meta), path('candidate_kmers.tsv'), path('candidate_kmer_occurrences.tsv'), emit: tables
    tuple val(meta), path('complexity_candidates.fasta'), optional: true, emit: complexity_candidates
    tuple val(meta), path('filtering_status.tsv'), optional: true, emit: terminal_status
    tuple val(meta), path('terminal_candidate_kmers.tsv'), optional: true, emit: terminal_manifest
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars
    tuple val(meta), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), emit: reported_biopython
    tuple val(meta), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: reported_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_candidate_tables.py \
        --design-id ${meta.id} \
        --source-sequences ${source_sequences} \
        --source-sequence-query-groups ${source_sequence_query_groups} \
        --meryl-source-counts ${meryl_source_counts} \
        --background-intersection-counts ${background_intersection_counts} \
        --kmer-size ${kmer_size} \
        --manifest-out candidate_kmers.tsv \
        --occurrences-out candidate_kmer_occurrences.tsv \
        --complexity-candidates-out complexity_candidates.fasta \
        --filtering-status-out filtering_status.tsv \
        --terminal-manifest-out terminal_candidate_kmers.tsv
    """
}
