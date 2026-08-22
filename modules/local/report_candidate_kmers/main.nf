process REPORT_CANDIDATE_KMERS {
    tag "${candidate_manifests.size()} designs"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple path(design_ids, stageAs: 'candidate_report_design_ids.txt'), path(candidate_manifests, stageAs: 'candidate_manifests/design????/*', arity: '1..*'), path(candidate_occurrences, stageAs: 'candidate_occurrences/design????/*', arity: '1..*')

    output:
    path('candidate_kmers.tsv'), emit: candidates
    path('candidate_kmer_occurrences.tsv'), emit: occurrences
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    report_candidate_kmers.py \
        --design-ids ${design_ids} \
        --candidate-manifests-dir candidate_manifests \
        --candidate-occurrences-dir candidate_occurrences \
        --candidates-out candidate_kmers.tsv \
        --occurrences-out candidate_kmer_occurrences.tsv
    """

    stub:
    """
    printf 'design_id\tkmer\tbackground_occurrences\toutcome\n' > candidate_kmers.tsv
    printf 'design_id\tkmer\tsource_sequence_id\tstart\tquery_group\n' > candidate_kmer_occurrences.tsv
    """
}
