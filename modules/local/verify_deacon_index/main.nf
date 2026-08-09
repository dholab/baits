process VERIFY_DEACON_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(baits, stageAs: 'verification_input/baits.fasta'), path(bait_roundtrip, stageAs: 'verification_input/bait_roundtrip.fasta'), path(background_roundtrip, stageAs: 'verification_input/background_roundtrip.fasta'), path(candidate_kmer_manifest, stageAs: 'verification_input/candidate_kmers.tsv'), path(bait_set_status, stageAs: 'verification_input/bait_set_status.tsv'), val(kmer_size), val(deacon_window)

    output:
    tuple val(meta), path('bait_set_status.tsv'), emit: bait_set_status
    tuple val(meta), path('verification_summary.tsv'), emit: summary
    tuple val(meta), path('verification_report.md'), emit: report
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    verify_deacon_index.py \
        --design-id ${meta.id} \
        --kmer-size ${kmer_size} \
        --deacon-window ${deacon_window} \
        --baits ${baits} \
        --bait-roundtrip ${bait_roundtrip} \
        --background-roundtrip ${background_roundtrip} \
        --manifest ${candidate_kmer_manifest} \
        --bait-set-status-in ${bait_set_status} \
        --bait-set-status-out bait_set_status.tsv \
        --summary-out verification_summary.tsv \
        --report-out verification_report.md
    """
}
