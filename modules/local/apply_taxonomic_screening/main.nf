process APPLY_TAXONOMIC_SCREENING {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), path(locally_filtered_baits), path(blast_hits), path(candidate_kmer_manifest, stageAs: 'screening_input/candidate_kmers.tsv'), path(bait_set_status, stageAs: 'screening_input/bait_set_status.tsv'), val(target_taxid), val(kmer_size)

    output:
    tuple val(meta), path('candidate_kmers.tsv'), emit: manifest
    tuple val(meta), path('screening_decisions.tsv'), emit: decisions
    tuple val(meta), path('screening_status.tsv'), emit: screening_status
    tuple val(meta), path('bait_set_status.tsv'), emit: bait_set_status
    tuple val(meta), path('taxonomically_screened_baits.fasta'), path('candidate_kmers.tsv'), path('bait_set_status.tsv'), optional: true, emit: baits
    tuple val(meta), path('terminal_bait_set_status.tsv'), optional: true, emit: terminal_bait_set_status
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars
    tuple val(meta), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), emit: reported_biopython
    tuple val(meta), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: reported_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    apply_taxonomic_screening.py \
        --design-id ${meta.id} \
        --target-taxid ${target_taxid} \
        --kmer-size ${kmer_size} \
        --baits ${locally_filtered_baits} \
        --blast-hits ${blast_hits} \
        --manifest-in ${candidate_kmer_manifest} \
        --bait-set-status-in ${bait_set_status} \
        --manifest-out candidate_kmers.tsv \
        --baits-out taxonomically_screened_baits.fasta \
        --decisions-out screening_decisions.tsv \
        --screening-status-out screening_status.tsv \
        --bait-set-status-out bait_set_status.tsv \
        --terminal-bait-set-status-out terminal_bait_set_status.tsv
    """
}
