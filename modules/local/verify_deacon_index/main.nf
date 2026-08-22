process VERIFY_DEACON_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_python:43cc0a959225cd71@sha256:abd18a716d4d3c5b158555116181c91dda4f1fa8edbcfb07cad291364ed436e4'

    input:
    tuple val(meta), path(baits, stageAs: 'verification_input/baits.fasta'), path(bait_roundtrip, stageAs: 'verification_input/bait_roundtrip.fasta'), path(background_roundtrip, stageAs: 'verification_input/background_roundtrip.fasta'), path(candidate_kmer_manifest, stageAs: 'verification_input/candidate_kmers.tsv'), path(deacon_index, stageAs: 'verification_input/deacon.idx'), val(kmer_size)

    output:
    tuple val(meta), path('verification_input/deacon.idx', includeInputs: true), emit: index
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    verify_deacon_index.py \
        --kmer-size ${kmer_size} \
        --baits ${baits} \
        --bait-roundtrip ${bait_roundtrip} \
        --background-roundtrip ${background_roundtrip} \
        --manifest ${candidate_kmer_manifest}
    """
}
