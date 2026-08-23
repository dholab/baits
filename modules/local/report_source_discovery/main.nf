process REPORT_SOURCE_DISCOVERY {
    tag "${candidate_loci.size()} designs"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple path(design_ids, stageAs: 'source_discovery_design_ids.txt'), path(candidate_loci, stageAs: 'candidate_loci/design????/*', arity: '0..*'), path(source_matches, stageAs: 'source_matches/design????/*', arity: '0..*')

    output:
    path('candidate_loci.tsv'), emit: candidate_loci
    path('source_matches.tsv.gz'), emit: source_matches
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    report_source_discovery.py \
        --design-ids ${design_ids} \
        --candidate-loci-dir candidate_loci \
        --source-matches-dir source_matches \
        --candidate-loci-out candidate_loci.tsv \
        --source-matches-out source_matches.tsv.gz
    """

    stub:
    """
    printf 'design_id\tcandidate_locus_id\tsource_sequence_id\tquery_id\tquery_group\tassembly_sequence_id\tassembly_start\tassembly_end\tstrand\tpercent_identity\tquery_coverage\n' > candidate_loci.tsv
    python -c "import gzip; open('source_matches.tsv.gz', 'wb').write(gzip.compress(b'design_id\\tquery_id\\tquery_length\\tquery_start\\tquery_end\\tassembly_sequence_id\\tassembly_start\\tassembly_end\\talignment_length\\tpercent_identity\\n', mtime=0))"
    """
}
