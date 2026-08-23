process REPORT_TAXONOMIC_MATCHES {
    tag "${candidate_manifests.size()} designs"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple path(design_ids, stageAs: 'taxonomic_match_design_ids.txt'), path(candidate_manifests, stageAs: 'candidate_manifests/design????/*', arity: '0..*'), path(taxonomic_matches, stageAs: 'taxonomic_matches/design????/*', arity: '0..*')

    output:
    path('taxonomic_matches.tsv.gz'), emit: matches
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    report_taxonomic_matches.py \
        --design-ids ${design_ids} \
        --candidate-manifests-dir candidate_manifests \
        --taxonomic-matches-dir taxonomic_matches \
        --matches-out taxonomic_matches.tsv.gz
    """

    stub:
    """
    python -c "import gzip; open('taxonomic_matches.tsv.gz', 'wb').write(gzip.compress(b'design_id\\tkmer\\tsubject_accession\\ttaxids\\tscientific_names\\tpercent_identity\\talignment_length\\tmismatches\\tgap_opens\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tsubject_title\\n', mtime=0))"
    """
}
