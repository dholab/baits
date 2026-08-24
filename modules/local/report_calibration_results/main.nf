process REPORT_CALIBRATION_RESULTS {
    tag "${calibration_reads.size()} designs"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple path(design_ids, stageAs: 'calibration_design_ids.txt'), path(calibration_reads, stageAs: 'calibration_reads/design????/*', arity: '0..*'), path(calibration_matches, stageAs: 'calibration_matches/design????/*', arity: '0..*')

    output:
    path('calibration_reads.tsv'), emit: reads
    path('calibration_matches.tsv.gz'), emit: matches
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), topic: versions, emit: versions_polars

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    report_calibration_results.py \
        --design-ids ${design_ids} \
        --calibration-reads-dir calibration_reads \
        --calibration-matches-dir calibration_matches \
        --reads-out calibration_reads.tsv \
        --matches-out calibration_matches.tsv.gz
    """

    stub:
    """
    printf 'design_id\tmetagenome_id\tread_id\tread_length\tbait_count\tquery_id\tbest_target_bit_score\tbest_other_bit_score\n' > calibration_reads.tsv
    python -c "import gzip; open('calibration_matches.tsv.gz', 'wb').write(gzip.compress(b'design_id\\tquery_id\\tquery_length\\tsubject_accession\\ttaxids\\tpercent_identity\\talignment_length\\tmismatches\\tgap_opens\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tevalue\\tbit_score\\tquery_coverage\\tsubject_title\\n', mtime=0))"
    """
}
