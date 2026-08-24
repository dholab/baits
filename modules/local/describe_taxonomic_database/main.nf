process DESCRIBE_TAXONOMIC_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd@sha256:f327d58dbb16c7ce228ce7079262a5e42d388ce47161397baa0d7e4313a452b7'

    input:
    path(taxonomic_reference_db, stageAs: 'reference_input/*')

    output:
    path('taxonomic_database.txt'), emit: description
    tuple val("${task.process}"), val('blast'), eval('blastdbcmd -version 2>&1 | sed "s/^.*blastdbcmd: //; s/ .*$//; s/+.*$//"'), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    db_dir="${taxonomic_reference_db}"
    export BLASTDB="\$db_dir"
    db_prefix=''
    for db_file in "\$db_dir"/*.nal; do
        if [ -e "\$db_file" ]; then
            db_prefix="\${db_file%.*}"
            break
        fi
    done
    if [ -z "\$db_prefix" ]; then
        for db_file in "\$db_dir"/*.ndb "\$db_dir"/*.nin; do
            if [ -e "\$db_file" ]; then
                db_prefix="\${db_file%.*}"
                break
            fi
        done
    fi
    if [ -z "\$db_prefix" ]; then
        printf 'Taxonomic reference database has no .nal, .nin, or .ndb entry.\n' >&2
        exit 1
    fi

    (
        cd "\$db_dir"
        unset BLASTDB
        blastdbcmd -db "./\$(basename "\$db_prefix")" -info
    ) > taxonomic_database.txt
    """

    stub:
    """
    printf 'Taxonomic reference database stub\n' > taxonomic_database.txt
    """
}
