process BLASTN_TAXONOMIC_SCREENING {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd@sha256:f327d58dbb16c7ce228ce7079262a5e42d388ce47161397baa0d7e4313a452b7'

    input:
    tuple val(meta), path(baits), path(taxonomic_reference_db, stageAs: 'reference_input/*'), val(kmer_size)

    output:
    tuple val(meta), path('taxonomic_blast_hits.tsv'), emit: hits
    tuple val(meta), path('taxonomic_reference_database.txt'), emit: reference_database_report
    tuple val("${task.process}"), val('blast'), eval('blastn -version 2>&1 | sed "s/^.*blastn: //; s/ .*$//; s/+.*$//"'), topic: versions, emit: versions_blast
    tuple val(meta), val('blast'), eval('blastn -version 2>&1 | sed "s/^.*blastn: //; s/ .*$//; s/+.*$//"'), emit: reported_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Resolve the first BLAST database entry point in the staged directory.
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
        printf 'Taxonomic reference database has no .nal, .nin, or .ndb entry.\\n' >&2
        exit 1
    fi

    (
        cd "\$db_dir"
        unset BLASTDB
        blastdbcmd -db "./\$(basename "\$db_prefix")" -info
    ) > taxonomic_reference_database.txt

    printf 'qseqid\\tsaccver\\tstaxids\\tsscinames\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tstitle\\n' > taxonomic_blast_hits.tsv
    blastn \\
        -query ${baits} \\
        -db "\$db_prefix" \\
        -task blastn \\
        -word_size ${kmer_size} \\
        -ungapped \\
        -perc_identity 100 \\
        -qcov_hsp_perc 100 \\
        -dust no \\
        -outfmt '6 qseqid saccver staxids sscinames pident length mismatch gapopen qstart qend sstart send stitle' \\
        >> taxonomic_blast_hits.tsv
    """

    stub:
    """
    printf 'qseqid\\tsaccver\\tstaxids\\tsscinames\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tstitle\\n' > taxonomic_blast_hits.tsv
    printf 'Taxonomic reference database stub\\n' > taxonomic_reference_database.txt
    """
}
