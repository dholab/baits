process BLASTN_WHOLE_READ_CLASSIFICATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd@sha256:f327d58dbb16c7ce228ce7079262a5e42d388ce47161397baa0d7e4313a452b7'

    input:
    tuple val(meta), path(queries), path(taxonomic_reference_db, stageAs: 'reference_input/*')

    output:
    tuple val(meta), path('whole_read_blast_hits.tsv'), emit: hits
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
        printf 'Taxonomic Reference Database has no .nal, .nin, or .ndb entry.\n' >&2
        exit 1
    fi

    printf 'qseqid\tqlen\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n' > whole_read_blast_hits.tsv
    blastn \
        -query ${queries} \
        -db "\$db_prefix" \
        -task blastn \
        -evalue 1e-10 \
        -max_target_seqs 25 \
        -dust no \
        -outfmt '6 qseqid qlen saccver staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle' \
        >> whole_read_blast_hits.tsv
    """

    stub:
    """
    printf 'qseqid\tqlen\tsaccver\tstaxids\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tstitle\n' > whole_read_blast_hits.tsv
    """
}
