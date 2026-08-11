#!/usr/bin/env nextflow

include { UTILS_NFSCHEMA_PLUGIN } from './subworkflows/nf-core/utils_nfschema_plugin'
include { BAITS_MAIN } from './workflows/baits'
include { RESOLVE_OPTIMIZATION_READS } from './modules/local/resolve_optimization_reads/main'

workflow {
    main:

    if (params.version) {
        println "${workflow.manifest.name} v${workflow.manifest.version}"
        exit 0
    }

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        params.validate_params,
        null,
        params.help,
        params.help_full,
        params.show_hidden,
        '',
        '',
        'nextflow run dholab/baits --help',
        false,
    )

    // Input normalization
    ch_direct_rows = params.source_sequences
        ? Channel.value([
            id: params.id ?: file(params.source_sequences).baseName,
            target_taxid: params.target_taxid,
            curated_source_sequences: params.source_sequences,
            representative_queries: '',
            query_rules: '',
            target_assembly: '',
            interference_background: params.interference_background,
            optimization_read_set: params.optimization_read_set,
        ])
        : Channel.empty()

    ch_samplesheet_rows = params.input
        ? Channel.fromPath(params.input).splitCsv(header: true)
        : Channel.empty()

    ch_rows = ch_direct_rows.mix(ch_samplesheet_rows)

    // Each row becomes [design metadata, normalized input row].
    ch_normalized_rows = ch_rows.map { row ->
        tuple(
            [id: row.id, source_sequence_origin: row.curated_source_sequences ? 'curated_input' : 'query_guided_discovery'],
            row,
        )
    }

    ch_optimization_read_sets = ch_normalized_rows
        .filter { meta, row -> row.optimization_read_set != null && row.optimization_read_set != '' }
        .map { meta, row -> tuple(meta, file(row.optimization_read_set)) }
    ch_without_optimization_keys = ch_normalized_rows
        .filter { meta, row -> row.optimization_read_set == null || row.optimization_read_set == '' }
        .map { meta, row -> tuple(meta, true) }
    RESOLVE_OPTIMIZATION_READS(ch_optimization_read_sets)
    ch_optimization_read_records = RESOLVE_OPTIMIZATION_READS.out.manifest
        .splitCsv(header: true, sep: '\t', quote: '"', elem: 1)
    ch_optimization_reads = ch_optimization_read_records
        .map { design_meta, row, optimization_read_set ->
            def read_meta = [
                id: row.id,
                design_id: row.design_id,
                metagenome_id: row.metagenome_id,
            ]
            def reads = [row.read_1, row.read_2]
                .findAll { read_name -> read_name }
                .collect { read_name -> file(optimization_read_set.resolve(read_name)) }
            tuple(design_meta, read_meta, reads)
        }
    // Aggregate [facts, roles, IDs, kinds, files] for each design's read files.
    ch_optimization_read_file_facts = ch_optimization_read_records
        .map { design_meta, row, optimization_read_set ->
            def readNames = [row.read_1, row.read_2].findAll { readName -> readName }
            def mates = readNames.size() == 1 ? ['single'] : ['R1', 'R2']
            def inputIds = mates.collect { mate ->
                mate == 'single' ? row.id : "${row.id}:${mate}"
            }
            def inputFacts = inputIds.indices.collectMany { index ->
                [
                    [input_role: 'optimization_read', input_id: inputIds[index], attribute: 'metagenome_id', value: row.metagenome_id],
                    [input_role: 'optimization_read', input_id: inputIds[index], attribute: 'mate', value: mates[index]],
                ]
            }
            def reads = readNames.collect { readName -> file(optimization_read_set.resolve(readName)) }
            tuple(
                design_meta,
                inputFacts,
                ['optimization_read'] * reads.size(),
                inputIds,
                ['file'] * reads.size(),
                reads,
            )
        }
        .groupTuple(by: 0)
        .map { design_meta, input_fact_groups, role_groups, id_groups, kind_groups, read_groups ->
            tuple(
                design_meta,
                input_fact_groups.collectMany { facts -> facts },
                role_groups.collectMany { roles -> roles },
                id_groups.collectMany { ids -> ids },
                kind_groups.collectMany { kinds -> kinds },
                read_groups.collectMany { reads -> reads },
            )
        }
    ch_optimization_resolver_versions = RESOLVE_OPTIMIZATION_READS.out.reported_python
        .join(RESOLVE_OPTIMIZATION_READS.out.reported_polars)
        .map { meta, python_component, python_version, polars_component, polars_version ->
            tuple(
                meta,
                [
                    [component: python_component, version: python_version],
                    [component: polars_component, version: polars_version],
                ],
            )
        }
    ch_optimization_provenance_facts = ch_optimization_read_file_facts
        .join(ch_optimization_resolver_versions)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, software_versions ->
            tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, software_versions)
        }

    ch_kmer_size = Channel.value(params.kmer_size)
    ch_deacon_window = Channel.value(params.deacon_window)
    ch_entropy_threshold = Channel.value(params.entropy_threshold)
    ch_taxonomic_reference_db = params.taxonomic_reference_db
        ? Channel.value(file(params.taxonomic_reference_db))
        : Channel.empty()
    ch_taxonomic_screening_not_run = params.taxonomic_reference_db
        ? Channel.empty()
        : Channel.value(true)

    // Shared design context: [metadata, target taxid, interference background].
    ch_design_context = ch_normalized_rows.map { meta, row ->
        tuple(meta, row.target_taxid as String, file(row.interference_background))
    }

    // Source Sequence acquisition
    source_inputs = ch_normalized_rows.branch { meta, row ->
        curated: row.curated_source_sequences
        discovery: true
    }

    ch_curated_inputs = source_inputs.curated.map { meta, row ->
        tuple(meta, file(row.curated_source_sequences))
    }

    ch_discovery_inputs = source_inputs.discovery.map { meta, row ->
        tuple(
            meta,
            file(row.representative_queries),
            file(row.query_rules),
            file(row.target_assembly),
        )
    }

    ch_curated_provenance_source_inputs = source_inputs.curated.map { meta, row ->
        tuple(
            meta,
            ['curated_source_sequences'],
            [file(row.curated_source_sequences)],
        )
    }

    ch_discovery_provenance_source_inputs = source_inputs.discovery.map { meta, row ->
        tuple(
            meta,
            ['representative_queries', 'query_rules', 'target_assembly'],
            [file(row.representative_queries), file(row.query_rules), file(row.target_assembly)],
        )
    }

    BAITS_MAIN(
        ch_design_context,
        ch_curated_inputs,
        ch_discovery_inputs,
        ch_curated_provenance_source_inputs,
        ch_discovery_provenance_source_inputs,
        ch_kmer_size,
        ch_deacon_window,
        ch_entropy_threshold,
        ch_taxonomic_reference_db,
        ch_taxonomic_screening_not_run,
        ch_optimization_reads,
        ch_optimization_provenance_facts,
        ch_without_optimization_keys,
    )

    // Publication
    publish:
    source_sequences = BAITS_MAIN.out.source_sequences
    candidate_loci = BAITS_MAIN.out.candidate_loci
    query_blast_hits = BAITS_MAIN.out.query_blast_hits
    discovery_status = BAITS_MAIN.out.discovery_status
    provenance = BAITS_MAIN.out.provenance
    candidate_kmers = BAITS_MAIN.out.candidate_kmers
    candidate_kmer_occurrences = BAITS_MAIN.out.candidate_kmer_occurrences
    filtering_status = BAITS_MAIN.out.filtering_status
    locally_filtered_baits = BAITS_MAIN.out.locally_filtered_baits
    taxonomic_blast_hits = BAITS_MAIN.out.taxonomic_blast_hits
    screening_decisions = BAITS_MAIN.out.screening_decisions
    screening_status = BAITS_MAIN.out.screening_status
    taxonomically_screened_baits = BAITS_MAIN.out.taxonomically_screened_baits
    bait_set_status = BAITS_MAIN.out.bait_set_status
    locally_filtered_deacon_index = BAITS_MAIN.out.locally_filtered_deacon_index
    taxonomically_screened_deacon_index = BAITS_MAIN.out.taxonomically_screened_deacon_index
    index_verification_summary = BAITS_MAIN.out.index_verification_summary
    index_verification_report = BAITS_MAIN.out.index_verification_report
    taxonomic_reference_database = BAITS_MAIN.out.taxonomic_reference_database
    candidate_read_counts = BAITS_MAIN.out.candidate_read_counts
    whole_read_blast_hits = BAITS_MAIN.out.whole_read_blast_hits
    classified_reads = BAITS_MAIN.out.classified_reads
    threshold_read_counts = BAITS_MAIN.out.threshold_read_counts
    threshold_curve = BAITS_MAIN.out.threshold_curve
    threshold_summary = BAITS_MAIN.out.threshold_summary
}

output {
    source_sequences {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/01_source_sequences/source_sequences.fasta" }
    }
    candidate_loci {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/01_source_sequences/candidate_loci.tsv" }
    }
    query_blast_hits {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/01_source_sequences/query_blast_hits.tsv" }
    }
    discovery_status {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/01_source_sequences/discovery_status.tsv" }
    }
    provenance {
        mode 'copy'
        path { record ->
            record[1] >> "${record[0].id}/07_provenance/inputs.tsv"
            record[2] >> "${record[0].id}/07_provenance/parameters.tsv"
            record[3] >> "${record[0].id}/07_provenance/software_versions.tsv"
        }
    }
    candidate_kmers {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/02_candidate_kmers/candidate_kmers.tsv" }
    }
    candidate_kmer_occurrences {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/02_candidate_kmers/candidate_kmer_occurrences.tsv" }
    }
    filtering_status {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/02_candidate_kmers/filtering_status.tsv" }
    }
    locally_filtered_baits {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/04_bait_sets/locally_filtered_baits.fasta" }
    }
    taxonomic_blast_hits {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/03_taxonomic_screening/blast_hits.tsv" }
    }
    screening_decisions {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/03_taxonomic_screening/screening_decisions.tsv" }
    }
    screening_status {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/03_taxonomic_screening/screening_status.tsv" }
    }
    taxonomically_screened_baits {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/04_bait_sets/taxonomically_screened_baits.fasta" }
    }
    bait_set_status {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/04_bait_sets/bait_set_status.tsv" }
    }
    locally_filtered_deacon_index {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/05_deacon_index/locally_filtered.idx" }
    }
    taxonomically_screened_deacon_index {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/05_deacon_index/taxonomically_screened.idx" }
    }
    index_verification_summary {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/05_deacon_index/verification_summary.tsv" }
    }
    index_verification_report {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/05_deacon_index/verification_report.md" }
    }
    taxonomic_reference_database {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/07_provenance/taxonomic_reference_database.txt" }
    }
    candidate_read_counts {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/06_calibration/candidate_read_counts.tsv" }
    }
    whole_read_blast_hits {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/06_calibration/whole_read_blast_hits.tsv" }
    }
    classified_reads {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/06_calibration/classified_reads.tsv" }
    }
    threshold_read_counts {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/06_calibration/threshold_read_counts.tsv" }
    }
    threshold_curve {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/06_calibration/threshold_curve.tsv" }
    }
    threshold_summary {
        mode 'copy'
        path { record -> record[1] >> "${record[0].id}/06_calibration/threshold_summary.tsv" }
    }
}
