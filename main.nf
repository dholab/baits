#!/usr/bin/env nextflow

include { UTILS_NFSCHEMA_PLUGIN } from './subworkflows/nf-core/utils_nfschema_plugin'
include { BAITS_MAIN } from './workflows/baits'
include { RESOLVE_CALIBRATION_READS } from './modules/local/resolve_calibration_reads/main'

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
            background: params.background,
            calibration_reads: params.calibration_reads,
        ])
        : Channel.empty()

    ch_samplesheet_rows = params.input
        ? Channel.fromPath(params.input).splitCsv(header: true)
        : Channel.empty()

    ch_rows = ch_direct_rows.mix(ch_samplesheet_rows)

    // Each row becomes [design metadata, normalized input row].
    ch_normalized_rows = ch_rows.map { row ->
        if (row.calibration_reads && !params.taxon_ref_db) {
            error "calibration_reads requires taxon_ref_db"
        }
        tuple(
            [id: row.id, source_sequence_origin: row.curated_source_sequences ? 'curated_input' : 'query_guided_discovery'],
            row,
        )
    }

    ch_calibration_read_sets = ch_normalized_rows
        .filter { meta, row -> row.calibration_reads != null && row.calibration_reads != '' }
        .map { meta, row -> tuple(meta, file(row.calibration_reads)) }
    ch_without_calibration_keys = ch_normalized_rows
        .filter { meta, row -> row.calibration_reads == null || row.calibration_reads == '' }
        .map { meta, row -> tuple(meta, true) }
    RESOLVE_CALIBRATION_READS(ch_calibration_read_sets)
    ch_calibration_read_records = RESOLVE_CALIBRATION_READS.out.manifest
        .splitCsv(header: true, sep: '\t', quote: '"', elem: 1)
    ch_calibration_reads = ch_calibration_read_records
        .map { design_meta, row, calibration_read_set ->
            def read_meta = [
                id: row.id,
                design_id: row.design_id,
                metagenome_id: row.metagenome_id,
            ]
            def reads = [row.read_1, row.read_2]
                .findAll { read_name -> read_name }
                .collect { read_name -> file(calibration_read_set.resolve(read_name)) }
            tuple(design_meta, read_meta, reads)
        }
    // Aggregate [facts, roles, IDs, kinds, files] for each design's read files.
    ch_calibration_read_file_facts = ch_calibration_read_records
        .map { design_meta, row, calibration_read_set ->
            def readNames = [row.read_1, row.read_2].findAll { readName -> readName }
            def mates = readNames.size() == 1 ? ['single'] : ['R1', 'R2']
            def inputIds = mates.collect { mate ->
                mate == 'single' ? row.id : "${row.id}:${mate}"
            }
            def inputFacts = inputIds.indices.collectMany { index ->
                [
                    [input_role: 'calibration_read', input_id: inputIds[index], attribute: 'metagenome_id', value: row.metagenome_id],
                    [input_role: 'calibration_read', input_id: inputIds[index], attribute: 'mate', value: mates[index]],
                ]
            }
            def reads = readNames.collect { readName -> file(calibration_read_set.resolve(readName)) }
            tuple(
                design_meta,
                inputFacts,
                ['calibration_read'] * reads.size(),
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
    ch_calibration_resolver_versions = RESOLVE_CALIBRATION_READS.out.reported_python
        .join(RESOLVE_CALIBRATION_READS.out.reported_polars)
        .map { meta, python_component, python_version, polars_component, polars_version ->
            tuple(
                meta,
                [
                    [component: python_component, version: python_version],
                    [component: polars_component, version: polars_version],
                ],
            )
        }
    ch_calibration_provenance_facts = ch_calibration_read_file_facts
        .join(ch_calibration_resolver_versions)
        .map { meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, software_versions ->
            tuple(meta, input_facts, input_file_roles, input_file_ids, input_file_kinds, input_files, software_versions)
        }

    ch_kmer_size = Channel.value(params.kmer_size)
    ch_deacon_window = Channel.value(params.deacon_window)
    ch_entropy_threshold = Channel.value(params.entropy_threshold)
    ch_taxonomic_reference_db = params.taxon_ref_db
        ? Channel.value(file(params.taxon_ref_db))
        : Channel.empty()
    ch_taxonomic_screening_not_run = params.taxon_ref_db
        ? Channel.empty()
        : Channel.value(true)

    // Shared design context: [metadata, target taxid, interference background].
    ch_design_context = ch_normalized_rows.map { meta, row ->
        tuple(meta, row.target_taxid as String, file(row.background))
    }

    // Source sequence acquisition
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
        ch_calibration_reads,
        ch_calibration_provenance_facts,
        ch_without_calibration_keys,
    )

}
