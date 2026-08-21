#!/usr/bin/env nextflow

def validateInputMode(pipelineParams) {
    if (pipelineParams.help || pipelineParams.help_full || pipelineParams.version || !pipelineParams.input) {
        return
    }

    def directOnlyParameters = [
        'id',
        'source_sequences',
        'target_taxid',
        'background',
        'calibration_reads',
        'calibration_target_taxids',
    ]
    def incompatible = directOnlyParameters.findAll { parameter ->
        def value = pipelineParams[parameter]
        value != null && value != '' && value != false
    }
    if (incompatible) {
        def flags = incompatible.collect { parameter -> "--${parameter}" }.join(', ')
        error """\
--input selects design-CSV mode and cannot be combined with: ${flags}
Put design-specific values in each applicable CSV row; --taxon_ref_db remains a run-level parameter.
""".stripIndent().trim()
    }
}

def resolveCalibrationTargetScope(row) {
    def targetTaxid = row.target_taxid as String
    if (!row.calibration_target_taxids) {
        return [[], [targetTaxid]]
    }

    def scopeFile = file(row.calibration_target_taxids)
    def lines = scopeFile.readLines()
    if (!lines || lines.first() != 'taxid') {
        error "calibration_target_taxids must have exactly one header: taxid"
    }
    def taxids = lines.drop(1)
    if (!taxids || taxids.any { taxid -> !(taxid ==~ /[1-9][0-9]*/) }) {
        error "calibration_target_taxids must contain canonical positive taxids"
    }
    if (taxids.toSet().size() != taxids.size()) {
        error "calibration_target_taxids must contain unique taxids"
    }
    if (!(targetTaxid in taxids)) {
        error "calibration_target_taxids must include target_taxid ${targetTaxid}"
    }
    def canonicalTaxids = taxids.sort { left, right ->
        left.size() <=> right.size() ?: left <=> right
    }
    return [scopeFile, canonicalTaxids]
}

include { UTILS_NFSCHEMA_PLUGIN } from './subworkflows/nf-core/utils_nfschema_plugin'
include { BAITS_MAIN } from './workflows/baits'
include { RESOLVE_CALIBRATION_READS } from './modules/local/resolve_calibration_reads/main'

workflow {
    main:

    validateInputMode(params)

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
            calibration_target_taxids: params.calibration_target_taxids,
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
        if (row.calibration_target_taxids && !row.calibration_reads) {
            error "calibration_target_taxids requires calibration_reads"
        }
        tuple(
            [id: row.id, source_sequence_origin: row.curated_source_sequences ? 'curated_input' : 'query_guided_discovery'],
            row,
        )
    }

    ch_calibration_scope_context = ch_normalized_rows
        .filter { meta, row -> row.calibration_reads != null && row.calibration_reads != '' }
        .map { meta, row ->
            def scope = resolveCalibrationTargetScope(row)
            tuple(meta, scope[0], scope[1])
        }
    ch_calibration_target_scopes = ch_calibration_scope_context
        .map { meta, calibration_target_taxids, canonical_taxids ->
            tuple(meta, calibration_target_taxids)
        }
    ch_calibration_scope_facts = ch_calibration_scope_context
        .map { meta, calibration_target_taxids, canonical_taxids ->
            def hasScopeFile = calibration_target_taxids as boolean
            tuple(
                meta,
                [[
                    input_role: 'calibration_target_scope',
                    input_id: meta.id,
                    attribute: 'taxids',
                    value: canonical_taxids.join(';'),
                ]],
                hasScopeFile ? ['calibration_target_scope'] : [],
                hasScopeFile ? [meta.id] : [],
                hasScopeFile ? ['file'] : [],
                hasScopeFile ? [calibration_target_taxids] : [],
            )
        }

    // Expand each calibration directory into its direct-child FASTQs. Keeping
    // the original paths preserves per-file cache identity when membership changes.
    ch_calibration_read_sets = ch_normalized_rows
        .filter { meta, row -> row.calibration_reads != null && row.calibration_reads != '' }
        .map { meta, row -> tuple(meta, file(row.calibration_reads)) }

    ch_calibration_read_paths = ch_calibration_read_sets
        .flatMap { meta, calibration_read_set ->
            files(
                "${calibration_read_set}/*.{fastq,fq,fastq.gz,fq.gz}",
                checkIfExists: true,
                hidden: true,
            ).collect { read -> tuple(meta, read.name, read) }
        }

    // Python still owns filename validation and stable source IDs. It receives
    // one small, deterministic filename list per design rather than the FASTQs.
    ch_calibration_read_name_files = ch_calibration_read_paths
        .map { meta, read_name, _read -> tuple(meta.id, read_name) }
        .collectFile(newLine: true, sort: true) { design_id, read_name ->
            ["${design_id}.calibration_read_names.txt", read_name]
        }

    // collectFile emits only file paths, so reattach each list to its design.
    ch_calibration_designs_by_name_file = ch_calibration_read_sets
        .map { meta, _calibration_read_set ->
            tuple("${meta.id}.calibration_read_names.txt", meta)
        }

    ch_calibration_resolver_inputs = ch_calibration_read_name_files
        .map { names_file -> tuple(names_file.name, names_file) }
        .join(
            ch_calibration_designs_by_name_file,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _key, names_file, meta -> tuple(meta, names_file) }

    ch_without_calibration_keys = ch_normalized_rows
        .filter { meta, row -> row.calibration_reads == null || row.calibration_reads == '' }
        .map { meta, row -> tuple(meta, true) }

    RESOLVE_CALIBRATION_READS(ch_calibration_resolver_inputs)

    // Join normalized metadata back to the corresponding original FASTQ path.
    // Strict joins turn any duplicate or missing filename into an explicit error.
    ch_resolved_calibration_read_records = RESOLVE_CALIBRATION_READS.out.manifest
        .splitCsv(header: true, sep: '\t', quote: '"', elem: 1)
        .map { design_meta, row -> tuple(design_meta.id, row.read, design_meta, row) }

    ch_original_calibration_read_paths = ch_calibration_read_paths
        .map { design_meta, read_name, read -> tuple(design_meta.id, read_name, read) }

    ch_calibration_read_records = ch_resolved_calibration_read_records
        .join(
            ch_original_calibration_read_paths,
            by: [0, 1],
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _design_id, _read_name, design_meta, row, read -> tuple(design_meta, row, read) }

    ch_calibration_reads = ch_calibration_read_records
        .map { design_meta, row, read ->
            def read_meta = [
                id: row.id,
                design_id: row.design_id,
                metagenome_id: row.metagenome_id,
            ]
            tuple(design_meta, read_meta, read)
        }

    // Aggregate [facts, roles, IDs, kinds, files] for each design's read files.
    ch_calibration_read_file_facts = ch_calibration_read_records
        .map { design_meta, row, read ->
            tuple(
                design_meta,
                [[input_role: 'calibration_read', input_id: row.id, attribute: 'metagenome_id', value: row.metagenome_id]],
                ['calibration_read'],
                [row.id],
                ['file'],
                [read],
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
        .join(ch_calibration_scope_facts)
        .join(ch_calibration_resolver_versions)
        .map { meta, read_facts, read_file_roles, read_file_ids, read_file_kinds, read_files, scope_facts, scope_file_roles, scope_file_ids, scope_file_kinds, scope_files, software_versions ->
            tuple(
                meta,
                read_facts + scope_facts,
                read_file_roles + scope_file_roles,
                read_file_ids + scope_file_ids,
                read_file_kinds + scope_file_kinds,
                read_files + scope_files,
                software_versions,
            )
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
        ch_calibration_target_scopes,
        ch_calibration_provenance_facts,
        ch_without_calibration_keys,
    )

}
