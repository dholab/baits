include { REPORT_SOURCE_DISCOVERY } from '../../../modules/local/report_source_discovery/main'

workflow REPORT_SOURCE_DISCOVERY_RESULTS {
    take:
    ch_candidate_loci
    ch_source_matches

    main:
    ch_report_rows = ch_candidate_loci
        .join(
            ch_source_matches,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, candidate_loci, source_matches ->
            tuple(meta.id as String, candidate_loci, source_matches)
        }
    ch_sorted_report_rows = ch_report_rows.toSortedList { left, right ->
        left[0] <=> right[0]
    }
    ch_design_ids = ch_sorted_report_rows
        .map { rows ->
            tuple(
                'source_discovery_design_ids.txt',
                rows.collect { design_id, _candidate_loci, _source_matches -> design_id }.join('\n') + (rows ? '\n' : ''),
            )
        }
        .collectFile(sort: false) { filename, contents -> [filename, contents] }
    ch_report_tables = ch_sorted_report_rows.map { rows ->
        tuple(
            rows.collect { _design_id, candidate_loci, _source_matches -> candidate_loci },
            rows.collect { _design_id, _candidate_loci, source_matches -> source_matches },
        )
    }
    ch_report_inputs = ch_report_tables
        .combine(ch_design_ids)
        .map { candidate_loci, source_matches, design_ids ->
            tuple(design_ids, candidate_loci, source_matches)
        }

    REPORT_SOURCE_DISCOVERY(ch_report_inputs)

    emit:
    candidate_loci = REPORT_SOURCE_DISCOVERY.out.candidate_loci
    source_matches = REPORT_SOURCE_DISCOVERY.out.source_matches
    versions = REPORT_SOURCE_DISCOVERY.out.versions_python.mix(
        REPORT_SOURCE_DISCOVERY.out.versions_polars,
    )
}
