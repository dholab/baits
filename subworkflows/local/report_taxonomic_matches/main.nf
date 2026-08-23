include { REPORT_TAXONOMIC_MATCHES } from '../../../modules/local/report_taxonomic_matches/main'

workflow REPORT_TAXONOMIC_MATCH_RESULTS {
    take:
    ch_taxonomic_matches
    ch_candidate_manifests

    main:
    ch_report_rows = ch_taxonomic_matches
        .join(
            ch_candidate_manifests,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, taxonomic_matches, candidate_manifest ->
            tuple(meta.id as String, candidate_manifest, taxonomic_matches)
        }
    ch_sorted_report_rows = ch_report_rows.toSortedList { left, right ->
        left[0] <=> right[0]
    }
    ch_design_ids = ch_sorted_report_rows
        .map { rows ->
            tuple(
                'taxonomic_match_design_ids.txt',
                rows.collect { design_id, _manifest, _matches -> design_id }.join('\n') + (rows ? '\n' : ''),
            )
        }
        .collectFile(sort: false) { filename, contents -> [filename, contents] }
    ch_report_tables = ch_sorted_report_rows.map { rows ->
        tuple(
            rows.collect { _design_id, manifest, _matches -> manifest },
            rows.collect { _design_id, _manifest, matches -> matches },
        )
    }
    ch_report_inputs = ch_report_tables
        .combine(ch_design_ids)
        .map { manifests, matches, design_ids -> tuple(design_ids, manifests, matches) }

    REPORT_TAXONOMIC_MATCHES(ch_report_inputs)

    emit:
    matches = REPORT_TAXONOMIC_MATCHES.out.matches
    versions = REPORT_TAXONOMIC_MATCHES.out.versions_python.mix(
        REPORT_TAXONOMIC_MATCHES.out.versions_polars,
    )
}
