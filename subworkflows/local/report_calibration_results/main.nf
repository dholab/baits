include { REPORT_CALIBRATION_RESULTS } from '../../../modules/local/report_calibration_results/main'

workflow REPORT_SHARED_CALIBRATION_RESULTS {
    take:
    ch_calibration_reads
    ch_calibration_matches

    main:
    ch_report_rows = ch_calibration_reads
        .join(
            ch_calibration_matches,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, calibration_reads, calibration_matches ->
            tuple(meta.id as String, calibration_reads, calibration_matches)
        }
    ch_sorted_report_rows = ch_report_rows.toSortedList { left, right ->
        left[0] <=> right[0]
    }
    ch_design_ids = ch_sorted_report_rows
        .map { rows ->
            tuple(
                'calibration_design_ids.txt',
                rows.collect { design_id, _reads, _matches -> design_id }.join('\n') + (rows ? '\n' : ''),
            )
        }
        .collectFile(sort: false) { filename, contents -> [filename, contents] }
    ch_report_tables = ch_sorted_report_rows.map { rows ->
        tuple(
            rows.collect { _design_id, reads, _matches -> reads },
            rows.collect { _design_id, _reads, matches -> matches },
        )
    }
    ch_report_inputs = ch_report_tables
        .combine(ch_design_ids)
        .map { reads, matches, design_ids -> tuple(design_ids, reads, matches) }

    REPORT_CALIBRATION_RESULTS(ch_report_inputs)

    emit:
    reads = REPORT_CALIBRATION_RESULTS.out.reads
    matches = REPORT_CALIBRATION_RESULTS.out.matches
    versions = REPORT_CALIBRATION_RESULTS.out.versions_python.mix(
        REPORT_CALIBRATION_RESULTS.out.versions_polars,
    )
}
