include { REPORT_CALIBRATION_RESULTS } from '../../../modules/local/report_calibration_results/main'
include { REPORT_CALIBRATION_THRESHOLDS } from '../../../modules/local/report_calibration_thresholds/main'

workflow REPORT_SHARED_CALIBRATION_RESULTS {
    take:
    ch_calibration_reads
    ch_calibration_matches
    ch_calibration_designs
    ch_recommended_thresholds

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

    ch_requested_design_ids = ch_calibration_designs
        .map { meta -> meta.id as String }
        .toSortedList()
        .map { design_ids ->
            tuple(
                'calibration_design_ids.txt',
                design_ids.join('\n') + (design_ids ? '\n' : ''),
            )
        }
        .collectFile(sort: false) { filename, contents -> [filename, contents] }
    ch_sorted_thresholds = ch_recommended_thresholds
        .map { meta, threshold -> tuple(meta.id as String, threshold) }
        .toSortedList { left, right -> left[0] <=> right[0] }
        .map { rows ->
            [thresholds: rows.collect { _design_id, threshold -> threshold }]
        }
    ch_threshold_inputs = ch_sorted_thresholds
        .combine(ch_requested_design_ids)
        .map { threshold_group, design_ids ->
            tuple(design_ids, threshold_group.thresholds)
        }
    REPORT_CALIBRATION_THRESHOLDS(ch_threshold_inputs)

    emit:
    reads = REPORT_CALIBRATION_RESULTS.out.reads
    matches = REPORT_CALIBRATION_RESULTS.out.matches
    thresholds = REPORT_CALIBRATION_THRESHOLDS.out.thresholds
    versions = REPORT_CALIBRATION_RESULTS.out.versions_python.mix(
        REPORT_CALIBRATION_RESULTS.out.versions_polars,
    ).mix(REPORT_CALIBRATION_THRESHOLDS.out.versions_python)
}
