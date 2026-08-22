include { REPORT_CANDIDATE_KMERS } from '../../../modules/local/report_candidate_kmers/main'

workflow REPORT_CANDIDATE_KMER_RESULTS {
    take:
    ch_candidate_manifests
    ch_candidate_occurrences

    main:
    ch_report_rows = ch_candidate_manifests
        .join(
            ch_candidate_occurrences,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, candidate_manifest, candidate_occurrences ->
            tuple(meta.id as String, candidate_manifest, candidate_occurrences)
        }

    ch_sorted_report_rows = ch_report_rows.toSortedList { left, right ->
        left[0] <=> right[0]
    }
    ch_design_ids = ch_sorted_report_rows
        .flatMap { rows -> rows.collect { design_id, _manifest, _occurrences -> design_id } }
        .collectFile(
            name: 'candidate_report_design_ids.txt',
            newLine: true,
            sort: false,
        )
    ch_report_tables = ch_sorted_report_rows.map { rows ->
        tuple(
            rows.collect { _design_id, manifest, _occurrences -> manifest },
            rows.collect { _design_id, _manifest, occurrences -> occurrences },
        )
    }
    ch_report_inputs = ch_report_tables
        .combine(ch_design_ids)
        .map { manifests, occurrences, design_ids -> tuple(design_ids, manifests, occurrences) }

    REPORT_CANDIDATE_KMERS(ch_report_inputs)

    emit:
    candidates = REPORT_CANDIDATE_KMERS.out.candidates
    occurrences = REPORT_CANDIDATE_KMERS.out.occurrences
    versions = REPORT_CANDIDATE_KMERS.out.versions_python.mix(
        REPORT_CANDIDATE_KMERS.out.versions_polars,
    )
}
