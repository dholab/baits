include { DEACON_INDEX_BUILD as DEACON_INDEX_BAIT_SET } from '../../../modules/local/deacon_index_build/main'
include { DEACON_FILTER_FASTA as DEACON_BAIT_ROUNDTRIP } from '../../../modules/local/deacon_filter_fasta/main'
include { DEACON_FILTER_FASTA as DEACON_BACKGROUND_ROUNDTRIP } from '../../../modules/local/deacon_filter_fasta/main'
include { VERIFY_DEACON_INDEX } from '../../../modules/local/verify_deacon_index/main'


workflow BUILD_VERIFY_DEACON_INDEX {
    take:
    ch_index_inputs
    ch_kmer_size
    ch_deacon_window

    main:

    // Build the operational Deacon Index from the selected Bait Set
    ch_index_build_inputs = ch_index_inputs
        .combine(ch_kmer_size)
        .combine(ch_deacon_window)
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, interference_background, kmer_size, deacon_window ->
            tuple(meta, baits, kmer_size, deacon_window, 0)
        }
    DEACON_INDEX_BAIT_SET(ch_index_build_inputs)

    // Round-trip the Bait Set and Interference Background at the permissive threshold
    ch_bait_roundtrip_inputs = DEACON_INDEX_BAIT_SET.out.index
        .join(ch_index_inputs)
        .map { meta, deacon_index, baits, candidate_kmer_manifest, bait_set_status, interference_background ->
            tuple(meta, deacon_index, baits)
        }
    ch_background_roundtrip_inputs = DEACON_INDEX_BAIT_SET.out.index
        .join(ch_index_inputs)
        .map { meta, deacon_index, baits, candidate_kmer_manifest, bait_set_status, interference_background ->
            tuple(meta, deacon_index, interference_background)
        }
    DEACON_BAIT_ROUNDTRIP(ch_bait_roundtrip_inputs)
    DEACON_BACKGROUND_ROUNDTRIP(ch_background_roundtrip_inputs)

    // Verify both round trips and finalize the Bait Set status
    ch_verification_inputs = ch_index_inputs
        .join(DEACON_BAIT_ROUNDTRIP.out.fasta)
        .join(DEACON_BACKGROUND_ROUNDTRIP.out.fasta)
        .combine(ch_kmer_size)
        .combine(ch_deacon_window)
        .map { meta, baits, candidate_kmer_manifest, bait_set_status, interference_background, bait_roundtrip, background_roundtrip, kmer_size, deacon_window ->
            tuple(meta, baits, bait_roundtrip, background_roundtrip, candidate_kmer_manifest, bait_set_status, kmer_size, deacon_window)
        }
    VERIFY_DEACON_INDEX(ch_verification_inputs)

    ch_verified_index = ch_index_inputs
        .join(DEACON_INDEX_BAIT_SET.out.index)
        .join(VERIFY_DEACON_INDEX.out.bait_set_status)
        .map { meta, baits, candidate_kmer_manifest, bait_set_status_draft, interference_background, deacon_index, bait_set_status ->
            tuple(meta, baits, candidate_kmer_manifest, bait_set_status, deacon_index)
        }

    emit:
    index = ch_verified_index
    bait_set_status = VERIFY_DEACON_INDEX.out.bait_set_status
    summary = VERIFY_DEACON_INDEX.out.summary
    report = VERIFY_DEACON_INDEX.out.report
    reported_deacon = DEACON_INDEX_BAIT_SET.out.reported_deacon
    versions = DEACON_INDEX_BAIT_SET.out.versions_deacon
        .mix(DEACON_BAIT_ROUNDTRIP.out.versions_deacon)
        .mix(DEACON_BACKGROUND_ROUNDTRIP.out.versions_deacon)
        .mix(VERIFY_DEACON_INDEX.out.versions_biopython)
        .mix(VERIFY_DEACON_INDEX.out.versions_polars)
        .mix(VERIFY_DEACON_INDEX.out.versions_python)
}
