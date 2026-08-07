#!/usr/bin/env nextflow

include { UTILS_NFSCHEMA_PLUGIN } from './subworkflows/nf-core/utils_nfschema_plugin'
include { BAITS_MAIN } from './workflows/baits'

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

    ch_kmer_size = Channel.value(params.kmer_size)
    ch_deacon_window = Channel.value(params.deacon_window)
    ch_entropy_threshold = Channel.value(params.entropy_threshold)

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
}
