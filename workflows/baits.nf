include { NORMALIZE_CURATED_SOURCE_SEQUENCES } from '../modules/local/normalize_curated_source_sequences/main'
include { WRITE_PROVENANCE } from '../modules/local/write_provenance/main'
include { DISCOVER_SOURCE_SEQUENCES } from '../subworkflows/local/discover_source_sequences/main'

workflow BAITS_MAIN {
    take:
    ch_design_context
    ch_curated_inputs
    ch_discovery_inputs
    ch_curated_provenance_source_inputs
    ch_discovery_provenance_source_inputs

    main:

    // Source Sequence acquisition
    NORMALIZE_CURATED_SOURCE_SEQUENCES(ch_curated_inputs)
    DISCOVER_SOURCE_SEQUENCES(ch_discovery_inputs)

    // Converge curated and discovered Source Sequences
    ch_source_sequences = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.source_sequences
        .mix(DISCOVER_SOURCE_SEQUENCES.out.source_sequences)

    // Provenance
    ch_curated_provenance_facts = ch_design_context
        .join(ch_curated_provenance_source_inputs)
        .join(NORMALIZE_CURATED_SOURCE_SEQUENCES.out.reported_biopython)
        .map { meta, target_taxid, interference_background, source_input_roles, source_input_files, component, version ->
            def inputFacts = [
                [input_role: 'design', input_id: meta.id, attribute: 'source_sequence_origin', value: meta.source_sequence_origin],
                [input_role: 'target_taxon', input_id: meta.id, attribute: 'target_taxid', value: target_taxid],
            ]
            def parameters = [[parameter: 'target_taxid', value: target_taxid]]
            def softwareVersions = [[component: component, version: version]]
            tuple(
                meta,
                inputFacts,
                ['interference_background'] + source_input_roles,
                ['interference_background'] + source_input_roles,
                ['file'] * (source_input_files.size() + 1),
                [interference_background] + source_input_files,
                parameters,
                softwareVersions,
            )
        }

    ch_discovery_versions = DISCOVER_SOURCE_SEQUENCES.out.reported_biopython
        .join(DISCOVER_SOURCE_SEQUENCES.out.reported_blast)

    ch_discovery_provenance_facts = ch_design_context
        .join(ch_discovery_provenance_source_inputs)
        .join(ch_discovery_versions)
        .map { meta, target_taxid, interference_background, source_input_roles, source_input_files, biopython_component, biopython_version, blast_component, blast_version ->
            def inputFacts = [
                [input_role: 'design', input_id: meta.id, attribute: 'source_sequence_origin', value: meta.source_sequence_origin],
                [input_role: 'target_taxon', input_id: meta.id, attribute: 'target_taxid', value: target_taxid],
            ]
            def parameters = [[parameter: 'target_taxid', value: target_taxid]]
            def softwareVersions = [
                [component: biopython_component, version: biopython_version],
                [component: blast_component, version: blast_version],
            ]
            tuple(
                meta,
                inputFacts,
                ['interference_background'] + source_input_roles,
                ['interference_background'] + source_input_roles,
                ['file'] * (source_input_files.size() + 1),
                [interference_background] + source_input_files,
                parameters,
                softwareVersions,
            )
        }

    WRITE_PROVENANCE(ch_curated_provenance_facts.mix(ch_discovery_provenance_facts))

    emit:
    source_sequences = ch_source_sequences
    candidate_loci = DISCOVER_SOURCE_SEQUENCES.out.candidate_loci
    query_blast_hits = DISCOVER_SOURCE_SEQUENCES.out.blast_hits
    discovery_status = DISCOVER_SOURCE_SEQUENCES.out.discovery_status
    provenance = WRITE_PROVENANCE.out.provenance
    versions = NORMALIZE_CURATED_SOURCE_SEQUENCES.out.versions_biopython
        .mix(DISCOVER_SOURCE_SEQUENCES.out.versions)
        .mix(WRITE_PROVENANCE.out.versions_python)
}
