process WRITE_PROVENANCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_polars_python:5f22dbb8160be6fa@sha256:20d95779744ec0d8a943473966c3403a9d4b09f3d53722ffd67169972ab1b2da'

    input:
    tuple val(meta), val(input_facts), val(input_file_roles), val(input_file_ids), val(input_file_kinds), path(input_files, stageAs: 'provenance_inputs/file????/*'), val(parameters), val(software_versions)

    output:
    tuple val(meta), path('inputs.tsv'), path('parameters.tsv'), path('software_versions.tsv'), emit: provenance
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def inputFactsBase64 = groovy.json.JsonOutput.toJson(input_facts).bytes.encodeBase64().toString()
    def inputFileRolesBase64 = groovy.json.JsonOutput.toJson(input_file_roles).bytes.encodeBase64().toString()
    def inputFileIdsBase64 = groovy.json.JsonOutput.toJson(input_file_ids).bytes.encodeBase64().toString()
    def inputFileKindsBase64 = groovy.json.JsonOutput.toJson(input_file_kinds).bytes.encodeBase64().toString()
    def parametersBase64 = groovy.json.JsonOutput.toJson(parameters).bytes.encodeBase64().toString()
    def softwareVersionsBase64 = groovy.json.JsonOutput.toJson(software_versions).bytes.encodeBase64().toString()
    """
    write_provenance.py \
        --input-facts-base64 ${inputFactsBase64} \
        --input-file-roles-base64 ${inputFileRolesBase64} \
        --input-file-ids-base64 ${inputFileIdsBase64} \
        --input-file-kinds-base64 ${inputFileKindsBase64} \
        --input-root provenance_inputs \
        --parameters-base64 ${parametersBase64} \
        --software-versions-base64 ${softwareVersionsBase64} \
        --inputs-out inputs.tsv \
        --parameters-out parameters.tsv \
        --software-versions-out software_versions.tsv
    """
}
