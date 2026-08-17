process MERYL_INTERSECT {
  tag "${meta.id}"
  label 'process_medium'

  conda "${moduleDir}/environment.yml"
  container 'quay.io/biocontainers/meryl:1.4.1--h4ac6f70_0@sha256:60ba02cde408b606fc1834ef3261c5abc33796d39bf9640dcee307c256501093'

  input:
  tuple val(meta), path(background_db, stageAs: 'background.meryl'), path(source_db, stageAs: 'source.meryl')

  output:
  tuple val(meta), path('intersection.meryl'), emit: db
  tuple val("${task.process}"), val('meryl'), eval("meryl --version |& sed 's/meryl //'"), topic: versions, emit: versions_meryl
  tuple val(meta), val('meryl'), eval("meryl --version |& sed 's/meryl //'"), emit: reported_meryl

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    meryl intersect \
        ${background_db} \
        ${source_db} \
        output intersection.meryl
    """

  stub:
  """
    mkdir intersection.meryl
    """
}
