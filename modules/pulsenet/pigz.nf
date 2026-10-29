process PIGZ_COMPRESS {
  errorStrategy 'terminate'
  tag {"${meta.project}:${meta.id}"}
  cpus 1

  input:
  tuple val(meta), path(raw_file)

  output:
  tuple val(meta), path("$gz_file"), emit: compressed

  script:
  def args = task.ext.args ?: ''
  gz_file = raw_file.toString() + ".gz"
  """
  pigz --processes $task.cpus --stdout --force ${args} ${raw_file} > ${gz_file}
  """
}
