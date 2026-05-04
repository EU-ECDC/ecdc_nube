process AMRFINDER {
  container "${params.containerRepository}/ncbi/amr:4.0.3-2024-10-22.1"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(meta), path(assembly)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.project}/amr/", overwrite: true

  output:
  tuple val(meta), path("${meta.id}_amrfinder.tsv")
 
  shell:
  """
  amrfinder -n ${assembly} -q --plus --organism Escherichia --output ${meta.id}_amrfinder.tsv
  """
}


workflow ECOLIISO {
take:
  data

main:
  AMRFINDER(data.filter{meta, assembly -> meta.experiment_list.contains("amrfinder")})
}
