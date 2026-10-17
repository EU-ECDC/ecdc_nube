process AMRFINDER {
  errorStrategy 'ignore'
  time '10m'
  tag {"${meta.project}:${meta.id}"}
  publishDir {"${params.output}/${meta.project}/amr/"}, overwrite: true

  input:
  tuple val(meta), path(assembly)

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
