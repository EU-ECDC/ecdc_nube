process AMRFINDER {
  container "${params.containerRepository}/ncbi/amr:4.0.3-2024-10-22.1"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(payLoad), path(assembly)

  tag {"${payLoad.project}:${payLoad.id}"}

  publishDir "${params.output}/${payLoad.project}/amr/", overwrite: true

  output:
  tuple val(payLoad), path("${payLoad.id}_amrfinder.tsv")
 
  shell:
  """
  amrfinder -n ${assembly} -q --plus --organism Salmonella --output ${payLoad.id}_amrfinder.tsv
  """
}

workflow SALMISO {
take:
  data

main:
  AMRFINDER(data.filter{payLoad, assembly -> payLoad.experiment_list.contains("amrfinder")})

}
