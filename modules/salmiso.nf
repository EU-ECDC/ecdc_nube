process AMRFINDER {
  container "${params.containerRepository}/ncbi/amr:4.0.3-2024-10-22.1"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(project), val(accession), val(sequencing_technology), path(assembly), val(organism), val(experiment_list)

  tag {"${project}:${accession}"}

  publishDir "az://iob/${project}/amr/", overwrite: true

  output:
  tuple val(project), val(accession), path("${accession}_amrfinder.tsv"), val(organism), val(experiment_list)
 
  shell:
  """
  amrfinder -n ${assembly} -q --plus --organism Salmonella --output ${accession}_amrfinder.tsv 
  """
}

workflow SALMISO {
take:
  data

main:
  AMRFINDER(data.filter{it -> it[5].contains("amrfinder")})

}
