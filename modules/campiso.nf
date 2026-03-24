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
  amrfinder -n ${assembly} -q --plus --organism Campylobacter --output ${accession}_amrfinder.tsv 
  """
}

process MLST_CGE {
  container "${params.containerRepository}/ejfresch/typing:1.0"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(project), val(accession), val(sequencing_technology), path(assembly), val(organism), val(experiment_list), path(path_to_mlst_schemes)

  tag {"${project}:${accession}"}

  publishDir "az://iob/${project}/MLST/", overwrite: true

  output:
  tuple val(project), val(accession), path("${accession}_mlst.json"), val(organism), val(experiment_list)
 
  shell:
  '''
  mlst.py -i !{assembly} -s cjejuni -p !{path_to_mlst_schemes} -o .
  ln data.json !{accession}_mlst.json
  '''
}


workflow CAMPISO {
take:
  data

main:
  AMRFINDER(data.filter{it -> it[5].contains("amrfinder")})
  MLST_CGE(data.filter{it -> it[5].contains("mlst")}.map{
    project, accession, technology, assembly, organism, experiment_list ->
    [project, accession, technology, assembly, organism, experiment_list, "az://iob/schemas/MLST/"]
  })

}
