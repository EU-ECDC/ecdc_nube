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
  amrfinder -n ${assembly} -q --plus --organism Campylobacter --output ${payLoad.id}_amrfinder.tsv 
  """
}

process MLST_CGE {
  container "${params.containerRepository}/ejfresch/typing:1.0"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(payLoad), path(assembly), path(path_to_mlst_schemes)

  tag {"${payLoad.project}:${payLoad.id}"}

  publishDir "${params.output}/${payLoad.project}/MLST/", overwrite: true

  output:
  tuple val(payLoad), path("${payLoad.id}_mlst.json")
 
  shell:
  '''
  mlst.py -i !{assembly} -s cjejuni -p !{path_to_mlst_schemes} -o .
  ln data.json !{payLoad.id}_mlst.json
  '''
}


workflow CAMPISO {
take:
  data

main:
  AMRFINDER(data.filter{payLoad, assembly -> payLoad.experiment_list.contains("amrfinder")})
  MLST_CGE(data.filter{payLoad, assembly -> payLoad.experiment_list.contains("mlst")}.map{
    payLoad, assembly ->
    [payLoad, assembly, "az://iob/schemas/MLST/"]
  })

}
