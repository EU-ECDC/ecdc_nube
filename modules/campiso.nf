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
  amrfinder -n ${assembly} -q --plus --organism Campylobacter --output ${meta.id}_amrfinder.tsv 
  """
}

process MLST_CGE {
  container "${params.containerRepository}/ejfresch/typing:1.0"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(meta), path(assembly), path(path_to_mlst_schemes)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.project}/MLST/", overwrite: true

  output:
  tuple val(meta), path("${meta.id}_mlst.json")
 
  shell:
  '''
  mlst.py -i !{assembly} -s cjejuni -p !{path_to_mlst_schemes} -o .
  ln data.json !{meta.id}_mlst.json
  '''
}


workflow CAMPISO {
take:
  data

main:
  AMRFINDER(data.filter{meta, assembly -> meta.experiment_list.contains("amrfinder")})
  MLST_CGE(data.filter{meta, assembly -> meta.experiment_list.contains("mlst")}.map{
    meta, assembly ->
    [meta, assembly, params.mlstSchemas]
  })

}
