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
  amrfinder -n ${assembly} -q --plus --organism Campylobacter --output ${meta.id}_amrfinder.tsv 
  """
}

process MLST_CGE {
  errorStrategy 'ignore'
  time '10m'
  tag {"${meta.project}:${meta.id}"}
  publishDir {"${params.output}/${meta.project}/MLST/"}, overwrite: true

  input:
  tuple val(meta), path(assembly), path(path_to_mlst_schemes)

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
