
process QUAST {
  container "${params.containerRepository}/ejfresch/quast:5.2.0"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(project), val(accession), val(sequencing_technology), path(assembly), val(organism), val(experiment_list)

  tag {"${project}:${accession}"}

  publishDir "az://iob/${project}/qc/", overwrite: true

  output:
  tuple val(project), val(accession), path("${accession}_quast.tsv", arity: '1..*'), val(organism), val(experiment_list)
 
  shell:
  """
  quast.py -o . --no-plots --no-html ${assembly}
  ln transposed_report.tsv ${accession}_quast.tsv
  """
}


workflow QC {
take:
  data

main:
  QUAST(data)

emit:
  QUAST.out
}
