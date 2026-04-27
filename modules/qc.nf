process QUAST {
  container "${params.containerRepository}/ejfresch/quast:5.2.0"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(payLoad), path(assembly)

  tag {"${payLoad.project}:${payLoad.id}"}

  publishDir "${params.output}/${payLoad.project}/qc/", overwrite: true

  output:
  tuple val(payLoad), path("${payLoad.id}_quast.tsv", arity: '1..*')

  shell:
  """
  quast.py -o . --no-plots --no-html ${assembly}
  ln transposed_report.tsv ${payLoad.id}_quast.tsv
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
