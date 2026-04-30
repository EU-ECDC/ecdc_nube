process QUAST {
  container "${params.containerRepository}/ejfresch/quast:5.2.0"
  errorStrategy 'ignore'
  time '10m'

  input:
  tuple val(meta), path(assembly)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.project}/qc/", overwrite: true

  output:
  tuple val(meta), path("${meta.id}_quast.tsv", arity: '1..*')

  shell:
  """
  quast.py -o . --no-plots --no-html ${assembly}
  ln transposed_report.tsv ${meta.id}_quast.tsv
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
