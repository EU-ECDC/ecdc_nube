process MIST {
  container "${params.containerRepository}/ejfresch/mist:1.1.0"
  errorStrategy 'ignore'
  time '30m'

  input:
  tuple val(meta), path(assembly), path(schema_path)

  tag {"${meta.project}:${meta.id}:${meta.schema}"}

  publishDir "${params.output}/${meta.project}/allele_call/${meta.schema}/", overwrite: true

  output:
  path "*.tsv", emit: tsv
  path "*.json", emit: json
  path "*.log", emit: log
  path "*_novel_alleles/*", emit: novel_alleles

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call-mist_${meta.schema}"
  """
  mist call \
  --db ${schema_path} \
  --fasta ${assembly} \
  --out-dir . \
  --out-json ${prefix}.json \
  --out-tsv ${prefix}.tsv \
  --threads ${task.cpus} \
  --log ${prefix}.log \
  --export-novel \
  ${args}

  mv novel_alleles/ ${prefix}_novel_alleles/
  """
}
