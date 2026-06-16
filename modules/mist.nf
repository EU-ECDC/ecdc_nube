process MIST {
  container "docker.io/ejfresch/mist:1.1.0"
  errorStrategy 'ignore'
  time '30m'
  tag {"${meta.project}:${meta.id}:${meta.schema}"}
  publishDir "${params.output}/${meta.project}/allele_call/${meta.schema}/", overwrite: true

  input:
  tuple val(meta), path(assembly), path(schema_path)

  output:
  path "*.tsv", emit: tsv_mist
  path "*.json", emit: json
  path "*.log", emit: log
  path "*_novel_alleles/*", emit: novel_alleles

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call-mist_${meta.schema}"
  def prefix2 = task.ext.prefix2 ?: "${meta.id}_allele-call-mist_${meta.schema}"
  """
  mist call \
  --db ${schema_path} \
  --fasta ${assembly} \
  --out-dir . \
  --out-json ${prefix}.json \
  --threads ${task.cpus} \
  --log ${prefix}.log \
  --export-novel \
  ${args}
  mv novel_alleles/ ${prefix}_novel_alleles/

  hash_mist.py \
  --mist-json ${prefix}.json \
  --schema-dir ${schema_path} \
  --output . \
  --prefix ${prefix2} \
  --sample ${meta.id}
  """
}
