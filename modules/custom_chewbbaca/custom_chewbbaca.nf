process CUSTOM_CHEWBBACA_SSI {
  container "docker.io/ejfresch/custom_chewbbaca_ssi:1.0"
  errorStrategy 'terminate'
  time '30m'
  
  input:
    tuple val(meta), path(assembly), path(schema_path)

  tag {"${meta.project}:${meta.id}:${meta.schema}"}

  publishDir "${params.output}/${meta.project}/allele_call/${meta.schema}/", overwrite: true

  output:
    path "*.tsv", emit: tsv_chewbbaca_custom

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call-SSI_${meta.schema}"
  """
  mv ${assembly} sample.fa
  chewBBACA.py AlleleCall -i . -g ${schema_path} -o . --no-inferred --hash-profiles crc32 --cds ${args}
  file=\$(ls results_*/results_alleles_hashed.tsv)
  mv "\$file" ${prefix}.tsv
  """
}
