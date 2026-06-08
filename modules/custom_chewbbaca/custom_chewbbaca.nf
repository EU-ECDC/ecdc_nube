process CUSTOM_CHEWBBACA_SSI {
  container "docker.io/ejfresch/custom_chewbbaca_ssi:1.0"
  errorStrategy 'ignore'
  time '30m'

  input:
    tuple val(meta), path(assembly), path(schema_path), path(gene_list)

  tag {"${meta.project}:${meta.id}:${meta.schema}"}

  publishDir "${params.output}/${meta.project}/allele_call/${meta.schema}/", overwrite: true

  output:
    path "*.tsv", emit: tsv_chewbbaca_custom

  script:
  def args = task.ext.args ?: '--cds'
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call-SSI_${meta.schema}"
  """
  generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -gl ${gene_list} -f ${prefix} -ao="${args}"
  """
}
