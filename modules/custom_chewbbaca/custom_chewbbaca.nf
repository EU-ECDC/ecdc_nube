process CUSTOM_CHEWBBACA_SSI {
  container "${params.containerRepository}/ejfresch/custom_chewbbaca_ssi:1.0"
  errorStrategy 'ignore'
  tag {"${meta.project}:${meta.id}:${meta.schema}"}
  publishDir "${params.output}/${meta.project}/allele_call/${meta.schema}/", overwrite: true
  time '30m'
  cpus 3

  input:
    tuple val(meta), path(assembly), path(schema_path), path(gene_list), val(advOptions)

  output:
    path "*.tsv", emit: tsv_chewbbaca_custom

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call-SSI_${meta.schema}"
  def advOpt_chewBBACA = advOptions.get(meta.sequencing_technology, [:]).chewBBACA ?: ''
  """
  generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -gl ${gene_list} -f ${prefix} -ao="--cds --cpu-cores ${task.cpus} ${advOpt_chewBBACA} ${args}"
  """
}
