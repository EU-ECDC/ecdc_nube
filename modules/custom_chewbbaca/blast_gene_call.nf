process BLAST_GENE_CALL_SSI {
  container "${params.containerRepository}/ejfresch/allele_blast_ssi:1.0"
  errorStrategy 'ignore'
  tag {"${meta.project}:${meta.id}:${meta.schema}"}
  publishDir "${params.output}/${meta.project}/blast/${meta.schema}/", overwrite: true
  time '2h'
  cpus 3

  input:
    tuple val(meta), path(assembly), path(schema_path)

  output:
    tuple val(meta), path("*.fa")

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  blast_gene_call_per_locus_moreChecks_optimized2_skipNonCDS.py --scheme ${schema_path} --fa ${assembly} --out . --max_workers ${task.cpus} ${args}
  """
}
