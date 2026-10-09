process PULSENET_ALLELECALLING {
  container "${params.containerRepository}/ejfresch/pulsenet_allelecall:1.0.0"
  errorStrategy 'ignore'
  time '2h'
  cpus 3
  tag {"${meta.project}:${meta.id}:${meta.schema}"}
  publishDir {"${params.output}/${meta.project}/allele_call_pulsenet/${meta.schema}/"}, overwrite: true

  input:
    tuple val(meta), path(assembly), path(blast_kb), path(qc_kb)

  output:
    tuple val(meta), path("*_outputs.json"), emit: outputs
    tuple val(meta), path("*_stats_calls.json.gz"), emit: stats
    tuple val(meta), path("*_allele_calls.json.gz"), emit: allele_calls_json
    tuple val(meta), path("*_messages.txt"), emit: log
    tuple val(meta), path("*.tsv"), emit: hashed_tsv

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_pulsenet_${meta.schema}"
  """
  ngs-run AlleleCalling \
  --sample-id ${meta.id} \
  --publish-dir output/ \
  --assembly ${assembly} \
  --blast-kb.path ${blast_kb} \
  --qc-kb.path ${qc_kb} \
  --n-threads ${task.cpus} \
  --organism.genus ${meta.organism} \
  ${args}

  for FILE in \
    outputs.json \
    stats_calls.json.gz \
    allele_calls.json.gz
  do
    mv "\$FILE" "${prefix}_\$FILE"
  done

  mv work/logs/messages.txt ${prefix}_messages.txt

  gzip -dc ${prefix}_allele_calls.json.gz > ${prefix}_allele_calls.json
  hash_pulsenet.py \
  --pulsenet-json ${prefix}_allele_calls.json \
  --output-dir . \
  --prefix ${prefix} \
  --sample ${meta.id}
  """
}
