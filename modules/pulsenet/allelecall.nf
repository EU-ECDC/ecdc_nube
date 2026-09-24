process PULSENET_ALLELECALLING {
  errorStrategy 'ignore'
  time '30m'
  tag {"${meta.project}:${meta.id}:${meta.schema}"}
  publishDir {"${params.output}/${meta.project}/allele_call_pulsenet/${meta.schema}/"}, overwrite: true

  input:
    tuple val(meta), path(assembly), path(blast_kb), path(qc_kb)

  output:
    tuple val(meta), path("*_outputs.json"), emit: outputs
    tuple val(meta), path("*_stats_calls.json.gz"), emit: stats
    tuple val(meta), path("*_allele_calls.xml.gz"), emit: allele_calls_xml
    tuple val(meta), path("*_allele_calls.json.gz"), emit: allele_calls_json
    tuple val(meta), path("*_allele_calls.bam"), emit: allele_calls_bam
    tuple val(meta), path("*_allele_calls.bam.bai"), emit: allele_calls_bai
    tuple val(meta), path("*_calls_standard.json.gz"), emit: standard_calls
    tuple val(meta), path("*_calls_core_standard.csv.gz"), path("*_calls_core_pcr.csv.gz"), emit: csv_core
    tuple val(meta), path("*_calls_accessory_standard.csv.gz"), path("*_calls_accessory_pcr.csv.gz"), emit: csv_accessory
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
    allele_calls.xml.gz \
    allele_calls.json.gz \
    allele_calls.bam \
    allele_calls.bam.bai \
    calls_standard.json.gz \
    calls_core_standard.csv.gz \
    calls_core_pcr.csv.gz \
    calls_accessory_standard.csv.gz \
    calls_accessory_pcr.csv.gz
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
