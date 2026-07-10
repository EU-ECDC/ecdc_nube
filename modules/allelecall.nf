process ALLELE_CALL {
  container "${params.containerRepository}/ejfresch/chewbbaca:3.3.10"
  errorStrategy 'ignore'
  time '30m'
  tag {"${meta.project}:${meta.id}:${meta.schema}"}
  publishDir {"${params.output}/${meta.project}/allele_call/${meta.schema}/"}, overwrite: true

  input:
    tuple val(meta), path(assembly), path(schema_path), path(trn_file), path(gene_list), val(advOptions)

  output:
    path "*.tsv", emit: tsv_chewbbaca
    path "additional_results/*_novel_alleles.fasta", emit: novel_alleles
    path "additional_results/*_results_alleles.tsv", emit: raw_results
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call_${meta.schema}"
  def advOpt = advOptions[meta.sequencing_technology]?.chewBBACA ?: "None"
  """
  generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt}" -f ${prefix}

  mkdir additional_results
  mv results_*/novel_alleles.fasta additional_results/${prefix}_novel_alleles.fasta
  mv results_*/results_alleles.tsv additional_results/${prefix}_results_alleles.tsv
  """
}


process TARANYS {
  container "${params.containerRepository}/ejfresch/taranys:3.0.1-d"
  errorStrategy 'ignore'
  time '2h'
  cpus 4
  tag {"${meta.project}:${meta.id}:${schema}"}
  publishDir {"${params.output}/${meta.project}/taranys/${schema}/"}, overwrite: true

  input:
    tuple val(meta),
          path(assembly),
          path(schema_path),
          path(reference_alleles_path),
          path(annotation_file),
          val(advOptions),
          val(schema)

  output:
    tuple val(meta), val(schema), path("${meta.id}_taranys_${schema}.tsv"), emit: hashed_tsv
    tuple val(meta), val(schema), path("additional_results/*_contig_alignment_info.csv"), emit: contig_alignment_info
    tuple val(meta), val(schema), path("additional_results/*_allele_calling_match.csv"), emit: allele_calling_match

  script:
  """
  # this is to ensure naming consistency in the hash step
  if [ "\$(basename ${assembly})" != "${meta.id}.fasta" ]; then
    ln -sf ${assembly} ${meta.id}.fasta
  fi

  taranys allele-calling \\
    -s ${schema_path} \\
    -r ${reference_alleles_path} \\
    -a ${annotation_file} \\
    --cpus ${task.cpus} \\
    -o out \\
    ${meta.id}.fasta

  hash_taranys_calls.py \\
    --schema-dir ${schema_path} \\
    --contig-alignment out/contig_alignment_info.csv \\
    --match-csv out/allele_calling_match.csv \\
    --sample ${meta.id} \\
    --output ${meta.id}_taranys_${schema}.tsv

  mkdir additional_results
  mv out/contig_alignment_info.csv additional_results/${meta.id}_contig_alignment_info.csv
  mv out/allele_calling_match.csv additional_results/${meta.id}_allele_calling_match.csv
  """
}
