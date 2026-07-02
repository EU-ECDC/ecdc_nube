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

  script:
  def prefix = task.ext.prefix ?: "${meta.id}_allele-call_${meta.schema}"
  if(meta.sequencing_technology == "ILLUMINA"){
    def advOpt_illumina = advOptions.containsKey("ILLUMINA") ? advOptions.ILLUMINA : "None"
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_illumina}" -f ${prefix}
    """
  } else if(meta.sequencing_technology == "IONTORRENT"){
    def advOpt_iontorrent = advOptions.containsKey("IONTORRENT") ? advOptions.IONTORRENT : "None"
    def advOpt_chewBBACA = advOpt_iontorrent.chewBBACA
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_chewBBACA}" -f ${prefix}
    """
  } else if(meta.sequencing_technology == "ONT"){
    def advOpt_ont = advOptions.containsKey("ONT") ? advOptions.ONT : "None"
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_ont}" -f ${prefix}
    """
  } else if(meta.sequencing_technology == "PACBIO"){
    def advOpt_pacbio = advOptions.containsKey("PACBIO") ? advOptions.PACBIO : "None"
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_pacbio}" -f ${prefix}
    """
  }
}


process TARANYS {
  container "${params.containerRepository}/ejfresch/taranys:3.0.1-d"
  errorStrategy 'ignore'
  time '2h'
  cpus 4
  tag {"${meta.project}:${meta.id}:${schema}"}
  publishDir {"${params.output}/${meta.project}/taranys/${schema}/"}, overwrite: true, pattern: "*.tsv"

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
  """
}