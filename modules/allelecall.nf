process ALLELE_CALL {
  container "${params.containerRepository}/ejfresch/chewbbaca:3.3.10"
  errorStrategy 'ignore'
  time '30m'

  input:
    tuple val(meta), path(assembly), path(schema_path), path(trn_file), path(gene_list), val(advOptions)

  tag {"${meta.project}:${meta.id}:${meta.schema}"}

  publishDir "${params.output}/${meta.project}/allele_call/${meta.schema}/", overwrite: true

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
