process ALLELE_CALL {
  container "${params.containerRepository}/ejfresch/chewbbaca:3.3.10"
  errorStrategy 'ignore'
  time '30m'

  input:
    tuple val(meta), path(assembly), path(schema_path), path(trn_file), path(gene_list), path(fasta_ref_seqs_alleles), val(advOptions), val(prefix), val(schema)

  tag {"${meta.project}:${meta.id}:${schema}"}

  publishDir "${params.output}/${meta.project}/allele_call/${schema}/", overwrite: true

  output:
    path "*.tsv", emit: tsv_chewbbaca

  script:
  if(meta.sequencing_technology == "ILLUMINA"){
    def advOpt_illumina = advOptions.containsKey("ILLUMINA") ? advOptions.ILLUMINA : "None"
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_illumina}" -f ${prefix}
    """
  } else if(meta.sequencing_technology == "IONTORRENT"){
    def advOpt_iontorrent = advOptions.containsKey("IONTORRENT") ? advOptions.IONTORRENT : "None"
    def advOpt_chewBBACA = advOpt_iontorrent.chewBBACA
    """
    iontorrent_error_correction.py ${assembly} ${fasta_ref_seqs_alleles} corrected_assembly.fasta
    generate_cgMLST_profiles.py -i corrected_assembly.fasta -a ${meta.id} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_chewBBACA}" -f ${prefix}
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
