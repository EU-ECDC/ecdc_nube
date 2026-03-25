process CG_MLST {
  container "${params.containerRepository}/ejfresch/chewbbaca:3.3.10"
  errorStrategy 'terminate'
  time '30m'

  input:
    tuple val(project), val(accession), val(technology), path(assembly), val(organism), val(experiment_list), path(schema_path), path(trn_file), path(gene_list), path(fasta_ref_seqs_alleles), val(advOptions), val(prefix)

  tag {"${project}:${accession}"}

  publishDir "${params.output}/${project}/cgMLST/", overwrite: true

  output:
    path "*_cgMLST.tsv", emit: tsv_chewbbaca

  script:
  if(technology == "ILLUMINA"){
    def advOpt_illumina = advOptions.containsKey("ILLUMINA") ? advOptions.ILLUMINA : "None" 
    """
    echo generate_cgMLST_profiles.py -i ${assembly} -a ${accession} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_illumina}" -f ${prefix}
    touch ${prefix}_cgMLST.tsv
    """
  } else if(technology == "IONTORRENT"){
    def advOpt_iontorrent = advOptions.containsKey("IONTORRENT") ? advOptions.IONTORRENT : "None"
    def advOpt_chewBBACA = advOpt_iontorrent.chewBBACA
    """
    iontorrent_error_correction.py ${assembly} ${fasta_ref_seqs_alleles} corrected_assembly.fasta
    generate_cgMLST_profiles.py -i corrected_assembly.fasta -a ${accession} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_chewBBACA}" -f ${prefix}    
    """
  } else if(technology == "ONT"){
    def advOpt_ont = advOptions.containsKey("ONT") ? advOptions.ONT : "None"
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${accession} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_ont}" -f ${prefix}
    """
  } else if(technology == "PACBIO"){
    def advOpt_pacbio = advOptions.containsKey("PACBIO") ? advOptions.PACBIO : "None"
    """
    generate_cgMLST_profiles.py -i ${assembly} -a ${accession} -g ${schema_path} -p ${trn_file} -gl ${gene_list} -ao "${advOpt_pacbio}" -f ${prefix}
    """
  }
}