process IONTORRENT_ERROR_CORRECTION {
  container "${params.containerRepository}/ejfresch/chewbbaca:3.3.10"
  errorStrategy 'ignore'
  time '30m'

  input:
    tuple val(meta), path(assembly), path(fasta_ref_seqs_alleles)

  tag {"${meta.project}:${meta.id}"}

  output:
    tuple val(meta), path("${meta.id}_corrected_assembly.fasta")

  script:
    """
    iontorrent_error_correction.py ${assembly} ${fasta_ref_seqs_alleles} ${meta.id}_corrected_assembly.fasta
    """
}
