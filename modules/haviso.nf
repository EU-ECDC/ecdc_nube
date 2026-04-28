process SEQ_ALIGN_PHYLOGENY{
  container "docker.io/ejfresch/augur:D"
  errorStrategy 'ignore'

  input:
  tuple val(payLoad), path(fasta_input_sequences), path(current_data, arity: '1..*')

  output:
  tuple path("alignments/*_aligned.fasta", arity: '1..*'), path("dist_matrices/*.tsv", arity: '1..*'), path("trees/*.treefile", arity: '0..*') 

  publishDir "${params.output}/${payLoad.project}/", overwrite: true

  tag {"${payLoad.project}:${payLoad.id}"}

  script:
  """
  mkdir -p data/
  # all already existing alignments and distance
  # matrices will be moved to data/
  shopt -s nullglob
  files=( *_references.fasta *_aligned.fasta *_dist.tsv )
  if (( \${#files[@]} )); then
    mv "\${files[@]}" data/
  fi
  

  calc_alignments.py -in ${fasta_input_sequences} -r data/${payLoad.organism}_references.fasta -a data/ -o alignments/ -f ${payLoad.flags}
  calc_dist_matrices.py -a alignments/ -d data/ -o dist_matrices/
  calc_trees.py -a alignments/ -o trees/
  """
}


workflow HAVISO {

take:
  data

main:
  ch_current_alignment = Channel.fromPath("${params.phylogenyAlignments}/HAVISO/**", type: 'file').collect()
  SEQ_ALIGN_PHYLOGENY(data
    .filter{payLoad, sequences -> payLoad.experiment_list.contains("add_to_alignment")}
    .combine(ch_current_alignment).map{i ->
       def files_current_alignment = i[2..-1]
       return([i[0], i[1], files_current_alignment])
    }
  )
}
