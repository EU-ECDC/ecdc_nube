
process SEQ_ALIGN_PHYLOGENY{
  container "docker.io/ejfresch/augur:D"
  errorStrategy 'ignore'

  input:
  tuple path(fasta_input_sequences), val(flags), val(project), val(organism), val(batch_id), val(experiment_list), path(current_data, arity: '1..*') 

  output:
  tuple path("alignments/*_aligned.fasta", arity: '1..*'), path("dist_matrices/*.tsv", arity: '1..*'), path("trees/*.treefile", arity: '0..*') 

  publishDir "${params.output}/${project}/", overwrite: true

  tag {"${project}:${batch_id}"}

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
  

  calc_alignments.py -in ${fasta_input_sequences} -r data/${organism}_references.fasta -a data/ -o alignments/
  calc_dist_matrices.py -a alignments/ -d data/ -o dist_matrices/
  calc_trees.py -a alignments/ -o trees/
  """
}


workflow HAVISO {

take:
  data

main:

  ch_current_data = Channel.fromPath('az://iob/HAVISO/**', type: 'file').collect()
  k = data.combine(ch_current_data).map{i ->
       def block_current_data = i[6..-1]
       return([i[0], i[1],
               i[2], i[3],
               i[4], i[5], block_current_data])
     }
  
  k.view()
  SEQ_ALIGN_PHYLOGENY(k.filter{it -> it[5].contains("add_to_alignment")})
}
