process correctAssembly {
  container '<your_image_repository>:tag'

  memory 8.GB
  cpus 4 

  publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/correctAssembly", mode: params.output_mode_in_use

  input:  
      tuple val(id), val(sample_id), val(organism), path(read1), path(read2), path(assembly)
      path(cleaning_kb)
      path(qc_kb)

  output:
      tuple val(id), path("outputs.json"), emit: outputs
      tuple val(id), path("corrected_assembly.fasta.gz"), emit: assembly
      tuple val(id), path("corrected_alignment.cram"), emit: alignment
      tuple val(id), path("depth_contigs.tsv"), emit: stats

  script:
  """
  ngs-run \
  --sample-id $sample_id \
  --publish-dir $task.publishDir.path \
  --read1 $read1 \
  --read2 $read2 \
  --assembly $assembly \
  --cleaning-kb.path $cleaning_kb \
  --qc-kb.path $qc_kb \
  --organism.genus $organism_genus \
  ${species ? '--organism.species ' + $species : ''}
  """

  stub:
  """
  ngs-run \
  --sample-id $sample_id \
  --publish-dir $task.publishDir.path \
  --read1 $read1 \
  --read2 $read2 \
  --assembly $assembly \
  --cleaning-kb.path $cleaning_kb \
  --qc-kb.path $qc_kb \
  --organism.genus $organism_genus \
  ${species ? '--organism.species ' + $species : ''}
  --stub
  """
}