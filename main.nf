#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import java.text.SimpleDateFormat
import groovy.json.JsonSlurper

include { QC } from './modules/qc.nf'
include { SALMISO } from './modules/salmiso.nf'
include { ECOLIISO } from './modules/ecoliiso.nf'
include { CAMPISO } from './modules/campiso.nf'
include { HAVISO } from './modules/haviso.nf'
include { POLIISO } from './modules/poliiso.nf'
include { ALLELE_CALL } from './modules/allelecall.nf'

def parseJson(input_file){
  def jsonSlurper = new JsonSlurper()
  String currentJSONContent = new File(input_file).text
  def jsonContent = jsonSlurper.parseText(currentJSONContent)
  return(jsonContent)
}

process PREFETCH {
  container "${params.containerRepository}/ejfresch/ncbi-tools:2.3"
  errorStrategy { 
    //TODO would be nice to have 'ignore' right away if the exit code is 80. 
    if (task.attempt == 1) {
      'retry'  // retry only once
    } else {
      'ignore'
    }
  }
  
  input:
    tuple val(meta), val(run_id)

  tag {"${meta.project}:${meta.id}"}

  output:
    tuple val(meta), path('*.fastq.gz', arity: '1..*')

  shell:
  """
  get_sequencing_data.py ${run_id}
  """
}

process COMBINE_FASTQ {
  container "${params.containerRepository}/ejfresch/ncbi-tools:2.3"
  errorStrategy 'ignore'
  time '30m'

  input:
  tuple val(meta), path(R1, arity: '1..*'), path(R2, arity: '1..*')

  tag {"${meta.project}:${meta.id}"}

  output:
  tuple val(meta), path("*-combined_?.fastq.gz", arity: '1..*')
 
  script:
  """
  cat ${R1.join(" ")} > ${meta.id}-combined_1.fastq.gz
  cat ${R2.join(" ")} > ${meta.id}-combined_2.fastq.gz
  """
}

process TRIM {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.3"
  errorStrategy 'ignore'
  time '30m'
  memory '12 GB'

  input:
  tuple val(meta), path(reads, arity: '1..*')

  tag {"${meta.project}:${meta.id}"}

  output:
  tuple val(meta), path("*-trimmed_?.fastq.gz", arity: '1..*')
 
  shell:
  if(meta.single_end){
    """
    fastp -i ${reads[0]} -o !{meta.id}-trimmed_1.fastq.gz
    """
  } else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o !{meta.id}-trimmed_1.fastq.gz -O !{meta.id}-trimmed_2.fastq.gz
    """
  }

}

process DOWNSAMPLE {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.3"
  errorStrategy 'ignore'
  time '30m'

  input:
  tuple val(meta), path(reads, arity: '1..*'), val(genome_size)

  tag {"${meta.project}:${meta.id}"}

  output:
  tuple val(meta), path("*-downsampled_?.fastq.gz", arity: '1..*')
 
  shell:
  
  if (meta.single_end){
    """
    rasusa reads -o !{meta.id}-downsampled_1.fastq.gz --coverage 150 --genome-size ${genome_size} ${reads[0]}
    """
  } else {
    """
    rasusa reads -o !{meta.id}-downsampled_1.fastq.gz -o !{meta.id}-downsampled_2.fastq.gz --coverage 150 --genome-size ${genome_size} ${reads[0]} ${reads[1]}
    """
  }

}

process ASSEMBLE {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.2"
  errorStrategy 'ignore'
  //time '20h'

  input:
  tuple val(meta), path(reads)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.id}/assemblies/", overwrite: true

  output:
  tuple val(meta), path("${meta.id}.fasta")

  shell:
  if(meta.sequencing_technology == "ILLUMINA"){
    if(meta.single_end){
      """
      spades.py -s ${reads[0]} --careful --only-assembler --cov-cutoff 10 -o ${meta.id}_assembly
      ln ${meta.id}_assembly/contigs.fasta ${meta.id}.fasta
      """
    } else {
      """
      spades.py -1 ${reads[0]} -2 ${reads[1]} --careful --only-assembler --cov-cutoff 10 -o ${meta.id}_assembly
      ln ${meta.id}_assembly/contigs.fasta ${meta.id}.fasta
      """
    }
  } else if(meta.sequencing_technology == "IONTORRENT"){
    if(meta.single_end){
      """
      spades.py -s ${reads[0]} -k 21,33,55,77,99,127 --iontorrent --careful -o ${meta.id}_assembly
      ln ${meta.id}_assembly/contigs.fasta ${meta.id}.fasta
      """
    } else {
      """
      spades.py -1 ${reads[0]} -2 ${reads[1]} -k 21,33,55,77,99,127 --iontorrent --careful -o ${meta.id}_assembly
      ln ${meta.id}_assembly/contigs.fasta ${meta.id}.fasta
      """
    }
  } else if(meta.sequencing_technology == "ONT"){
    if(meta.single_end){
      """
      flye --nano-hq ${reads[0]} --out-dir ${meta.id}_assembly
      ln ${meta.id}_assembly/assembly.fasta ${meta.id}.fasta
      """
    }
  } else if(meta.sequencing_technology == "PACBIO"){
    if(meta.single_end){
      """
      flye --pacbio-hifi ${reads[0]} --out-dir ${meta.id}_assembly
      ln ${meta.id}_assembly/assembly.fasta ${meta.id}.fasta
      """
    }
  }
}

process KLEBORATE {
  container "${params.containerRepository}/ejfresch/kleborate:1.0"
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(assembly)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.project}/amr/", overwrite: true

  output:
  path "*.txt", emit: txt_kleborate

  shell:
  """
  kleborate --assemblies !{assembly} --outfile !{meta.id}_kleborate.txt --resistance --kaptive
  """
}

process RESFINDER {
  container 'docker.io/genomicepidemiology/resfinder'
  containerOptions '--volume $(pwd):/app --user root'
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(assembly)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.project}/amr/", overwrite: true

  output:
  path "*.json", emit: json_resfinder

  shell:
  """
  -ifa !{assembly} -o . -s other --acquired
  """
}

process SPECIES_VERIFICATION {
  container "docker.io/ejfresch/fastani:1.34"
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(assembly), path(references_path)

  tag {"${meta.project}:${meta.id}"}

  publishDir "${params.output}/${meta.project}/species_verification/", overwrite: true

  output:
  tuple val(meta), path("${meta.id}_species.tsv")

  shell:
  """
  ls -1 ${meta.organism}/*.fna > list_ref_genomes.txt
  fastANI -q ${meta.id}.fasta --rl list_ref_genomes.txt -o fastani_out.txt
  parse_fastani.py -in ${meta.id}.fasta --data_summary ${meta.organism}/data_summary.tsv
  """
}

workflow {

  target="data/signals/*"

  // Loading data on settings
  settings=parseJson("./settings.json")

  // Creating a channel with the signals 
  ch_signals=Channel.watchPath(target, 'create, modify').until{ it -> it.baseName == 'STOP' }

  // Filtering (looking for .json files) and renaming the files of the signals
  ch_json=ch_signals.filter(~/.+\.json$/).map{
    def orig_file = new File("${it}")
    def new_file = new File("${it}.lock")
    orig_file.renameTo(new_file)
    return(new_file)
  }
  ch_json.view()
  
  // Reshaping the signals
  ch_input = ch_json.splitJson().map{it ->
    def data = []
    def meta = [:]
    meta.id                      = it.key
    meta.project                 = it.value.project
    meta.organism                = it.value.organism
    meta.experiment_list         = it.value.experiment_list
    meta.sequencing_technology   = it.value.sequencing_technology?.getAt(0) ?: ""
    meta.num_seq_tech            = (it.value.sequencing_technology ?: []).flatten().size()
    meta.schemas                 = it.value?.schemas ?: settings["organism"][it.value.organism].defaultSchemas
    meta.flags                   = it.value?.flags ?: ""

    if(it.value.accessions){
      meta.entrypoint = "accessions"
      data = it.value.accessions.flatten().collect{ it.replace("NCBI|","").replace("ENA|","") }[0]
    }
    else if(it.value.reads){
      meta.entrypoint = "reads"
      // Fix the s3 paths
      data = it.value.reads.collect{
        nested_list ->
        nested_list.collect { item ->
        item.replace("S3|", "s3://")
        }
      }
      // Get the number of read sets
      meta.num_read_sets = data.size()
      // Get the layouts
      def layouts = data.collect{
        k ->
        def cardinality_read_set = k.size() 
        def lout = "error"
        if ( cardinality_read_set == 1 ){
          lout = "single"
        }
        else if ( cardinality_read_set == 2 ){
          lout = "paired"
        }
        else if ( (cardinality_read_set % 2) == 0 ){ 
          lout = "inferred_paired"
        }
        else {
          lout = "error"
        }        
        lout
      }
      def layout = "not set"
      if ( layouts.toSet().size() == 1 ){
        layout = layouts[0]
      }
      else {
        layout = "mixed/error"
      }
      meta.layout = layout
      if ( meta.layout == "single" ){
          meta.single_end = true
        }
      else {
        meta.single_end = false
      }
      println data.flatten().collect({it.replace("S3|","s3://")})
    }
    else if(it.value.assembly){
      meta.entrypoint = "assemblies"
      data = it.value.assembly
    }
    else if(it.value.sequences){
      meta.entrypoint = "sequences"
      data = it.value.sequences
    }
    return([meta, data])
  }

  // Branching by entrypoint
  ch_data = ch_input.branch{meta, data ->
    accessions: meta.entrypoint == "accessions"
    reads: meta.entrypoint == "reads"
    assemblies: meta.entrypoint == "assemblies"
    sequences: meta.entrypoint == "sequences"
  }
  
  // Getting reads from accessions
  // Current limitations: it will process only the first accession / assumes that you will provide a single accession
  PREFETCH(ch_data.accessions)
  
  // Split the entries with reads based on the sequencing technologies
  ch_reads_by_seq_tech = ch_data.reads.branch{
    meta, reads ->
    one: meta.num_seq_tech == 1
    two: meta.num_seq_tech == 2
    other: true
  }

  // Split the entries with reads based on the num_read_sets and the layout
  ch_reads_combine = ch_reads_by_seq_tech.one.branch{
    meta, reads ->
    combine_fastqs: meta.num_read_sets == 1 && meta.layout == "inferred_paired"
      R1 = reads[0].sort().indexed(1).findAll { i, v -> i % 2 == 1 }.collect{it.value}
      R2 = reads[0].sort().indexed(1).findAll { i, v -> i % 2 == 0 }.collect{it.value}
      return tuple(meta, R1, R2)
    no_need: true
      return tuple(meta, reads[0])
  }
  
  COMBINE_FASTQ(ch_reads_combine.combine_fastqs)
  
  // Creating a channel with all samples with reads (including those for which we downloaded the reads from NCBI/ENA)
  ch_reads = PREFETCH.out.mix(ch_reads_combine.no_need).mix(COMBINE_FASTQ.out)
  
  // Trimming reads, downsampling and generating assemblies
  TRIM(ch_reads)
  DOWNSAMPLE(TRIM.out.map{
    meta, reads ->
      [meta, reads, settings["organism"][meta.organism].genomeSize]
    }
  )
  ASSEMBLE(DOWNSAMPLE.out)
  
  // Creating a channel with all samples with assemblies
  ch_assemblies = ASSEMBLE.out.mix(ch_data.assemblies)
  ch_assemblies.view()

  // Generating cgMLST profiles
  ALLELE_CALL(ch_assemblies
    .filter{ meta, assembly -> meta.experiment_list.contains("allele_call") }
    .flatMap{ meta, assembly ->
      meta.schemas.collect { schema ->
        [meta,
        assembly,
        "${params.allelecallSchemas}/${settings["schemas"][schema].schemaPath}",
        "${params.allelecallSchemas}/${settings["schemas"][schema].trnFile}",
        "${params.allelecallSchemas}/${settings["schemas"][schema].geneList}",
        settings["schemas"][schema].containsKey("advOptions") ? "${params.allelecallSchemas}/${settings["schemas"][schema].advOptions.IONTORRENT.refAllelesErrorCorrection}" : "",
        settings["schemas"][schema].containsKey("advOptions") ? settings["schemas"][schema].advOptions : [:],
        "${meta.id}_allele-call_${schema}",
        schema
        ]
      }
    }
  )
  
  // QC
  QC(ch_assemblies.filter{ meta, assembly -> meta.experiment_list.contains("qc")})

  // Species verification
  SPECIES_VERIFICATION(ch_assemblies.filter{ meta, assembly -> meta.experiment_list.contains("species_verification")}
    .map{meta, assembly ->
      [meta, assembly, "${params.speciesReferences}/${meta.organism}/"]
    }
  )
  
  // Pathogen-specific sub-workflows
  SALMISO(ch_assemblies.filter{meta, assembly -> meta.organism == "SALMISO"})
  ECOLIISO(ch_assemblies.filter{meta, assembly -> meta.organism == "ECOLIISO"})
  CAMPISO(ch_assemblies.filter{meta, assembly -> meta.organism == "CAMPISO"})
  HAVISO(ch_data.sequences.filter{meta, sequences -> meta.organism == "HAVISO"})
  POLIISO(ch_data.sequences.filter{meta, sequences -> meta.organism == "POLIISO"})

}

workflow.onComplete{
  if (workflow.success) {
    println "Workflow completed successfully."
    // determine when the onComplete event handler was triggered
    def date = new Date()
    def custom_date_format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
    def timestamp_onComplete=custom_date_format.format(date)

    /* remove the signals, given that all the expected
     results have been generated */
    def proc = ['./scripts/remove_signals.sh'].execute()
    proc.waitForProcessOutput()

    // write a NF_READY file
    def oncomplete_signal = "./data/signals/NF_READY"
    def content = "$timestamp_onComplete $workflow.runName $workflow.sessionId"
    file(oncomplete_signal).write(content)
  }   
}

workflow.onError{
  // determine when the onError event handler was triggered
  def date = new Date()
  def custom_date_format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
  def timestamp_onError=custom_date_format.format(date)

  // write a NF_ERROR file
  def filename = "./data/signals/NF_ERROR"
  def content = "$timestamp_onError\t$workflow.runName\t$workflow.sessionId\t$workflow.errorMessage"
  file(filename).write(content)
  println "[onError] errorMessage: $workflow.errorMessage"
  // remove the ".lock" extensions
  def proc = ['./scripts/remove_locks.sh'].execute()
  proc.waitForProcessOutput()

}
