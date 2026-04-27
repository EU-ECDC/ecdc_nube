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
    tuple val(payLoad), val(run_id)

  tag {"${payLoad.project}:${payLoad.id}"}

  output:
    tuple val(payLoad), path('*.fastq.gz', arity: '1..*')

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
  tuple val(payLoad), path(R1, arity: '1..*'), path(R2, arity: '1..*')

  tag {"${payLoad.project}:${payLoad.id}"}

  output:
  tuple val(payLoad), path("*-combined_?.fastq.gz", arity: '1..*')
 
  script:
  """
  cat ${R1.join(" ")} > ${payLoad.id}-combined_1.fastq.gz
  cat ${R2.join(" ")} > ${payLoad.id}-combined_2.fastq.gz
  """
}

process TRIM {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.3"
  errorStrategy 'ignore'
  time '30m'
  memory '12 GB'

  input:
  tuple val(payLoad), path(reads, arity: '1..*')

  tag {"${payLoad.project}:${payLoad.id}"}

  output:
  tuple val(payLoad), path("*-trimmed_?.fastq.gz", arity: '1..*')
 
  shell:
  def cardinality = reads.size()
  if(cardinality == 1){
    """
    fastp -i ${reads[0]} -o !{payLoad.id}-trimmed_1.fastq.gz
    """
  } else if(cardinality == 2){
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o !{payLoad.id}-trimmed_1.fastq.gz -O !{payLoad.id}-trimmed_2.fastq.gz
    """
  }

}

process DOWNSAMPLE {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.3"
  errorStrategy 'ignore'
  time '30m'

  input:
  tuple val(payLoad), path(reads, arity: '1..*'), val(genome_size)

  tag {"${payLoad.project}:${payLoad.id}"}

  output:
  tuple val(payLoad), path("*-downsampled_?.fastq.gz", arity: '1..*')
 
  shell:
  def cardinality = reads.size()
  if(cardinality == 1){
    """
    rasusa reads -o !{payLoad.id}-downsampled_1.fastq.gz --coverage 150 --genome-size ${genome_size} ${reads[0]}
    """
  } else if(cardinality == 2){
    """
    rasusa reads -o !{payLoad.id}-downsampled_1.fastq.gz -o !{payLoad.id}-downsampled_2.fastq.gz --coverage 150 --genome-size ${genome_size} ${reads[0]} ${reads[1]}
    """
  }

}

process ASSEMBLE {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.2"
  errorStrategy 'ignore'
  //time '20h'

  input:
  tuple val(payLoad), path(reads)

  tag {"${payLoad.project}:${payLoad.id}"}

  publishDir "${params.output}/${payLoad.id}/assemblies/", overwrite: true

  output:
  tuple val(payLoad), path("${payLoad.id}.fasta")

  shell:
  def cardinality = reads.size()
  if(payLoad.sequencing_technology == "ILLUMINA"){
    if(cardinality == 1){
      """
      spades.py -s ${reads[0]} --careful --only-assembler --cov-cutoff 10 -o ${payLoad.id}_assembly
      ln ${payLoad.id}_assembly/contigs.fasta ${payLoad.id}.fasta
      """
    } else if(cardinality == 2){
      """
      spades.py -1 ${reads[0]} -2 ${reads[1]} --careful --only-assembler --cov-cutoff 10 -o ${payLoad.id}_assembly
      ln ${payLoad.id}_assembly/contigs.fasta ${payLoad.id}.fasta
      """
    }
  } else if(payLoad.sequencing_technology == "IONTORRENT"){
    if(cardinality == 1){
      """
      spades.py -s ${reads[0]} -k 21,33,55,77,99,127 --iontorrent --careful -o ${payLoad.id}_assembly
      ln ${payLoad.id}_assembly/contigs.fasta ${payLoad.id}.fasta
      """
    } else if(cardinality == 2){
      """
      spades.py -1 ${reads[0]} -2 ${reads[1]} -k 21,33,55,77,99,127 --iontorrent --careful -o ${payLoad.id}_assembly
      ln ${payLoad.id}_assembly/contigs.fasta ${payLoad.id}.fasta
      """
    }
  } else if(payLoad.sequencing_technology == "ONT"){
    if(cardinality == 1){
      """
      flye --nano-hq ${reads[0]} --out-dir ${payLoad.id}_assembly
      ln ${payLoad.id}_assembly/assembly.fasta ${payLoad.id}.fasta
      """
    }
  } else if(payLoad.sequencing_technology == "PACBIO"){
    if(cardinality == 1){
      """
      flye --pacbio-hifi ${reads[0]} --out-dir ${payLoad.id}_assembly
      ln ${payLoad.id}_assembly/assembly.fasta ${payLoad.id}.fasta
      """
    } 
  }
}

process KLEBORATE {
  container "${params.containerRepository}/ejfresch/kleborate:1.0"
  errorStrategy 'ignore'

  input:
  tuple val(payLoad), path(assembly)

  tag {"${payLoad.project}:${payLoad.id}"}

  publishDir "${params.output}/${payLoad.project}/amr/", overwrite: true

  output:
  path "*.txt", emit: txt_kleborate

  shell:
  """
  kleborate --assemblies !{assembly} --outfile !{payLoad.id}_kleborate.txt --resistance --kaptive
  """
}

process RESFINDER {
  container 'docker.io/genomicepidemiology/resfinder'
  containerOptions '--volume $(pwd):/app --user root'
  errorStrategy 'ignore'

  input:
  tuple val(payLoad), path(assembly)

  tag {"${payLoad.project}:${payLoad.id}"}

  publishDir "${params.output}/${payLoad.project}/amr/", overwrite: true

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
  tuple val(project), val(accession), path(assembly), val(organism), val(experiment_list), path(references_path)

  tag {"${project}:${accession}"}

  publishDir "${params.output}/${project}/species_verification/", overwrite: true

  output:
  tuple val(project), val(accession), val(organism), val(experiment_list), path("${accession}_species.tsv")

  shell:
  """
  ls -1 ${organism}/*.fna > list_ref_genomes.txt
  fastANI -q ${accession}.fasta --rl list_ref_genomes.txt -o fastani_out.txt
  parse_fastani.py -in ${accession}.fasta --data_summary ${organism}/data_summary.tsv
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
    def payLoad = [:]
    payLoad.id                      = it.key
    payLoad.project                 = it.value.project
    payLoad.organism                = it.value.organism
    payLoad.sequencing_technology   = it.value.sequencing_technology[0]
    payLoad.experiment_list         = it.value.experiment_list
    payLoad.schemas                 = it.value?.schemas ?: settings["organism"][it.value.organism].defaultSchemas
    payLoad.num_seq_tech            = it.value.sequencing_technology.flatten().size()

    if(it.value.accessions){
      payLoad.entrypoint = "accessions"
      data = it.value.accessions.flatten().collect{ it.replace("NCBI|","").replace("ENA|","") }[0]
    }
    else if(it.value.reads){
      payLoad.entrypoint = "reads"
      // Fix the s3 paths
      data = it.value.reads.collect{
        nested_list ->
        nested_list.collect { item ->
        item.replace("S3|", "s3://")
        }
      }
      // Get the number of read sets
      payLoad.num_read_sets = data.size()
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
      payLoad.layout = layout
      println data.flatten().collect({it.replace("S3|","s3://")})
    }
    else if(it.value.assembly){
      payLoad.entrypoint = "assemblies"
      data = it.value.assembly
    }
    else if(it.value.sequences){
      payLoad.entrypoint = "sequences"
      payLoad.flags = it.value?.flags ?: ""
      data = it.value.sequences
    }
    return([payLoad, data])
  }

  // Branching by entrypoint
  ch_data = ch_input.branch{payLoad, data ->
    accessions: payLoad.entrypoint == "accessions"
    reads: payLoad.entrypoint == "reads"
    assemblies: payLoad.entrypoint == "assemblies"
    sequences: payLoad.entrypoint == "sequences"
  }
  
  // Getting reads from accessions
  // Current limitations: it will process only the first accession / assumes that you will provide a single accession
  PREFETCH(ch_data.accessions)
  
  // Split the entries with reads based on the sequencing technologies
  ch_reads_by_seq_tech = ch_data.reads.branch{
    payLoad, reads ->
    one: payLoad.num_seq_tech == 1
    two: payLoad.num_seq_tech == 2
    other: true
  }

  // Split the entries with reads based on the num_read_sets and the layout
  ch_reads_combine = ch_reads_by_seq_tech.one.branch{
    payLoad, reads ->
    combine_fastqs: payLoad.num_read_sets == 1 && payLoad.layout == "inferred_paired"
      R1 = reads[0].sort().indexed(1).findAll { i, v -> i % 2 == 1 }.collect{it.value}
      R2 = reads[0].sort().indexed(1).findAll { i, v -> i % 2 == 0 }.collect{it.value}
      return tuple(payLoad, R1, R2)
    no_need: true
      return tuple(payLoad, reads[0])
  }
  
  COMBINE_FASTQ(ch_reads_combine.combine_fastqs)
  
  // Creating a channel with all samples with reads (including those for which we downloaded the reads from NCBI/ENA)
  ch_reads = PREFETCH.out.mix(ch_reads_combine.no_need).mix(COMBINE_FASTQ.out)
  
  // Trimming reads, downsampling and generating assemblies
  TRIM(ch_reads)
  DOWNSAMPLE(TRIM.out.map{
    payLoad, reads ->
      [payLoad, reads, settings["organism"][payLoad.organism].genomeSize]
    }
  )
  ASSEMBLE(DOWNSAMPLE.out)
  
  // Creating a channel with all samples with assemblies
  ch_assemblies = ASSEMBLE.out.mix(ch_data.assemblies)
  ch_assemblies.view()

  // Generating cgMLST profiles
  ALLELE_CALL(ch_assemblies
    .filter{ payLoad, assembly -> payLoad.experiment_list.contains("allele_call") }
    .flatMap{ payLoad, assembly ->
      payLoad.schemas.collect { schema ->
        [payLoad,
        assembly,
        settings["schemas"][schema].schemaPath,
        settings["schemas"][schema].trnFile,
        settings["schemas"][schema].geneList,
        settings["schemas"][schema].containsKey("advOptions") ? settings["schemas"][schema].advOptions.IONTORRENT.refAllelesErrorCorrection : "",
        settings["schemas"][schema].containsKey("advOptions") ? settings["schemas"][schema].advOptions : [:],
        "${payLoad.id}_allele-call_${schema}",
        schema
        ]
      }
    }
  )
  
  // QC
  QC(ch_assemblies.filter{ payLoad, assembly -> payLoad.experiment_list.contains("qc")})

  // Species verification
  SPECIES_VERIFICATION(ch_assemblies.filter{ project, accession, technology, assembly, organism, experiment_list, schemas -> 
    experiment_list.contains("species_verification")
    }.map{project, accession, technology, assembly, organism, experiment_list, schemas ->
      [project, accession, assembly, organism, experiment_list, "${params.speciesReferences}/${organism}/"]
    }
  )
  
  // Pathogen-specific sub-workflows
  SALMISO(ch_assemblies.filter{payLoad, assembly -> payLoad.organism == "SALMISO"})
  ECOLIISO(ch_assemblies.filter{payLoad, assembly -> payLoad.organism == "ECOLIISO"})
  CAMPISO(ch_assemblies.filter{payLoad, assembly -> payLoad.organism == "CAMPISO"})
  HAVISO(ch_data.sequences.filter{payLoad, sequences -> payLoad.organism == "HAVISO"})
  POLIISO(ch_data.sequences.filter{payLoad, sequences -> payLoad.organism == "POLIISO"})

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
