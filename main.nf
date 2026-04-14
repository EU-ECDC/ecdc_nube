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
    tuple val(project), val(accession), val(run_id), val(technology), val(organism), val(experiment_list), val(schemas)

  tag {"${project}:${accession}"}

  output:
    tuple val(project), val(accession), val(technology), path('*.fastq.gz', arity: '1..*'), val(organism), val(experiment_list), val(schemas)

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
  tuple val(project), val(accession), val(technology), path(R1, arity: '1..*'), path(R2, arity: '1..*'), val(layout), val(organism), val(experiment_list), val(schemas)

  tag {"${project}:${accession}"}

  output:
  tuple val(project), val(accession), val(technology), path("*-combined_?.fastq.gz", arity: '1..*'), val(layout), val(organism), val(experiment_list), val(schemas)
 
  script:
  """
  cat ${R1.join(" ")} > ${accession}-combined_1.fastq.gz
  cat ${R2.join(" ")} > ${accession}-combined_2.fastq.gz
  """
}


process TRIM {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.3"
  errorStrategy 'ignore'
  time '30m'
  memory '12 GB'

  input:
  tuple val(project), val(accession), val(technology), path(reads, arity: '1..*'), val(cardinality), val(organism), val(experiment_list), val(genome_size), val(schemas)

  tag {"${project}:${accession}"}

  output:
  tuple val(project), val(accession), val(technology), path("*-trimmed_?.fastq.gz", arity: '1..*'), val(cardinality), val(organism), val(experiment_list), val(genome_size), val(schemas)
 
  shell:
  if(cardinality == 1){
    """
    fastp -i ${reads[0]} -o !{accession}-trimmed_1.fastq.gz
    """
  } else if(cardinality == 2){
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o !{accession}-trimmed_1.fastq.gz -O !{accession}-trimmed_2.fastq.gz
    """
  }

}


process DOWNSAMPLE {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.3"
  errorStrategy 'ignore'
  time '30m'

  input:
  tuple val(project), val(accession), val(technology), path(reads, arity: '1..*'), val(cardinality), val(organism), val(experiment_list), val(genome_size), val(schemas)

  tag {"${project}:${accession}"}

  output:
  tuple val(project), val(accession), val(technology), path("*-downsampled_?.fastq.gz", arity: '1..*'), val(cardinality), val(organism), val(experiment_list), val(schemas)
 
  shell:
  if(cardinality == 1){
    """
    rasusa reads -o !{accession}-downsampled_1.fastq.gz --coverage 150 --genome-size ${genome_size} ${reads[0]}
    """
  } else if(cardinality == 2){
    """
    rasusa reads -o !{accession}-downsampled_1.fastq.gz -o !{accession}-downsampled_2.fastq.gz --coverage 150 --genome-size ${genome_size} ${reads[0]} ${reads[1]}
    """
  }

}

process ASSEMBLE {
  container "${params.containerRepository}/ejfresch/denovo_assembly-tools:1.2"
  errorStrategy 'ignore'
  //time '20h'

  input:
  tuple val(project), val(accession), val(technology), path(read_set), val(cardinality), val(organism), val(experiment_list), val(schemas)

  tag {"${project}:${accession}"}

  publishDir "${params.output}/${project}/assemblies/", overwrite: true

  output:
  tuple val(project), val(accession), val(technology), path("${accession}.fasta"), val(organism), val(experiment_list), val(schemas)

  shell:
  if(technology == "ILLUMINA"){
    if(cardinality == 1){
      """
      spades.py -s ${read_set[0]} --careful --only-assembler --cov-cutoff 10 -o ${accession}_assembly
      ln ${accession}_assembly/contigs.fasta ${accession}.fasta
      """
    } else if(cardinality == 2){
      """
      spades.py -1 ${read_set[0]} -2 ${read_set[1]} --careful --only-assembler --cov-cutoff 10 -o ${accession}_assembly
      ln ${accession}_assembly/contigs.fasta ${accession}.fasta
      """
    }
  } else if(technology == "IONTORRENT"){
    if(cardinality == 1){
      """
      spades.py -s ${read_set[0]} -k 21,33,55,77,99,127 --iontorrent --careful -o ${accession}_assembly
      ln ${accession}_assembly/contigs.fasta ${accession}.fasta
      """
    } else if(cardinality == 2){
      """
      spades.py -1 ${read_set[0]} -2 ${read_set[1]} -k 21,33,55,77,99,127 --iontorrent --careful -o ${accession}_assembly
      ln ${accession}_assembly/contigs.fasta ${accession}.fasta
      """
    }
  } else if(technology == "ONT"){
    if(cardinality == 1){
      """
      flye --nano-hq ${read_set[0]} --out-dir ${accession}_assembly
      ln ${accession}_assembly/assembly.fasta ${accession}.fasta
      """
    }
  } else if(technology == "PACBIO"){
    if(cardinality == 1){
      """
      flye --pacbio-hifi ${read_set[0]} --out-dir ${accession}_assembly 
      ln ${accession}_assembly/assembly.fasta ${accession}.fasta
      """
    } 
  }
}

process KLEBORATE {
  container "${params.containerRepository}/ejfresch/kleborate:1.0"
  errorStrategy 'ignore'

  input:
  tuple val(project), val(accession), path(assembly), val(organism)

  tag {"${project}:${accession}"}

  publishDir "${params.output}/${project}/amr/", overwrite: true

  output:
  path "*.txt", emit: txt_kleborate

  shell:
  """
  kleborate --assemblies !{assembly} --outfile !{accession}_kleborate.txt --resistance --kaptive
  """
}

process RESFINDER {
  container 'docker.io/genomicepidemiology/resfinder'
  containerOptions '--volume $(pwd):/app --user root'
  errorStrategy 'ignore'

  input:
  tuple val(project), val(accession), path(assembly), val(organism)

  tag {"${project}:${accession}"}

  publishDir "${params.output}/${project}/amr/", overwrite: true

  output:
  path "*.json", emit: json_resfinder

  shell:
  """
  -ifa !{assembly} -o . -s other --acquired
  """
}

workflow {

  target="data/signals/*"

  // Loading data on settings
  settings=parseJson("./settings.json")

  // Creating a channel with the signals 
  ch_signals=Channel.watchPath(target, 'create, modify').until{ it -> it.baseName == 'STOP' }
  //ch_signals=Channel.fromPath(target) 

  // Filtering (looking for .json files) and renaming the files of the signals
  a01=ch_signals.filter(~/.+\.json$/).map{
    def orig_file = new File("${it}")
    def new_file = new File("${it}.lock")
    orig_file.renameTo(new_file)
    return(new_file)
  }
  a01.view()
  
  // Reshaping the signals
  a=a01.splitJson().map{it -> 
  def data=[:]
  data["tag"]="${it.value.project}:${it.key}"
  data["payLoad"] = it.value
  data.payLoad.key = it.key
  data.payLoad.schemas = it.value?.schemas ?: settings["organism"][data.payLoad.organism].defaultSchemas
  data.payLoad.num_seq_tech = data.payLoad.sequencing_technology.flatten().size()
  if(data.payLoad.accessions){
      data.payLoad.entrypoint = "accessions"
      data.payLoad.accessions = data.payLoad.accessions.flatten()
      data.payLoad.accessions = data.payLoad.accessions.collect{ it.replace("NCBI|","").replace("ENA|","") }
    }
  else if(data.payLoad.reads){
      data.payLoad.entrypoint = "reads"
      // I fix the s3 paths
      data.payLoad.reads = data.payLoad.reads.collect{
        nested_list ->
          nested_list.collect { item ->
          item.replace("S3|", "s3://")
          }
      }
      // I get the number of read sets
      data.payLoad.num_read_sets = data.payLoad.reads.size()
      // I get the layouts
      def layouts = data.payLoad.reads.collect{ 
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

      data.payLoad.layout = layout

      println data.payLoad.reads.flatten().collect({it.replace("S3|","s3://")})
    }
    else if(data.payLoad.assembly){
      data.payLoad.entrypoint = "assemblies"
    }
    else if(data.payLoad.sequences){
      data.payLoad.entrypoint = "sequences"
      data.payLoad.flags = it.value?.flags ?: ""
    }
  return(data)
}


// Branching by entrypoint
b = a.branch{i -> 
    accessions: i.payLoad.entrypoint == "accessions" 
      return i 
    reads: i.payLoad.entrypoint == "reads" 
      return i 
    assemblies: i.payLoad.entrypoint == "assemblies" 
      return i
    sequences: i.payLoad.entrypoint == "sequences"
      return i
  } 

b.sequences.view()


  // Getting reads from accessions
  // Current limitations: it will process only the first accession / assumes that you will provide a single accession
  PREFETCH(b.accessions.map{
    it -> 
    [it.payLoad.project, it.payLoad.key, it.payLoad.accessions[0], it.payLoad.sequencing_technology[0], it.payLoad.organism, it.payLoad.experiment_list, it.payLoad.schemas]
    }
  )

  // I split the entries with reads based on the sequencing technologies 
  ch_reads_by_seq_tech = b.reads.branch{
    it ->
    one: it.payLoad.num_seq_tech == 1
      return it
    two: it.payLoad.num_seq_tech == 2
      return it
    other: true
      return it
  }
  
  // I split the entries with reads based on the num_read_sets and the layout
  ch_combine = ch_reads_by_seq_tech.one.branch{
    it -> 
    combine_fastq: it.payLoad.num_read_sets == 1 && it.payLoad.layout == "inferred_paired"
      it.payLoad.R1 = it.payLoad.reads[0].sort().indexed(1).findAll { i, v -> i % 2 == 1 }.collect{it.value}
      it.payLoad.R2 = it.payLoad.reads[0].sort().indexed(1).findAll { i, v -> i % 2 == 0 }.collect{it.value}
      return it
    no_need: true
      return it
  }

  ch_combine.combine_fastq.view()

  COMBINE_FASTQ(ch_combine.combine_fastq.map{
    it ->
    [it.payLoad.project, 
    it.payLoad.key, 
    it.payLoad.sequencing_technology[0], 
    it.payLoad.R1,
    it.payLoad.R2,
    it.payLoad.layout,
    it.payLoad.organism,
    it.payLoad.experiment_list,
    it.payLoad.schemas
    ]
  })

  COMBINE_FASTQ.out.view()





  // Creating a channel with all samples with reads (including those for which we downloaded the reads from NCBI/ENA)
  ch_reads = PREFETCH.out.map{
    project, accession, technology, reads, organism, experiment_list, schemas ->
    [project, accession, technology, reads, reads.size(), organism, experiment_list, settings["organism"][organism].genomeSize, schemas]
    }.mix(
      ch_combine.no_need.map{
        it -> 
        def read_set = it.payLoad.reads[0]
        [it.payLoad.project, it.payLoad.key, it.payLoad.sequencing_technology[0], read_set, read_set.size(), it.payLoad.organism, it.payLoad.experiment_list, settings["organism"][it.payLoad.organism].genomeSize, it.payLoad.schemas]
      }
    ).mix(
      COMBINE_FASTQ.out.map{
        project, accession, technology, reads, layout, organism, experiment_list, schemas ->
        [project, accession, technology, reads, reads.size(), organism, experiment_list, settings["organism"][organism].genomeSize, schemas]
      }
    )
  
  ch_reads.view()
  // Trimming reads, downsampling and generating assemblies
  TRIM(ch_reads) | DOWNSAMPLE | ASSEMBLE

  // Creating a channel with all samples with assemblies
  ch_assemblies = ASSEMBLE.out.map{
    project, accession, technology, assembly, organism, experiment_list, schemas ->
      [project, accession, technology, assembly, organism, experiment_list, schemas]
  }.mix(
    b.assemblies.map{
      it ->
        [it.payLoad.project, 
        it.payLoad.key, 
        it.payLoad.sequencing_technology[0], 
        it.payLoad.assembly, 
        it.payLoad.organism, 
        it.payLoad.experiment_list,
        it.payLoad.schemas
        ]
    }
  )
  ch_assemblies.view()

  // Generating cgMLST profiles
  ALLELE_CALL(ch_assemblies
    .filter{ it -> it[5].contains("allele_call") }
    .flatMap{ project, accession, technology, assembly, organism, experiment_list, schemas ->
      schemas.collect { schema ->
        [project, accession, technology, assembly, organism, experiment_list,
        settings["schemas"][schema].schemaPath,
        settings["schemas"][schema].trnFile,
        settings["schemas"][schema].geneList,
        settings["schemas"][schema].containsKey("advOptions") ? settings["schemas"][schema].advOptions.IONTORRENT.refAllelesErrorCorrection : "",
        settings["schemas"][schema].containsKey("advOptions") ? settings["schemas"][schema].advOptions : [:],
        "${accession}_allele-call_${schema}",
        schema
        ]
      }
    }
  )
  
  // QC
  QC(ch_assemblies.filter{it -> it[5].contains("qc")})

  // Pathogen-specific sub-workflows
  SALMISO(ch_assemblies.filter{it -> it[4] == "SALMISO"})
  ECOLIISO(ch_assemblies.filter{it -> it[4] == "ECOLIISO"})
  CAMPISO(ch_assemblies.filter{it -> it[4] == "CAMPISO"})
  HAVISO(b.sequences.filter{it -> it.payLoad.organism == "HAVISO"}.map{it -> [
      it.payLoad.sequences,
      it.payLoad.flags,
      it.payLoad.project,
      it.payLoad.organism,
      it.payLoad.key, 
      it.payLoad.experiment_list]
    }
)
POLIISO(b.sequences.filter{it -> it.payLoad.organism == "POLIISO"}.map{it -> [
      it.payLoad.sequences,
      it.payLoad.flags,
      it.payLoad.project,
      it.payLoad.organism,
      it.payLoad.key, 
      it.payLoad.experiment_list]
    }
)




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
