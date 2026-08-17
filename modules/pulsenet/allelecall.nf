process alleleCalling {

    publishDir "${params.output_dir}/${params.run_directory}/${(params.use_sample_id) ? id : hashed_id}/alleleCalling", mode: params.output_mode_in_use

    input:
        tuple val(hashed_id), val(id), val(genus), val(species), path(blast_kb), path(assembly)
        path(qc_kb), stageAs: "qc_kb"

    output:
        tuple val(hashed_id), path("outputs.json"), emit: outputs
        tuple val(hashed_id), path("stats_calls.json.gz"), emit: stats
        tuple val(hashed_id), path("allele_calls.xml.gz"), path("allele_calls.bam"), path("allele_calls.json.gz"), emit: allele_calls
        tuple val(hashed_id), path("calls_standard.json.gz"), emit: standard_calls
        tuple val(hashed_id), path("calls_core_standard.csv.gz"), path("calls_core_pcr.csv.gz"), emit: csv_core
        tuple val(hashed_id), path("calls_accessory_standard.csv.gz"), path("calls_accessory_pcr.csv.gz"), emit: csv_accessory


    script:
    """
    ngs-run \
    --sample-id $id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --blast-kb.path $blast_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus}
    """

    stub:
    """
    ngs-run \
    --sample-id $id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --blast-kb.path $blast_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus} \
    --stub
    """
}
