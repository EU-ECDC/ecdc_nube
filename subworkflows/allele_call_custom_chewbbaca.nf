include { BLAST_GENE_CALL_SSI } from '../modules/custom_chewbbaca/blast_gene_call.nf'
include { CUSTOM_CHEWBBACA_SSI } from '../modules/custom_chewbbaca/custom_chewbbaca.nf'

workflow CUSTOM_ALLELE_CALL_SSI {
take:
  data
  settings

main:
  BLAST_GENE_CALL_SSI(data
    .map{ meta, assembly ->
      [meta,
      assembly,
      "${params.allelecallSchemas}/${settings["schemas"][meta.schema].schemaPath}"
      ]
    }
  )
  CUSTOM_CHEWBBACA_SSI(BLAST_GENE_CALL_SSI.out)
}
