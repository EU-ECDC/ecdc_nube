include { PULSENET_ALLELECALLING } from '../modules/pulsenet/allelecall.nf'
include { PIGZ_COMPRESS } from '../modules/pulsenet/pigz.nf'

workflow ALLELE_CALL_PULSENET {
take:
  data
  settings

main:
  compressed_data = PIGZ_COMPRESS(data).compressed
  PULSENET_ALLELECALLING(compressed_data
    .map{ meta, assembly ->
      [meta,
      assembly,
      "${params.allelecallSchemas}/${settings["schemas"][meta.schema].schemaPath}",
      "${params.pulsenetQCKB}"
      ]
    }
  )
}