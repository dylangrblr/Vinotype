include { SNPEFF_DOWNLOAD      } from '../../../modules/nf-core/snpeff/download/main'
include { SNPEFF_SNPEFF     } from '../../../modules/nf-core/snpeff/snpeff/main'

workflow VARIANT_EFFECT_PREDICTION {

    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
