include { SNPEFF_SNPEFF     } from '../../../modules/nf-core/snpeff/snpeff/main'

workflow VARIANT_EFFECT_PREDICTION {

    take:
    ch_gvcf_tbi // channel: [ val(meta), gvcf, tbi ]
    snpeff_db      // path to snpeff database
    snpeff_config  // path to snpeff config file

    main:

    ch_versions = Channel.empty()

    SNPEFF_SNPEFF( ch_gvcf_tbi, snpeff_db, snpeff_config, [[],[]])

    // SNPSIFT_FILTER (
    //     SNPEFF_SNPEFF.out.vcf,
    //     [[],[]]
    // )

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
