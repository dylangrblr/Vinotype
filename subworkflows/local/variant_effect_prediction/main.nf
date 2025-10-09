include { SNPEFF_SNPEFF     } from '../../../modules/nf-core/snpeff/snpeff/main'
include { SNPSIFT_FILTER    } from '../../../modules/local/snpsift/filter/main'

workflow VARIANT_EFFECT_PREDICTION {

    take:
    ch_gvcf_tbi // channel: [ val(meta), gvcf, tbi ]
    snpeff_db_name
    snpeff_db
    snpeff_config

    main:

    ch_versions = Channel.empty()

    ch_snpeff_in = Channel.of([snpeff_db_name, snpeff_db, snpeff_config])

    SNPEFF_SNPEFF( ch_gvcf_tbi, ch_snpeff_in, [[],[]])

    SNPSIFT_FILTER ( SNPEFF_SNPEFF.out.vcf )

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
