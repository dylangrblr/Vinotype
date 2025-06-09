// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow JOINT_GENOTYPING {

    take:

    // Vitis_DB="${base_dir}/concat/Vitis_DB"
    // sample_map="${vcf_input}/Vv_sample.map"
    // output_norm="${vcf_input}/processed/normalised"
    // output_reblock="${vcf_input}/processed/reblocked"
    // output_variants="${vcf_input}/processed/raw_variants"
    // output_gvcf="${vcf_input}/processed/genotype_GVCF"

    main:

    ch_versions = Channel.empty()

    GATK4_GENOMICSDBIMPORT

    GATK4_GENOTYPEGVCFS


    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
