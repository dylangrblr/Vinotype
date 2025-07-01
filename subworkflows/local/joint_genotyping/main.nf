/*
 * Subworkflow: joint_genotyping
 * Purpose: Combine per-chromosome GVCFs, normalize, reblock, import to GenomicsDB, genotype, and select variants.
 * Inputs:
 *   - ch_gvcf_bins: Channel emitting tuples [ sample_id, path(gvcf) ] for all bins
 *   - ch_fasta: Channel emitting path to reference FASTA
 *   - ch_intervals: Channel emitting path to intervals.list (optional, can be generated)
 * Outputs:
 *   - snp_vcf: Path to final SNP VCF
 *   - indel_vcf: Path to final INDEL VCF
 *   - versions: Channel with version information from all modules
 */

include { BCFTOOLS_INDEX         } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_CONCAT        } from '../../../modules/nf-core/bcftools/concat/main.nf'
include { BCFTOOLS_NORM          } from '../../../modules/nf-core/bcftools/norm/main.nf'
include { GATK4_REBLOCKGVCF      } from '../../../modules/nf-core/gatk4/reblockgvcf/main.nf'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main.nf'
include { GATK4_GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main.nf'
include { GATK4_SELECTVARIANTS   } from '../../../modules/nf-core/gatk4/selectvariants/main.nf'

workflow JOINT_GENOTYPING {

    take:
        ch_gvcf_tbi
        ch_intervals_list   // path to intervals.list (optional)
        ch_fasta       
        ch_faidx_index
        ch_dict_index
        ch_gvcf
        ch_tbi
        //ch_gvcf_bins   // [ sample_id, path(gvcf) ] for all bins

    main:
    //NOTE: @abhi18av these are necessary only if we split the BAM files per chromosome
        // // 1. Index all GVCFs
        // BCFTOOLS_INDEX(ch_gvcf_bins.map { sample, gvcf -> tuple(sample, gvcf) })

        // // 2. Group GVCFs by sample for concatenation
        // ch_gvcf_by_sample = BCFTOOLS_INDEX.out.tbi
        //     .groupTuple(0) // group by sample_id

        // // 3. Concatenate per-sample GVCFs
        // BCFTOOLS_CONCAT(
        //     ch_gvcf_by_sample.map { sample, gvcf_list -> tuple(sample, gvcf_list) }
        // )

        // // 4. Normalize concatenated GVCFs
        // BCFTOOLS_NORM(
        //     BCFTOOLS_CONCAT.out.vcf.map { sample, vcf -> tuple(sample, vcf, ch_fasta.first()) },
        //     ch_fasta
        // )

        // // 5. Reblock GVCFs
        // GATK4_REBLOCKGVCF(
        //     BCFTOOLS_NORM.out.vcf.map { sample, norm_vcf -> tuple(sample, norm_vcf, ch_fasta.first()) },
        //     ch_fasta,
        //     ch_samtools_faidx,
        //     ch_samtools_dict,
        //     [],
        //     []
        // )

        // // 6. Create sample map for GenomicsDBImport
        // ch_sample_map = GATK4_REBLOCKGVCF.out.vcf
        //     .map { sample, reblocked_vcf -> "${sample}\t${reblocked_vcf}" }
        //     .collectFile(name: 'Vv_sample.map', sep: '\n')

        ch_gvcf_tbi_extended = ch_gvcf_tbi
                                .combine(ch_intervals_list.map{it[1]})
                                .map { meta, vcf, tbi, interval -> tuple( meta, vcf, tbi, interval, [], []) }
                                //.dump(tag: 'GVCFs and TBIs with intervals', pretty: true)


        // Optionally use an existing GenomicsDB workspace or create a new one
        if (params.base_genomicsdb) {

            ch_genomicsdb = Channel.value([[id:file(params.base_genomicsdb).baseName], params.base_genomicsdb])

        } else {

            ch_genomicsdb = Channel.of( [ [id:"Vitis_DB.new"], []] ) 
        }

        // Update existing workspace
        GATK4_GENOMICSDBIMPORT(
            ch_gvcf.collect{it[1]},
            ch_tbi.collect{it[1]},
            ch_genomicsdb,
            ch_fasta,
            ch_faidx_index,
            ch_dict_index,
            ch_intervals_list
        )


        // 8. GenotypeGVCFs
        GATK4_GENOTYPEGVCFS(
                  GATK4_GENOMICSDBIMPORT.out.updatedb,
                  ch_intervals_list,
                  ch_fasta, 
                  ch_faidx_index, 
                  ch_dict_index)


//GATK4_GENOTYPEGVCFS.out.vcf.dump(tag: 'GenotypeGVCFs output', pretty: true)

        // 9. SelectVariants for SNPs and INDELs
        GATK4_SELECTVARIANTS(
            GATK4_GENOTYPEGVCFS.out.vcf_tbi,
            ch_intervals_list,
            'snp'
        )
        // ch_snp_vcf = GATK4_SELECTVARIANTS.out.vcf

        // GATK4_SELECTVARIANTS(
        //     GATK4_GENOTYPEGVCFS.out.vcf.map { vcf -> tuple(vcf, ch_fasta.first(), 'INDEL') }
        // )
        // ch_indel_vcf = GATK4_SELECTVARIANTS.out.vcf


        // Collect version info
        ch_versions = Channel.empty()
        // ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
        // ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        // ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
        // ch_versions = ch_versions.mix(GATK4_REBLOCKGVCF.out.versions)
        // ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
        // ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
        // ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)

    emit:
        // snp_vcf   = ch_snp_vcf
        // indel_vcf = ch_indel_vcf
        versions  = ch_versions
}
