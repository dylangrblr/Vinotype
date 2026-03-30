/*
 * Variant calling subworkflow: BAM → per-chromosome BAM → GVCF
 * Inputs:
 *   - ch_bam: tuple(meta, bam file)
 *   - ch_fasta: tuple(meta, fasta file)
 * Outputs:
 *   - gvcf_files: tuple(meta, chr, gvcf file)
 *   - versions: version info from all modules
 */

include { GENOZIP_COMPRESS as GENOZIP_COMPRESS__GVCF      } from '../../../modules/local/genozip/compress/main'
include { GENOZIP_DECOMPRESS as GENOZIP_DECOMPRESS__GVCF  } from '../../../modules/local/genozip/decompress/main'
include { BCFTOOLS_SPLITGVCFBYCHROMOSOME } from '../../../modules/local/bcftools/splitgvcfbychromosome/main'
include { SAMTOOLS_INDEX      } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_IDXSTATS  } from '../../../modules/nf-core/samtools/idxstats/main.nf'
include { SAMTOOLS_VIEW      } from '../../../modules/nf-core/samtools/view/main.nf'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main.nf'
include { GATK4SPARK_HAPLOTYPECALLER } from '../../../modules/local/gatk4spark/haplotypecaller/main.nf'
include { GATK4_INDEXFEATUREFILE } from '../../../modules/nf-core/gatk4/indexfeaturefile/main.nf'


workflow VARIANT_CALLING {

    take:
        ch_samplesheet_gvcf // channel: samplesheet read in from --input
        ch_bam_bai      // tuple(meta, bam, bai)
        ch_fasta    // tuple(meta, fasta)
        ch_faidx_index  // tuple(meta, faidx_index)
        ch_dict_index  // tuple(meta, dict_index)
        ch_fasta_genozip  // tuple(meta, fasta_genozip)


    main:

        // Collect version info
        ch_versions = Channel.empty()
        ch_gvcf_tbi = Channel.empty()
        ch_gvcf = Channel.empty()
        ch_tbi = Channel.empty()

        // Assume ch_bam_bai is a channel emitting (meta, bam, bai)
        ch_bam_bai_extended = ch_bam_bai.map { meta, bam, bai -> tuple(meta, bam, bai, [], []) }\

        if(params.skip_gatk4_haplotypecaller_spark) {

            GATK4_HAPLOTYPECALLER(
                ch_bam_bai_extended,
                ch_fasta,
                ch_faidx_index,
                ch_dict_index,
                [[],[]],
                [[],[]],
            )
            ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
            ch_gvcf_tbi = GATK4_HAPLOTYPECALLER.out.vcf_tbi
            ch_tbi = GATK4_HAPLOTYPECALLER.out.tbi
            ch_gvcf = GATK4_HAPLOTYPECALLER.out.vcf

        } else {

            GATK4SPARK_HAPLOTYPECALLER(
                ch_bam_bai_extended,
                ch_fasta,
                ch_faidx_index,
                ch_dict_index,
                [[],[]],
                [[],[]],
            )
            ch_versions = ch_versions.mix(GATK4SPARK_HAPLOTYPECALLER.out.versions)
            ch_gvcf_tbi = GATK4SPARK_HAPLOTYPECALLER.out.vcf_tbi
            ch_tbi = GATK4SPARK_HAPLOTYPECALLER.out.tbi
            ch_gvcf = GATK4SPARK_HAPLOTYPECALLER.out.vcf

        }

        if(params.use_genozip && params.genozip_license_file) {

            GENOZIP_COMPRESS__GVCF(
                ch_gvcf_tbi,
                params.genozip_license_file,
                ch_fasta_genozip
            )
            ch_versions = ch_versions.mix(GENOZIP_COMPRESS__GVCF.out.versions.first())

            // GENOZIP_COMPRESS__GVCF.out.genozipped.dump(tag: 'GENOZIP_COMPRESS__GVCF')


            // If samplesheet has GENOZIP VCF files, decompress them
            GENOZIP_DECOMPRESS__GVCF(
                ch_samplesheet_gvcf,
                params.genozip_license_file,
                ch_fasta_genozip
            )
            ch_versions = ch_versions.mix(GENOZIP_DECOMPRESS__GVCF.out.versions.first())

            // GENOZIP_DECOMPRESS__GVCF.out.gvcf.dump(tag: 'GENOZIP_DECOMPRESS__GVCF_gvcf')

            // GENOZIP_DECOMPRESS__GVCF.out.gvcf.view()


            GATK4_INDEXFEATUREFILE( GENOZIP_DECOMPRESS__GVCF.out.gvcf )

            GATK4_INDEXFEATUREFILE.out.index.dump(tag: "GATK4_INDEXFEATUREFILE_out_index")


            ch_gvcf_tbi = ch_gvcf_tbi.mix((GENOZIP_DECOMPRESS__GVCF.out.gvcf).join(GATK4_INDEXFEATUREFILE.out.index))
            ch_gvcf = ch_gvcf.mix(GENOZIP_DECOMPRESS__GVCF.out.gvcf)
            ch_tbi = ch_tbi.mix(GATK4_INDEXFEATUREFILE.out.index)

            // ch_gvcf_tbi.dump(tag:'ch_gvcf_tbi')

            // ch_bam_bai = SAMTOOLS_VIEW.out.bam_bai.mix(GENOZIP_DECOMPRESS__GVCF.out.bam_bai)
        }

        // Create one task per (sample gVCF, chromosome) and split with bcftools view.
        ch_chromosomes = ch_faidx_index
            .map { _meta, fai ->
                file(fai).text.readLines()
                    .findAll { line -> line?.trim() }
                    .collect { line -> line.split('\t')[0] }
            }
            .flatten()
            .unique()

        ch_gvcf_tbi_per_chromosome = ch_gvcf_tbi
            .cross(ch_chromosomes)
            .map { gvcf_tuple, chromosome ->
                def (meta, gvcf, tbi) = gvcf_tuple
                tuple(meta + [chromosome: chromosome], gvcf, tbi, chromosome)
            }

        BCFTOOLS_SPLITGVCFBYCHROMOSOME(ch_gvcf_tbi_per_chromosome)
        ch_versions = ch_versions.mix(BCFTOOLS_SPLITGVCFBYCHROMOSOME.out.versions)

        ch_gvcf_tbi = BCFTOOLS_SPLITGVCFBYCHROMOSOME.out.vcf_tbi
        ch_gvcf = BCFTOOLS_SPLITGVCFBYCHROMOSOME.out.vcf
        ch_tbi = BCFTOOLS_SPLITGVCFBYCHROMOSOME.out.tbi



        // // Get chromosome names for each BAM
        // SAMTOOLS_IDXSTATS(ch_bam_bai)
        // ch_chr_info = SAMTOOLS_IDXSTATS.out.idxstats.map { meta, idxstats_file ->
        //     def chromosomes = file(idxstats_file).text.readLines()
        //         .findAll { !it.startsWith('*') }
        //         .collect { it.split('\t')[0] }
        //     chromosomes.collect { chr -> tuple(meta, chr) }
        // }.flatten()

        // // Subset BAMs by chromosome
        // ch_bam_chr = ch_chr_info.map { meta, chr ->
        //     def bam = ch_bam_bai.filter { it[0] == meta }[0][1]
        //     tuple(meta, bam, chr)
        // }
        // SAMTOOLS_VIEW(ch_bam_chr)

        // // Run HaplotypeCaller per chromosome
        // ch_hc_input = SAMTOOLS_VIEW.out.bam.map { meta, chr_bam, chr ->
        //     tuple(meta, chr_bam, chr, ch_fasta.filter { it[0] == meta }[0][1])
        // }
        // GATK4_HAPLOTYPECALLER(ch_hc_input)


        // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        // ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
        // ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)


        ch_gvcf_tbi.dump(tag:'ch_gvcf_tbi')
        ch_gvcf.dump(tag:'ch_gvcf')
        ch_tbi.dump(tag:'tbi_files')

    emit:
        gvcf_tbi_files = ch_gvcf_tbi
        gvcf_files = ch_gvcf
        tbi_files = ch_tbi
        versions   = ch_versions
}
