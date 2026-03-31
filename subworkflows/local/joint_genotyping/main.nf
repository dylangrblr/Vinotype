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
include { GATK4_MERGEVCFS        } from '../../../modules/nf-core/gatk4/mergevcfs/main.nf'
workflow JOINT_GENOTYPING {

    take:
        ch_gvcf_tbi
        ch_intervals_list   // path to intervals.list (optional)
        ch_fasta
        ch_faidx_index
        ch_dict_index
        _ch_gvcf
        _ch_tbi
        //ch_gvcf_bins   // [ sample_id, path(gvcf) ] for all bins

    main:

        ch_versions = Channel.empty()

        ch_reblocked_by_chromosome = ch_gvcf_tbi
            .map { meta, gvcf, tbi -> tuple(meta.chromosome ?: 'unknown', gvcf, tbi) }
            .groupTuple(by: 0)

        ch_import_input = ch_reblocked_by_chromosome.map { chromosome_name, gvcfs, tbis ->
            def db_meta = [id: "Vitis_DB.${chromosome_name}", chromosome: chromosome_name]
            def existing_db = params.base_genomicsdb ? file(params.base_genomicsdb) : []
            tuple(db_meta, gvcfs, tbis, chromosome_name, existing_db)
        }

        // One GenomicsDBImport instance per chromosome using reblocked per-chromosome gVCFs.
        GATK4_GENOMICSDBIMPORT(
            ch_import_input.map { row -> row[1] },
            ch_import_input.map { row -> row[2] },
            ch_import_input.map { row -> tuple(row[0], row[4]) },
            ch_fasta,
            ch_faidx_index,
            ch_dict_index,
            ch_import_input.map { row -> row[3] }
        )
        ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())


    //NOTE: @abhi18av - now 20 - instances of GENOMICSDBIMPORT have run and they
    //need to be passed individually to GENOTYPEGVCFS


        // One GenotypeGVCFs instance per imported chromosome workspace.
        GATK4_GENOTYPEGVCFS(
                  GATK4_GENOMICSDBIMPORT.out.updatedb,
              ch_intervals_list.first(),
                  ch_fasta,
                  ch_faidx_index,
                  ch_dict_index)

        ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

        // Merge chromosome-level GenotypeGVCFs outputs before variant selection/filtering.
        ch_merge_input = GATK4_GENOTYPEGVCFS.out.vcf_tbi
            .map { meta, vcf, _tbi ->
                def merge_id = (meta.id as String).replaceFirst(/\.[^.]+$/, '')
                tuple(merge_id, vcf)
            }
            .groupTuple(by: 0)
            .map { merge_id, vcfs -> tuple([id: merge_id], vcfs) }

        GATK4_MERGEVCFS(
            ch_merge_input,
            ch_dict_index.first()
        )

        ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())

    //GATK4_GENOTYPEGVCFS.out.vcf.dump(tag: 'GenotypeGVCFs output', pretty: true)

        // Collect version info
        // ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
        // ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        // ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
        // ch_versions = ch_versions.mix(GATK4_REBLOCKGVCF.out.versions)
        // ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
        // ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
        // ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)


    //NOTE: @abhi18av, then the output VCF files need to be reconciled using
    //gatk-MergeVCF


    // Then another gatk-MergeVCF to combine all the vcf files into a single
    // multi-sample vcf

    emit:
        gvcf_tbi_files = GATK4_MERGEVCFS.out.vcf.join(GATK4_MERGEVCFS.out.tbi)
        versions  = ch_versions
}
