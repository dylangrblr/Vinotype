/*
 * Subworkflow: joint_genotyping
 * Purpose: Import per-chromosome GVCFs to GenomicsDB and emit chromosome-level multi-sample VCFs.
 * Inputs:
 *   - ch_gvcf_bins: Channel emitting tuples [ sample_id, path(gvcf) ] for all bins
 *   - ch_fasta: Channel emitting path to reference FASTA
 *   - ch_intervals: Channel emitting path to intervals.list (optional, can be generated)
 * Outputs:
 *   - gvcf_tbi_files: Channel of chromosome-level tuple(meta, vcf, tbi)
 *   - versions: Channel with version information from executed modules
 */

include { BCFTOOLS_INDEX         } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_CONCAT        } from '../../../modules/nf-core/bcftools/concat/main.nf'
include { BCFTOOLS_NORM          } from '../../../modules/nf-core/bcftools/norm/main.nf'
include { GATK4_REBLOCKGVCF      } from '../../../modules/nf-core/gatk4/reblockgvcf/main.nf'
include { GATK4_INDEXFEATUREFILE } from '../../../modules/nf-core/gatk4/indexfeaturefile/main.nf'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main.nf'
include { GATK4_GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main.nf'
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

        ch_intervals_path = ch_intervals_list.map { _meta, intervals -> intervals }

        ch_reblock_input = ch_gvcf_tbi
            .combine(ch_intervals_path)
            .map { meta, gvcf, tbi, intervals -> tuple(meta, gvcf, tbi, intervals) }

        ch_reblock_chromosome = ch_gvcf_tbi
            .map { meta, _gvcf_unused, _tbi_unused -> meta.chromosome ?: 'unknown' }

        GATK4_REBLOCKGVCF(
            ch_reblock_input,
            ch_reblock_chromosome,
            ch_fasta,
            ch_faidx_index,
            ch_dict_index,
            [],
            []
        )
        ch_versions = ch_versions.mix(GATK4_REBLOCKGVCF.out.versions.first())

        ch_gvcf_tbi_for_import = GATK4_REBLOCKGVCF.out.vcf

        ch_reblocked_by_chromosome = ch_gvcf_tbi_for_import
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
              ch_intervals_list,
                  ch_fasta,
                  ch_faidx_index,
                  ch_dict_index)

        ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

    //GATK4_GENOTYPEGVCFS.out.vcf.dump(tag: 'GenotypeGVCFs output', pretty: true)

        // Collect version info
        // ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
        // ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        // ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
        // ch_versions = ch_versions.mix(GATK4_REBLOCKGVCF.out.versions)
        // ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
        // ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
        // ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)


    emit:
        gvcf_tbi_files = GATK4_GENOTYPEGVCFS.out.vcf_tbi
        versions  = ch_versions
}
