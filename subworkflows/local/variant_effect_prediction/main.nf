include { SNPEFF_SNPEFF     } from '../../../modules/nf-core/snpeff/snpeff/main'
include { SNPSIFT_FILTER    } from '../../../modules/local/snpsift/filter/main'
include { GATK4_MERGEVCFS   } from '../../../modules/nf-core/gatk4/mergevcfs/main.nf'
include { ESM_PREPARATION   } from '../esm_preparation/main'
include { ESM_SCORING       } from '../esm_scoring/main'

workflow VARIANT_EFFECT_PREDICTION {

    take:
    ch_gvcf_tbi        // channel: [ val(meta), gvcf, tbi ]
    ch_dict_index      // channel: [ val(meta), dict ]
    snpeff_db_name
    snpeff_db
    snpeff_config
    protein_ref_fasta  // string path or null — enables ESM scoring when set
    protein_ref_fai    // string path or null — .fai index for protein_ref_fasta
    esm_model_cache    // string path or null — pre-downloaded HuggingFace cache

    main:

    // Keep chromosome files split through filtering, then merge once before annotation.

    ch_versions = Channel.empty()

    ch_merge_input = ch_gvcf_tbi
        .map { meta, vcf, _tbi ->
            def merge_id = (meta.id as String).replaceFirst(/\.[^.]+$/, '')
            tuple(merge_id, vcf)
        }
        .groupTuple(by: 0)
        .map { merge_id, vcfs -> tuple([id: merge_id], vcfs) }

    GATK4_MERGEVCFS(
        ch_merge_input,
        ch_dict_index
    )
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())

    ch_merged_vcf_tbi = GATK4_MERGEVCFS.out.vcf.join(GATK4_MERGEVCFS.out.tbi)

    ch_snpeff_in = Channel.of([snpeff_db_name, snpeff_db, snpeff_config])

    SNPEFF_SNPEFF( ch_merged_vcf_tbi, ch_snpeff_in, [[],[]])
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions.first())

    SNPSIFT_FILTER ( SNPEFF_SNPEFF.out.vcf )
    ch_versions = ch_versions.mix(SNPSIFT_FILTER.out.versions.first())

    //
    // Optional: ESM2 variant-effect scoring
    // This is only enabled if users supply a protein reference FASTA file
    //
    if (protein_ref_fasta) {
        ch_protein_ref = Channel.of([
            file(protein_ref_fasta, checkIfExists: true),
            protein_ref_fai ? file(protein_ref_fai, checkIfExists: true) : []
        ])

        ESM_PREPARATION ( SNPSIFT_FILTER.out.vcf, ch_protein_ref )
        ch_versions = ch_versions.mix(ESM_PREPARATION.out.versions)

        ch_mc = Channel.value(
            esm_model_cache ? file(esm_model_cache, checkIfExists: true) : []
        )

        ESM_SCORING (
            ESM_PREPARATION.out.tsv,
            ESM_PREPARATION.out.transcripts,
            ch_mc
        )
        ch_versions = ch_versions.mix(ESM_SCORING.out.versions)
    }

    emit:
    vcf           = SNPSIFT_FILTER.out.vcf     // channel: [ val(meta), vcf ]
    snpeff_report = SNPEFF_SNPEFF.out.report   // channel: [ val(meta), csv ]
    versions      = ch_versions                // channel: [ versions.yml ]
}
