include { SNPSIFT_VCFEFFONEPERLLINE } from '../../../modules/local/snpsift/vcfEffOnePerLine/main'
include { SNPSIFT_EXTRACTFIELDS     } from '../../../modules/local/snpsift/extractFields/main'
include { SAMTOOLS_FAIDX_EXTRACT    } from '../../../modules/local/samtools/faidx_extract/main'

workflow ESM_PREPARATION {

    take:
    ch_vcf              // channel: [ val(meta), vcf.gz ] — filtered+annotated VCF from VARIANT_EFFECT_PREDICTION
    ch_transcriptome    // channel: [ path(fasta), path(fai) ] — indexed reference proteins FASTA

    main:

    ch_versions = Channel.empty()

    //
    // Split multi-effect ANN entries so that each variant effect occupies its own VCF line
    //
    SNPSIFT_VCFEFFONEPERLLINE ( ch_vcf )
    ch_versions = ch_versions.mix(SNPSIFT_VCFEFFONEPERLLINE.out.versions.first())

    //
    // Extract annotation fields and filter to HIGH/MODERATE impact variants
    // with HGVS_P annotations
    //
    SNPSIFT_EXTRACTFIELDS ( SNPSIFT_VCFEFFONEPERLLINE.out.vcf )
    ch_versions = ch_versions.mix(SNPSIFT_EXTRACTFIELDS.out.versions.first())

    //
    // Retrieve the reference protein sequence for each variant
    //
    ch_transcriptome_fasta = ch_transcriptome.map { fasta, fai -> fasta }
    ch_transcriptome_fai   = ch_transcriptome.map { fasta, fai -> fai   }

    SAMTOOLS_FAIDX_EXTRACT (
        SNPSIFT_EXTRACTFIELDS.out.tsv,
        ch_transcriptome_fasta,
        ch_transcriptome_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_EXTRACT.out.versions.first())

    emit:
    tsv          = SNPSIFT_EXTRACTFIELDS.out.tsv       // channel: [ val(meta), tsv ]
    transcripts  = SAMTOOLS_FAIDX_EXTRACT.out.fasta    // channel: [ val(meta), fasta ]
    versions     = ch_versions                         // channel: [ versions.yml ]
}
