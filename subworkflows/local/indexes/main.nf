
include { BWA_INDEX      } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_DICT  } from '../../../modules/nf-core/samtools/dict/main'

//TODO: @abhi18av, need to actually make use of these
include { UNTAR as UNTAR_BWA_INDEX } from '../../../modules/nf-core/untar/main'
include { UNTARFILES           } from '../../../modules/nf-core/untarfiles/main'

workflow INDEXES {

    /*
     * Input:
     *   ch_fasta: Channel emitting tuples [ val(meta), path(fasta) ]
     * Output:
     *   versions: Channel emitting all versions.yml files from the modules
     */

    take:
        ch_fasta

    main:

        //ref_bwa_index = Channel.empty()
        // Index reference with BWA
        BWA_INDEX(ch_fasta)

        //ref_samtools_faidx = Channel.empty()
        // Index reference with samtools faidx
        SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta_ -> [meta, fasta_, []] }, false)
        
        //ref_samtools_dict = Channel.empty()
        // Create sequence dictionary with samtools dict
        SAMTOOLS_DICT(ch_fasta)

        // Collect all index files for downstream use
        ch_bwa_index   = BWA_INDEX.out.index
        ch_faidx_index = SAMTOOLS_FAIDX.out.fai
        ch_dict_index  = SAMTOOLS_DICT.out.dict

        // Collect version info from all modules
        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)

    emit:
        bwa_index   = ch_bwa_index
        faidx_index = ch_faidx_index
        dict_index  = ch_dict_index
        versions    = ch_versions
}
