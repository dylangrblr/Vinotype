include { GUNZIP             } from '../../../modules/nf-core/gunzip/main'
include { BWA_INDEX          } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX      } from '../../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX     } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_DICT      } from '../../../modules/nf-core/samtools/dict/main'
include { GENOZIP_REFERENCE  } from '../../../modules/local/genozip/reference/main'

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

        ch_versions = Channel.empty()
        ch_fasta_genozip = Channel.empty()
        ch_bwa_index = Channel.empty()
        ch_faidx_index = Channel.empty()
        ch_dict_index = Channel.empty()
        ch_versions = Channel.empty()

        // -----------------------------------------------------------------------
        // Decompress the reference when it is plain-gzip (.gz).
        // samtools faidx / dict and bwa-mem2 index require an uncompressed FASTA
        // (or a bgzip-indexed one).  NCBI .fna.gz files are plain gzip and
        // cannot be indexed directly — GUNZIP is run only when needed.
        // GENOZIP_REFERENCE uses the original ch_fasta because genozip accepts .gz.
        // -----------------------------------------------------------------------
        ch_fasta_gz = ch_fasta.branch {
            gz:    it[1].toString().endsWith('.gz')
            plain: true
        }
        GUNZIP(ch_fasta_gz.gz)
        ch_versions = ch_versions.mix(GUNZIP.out.versions)

        // ch_fasta_ready: always plain FASTA — used by all three indexing tools
        ch_fasta_ready = ch_fasta_gz.plain.mix(GUNZIP.out.gunzip)

        // BWA index
        if (params.ref_bwa_index) {
            ch_bwa_index = Channel.value([[id:file(params.ref_fasta).baseName], params.ref_bwa_index])
        } else {
            BWAMEM2_INDEX(ch_fasta_ready)
            ch_bwa_index = BWAMEM2_INDEX.out.index.first()
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }

        // samtools faidx
        if (params.ref_samtools_faidx) {
            ch_faidx_index = Channel.value([[id:file(params.ref_fasta).baseName], params.ref_samtools_faidx])
        } else {
            SAMTOOLS_FAIDX(ch_fasta_ready.map { meta, fasta_ -> [meta, fasta_, []] }, false)
            ch_faidx_index = SAMTOOLS_FAIDX.out.fai.first()
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }

        // samtools dict
        if (params.ref_samtools_dict) {
            ch_dict_index = Channel.value([[id:file(params.ref_fasta).baseName], params.ref_samtools_dict])
        } else {
            SAMTOOLS_DICT(ch_fasta_ready)
            ch_dict_index = SAMTOOLS_DICT.out.dict.first()
            ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
        }

        // genozip reference — uses the original ch_fasta (genozip handles .gz natively)
        if(params.use_genozip && params.genozip_license_file) {
            if (params.ref_fasta_genozip ) {
                ch_fasta_genozip = Channel.value([[id:file(params.ref_fasta).baseName], params.ref_fasta_genozip])
            } else {
                GENOZIP_REFERENCE(ch_fasta, params.genozip_license_file)
                ch_fasta_genozip = GENOZIP_REFERENCE.out.fasta_genozip
                ch_versions = ch_versions.mix(GENOZIP_REFERENCE.out.versions)
            }
        }

        //TODO @abhi18av
        // // Intervals list
        // if (params.ref_intervals) {
        //     ch_bwa_index = Channel.value([[id:file(params.ref_fasta).baseName], params.ref_intervals])
        // } else {
        //     UTILS_INTERVALS_LIST(ch_fasta)
        //     ch_intervals_list = UTILS_INTERVALS_LIST.out.intervals
        // }


    emit:
        fasta_genozip = ch_fasta_genozip
        bwa_index   = ch_bwa_index
        faidx_index = ch_faidx_index
        dict_index  = ch_dict_index
        versions    = ch_versions
}
