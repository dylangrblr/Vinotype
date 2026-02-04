process GATK4SPARK_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"


    // FIXME for singularity container
    // container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    //     ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/498aea9c9bcaf736b9fb2a01366c1b7b38ccc0d38143178afc325d6a93241447/data'
    //     : 'community.wave.seqera.io/library/gatk4-spark:4.6.2.0--8b5cd67ee60a714e'}"

    container "broadinstitute/gatk:4.6.2.0"

    input:
    tuple val(meta),  path(input), path(input_index), path(intervals), path(dragstr_model)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(dbsnp)
    tuple val(meta6), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.tbi")       , emit: vcf_tbi
    tuple val(meta), path("*.vcf.gz")                      , emit: vcf
    tuple val(meta), path("*.tbi")                         , optional:true, emit: tbi
    path "versions.yml"                                    , emit: versions
    //NOTE: The BAM_OUTPUT is not supported as of 05-FEB-2026 https://gatk.broadinstitute.org/hc/en-us/articles/360037433931-HaplotypeCallerSpark-BETA
    // tuple val(meta), path("*.realigned.bam")               , optional:true, emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    def interval_command = intervals ? "--intervals $intervals" : ""
    def dragstr_command = dragstr_model ? "--dragstr-params-path $dragstr_model" : ""
    def bamout_command = "" //args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    cp ${dict} ${meta4.id}.dict

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCallerSpark \\
        --input $input \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        --spark-master local[${task.cpus}] \\
        $dbsnp_command \\
        $interval_command \\
        $dragstr_command \\
        $bamout_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """


    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = "" //args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${stub_realigned_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
