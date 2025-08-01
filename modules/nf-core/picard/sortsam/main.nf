process PICARD_SORTSAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    //TODO: @abhi18av update all picard containers to 2.20.3 in the environment files, as per the shell scripts
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.20.3--0' :
        'biocontainers/picard:2.20.3--0' }"

    input:
    tuple val(meta), path(bam)
    val sort_order

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard SortSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        SortSam \\
        -Xmx${avail_mem}M \\
        I=$bam \\
        O=${prefix}.bam \\
        SORT_ORDER=$sort_order

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SortSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
