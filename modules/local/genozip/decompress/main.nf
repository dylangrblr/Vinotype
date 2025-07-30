process GENOZIP_DECOMPRESS {
    tag "$meta.id"
    label 'process_medium'

//FIXME: Use the latest containers for genozip. 
//NOTE: Also needs samtools in the container for generating index, check the environment.yml file

//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'oras://community.wave.seqera.io/library/genozip:15.0.73--5e7a11b582e9aebe':
//        'community.wave.seqera.io/library/genozip:15.0.73--b885f85775e92bfa' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/samtools_genozip:f15b7c96d5843aea':
        'community.wave.seqera.io/library/samtools_genozip:f15b7c96d5843aea' }"


    input:
    tuple val(meta), path(bam)
    path genozip_license_file
    tuple val(meta2), path(fasta_genozip) // genozip reference fasta file

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bam"), path("*bai"), emit: bam_bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    genounzip \\
        --threads $task.cpus \\
        $args \\
        $bam \\
        --no-cache \\
        --index \\
        --reference ${fasta_genozip} \\
        --output ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genounzip: \$(genounzip --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genounzip: \$(genounzip --version)
    END_VERSIONS
    """
}

