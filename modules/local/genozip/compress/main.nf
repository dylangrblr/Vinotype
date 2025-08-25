process GENOZIP_COMPRESS {
    tag "$meta.id"
    label 'process_medium'

//FIXME: Use the latest containers for genozip
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'oras://community.wave.seqera.io/library/genozip:15.0.73--5e7a11b582e9aebe':
//        'community.wave.seqera.io/library/genozip:15.0.73--b885f85775e92bfa' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bcftools_samtools_genozip:6edb384e04414b50':
        'community.wave.seqera.io/library/bcftools_samtools_genozip:6edb384e04414b50' }"


    input:
    tuple val(meta), path(bam_or_gvcf)
    path genozip_license_file
    tuple val(meta2), path(fasta_genozip) // genozip reference fasta file

    output:
    tuple val(meta), path("*.genozip"), emit: genozipped
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    genozip \\
        --licfile ${genozip_license_file} \\
        --threads $task.cpus \\
        $args \\
        $bam_or_gvcf \\
        --reference ${fasta_genozip} \\
        --output ${prefix}.genozip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genozip: \$(genozip --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.genozip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genozip: \$(genozip --version)
    END_VERSIONS
    """
}
