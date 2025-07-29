process GENOZIP_COMPRESS {
    tag "$meta.id"
    label 'process_medium'

//FIXME: Use the latest containers for genozip
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'oras://community.wave.seqera.io/library/genozip:15.0.73--5e7a11b582e9aebe':
//        'community.wave.seqera.io/library/genozip:15.0.73--b885f85775e92bfa' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/genozip:15.0.16--205028cd2b169348':
        'community.wave.seqera.io/library/genozip:15.0.16--205028cd2b169348' }"


    input:
    tuple val(meta), path(bam)
    path genozip_license_file
    tuple val(meta2), path(fasta_genozip) // genozip reference fasta file

    output:
    tuple val(meta), path("*.genozip"), emit: bam_genozip
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
        $bam \\
        --reference ${fasta_genozip} \\
        --output ${prefix}.bam.genozip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genozip: \$(genozip --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genozip: \$(genozip --version)
    END_VERSIONS
    """
}
