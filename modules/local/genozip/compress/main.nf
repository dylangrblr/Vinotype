process GENOZIP_COMPRESS {
    tag "$meta.id"
    label 'process_medium'

    //FIXME: Use the latest containers for genozip
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'docker://community.wave.seqera.io/library/bcftools_samtools_genozip:600b21e770b9aafc':
    //     'community.wave.seqera.io/library/bcftools_samtools_genozip:600b21e770b9aafc' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'docker://community.wave.seqera.io/library/genozip:15.0.16--205028cd2b169348':
    //     'community.wave.seqera.io/library/genozip:15.0.16--205028cd2b169348' }"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bcftools_samtools_tabix_genozip:77b4ccef5d13e48b':
        'community.wave.seqera.io/library/bcftools_samtools_tabix_genozip:77b4ccef5d13e48b' }"




    input:
    tuple val(meta), path(bam_or_gvcf), path(bai_or_tbi)
    path genozip_license_file
    tuple val(meta2), path(fasta_genozip) // genozip reference fasta file

    output:
    tuple val(meta), path("*.genozip"), emit: genozipped
    tuple val(meta), path("*.genozip"), path(bai_or_tbi), emit: genozipped_index
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
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
