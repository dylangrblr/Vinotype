process GENOZIP_DECOMPRESS {
    tag "$meta.id"
    label 'process_medium'

    //FIXME: Use the latest containers for genozip
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'docker://community.wave.seqera.io/library/bcftools_samtools_genozip:600b21e770b9aafc':
    //     'community.wave.seqera.io/library/bcftools_samtools_genozip:600b21e770b9aafc' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bcftools_samtools_tabix_genozip:77b4ccef5d13e48b':
        'community.wave.seqera.io/library/bcftools_samtools_tabix_genozip:77b4ccef5d13e48b' }"


    input:
    tuple val(meta), path(bam_or_gvcf)
    path genozip_license_file
    tuple val(meta2), path(fasta_genozip) // genozip reference fasta file

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true

    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf, optional: true

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def isBAM = file(bam_or_gvcf).getExtension() == 'bam' ? true : false

    """
    genounzip \\
        --threads $task.cpus \\
        $args \\
        $bam_or_gvcf \\
        --no-cache \\
        --reference ${fasta_genozip} \\
        --output ${prefix}


    ${isBAM ? "" : "bgzip -c ${prefix} > ${prefix}.gz"}

    ${isBAM ? "" : "rm ${prefix}"}

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
    touch ${prefix}.bam.bai
    touch ${prefix}.g.vcf
    touch ${prefix}.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genounzip: \$(genounzip --version)
    END_VERSIONS
    """
}
