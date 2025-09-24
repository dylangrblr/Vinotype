process SNPSIFT_FILTER {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:5.2--hdfd78af_0':
        'biocontainers/snpsift:5.2--hdfd78af_0' }"


    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*ann.vcf.gz"), emit: vcf
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat ${vcf} \\
        | java -jar SnpSift.jar filter "${args}" \\
        | gzip > ${prefix}.vcf.gz || { echo "SnpSift filter failed"; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        snpsift: \$(SnpSift.jar 2>&1 | grep SnpSift | head -n1)
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$(snpsift --version)
    END_VERSIONS
    """
}
