process SNPSIFT_VCFEFFONEPERLLINE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:5.2--hdfd78af_0':
        'biocontainers/snpsift:5.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.split.vcf.gz"), emit: vcf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PERL_SCRIPT="${moduleDir}/bin/vcfEffOnePerLine.pl"

    zcat ${vcf} \\
        | perl "\$PERL_SCRIPT" ${args} \\
        | gzip > ${prefix}.split.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$(echo \$(SnpSift 2>&1) | grep SnpSift | head -n1 | cut -d ' ' -f3)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.split.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: 5.2
    END_VERSIONS
    """
}
