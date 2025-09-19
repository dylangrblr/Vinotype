process UTILS_FINALFILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bcftools_vcftools:0c2ee0aa2e344d87':
        'community.wave.seqera.io/library/bcftools_vcftools:0c2ee0aa2e344d87' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def optional_indel_filter = task.ext.args_bcftools_filter_indels ? "-e '${task.ext.args_bcftools_filter_indels}'" : ""

    """
    # Complex VCF filtering pipeline
    bcftools norm \\
        ${task.ext.args_bcftools_norm_1} \\
        -f ${fasta} \\
        ${vcf} | \\
    bcftools filter \\
        -e '${task.ext.args_bcftools_filter_snps}' | \\
    ${optional_indel_filter} \\
    bcftools norm \\
        ${task.ext.args_bcftools_norm_2} | \\
    bcftools view \\
        ${task.ext.args_bcftools_view} | \\
    vcftools \\
        --vcf - \\
        --stdout \\
        ${task.ext.args_vcftools} | \\
    bgzip > ${prefix}_filtered.vcf.gz

    # Index the output VCF
    tabix -p vcf ${prefix}_filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        vcftools: \$(vcftools --version 2>&1 | head -n1 | sed 's/^.*VCFtools (//; s/).*\$//')
        bgzip: \$(bgzip --version 2>&1 | head -n1 | sed 's/^.*bgzip //; s/ .*\$//')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_filtered.vcf.gz
    touch ${prefix}_filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        vcftools: \$(vcftools --version 2>&1 | head -n1 | sed 's/^.*VCFtools (//; s/).*\$//')
        bgzip: \$(bgzip --version 2>&1 | head -n1 | sed 's/^.*bgzip //; s/ .*\$//')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*\$//')
    END_VERSIONS
    """
}
