// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

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

    """
    # Complex VCF filtering pipeline
    bcftools norm \\
        ${task.ext.args_bcftools_norm_1} \\
        -f ${fasta} \\
        ${vcf} | \\
    bcftools filter \\
        -e '${task.ext.args_bcftools_filter}' | \\
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
