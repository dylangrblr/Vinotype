process SNPSIFT_EXTRACTFIELDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:5.2--hdfd78af_0':
        'biocontainers/snpsift:5.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.ann_fields.tsv"), emit: tsv
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Keep header plus rows where IMPACT is HIGH/MODERATE and HGVS_P is present and well-formed (no '?')
    """
    SnpSift extractFields ${args} ${vcf} \\
        CHROM POS REF ALT "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].HGVS_P" \\
        | awk -F'\t' 'NR==1 || ((\$5=="HIGH"||\$5=="MODERATE") && \$9!="" && \$9!~/[?]/)' \
        > ${prefix}.ann_fields.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$(echo \$(SnpSift 2>&1) | grep SnpSift | head -n1 | cut -d ' ' -f3)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann_fields.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: 5.2
    END_VERSIONS
    """
}
