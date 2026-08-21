process ESM_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../../esm/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0' :
        'ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0' }"

    input:
    tuple val(meta), path(mapping), path(scores)

    output:
    tuple val(meta), path("*.esm_annotated.tsv"), emit: tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_esm_results.py \\
        --mapping ${mapping} \\
        --scores  ${scores}  \\
        --prefix  ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f2)
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.esm_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.2.2
    END_VERSIONS
    """
}
