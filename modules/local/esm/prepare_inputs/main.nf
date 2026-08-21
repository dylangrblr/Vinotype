process ESM_PREPARE_INPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../../esm/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0' :
        'ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0' }"

    input:
    tuple val(meta), path(tsv), path(fasta)

    output:
    tuple val(meta), path("fasta_dir"), path("var_dir"), emit: dirs
    tuple val(meta), path("*.mapping.tsv")             , emit: mapping
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prepare_esm_inputs.py \\
        --tsv       ${tsv}      \\
        --fasta     ${fasta}    \\
        --fasta_dir fasta_dir   \\
        --var_dir   var_dir     \\
        --prefix    ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f2)
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p fasta_dir var_dir
    touch ${prefix}.mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        biopython: 1.83
    END_VERSIONS
    """
}
