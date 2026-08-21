process ESM_PREDICT {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0' :
        'ghcr.io/dylangrblr/vinotype/esm_vep:1.0.0' }"

    input:
    tuple val(meta), path(fasta_dir), path(var_dir)
    path model_cache  // pre-downloaded HuggingFace cache dir; pass [] to download at runtime

    output:
    tuple val(meta), path("*.scores.csv"), emit: scores
    tuple val(meta), path("figures/")    , emit: figures, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def cache_env = model_cache     ? "export HF_HOME=${model_cache}" : ''
    """
    ${cache_env}

    esm_vep.py \\
        --fasta_dir   ${fasta_dir}         \\
        --var_dir     ${var_dir}           \\
        --output      ${prefix}.scores.csv \\
        --figures_dir figures              \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f2)
        transformers: \$(python -c "import transformers; print(transformers.__version__)")
        torch: \$(python -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.scores.csv
    mkdir -p figures

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        transformers: 4.41.2
        torch: 2.3.0
    END_VERSIONS
    """
}
