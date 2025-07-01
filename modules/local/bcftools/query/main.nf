process BCFTOOLS_QUERY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val annotation_types


    output:
    tuple val(meta), path("*.txt"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   
    // Process annotation types and create safe filenames
    def safe_annotation_types = annotation_types.collect { annotation ->
        def safe_name = annotation.replaceAll(/[\/\%\[\]\s]/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, '')
        [annotation: annotation, safe_name: safe_name]
    }

    """
    # Function to safely extract annotation with error handling
    extract_annotation_safe() {
        local vcf_input="\$1"
        local annotation="\$2"
        local safe_name="\$3"
        local output_file="${prefix}.annotation"."\${safe_name}.txt"
        
        echo "Extracting annotation: \${annotation}"
        
        # Try to extract the annotation, handle errors gracefully
        if bcftools query -f "\${annotation}\\n" "\${vcf_input}" 2>/dev/null > "\${output_file}"; then
            # Check if file has content
            if [[ -s "\${output_file}" ]]; then
                # Remove empty lines and spaces, convert to proper format
                sed '/^\$/d' "\${output_file}" | tr ' ' '\\n' | grep -v '^\$' > "\${output_file}.tmp" && mv "\${output_file}.tmp" "\${output_file}"
                echo "Successfully extracted \${annotation} to \${output_file}"
            else
                echo "Warning: No data found for annotation \${annotation}, creating empty file"
            fi
        else
            echo "Warning: Failed to extract annotation \${annotation}, annotation may not exist in VCF"
        fi
    }

    # Process each annotation type
    ${safe_annotation_types.collect { item ->
        "extract_annotation_safe \"${vcf}\" \"${item.annotation}\" \"${item.safe_name}\""
    }.join('\n    ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def safe_annotation_types = annotation_types.collect { annotation ->
        def safe_name = annotation.replaceAll(/[\/\%\[\]\s]/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, '')
        [annotation: annotation, safe_name: safe_name]
    }
    
    """
    ${safe_annotation_types.collect { item ->
        "touch ${item.safe_name}.txt"
    }.join('\n    ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}