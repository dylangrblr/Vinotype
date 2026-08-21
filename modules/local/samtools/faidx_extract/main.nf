process SAMTOOLS_FAIDX_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    //TODO: Update naming of transcriptome_fasta to protein_fasta
    //Double check that fasta headers in transcriptome_fasta match the transcript IDs in the FEATUREID column of the SNPSIFT_EXTRACTFIELDS output TSV
    input:
    tuple val(meta), path(tsv)      // filtered annotation TSV from SNPSIFT_EXTRACTFIELDS
    path fasta                       // indexed transcriptome FASTA
    path fai                         // samtools .fai index for the FASTA

    output:
    tuple val(meta), path("*.transcripts.fa"), emit: fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Extract unique transcript IDs from FEATUREID column (col 8), then pull sequences
    """
    tail -n +2 ${tsv} \\
        | cut -f8 \\
        | sort -u \\
        > transcript_ids.txt

    samtools faidx \\
        ${args} \\
        ${fasta} \\
        -r transcript_ids.txt \\
        > ${prefix}.transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
