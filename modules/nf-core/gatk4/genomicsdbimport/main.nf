process GATK4_GENOMICSDBIMPORT {
    tag "$chromosome_name"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
        path(gvcfs)
        path(tbis)
        tuple val(meta), path(existing_db)
        tuple val(meta1), path(reference_fasta)
        tuple val(meta2), path(reference_fasta_index)
        tuple val(meta3), path(reference_dict)
        val(chromosome_name)


    output:
        tuple val(meta), path("${updated_db}"), emit: updatedb
        tuple val(meta), path("${prefix}"), optional: true, emit: genomicsdb
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // Determine input command: sample map or list of --variant
    input_command = gvcfs.collect { gvcf -> "--variant ${gvcf}" }.join(' ')
    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    interval_command = chromosome_name ? "--intervals ${chromosome_name}" : ""
    updated_db = "${meta.id}"

    // // Interval list mode
    // if (run_intlist) {
    //     genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
    //     interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    // }

    // Update workspace mode
    if (existing_db) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${existing_db}"
        interval_command = ''
        updated_db = "${meta.id}"
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GenomicsDBImport] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    cp ${reference_dict} ${meta3.id}.dict

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenomicsDBImport \\
        $input_command \\
        $genomicsdb_command \\
        $interval_command \\
        --reference ${reference_fasta} \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    // stub:
    // prefix   = task.ext.prefix ?: "${meta.id}"

    // genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    // interval_command = "--intervals ${interval_file}"
    // updated_db = ""

    // // if (run_intlist) {
    // //     genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
    // //     interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    // // }

    // if (existing_db) {
    //     genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
    //     interval_command = ''
    //     updated_db = "${meta.id}.new"
    // }

    // def stub_genomicsdb = genomicsdb_command == "--genomicsdb-workspace-path ${prefix}" ? "touch ${prefix}" : ""
    // def stub_interval   = interval_command == "--output-interval-list-to-file ${prefix}.interval_list" ? "touch ${prefix}.interval_list" : ""
    // def stub_update     = updated_db != "" ? "touch ${updated_db}" : ""

    // """
    // ${stub_genomicsdb}
    // ${stub_interval}
    // ${stub_update}

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    // END_VERSIONS
    // """
}
