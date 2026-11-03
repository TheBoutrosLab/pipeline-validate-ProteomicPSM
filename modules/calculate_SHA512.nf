process calculate_SHA512 {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { (params.run_protloom || params.protloom_split_fdr) ?
        "${params.final_output_dir}/${sample_name}/ProtLoom-${params.protloom_version}/output" :
        "${params.final_output_dir}/${sample_name}/OpenMS-${params.openms_version}/output" },
        mode: 'copy', pattern: '*.idXML.sha512'

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy', pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple val(sample_name), path(input_idxml)

    output:
        tuple val(sample_name), path("${input_idxml.name}.sha512"), emit: checksum
        path '.command.*'

    script:
    """
    python3 ${projectDir}/src/WriteSHA512.py ${input_idxml}
    """
}
