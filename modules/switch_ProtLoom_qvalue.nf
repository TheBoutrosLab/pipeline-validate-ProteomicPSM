process switch_ProtLoom_qvalue {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/ProtLoom-${params.protloom_version}/intermediate/${task.process.replace(':', '/')}/" },
        mode: 'copy', enabled: params.save_intermediate_files, pattern: '*.idXML'

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy', pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple val(sample_name), path(input_idxml)

    output:
        tuple val(sample_name), path("${sample_name}_protloom_qvalue.idXML"), emit: switched
        path '.command.*'

    script:
    """
    IDScoreSwitcher \\
        -in ${input_idxml} \\
        -out ${sample_name}_protloom_qvalue.idXML \\
        -threads ${task.cpus} \\
        -new_score 'q_value' \\
        -new_score_type 'ProtLoom q-value' \\
        -new_score_orientation 'lower_better'
    """
}
