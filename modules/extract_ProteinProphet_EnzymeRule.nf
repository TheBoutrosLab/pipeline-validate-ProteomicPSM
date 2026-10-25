process extract_ProteinProphet_EnzymeRule {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/ProteinProphet/intermediate/${task.process}/" },
        mode: 'copy', enabled: params.save_intermediate_files, pattern: '*.json'

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine)

    output:
        tuple val(sample_name), path(input_idxml),
            path("${sample_name}_${search_engine}_EnzymeRule.json"), val(search_engine)

    script:
    """
    python3 ${projectDir}/src/ExtractEnzymeRule.py \\
        --input-path ${input_idxml} \\
        --output-path ${sample_name}_${search_engine}_EnzymeRule.json
    """
}
