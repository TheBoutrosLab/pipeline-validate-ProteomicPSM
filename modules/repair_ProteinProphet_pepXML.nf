process repair_ProteinProphet_pepXML {
    containerOptions params.project_script_mount_options
    container params.docker_image_fragpipe

    publishDir { "${params.intermediate_file_dir}/${sample_name}/ProteinProphet/intermediate/${task.process}/" },
        mode: 'copy', enabled: params.save_intermediate_files, pattern: '*.pepXML'

    input:
        tuple val(sample_name), path(input_pepxml), path(enzyme_rule_json), val(search_engine)

    output:
        tuple val(sample_name),
            path("${sample_name}_${search_engine}_ProteinProphet_terminiFixed.pepXML"),
            val(search_engine)

    script:
    """
    python3 ${projectDir}/src/FixPepXMLTermini.py \\
        --input-path ${input_pepxml} \\
        --enzyme-rule-path ${enzyme_rule_json} \\
        --output-path ${sample_name}_${search_engine}_ProteinProphet_terminiFixed.pepXML
    """
}
