process convert_ProteinProphet_pepXML {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/ProteinProphet/intermediate/${task.process}/" },
        mode: 'copy', enabled: params.save_intermediate_files, pattern: '*.pepXML'

    input:
        tuple val(sample_name), path(input_idxml), path(enzyme_rule_json), val(search_engine)

    output:
        tuple val(sample_name), path("${sample_name}_${search_engine}_ProteinProphet.pepXML"),
            path(enzyme_rule_json), val(search_engine)

    script:
    """
    IDFileConverter \\
        -in ${input_idxml} \\
        -out ${sample_name}_${search_engine}_ProteinProphet.pepXML
    """
}
