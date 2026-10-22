process prepare_ProteinProphet_idXML {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/ProteinProphet/intermediate/${task.process}/" },
        mode: 'copy', enabled: params.save_intermediate_files, pattern: '*.idXML'

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine)
        path(database_fasta)

    output:
        tuple val(sample_name), path("${sample_name}_${search_engine}_ProteinProphet.idXML"), val(search_engine)

    script:
    """
    python3 ${projectDir}/src/SetIdXMLDatabasePath.py \\
        --input ${input_idxml} \\
        --database ${database_fasta} \\
        --output ${sample_name}_${search_engine}_ProteinProphet.idXML
    """
}
