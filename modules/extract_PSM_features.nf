
process extract_PSM_features {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/OpenMS-${params.openms_version}/intermediate/${task.process.replace(':', '/')}/" },
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: "*.idXML"

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}/" },
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "${input_idxml.baseName}/log${file(it).name}" }

    input:
        tuple(
            val(sample_name),
            file(input_idxml),
            val(search_engine)
        )

    output:
        tuple(
            val(sample_name),
            file(output_idxml),
            val(search_engine)
        )
        file ".command.*"

    script:
    output_idxml = "${sample_name}_${search_engine}_featureExtracted.idXML"
    """
    PSMFeatureExtractor \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus}
    """
}
