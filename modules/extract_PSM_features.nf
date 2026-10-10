
process extract_PSM_features {
    container params.docker_image_openms

    publishDir "${params.intermediate_file_dir}/${task.process.replace(':', '/')}/",
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: "*.idXML"

    publishDir "${params.process_log_dir}/${task.process.replace(':', '/')}/",
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "${input_idxml.baseName}/log${file(it).name}" }

    input:
        tuple(
            file(input_idxml),
            val(search_engine)
        )

    output:
        tuple(
            file(output_idxml),
            val(search_engine)
        )
        file ".command.*"

    script:
    output_idxml = "${params.sample_name}_${search_engine}_featureExtracted.idXML"
    """
    PSMFeatureExtractor \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus}
    """
}
