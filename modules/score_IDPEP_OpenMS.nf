process score_IDPEP_OpenMS {
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
        file(input_idxml)

    output:
        file(output_idxml)
        file ".command.*"

    script:
    output_idxml = "${params.sample_name}_idpep.idXML"
    """
    IDPosteriorErrorProbability \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -fit_algorithm:outlier_handling ${params.outlier_handling} \\
        -threads ${task.cpus}
    """
}
