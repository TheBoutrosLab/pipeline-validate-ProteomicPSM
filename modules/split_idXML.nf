/**
* Split idXML file into separate fractions.
*/

process split_idXML {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/OpenMS-${params.openms_version}/intermediate/${task.process.replace(':','/')}" },
        enabled: params.save_intermediate_files,
        mode: 'copy',
        pattern: "${split_dir}/*.idXML"

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple(
            val(sample_name),
            file(input_idxml),
            val(search_engine)
        )

    output:
        tuple(
            val(sample_name),
            file("${split_dir}/*.idXML"),
            val(search_engine)
        )
        file ".command.*"

    script:
    split_dir = "split"
    """
    mkdir -p ${split_dir}
    IDRipper \\
        -in ${input_idxml} \\
        -out ${split_dir} \\
        -threads ${task.cpus}

    for file in \$(ls ${split_dir}); do
        source=${split_dir}/\$file
        target=${split_dir}/\$(basename \$file .idXML)_${search_engine}.idXML
        mv \$source \$target
    done
    """
}
