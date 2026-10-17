/**
* This nextflow module takes in an idXML file and replace the `spectra_data` UserParam attribute to
* a specific value to make it consistent in all search engine results. This is to make ConsensusID
* to be able to match PSMs.
*/

process fix_idXMLSpectraDataPath {
    container params.docker_image_pyopenms

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
        val fractionated
        path helper_script

    output:
        tuple(
            val(sample_name),
            file(output_idxml),
            val(search_engine)
        )
        file ".command.*"

    script:
    if (fractionated == true) {
        output_idxml = "${input_idxml.baseName}_FixSpectraDataPath.idXML"
    } else {
        output_idxml = "${sample_name}_${search_engine}_FixSpectraDataPath.idXML"
    }

    """
    python ${helper_script} \\
        --input-path ${input_idxml} \\
        --output-path ${output_idxml} \\
        --mzML ${sample_name}_merged.mzML
    """
}
