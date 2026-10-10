
process merge_PSM {
    container params.docker_image_openms

    publishDir "${params.intermediate_file_dir}/${task.process.replace(':', '/')}/",
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: "*.idXML"

    publishDir "${params.process_log_dir}/${task.process.replace(':', '/')}-${task.index}/",
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        file input_idxmls
        val pepxml_protxml
        val fractionated

    output:
        file(output_idxml)
        file ".command.*"

    script:
    if (fractionated == true) {
        // file name pattern: <original_file_basename>_<search_engine>_<FixSpectraDatapath>
        // e.g., CPT0001_reindex_1_MSFragger_FixSpectraDataPath.idXML
        output_idxml = "${input_idxmls[0].baseName.replaceFirst('_[a-zA-Z0-9]+_[a-zA-Z0-9]+$', '')}_merged.idXML"
    } else {
        output_idxml = "${params.sample_name}_merged.idXML"
    }
    
    arg_pepxml_protxml = pepxml_protxml ? ' -pepxml_protxml' : ''
    """
    IDMerger \\
        -in ${input_idxmls} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -debug 1
    """
}
