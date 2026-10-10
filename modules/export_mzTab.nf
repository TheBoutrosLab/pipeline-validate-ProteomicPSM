/**
*  NF module to export idXML to mzTab
*/

process export_mzTab {
    container params.docker_image_openms

    publishDir params.final_output_dir,
        mode: 'copy',
        pattern: '*.mzTab'

    publishDir "${params.process_log_dir}/${task.process.replace(':', '/')}-${task.index}/",
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        file idXML_file

    output:
        file mzTab_file
        file ".command.*"

    script:
    mzTab_file = "${idXML_file.baseName}.mzTab"
    arg_string = ''
    valid_args = [
        'export_all_psms'
    ]
    valid_args.each {
        if (params.containsKey(it)) {
            arg_string += " -${it} ${params[it]}"
        }
    }
    """
    MzTabExporter \\
        -in ${idXML_file} \\
        -out ${mzTab_file} \\
        ${arg_string}
    """
}
