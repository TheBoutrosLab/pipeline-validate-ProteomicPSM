/**
 * Common groovy functions
 */

/**
 * Get extra arguments from the params ConfigMap object to be parsed to processes
 */
def getCommandArgs(Map parameter_mapper, Map flag_mapper) {
    def args = ""
    parameter_mapper.each { key, val ->
        if (params.containsKey(key) && params[key] != null) {
            args += " ${val} ${params[key]}"
        }
    }
    flag_mapper.each { key, val ->
        if (params.containsKey(key) && params[key] == true) {
            args += " ${val}"
        }
    }
    return args
}

def printPrelogue(){
    def options = [:]
    exclude_args = [
        'input_csv', 'decoy_database', 'output_dir', 'openms_version',
        'docker_image_openms', 'docker_image_pyopenms',
        'docker_image_fragpipe', 'docker_image_protloom',
        'max_cpus', 'max_memory', 'philosopher_version', 'intermediate_file_dir', 'process_log_dir',
        'pipeline_log_dir'
    ]
    params.each { key, val ->
        if (key in exclude_args) {
            return
        }
        if ('-' in key) {
            return
        }
        val = val as String
        options[key] = val
    }
    def options_str = ''
    def max_len = options.keySet().collect{it.size()}.max()
    options.each { key, val ->
        options_str += '\n' + ' ' * 12 + key.padRight(max_len + 1, ' ') + ": ${val}"
    }
    def images = [OpenMS: params.docker_image_openms,
                  pyOpenMS: params.docker_image_pyopenms]
    if ((params.run_percolator && params.percolator_workflow == 'msbooster') ||
        (params.export_easypqp_reports && !params.run_protloom)) {
        images.FragPipe = params.docker_image_fragpipe
    }
    if (params.run_protloom) {
        images.ProtLoom = params.docker_image_protloom
    }
    def images_str = images.collect { name, image ->
        "${name.padRight(9)}: ${image}"
    }.join('\n            ')
    log.info """\
        =============================================================
        P I P E L I N E   V A L I D A T E   P R O T E O M I C   P S M
        =============================================================
        Boutros Lab

        Current Configuration:
        - pipeline:
            name   : ${workflow.manifest.name}
            version: ${workflow.manifest.version}

        - input:
            input_csv      : ${params.input_csv}
            decoy_database : ${params.decoy_database}

        - output:
            output_dir: ${params.output_dir}

        - options: ${options_str}

        Tools Used:
            ${images_str}

        ------------------------------------
        Starting workflow...
        ------------------------------------
    """
    .stripIndent()
}
