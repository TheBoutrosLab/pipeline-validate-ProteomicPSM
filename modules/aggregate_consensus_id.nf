include { getCommandArgs } from "${moduleDir}/common"

PARAMETER_MAPPER = [
    'consensusid_rt_delta': '-rt_delta',
    'consensusid_mz_delta': '-mz_delta',
    'consensusid_considered_hits': '-filter:considered_hits',
    'consensusid_min_support': '-filter:min_support'
]

FLAG_MAPPER = [
    'consensusid_keep_old_scores': '-filter:keep_old_scores',
    'consensusid_count_empty': '-filter:count_empty'
]

process aggregate_consensus_id {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/OpenMS-${params.openms_version}/intermediate/${task.process.replace(':', '/')}/" },
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: "*.idXML"

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple val(sample_name), file(input_idxml)
        val fractionated

    output:
        tuple val(sample_name), file(output_idxml)
        file ".command.*"

    script:
    if (fractionated == true) {
        output_idxml = "${input_idxml.baseName}_consensus.idXML"
    } else {
        output_idxml = "${sample_name}_consensus.idXML"
    }
    extra_args = getCommandArgs(PARAMETER_MAPPER, FLAG_MAPPER)
    algorithm = params.get('consensusid_algorithm') ?: (params.run_percolator ? 'PEPMatrix' : 'ranks')
    """
    ConsensusID \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -algorithm ${algorithm} \\
        ${extra_args}
    """
}
