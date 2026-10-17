include { getCommandArgs } from "${moduleDir}/common"

PARAMETER_MAPPER = [
    'fdr_algorithm_conservative': '-algorithm:conservative'
]

FLAG_MAPPER = [
    'fdr_algorithm_use_all_hits': '-algorithm:use_all_hits',
    'fdr_algorithm_treat_runs_separately': '-algorithm:treat_runs_separately',
    'fdr_algorithm_split_charge_variants': '-algorithm:split_charge_variants',
    'fdr_algorithm_no_qvalues': '-algorithm:no_qvalues'
]

process score_FDR_OpenMS {
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
    output_idxml = "${sample_name}_${search_engine}_FDR.idXML"
    extra_args = getCommandArgs(PARAMETER_MAPPER, FLAG_MAPPER)
    """
    FalseDiscoveryRate \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -PSM true \\
        -protein false \\
        -algorithm:add_decoy_peptides \\
        -algorithm:add_decoy_proteins \\
        ${extra_args}
    """
}
