process switch_score_2_pep {
    container params.docker_image_openms

    publishDir { "${params.final_output_dir}/${sample_name}/OpenMS-${params.openms_version}/output" },
        mode: 'copy',
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
    output_idxml = "${sample_name}_${search_engine}_idswitch_pep.idXML"
    old_score = params?['fdr_algorithm_no_qvalues'] == true ? 'FDR' : 'q-value'
    """
    IDScoreSwitcher \
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -old_score '${old_score}' \\
        -new_score 'Posterior Error Probability_score' \\
        -new_score_type 'Posterior Error Probability' \\
        -new_score_orientation 'lower_better'
    """
}
