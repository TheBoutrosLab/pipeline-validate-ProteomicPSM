process score_FDR_ProtLoomSplit {
    container params.docker_image_protloom

    publishDir { "${params.final_output_dir}/${sample_name}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: '*_psm_fdr.*'

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy', pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine)
        path(entries_tsv)
        val(markers)

    output:
        tuple val(sample_name), path("${sample_name}_${search_engine}_psm_fdr.idXML"),
            val(search_engine)
        path "${sample_name}_${search_engine}_psm_fdr.tsv", emit: table
        path '.command.*'

    script:
    output_prefix = "${sample_name}_${search_engine}"
    """
    protloom computeFDR \\
        --psm-score-source id-file \\
        --id-file ${input_idxml} \\
        --entries-tsv ${entries_tsv} \\
        --pt-order '${params.protloom_pt_order}' \\
        --pc-order '${params.protloom_pc_order}' \\
        --decoy-string '${markers.decoy_string}' \\
        --decoy-position '${markers.decoy_position}' \\
        --contam-string '${markers.contam_string}' \\
        --contam-position '${markers.contam_position}' \\
        --output-prefix ${output_prefix}
    """
}
