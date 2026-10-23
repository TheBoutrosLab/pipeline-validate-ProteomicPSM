process score_FDR_ProtLoom {
    container params.docker_image_protloom

    publishDir { "${params.final_output_dir}/${sample_name}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: "${sample_name}_*fdr.tsv"

    input:
        tuple val(sample_name), path(input_idxml), path(group_table), path(assignment_table)
        val(markers)

    output:
        tuple val(sample_name), path("${sample_name}_psm_fdr.tsv"),
            path("${sample_name}_protein_group_fdr.tsv")

    script:
    """
    protloom computeFDR \\
        --id-file ${input_idxml} \\
        --group-table ${group_table} \\
        --assignment-table ${assignment_table} \\
        --pt-order '${params.protloom_pt_order}' \\
        --pc-order '${params.protloom_pc_order}' \\
        --decoy-string '${markers.decoy_string}' \\
        --decoy-position '${markers.decoy_position}' \\
        --contam-string '${markers.contam_string}' \\
        --contam-position '${markers.contam_position}' \\
        --output-prefix ${sample_name}
    """
}
