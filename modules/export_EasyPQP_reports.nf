process export_EasyPQP_reports {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.final_output_dir}/${params.dataset_id}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: '*.tsv'

    input:
        tuple val(manifest_json), path(input_idxmls), path(assignment_tables),
            path(psm_fdr_tables), path(group_fdr_tables), path(input_mzmls)

    output:
        path("${params.dataset_id}-psm.tsv"), emit: psm
        path("${params.dataset_id}-peptides.tsv"), emit: peptide
        path("${params.dataset_id}-easypqp_export_summary.tsv"), emit: summary

    script:
    """
    cat > easypqp_export_manifest.json <<'JSON'
    ${manifest_json}
    JSON
    python3 ${projectDir}/src/ExportEasyPQPReports.py \\
        --manifest easypqp_export_manifest.json \\
        --psm-qvalue ${params.easypqp_psm_qvalue} \\
        --group-qvalue ${params.protloom_easypqp_group_qvalue} \\
        --output-prefix ${params.dataset_id}
    """
}
