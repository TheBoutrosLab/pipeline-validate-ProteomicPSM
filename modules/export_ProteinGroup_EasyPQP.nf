process export_ProteinGroup_EasyPQP {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.final_output_dir}/${params.dataset_id}/${params.protein_inference_algorithm}/output" },
        mode: 'copy', pattern: '*.tsv'

    input:
        tuple val(manifest_json), path(input_idxmls), path(input_mzmls)
        path(group_idxml)

    output:
        path("${params.dataset_id}-psm.tsv"), emit: psm
        path("${params.dataset_id}-peptides.tsv"), emit: peptide
        path("${params.dataset_id}-easypqp_export_summary.tsv"), emit: summary

    script:
    min_probability = params.protein_inference_algorithm == 'ProteinProphet' ?
        params.proteinprophet_min_probability :
        (params.protein_inference_algorithm == 'Fido' ?
            params.fido_min_probability : params.epifany_min_probability)
    """
    cat > easypqp_export_manifest.json <<'JSON'
    ${manifest_json}
    JSON
    python3 ${projectDir}/src/ExportProteinGroupEasyPQP.py \\
        --manifest easypqp_export_manifest.json \\
        --groups ${group_idxml} \\
        --psm-qvalue ${params.proteinprophet_easypqp_psm_qvalue} \\
        --min-protein-probability ${min_probability} \\
        --algorithm ${params.protein_inference_algorithm} \\
        --decoy-string '${params.decoy_string}' \\
        --output-prefix ${params.dataset_id}
    """
}
