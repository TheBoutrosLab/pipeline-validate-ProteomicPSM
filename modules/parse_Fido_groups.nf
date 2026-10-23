process parse_Fido_groups {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.final_output_dir}/${params.dataset_id}/Fido/output" },
        mode: 'copy', pattern: "${params.dataset_id}-Fido-groups.json"

    input:
        path(fido_results)
        path(accession_mapping)

    output:
        path("${params.dataset_id}-Fido-groups.json"), emit: groups

    script:
    """
    python3 ${projectDir}/src/FidoBridge.py parse --results ${fido_results} --mapping ${accession_mapping} --output ${params.dataset_id}-Fido-groups.json
    """
}
