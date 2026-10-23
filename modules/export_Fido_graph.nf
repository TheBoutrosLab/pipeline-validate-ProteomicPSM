process export_Fido_graph {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.intermediate_file_dir}/${params.dataset_id}/Fido/intermediate/${task.process}/" },
        mode: 'copy', enabled: params.save_intermediate_files

    input:
        path(input_idxmls)

    output:
        path("${params.dataset_id}-Fido-graph.txt"), emit: graph
        path("${params.dataset_id}-Fido-target-decoy.txt"), emit: proteins
        path("${params.dataset_id}-Fido-accessions.json"), emit: mapping

    script:
    """
    python3 ${projectDir}/src/FidoBridge.py prepare \\
        --idxml ${input_idxmls.join(' ')} \\
        --decoy-prefix '${params.decoy_string}' \\
        --graph ${params.dataset_id}-Fido-graph.txt \\
        --proteins ${params.dataset_id}-Fido-target-decoy.txt \\
        --mapping ${params.dataset_id}-Fido-accessions.json
    """
}
