process infer_ProteinGroups_Fido {
    container params.docker_image_fido

    publishDir { "${params.final_output_dir}/${params.dataset_id}/Fido/output" },
        mode: 'copy', pattern: "${params.dataset_id}-Fido.txt"

    input:
        path(input_graph)
        path(target_decoy)

    output:
        path("${params.dataset_id}-Fido.txt"), emit: results

    script:
    """
    set -euo pipefail
    command -v FidoChooseParameters >/dev/null
    FidoChooseParameters \\
        -g \\
        ${input_graph} \\
        ${target_decoy} \\
        > ${params.dataset_id}-Fido.txt
    test -s ${params.dataset_id}-Fido.txt
    """
}
