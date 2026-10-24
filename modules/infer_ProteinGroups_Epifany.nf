process infer_ProteinGroups_Epifany {
    container params.docker_image_openms

    publishDir { "${params.final_output_dir}/${params.dataset_id}/Epifany/output" },
        mode: 'copy', pattern: "${params.dataset_id}-Epifany.idXML"

    input:
        path(input_idxmls)

    output:
        path("${params.dataset_id}-Epifany.idXML"), emit: groups

    script:
    """
    Epifany \\
        -in ${input_idxmls.join(' ')} \\
        -out ${params.dataset_id}-Epifany.idXML \\
        -greedy_group_resolution ${params.epifany_greedy_group_resolution} \\
        -threads ${task.cpus}
    """
}
