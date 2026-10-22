process convert_ProteinProphet_idXML {
    container params.docker_image_openms

    publishDir { "${params.final_output_dir}/${params.dataset_id}/ProteinProphet/output" },
        mode: 'copy', pattern: '*.idXML'

    input:
        path(input_protxml)

    output:
        path("${params.dataset_id}-ProteinProphet.idXML"), emit: groups

    script:
    """
    IDFileConverter \\
        -in ${input_protxml} \\
        -out ${params.dataset_id}-ProteinProphet.idXML
    """
}
