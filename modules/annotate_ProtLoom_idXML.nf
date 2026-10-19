process annotate_ProtLoom_idXML {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    publishDir { "${params.final_output_dir}/${sample_name}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: "${sample_name}_protloom.idXML"

    input:
        tuple val(sample_name), path(input_idxml), path(assignment_table),
            path(psm_fdr_table), path(group_fdr_table)

    output:
        tuple val(sample_name), path("${sample_name}_protloom.idXML"),
            path("${sample_name}_protloom_mapped.idXML"), emit: annotated

    script:
    """
    python3 ${projectDir}/src/AnnotateIdXMLWithProtLoom.py \\
        --idxml ${input_idxml} \\
        --assignments ${assignment_table} \\
        --psms ${psm_fdr_table} \\
        --groups ${group_fdr_table} \\
        --primary-score ${params.protloom_idxml_primary_score} \\
        --unmatched-psms ${params.protloom_idxml_unmatched_psms} \\
        --matched-output ${sample_name}_protloom_mapped.idXML \\
        --output ${sample_name}_protloom.idXML
    """
}
