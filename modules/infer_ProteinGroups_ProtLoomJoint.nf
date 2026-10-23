process infer_ProteinGroups_ProtLoomJoint {
    container params.docker_image_protloom

    publishDir { "${params.final_output_dir}/${sample_name}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: "${sample_name}_*.tsv"
    publishDir { "${params.final_output_dir}/${sample_name}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: "${sample_name}_run.json"

    input:
        tuple val(sample_name), path(input_idxml), path(cohort_graph), val(annotation_id)
        path(database_fasta)
        path(entries_tsv)
        val(annotation_ids)
        path(annotation_gtfs)
        val(markers)

    output:
        tuple val(sample_name), path(input_idxml), path("${sample_name}_groups.tsv"),
            path("${sample_name}_assignments.tsv"), emit: inference
        path("${sample_name}_peptides.tsv")
        path("${sample_name}_run.json")

    script:
    manifest_json = groovy.json.JsonOutput.toJson([
        schema_version: 1,
        database: [output_fasta: database_fasta.name, output_entries_tsv: entries_tsv.name],
        annotations: (0..<annotation_ids.size()).collect { index ->
            [annotation_id: annotation_ids[index], gtf_path: annotation_gtfs[index].name]
        },
        samples: [[sample_id: sample_name, id_files: [input_idxml.name],
                   annotation_id: annotation_id]]
    ])
    """
    cat > project.json <<'JSON'
    ${manifest_json}
    JSON
    protloom jointInfer \\
        --graph-json ${cohort_graph} \\
        --manifest project.json \\
        --sample-id '${sample_name}' \\
        --score-kind pep \\
        --enzyme '${params.protloom_enzyme}' \\
        --pt-order '${params.protloom_pt_order}' \\
        --pc-order '${params.protloom_pc_order}' \\
        --decoy-string '${markers.decoy_string}' \\
        --decoy-position '${markers.decoy_position}' \\
        --contam-string '${markers.contam_string}' \\
        --contam-position '${markers.contam_position}' \\
        --output-prefix ${sample_name} \\
        --disable-ccs
    """
}
