process build_CohortGraph_ProtLoom {
    container params.docker_image_protloom

    publishDir { "${params.final_output_dir}/${params.dataset_id}/ProtLoom-${params.protloom_version}/output" },
        mode: 'copy', pattern: "${params.dataset_id}-cohort_*.json"

    input:
        tuple val(sample_names), path(input_idxmls)
        path(database_fasta)
        path(entries_tsv)
        val(markers)

    output:
        path("${params.dataset_id}-cohort_graph.json"), emit: graph
        path("${params.dataset_id}-cohort_manifest.json")

    script:
    manifest_json = groovy.json.JsonOutput.toJson([
        schema_version: 1,
        database: [output_fasta: database_fasta.name, output_entries_tsv: entries_tsv.name],
        samples: (0..<sample_names.size()).collect { index -> [
            sample_id: sample_names[index], id_files: [input_idxmls[index].name]
        ] }
    ])
    """
    cat > ${params.dataset_id}-cohort_manifest.json <<'JSON'
    ${manifest_json}
    JSON
    protloom jointGraph \\
        --manifest ${params.dataset_id}-cohort_manifest.json \\
        --score-kind pep \\
        --output-graph-json ${params.dataset_id}-cohort_graph.json \\
        --enzyme '${params.protloom_enzyme}' \\
        --decoy-string '${markers.decoy_string}' \\
        --decoy-position '${markers.decoy_position}' \\
        --contam-string '${markers.contam_string}' \\
        --contam-position '${markers.contam_position}'
    """
}
