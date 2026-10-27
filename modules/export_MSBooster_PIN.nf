process export_MSBooster_PIN {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(input_mzmls), val(fraction_labels)

    output:
        tuple val(sample_name), path(input_idxml), val(search_engine), path('run_*.pin'),
            path('mapping.tsv'), path('run_*.mzML'), val(fraction_labels)

    script:
    raw_files = input_mzmls instanceof List ? input_mzmls : [input_mzmls]
    """
    python3 ${projectDir}/src/IdXMLPinBridge.py export \\
        --idxml ${input_idxml} --mzml ${raw_files.join(' ')} \\
        --mapping mapping.tsv --enzyme '${params.enzyme_percolator}' \\
        --search-engine '${search_engine}'
    """
}
