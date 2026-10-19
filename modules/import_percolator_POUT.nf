process import_percolator_POUT {
    containerOptions params.project_script_mount_options
    container params.docker_image_pyopenms

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(mapping), path(target_pout), path(decoy_pout)

    output:
        tuple val(sample_name), path(output_idxml), val(search_engine)

    script:
    output_idxml = "${sample_name}_${search_engine}_percolator.idXML"
    """
    python3 ${projectDir}/src/IdXMLPinBridge.py integrate \\
        --idxml ${input_idxml} --mapping ${mapping} \\
        --target ${target_pout} --decoy ${decoy_pout} --output ${output_idxml}
    """
}
