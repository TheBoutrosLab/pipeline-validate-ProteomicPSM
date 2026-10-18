process import_percolator_POUT {
    container params.docker_image_pyopenms

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(mapping), path(target_pout), path(decoy_pout)
        path(helper_script)

    output:
        tuple val(sample_name), path(output_idxml), val(search_engine)

    script:
    output_idxml = "${sample_name}_${search_engine}_percolator.idXML"
    """
    python3 ${helper_script} integrate \\
        --idxml ${input_idxml} --mapping ${mapping} \\
        --target ${target_pout} --decoy ${decoy_pout} --output ${output_idxml}
    """
}
