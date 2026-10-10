include { getCommandArgs } from "${moduleDir}/common"

PARAMETER_MAPPER = [
    'percolator_decoy_pattern': '-decoy_pattern'
]

FLAG_MAPPER = [
    'percolator_peptide_level_fdrs': '-peptide_level_fdrs',
    'percolator_protein_level_fdrs': '-protein_level_fdrs',
    'percolator_post_processing_tdc': '-post_processing_tdc'
]

process score_PEP_percolator {
    container params.docker_image_openms

    publishDir "${params.intermediate_file_dir}/${task.process.replace(':', '/')}/",
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: "*.idXML"

    publishDir "${params.process_log_dir}/${task.process.replace(':', '/')}/",
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "${input_idxml.baseName}/log${file(it).name}" }

    input:
        tuple(
            file(input_idxml),
            val(search_engine)
        )

    output:
        tuple(
            file(output_idxml),
            val(search_engine)
        )
        file ".command.*"

    script:
    output_idxml = "${params.sample_name}_${search_engine}_percolator.idXML"
    extra_args = getCommandArgs(PARAMETER_MAPPER, FLAG_MAPPER)
    """
    OMP_NUM_THREADS=${task.cpus} PercolatorAdapter \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -score_type pep \\
        -enzyme ${params.enzyme_percolator} \\
        ${extra_args}
    """
}
