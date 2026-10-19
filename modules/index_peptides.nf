include { getCommandArgs } from "${moduleDir}/common"

PARAMETER_MAPPER = [
    'enzyme': '-enzyme:name',
    'decoy_string': '-decoy_string',
    'decoy_string_position': '-decoy_string_position',
    'index_peptides_enzyme_specifity': '-enzyme:specificity',
    'index_peptides_missing_decoy_action': '-missing_decoy_action',
    'index_peptides_aaa_max': '-aaa_max',
    'index_peptides_mismatches_max': '-mismatches_max',
    'index_peptides_unmatched_action': '-unmatched_action',
    'allow_nterm_protein_cleavage': '-allow_nterm_protein_cleavage'
]

FLAG_MAPPER = [
    'index_peptides_IL_equivalent': '-IL_equivalent',
    'index_peptides_write_protein_sequence': '-write_protein_sequence',
    'index_peptides_write_protein_description': '-write_protein_description',
    'index_peptides_keep_unreferenced_proteins': '-keep_unreferenced_proteins'
]


process index_peptides {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/OpenMS-${params.openms_version}/intermediate/${task.process.replace(':', '/')}/" },
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: "*.idXML"

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple(
            val(sample_name),
            file(input_idxml),
            val(search_engine)
        )
        file decoy_database

    output:
        tuple(
            val(sample_name),
            file(output_idxml),
            val(search_engine)
        )
        file ".command.*"

    script:
    output_idxml = "${sample_name}_${search_engine}_PeptideIndexer.idXML"
    extra_args = getCommandArgs(PARAMETER_MAPPER, FLAG_MAPPER)
    """
    PeptideIndexer \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -fasta ${decoy_database} \\
        -threads ${task.cpus} \\
        ${extra_args}
    """
}
