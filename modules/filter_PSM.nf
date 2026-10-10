include { getCommandArgs } from "${moduleDir}/common"

ADDITIONAL_MAPPER = [
    'filter_psm_whitelist_proteins': '-whitelist:proteins',
    'filter_psm_whitelist_peptides': '-whitelist:peptides',
    'filter_psm_blacklist_proteins': '-blacklist:proteins',
    'filter_psm_blacklist_peptides': '-blacklist:peptides',
    'filter_psm_in_silico_digestion': '-in_silico_digestion:fasta'
]

PARAMETER_MAPPER = [
    'filter_psm_precursor_rt': '-precursor:rt',
    'filter_psm_precursor_mz': '-precursor:mz',
    'filter_psm_precursor_length': '-precursor:length',
    'filter_psm_precursor_charge': '-precursor:charge',
    'filter_psm_score_peptide': '-score:psm',
    'filter_psm_whitelist_protein_accessions': '-whitelist:protein_accessions',
    'filter_psm_whitelist_modifications': '-whitelist:modifications',
    'filter_psm_whitelist_ignore_modifications': '-whitelist:ignore_modifications',
    'filter_psm_blacklist_protein_accessions': '-blacklist:protein_accessions',
    'filter_psm_blacklist_modifications': '-blacklist:modifications',
    'filter_psm_blacklist_ignore_modifications': '-blacklist:ignore_modifications',
    'filter_psm_number_of_missed_cleavages': '-missed_cleavages:number_of_missed_cleavages',
    'filter_psm_best_n_spectra': '-best:n_spectra',
    'filter_psm_best_n_peptide_hits': '-best:n_peptide_hits',
    'filter_psm_best_n_protein_hits': '-best:n_protein_hits'
]

FLAG_MAPPER = [
    'filter_psm_best_strict': '-best:strict',
    'filter_psm_var_mods': '-var_mods',
    'filter_psm_remove_shared_peptides': '-remove_shared_peptides',
    'filter_psm_keep_unreferenced_protein_hits': '-keep_unreferenced_protein_hits',
    'filter_psm_remove_decoys': '-remove_decoys',
    'filter_psm_delete_unreferenced_peptide_hits': '-delete_unreferenced_peptide_hits'
]

process filter_PSM {
    container params.docker_image_openms

    publishDir params.final_output_dir,
        mode: 'copy',
        pattern: "*.idXML"

    publishDir "${params.process_log_dir}/${task.process.replace(':', '/')}-${task.index}/",
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple(
            file(input_idxml),
            val(search_engine)
        )
        path whitelist_proteins, stageAs: 'whitelist_proteins/*'
        path whitelist_peptides, stageAs: 'whitelist_peptides/*'
        path blacklist_proteins, stageAs: 'blacklist_proteins/*'
        path blacklist_peptides, stageAs: 'blacklist_peptides/*'
        path in_silico_digestion_fasta, stageAs: 'in_silico_digestion_fasta/*'

    output:
        file output_idxml
        file ".command.*"

    script:
    output_idxml = "${params.sample_name}_${search_engine}_filter.idXML"
    extra_args = getCommandArgs(PARAMETER_MAPPER, FLAG_MAPPER)
    if (whitelist_proteins.name.tokenize('/').last() != 'NO_FILE') {
        extra_args += " -whitelist:proteins ${whitelist_proteins}"
    }
    if (whitelist_peptides.name.tokenize('/').last() != 'NO_FILE') {
        extra_args += " -whitelist:peptides ${whitelist_peptides}"
    }
    if (blacklist_proteins.name.tokenize('/').last() != 'NO_FILE') {
        extra_args += " -blacklist:proteins ${blacklist_proteins}"
    }
    if (blacklist_peptides.name.tokenize('/').last() != 'NO_FILE') {
        extra_args += " -blacklist:peptides ${blacklist_peptides}"
    }
    if (in_silico_digestion_fasta.name.tokenize('/').last() != 'NO_FILE') {
        extra_args += " -in_silico_digestion:fasta ${in_silico_digestion_fasta}"
    }
    """
    IDFilter \\
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -missed_cleavages:enzyme ${params.enzyme} \\
        ${extra_args}
    """
}
