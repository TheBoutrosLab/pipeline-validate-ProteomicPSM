/**
 * Nextflow module to call PeptideProphet
 */

include { getCommandArgs } from "${moduleDir}/common"

PARAMETER_MAPPER = [
    'enzyme': '--enzyme',
    'peptideprophet_clevel': '--clevel',
    'decoy_string': '--decoy',
    'peptideprophet_ignorechg': '--ignorechg',
    'peptideprophet_masswidth': '--masswidth',
    'peptideprophet_minpeplen': '--minpeplen',
    'peptideprophet_minprob': '--minprob'
]

FLAG_MAPPER = [
    'peptideprophet_accmass': '--accmass',
    'peptideprophet_combine': '--combine',
    'peptideprophet_decoyprobs': '--decoyprobs',
    'peptideprophet_expectscore': '--expectscore',
    'peptideprophet_glyc': '--glyc',
    'peptideprophet_nomass': '--nomass',
    'peptideprophet_nonmc': '--nonmc',
    'peptideprophet_nonparam': '--nonparam',
    'peptideprophet_nontt': '--nontt',
    'peptideprophet_phospho': '--phospho',
    'peptideprophet_ppm': '--ppm'
]

process run_PeptideProphet {
    container params.docker_image_philosopher

    publishDir { "${params.final_output_dir}/${sample_name}/OpenMS-${params.openms_version}/output" },
        mode: 'copy',
        enabled: params.save_peptideprophet_pepxml,
        pattern: '*_PeptideProphet.pepXML'

    publishDir { "${params.intermediate_file_dir}/${sample_name}/OpenMS-${params.openms_version}/intermediate/${task.process.replace(':', '/')}/" },
        mode: 'copy',
        enabled: params.save_intermediate_files && !(params.save_peptideprophet_pepxml),
        pattern: '*_PeptideProphet.pepXML'

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple(
            val(sample_name),
            file(input_pepxml),
            val(search_engine)
        )
        file decoy_database

    output:
        tuple(
            val(sample_name),
            file(output_pepxml),
            val(search_engine)
        )
        file ".command.*"

    script:
    output_pepxml = "${sample_name}_${search_engine}_PeptideProphet.pepXML"
    extra_args = getCommandArgs(PARAMETER_MAPPER, FLAG_MAPPER)
    """
    philosopher workspace --init
    philosopher database --annotate ${params.decoy_database} --prefix ${params.decoy_string}
    philosopher peptideprophet \\
        ${extra_args} \\
        --database ${params.decoy_database} \\
        --output 'PeptideProphet' \\
        ${input_pepxml}
    mv PeptideProphet-${input_pepxml.baseName}.pep.xml ${output_pepxml}
    """
}
