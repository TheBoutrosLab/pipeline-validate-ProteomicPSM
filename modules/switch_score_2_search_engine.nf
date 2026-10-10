
SEARCH_ENGINE_SCORE_MAP = [
    'comet': [
        'name': 'MS:1002257',
        'type': 'expect',
        'orientation': 'lower_better'
    ],
    'xtandem': [
        'name': 'hyperscore',
        'type': 'XTandem',
        'orientation': 'higher_better'
    ],
    'msfragger': [
        'name': 'MS:1001330',
        'type': 'expect',
        'orientation': 'lower_better'
    ]
]

process switch_score_2_search_engine {
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
    output_idxml = "${params.sample_name}_${search_engine}_idswitch_search_engine.idXML"
    search_engine_clean = search_engine.toLowerCase().replace('!', '')
    """
    IDScoreSwitcher \
        -in ${input_idxml} \\
        -out ${output_idxml} \\
        -threads ${task.cpus} \\
        -old_score 'Posterior Error Probability' \\
        -new_score ${SEARCH_ENGINE_SCORE_MAP[search_engine_clean]['name']} \\
        -new_score_type ${SEARCH_ENGINE_SCORE_MAP[search_engine_clean]['type']} \\
        -new_score_orientation ${SEARCH_ENGINE_SCORE_MAP[search_engine_clean]['orientation']}
    """
}
