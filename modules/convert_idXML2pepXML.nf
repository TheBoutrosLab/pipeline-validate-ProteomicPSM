/**
 * Module for convertion from idXML to pepXML
 *
 * Note: In the pepXML output by MSFragger, the `search_engeine` is set to 'X! Tandem' because of
 * compatibility issue with PeptideProphet according to the developers (see
 * https://github.com/Nesvilab/philosopher/discussions/333). OpenMS' IDFileConverter however sets
 * this attribute to MSFragger, so here it is manually converted to 'X! Tandem'. This work around
 * should be removed once MSFragger resolves the compatibility issue. We do not have an ETA on this.
 */

process convert_idXML2pepXML {
    container params.docker_image_openms

    publishDir "${params.intermediate_file_dir}/${task.process.replace(':', '/')}/",
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: '*.pepXML'

    publishDir "${params.process_log_dir}/${task.process.replace(':', '/')}-${task.index}/",
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple(
            file(idxml_file),
            val(search_engine)
        )
        file(mzml_file)

    output:
        tuple(
            file(pepxml_file),
            val(search_engine)
        )
        file '.command.*'

    script:
    pepxml_file = "${params.sample_name}_${search_engine}.pepXML"
    """
    set -euo pipefail

    IDFileConverter \\
        -in ${idxml_file} \\
        -out ${pepxml_file} \\
        -mz_file ${mzml_file} \\
        -peptideprophet_analyzed

    search_engine=\$(cat ${pepxml_file} | grep '<search_summary' | grep -P -o '(?<=search_engine=")[A-Za-z0-9]*')

    # For compatibility issue of MSFragger and PeptidePropphet
    if [[ "\${search_engine}" == 'MSFragger' ]]; then
        sed -i 's/search_engine="MSFragger"/search_engine="X! Tandem" actual_search_engine="MSFragger"/g' ${pepxml_file}
    fi
    """
}
