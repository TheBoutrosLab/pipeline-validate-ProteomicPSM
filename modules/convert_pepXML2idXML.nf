/**
 * Module for convertion from pepXML to idXML
 */

process convert_pepXML2idXML {
    container params.docker_image_openms

    publishDir { "${params.intermediate_file_dir}/${sample_name}/OpenMS-${params.openms_version}/intermediate/${task.process.replace(':', '/')}/" },
        mode: 'copy',
        enabled: params.save_intermediate_files,
        pattern: '*.idXML'

    publishDir { "${params.process_log_dir}/${sample_name}/${task.process.replace(':', '/')}-${task.index}/" },
        mode: 'copy',
        pattern: '.command.*',
        saveAs: { "log${file(it).name}" }

    input:
        tuple(
            val(sample_name),
            file(pepXML),
            val(search_engine)
        )
        file(mzML)

    output:
        tuple(
            val(sample_name),
            file(idXML),
            val(search_engine)
        )
        file '.command.*'

    script:
    pepXML_fixed = "${pepXML.baseName}_fixed.pepXML"
    idXML = "${sample_name}_${search_engine}.idXML"
    """
    set -euo pipefail

    # search_engine is set to 'X! Tandem' when converting to pepXML for compatibility issue between
    # MSFragger and PeptideProphet, and is converted back here.
    if grep -q 'search_engine="X! Tandem" actual_search_engine="MSFragger"' ${pepXML}; then
        cat ${pepXML} | sed 's/search_engine="X! Tandem" actual_search_engine="MSFragger"/search_engine="MSFragger"/g' > ${pepXML_fixed}
    else
        cp -L ${pepXML} ${pepXML_fixed}
    fi

    IDFileConverter \\
        -in ${pepXML_fixed} \\
        -out ${idXML} \\
        -mz_file ${mzML}
    if ! grep -Fq 'name="spectra_data"' ${idXML}; then
        sed -i '/<\\/ProteinIdentification>/ i \\${' ' * 24}<UserParam type="stringList" name="spectra_data" value="[${mzML.name}]"\\/>' ${idXML}
    fi
    sed -i 's/higher_score_better="true"/higher_score_better="false"/g' ${idXML}
    """
}
