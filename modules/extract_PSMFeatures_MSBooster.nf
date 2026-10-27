process extract_PSMFeatures_MSBooster {
    container params.docker_image_fragpipe

    publishDir { "${params.final_output_dir}/${sample_name}/MSBooster/output" },
        mode: 'copy', pattern: 'run_*_edited.pin',
        saveAs: { source ->
            def name = file(source).name
            def match = name =~ /^run_(\d+)_edited\.pin$/
            if (!match.matches()) {
                error "Unexpected MSBooster PIN name: ${name}"
            }
            def index = match[0][1].toInteger() - 1
            if (index < 0 || index >= fraction_labels.size()) {
                error "MSBooster PIN ${name} has no matching fraction in the mzML CSV."
            }
            def engine = java.net.URLEncoder.encode(search_engine, 'UTF-8')
            def fraction = java.net.URLEncoder.encode(fraction_labels[index], 'UTF-8')
            "engine_${engine}/run_${match[0][1]}-fraction_${fraction}_edited.pin"
        }

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(input_pins),
            path(mapping), path(input_mzmls), val(fraction_labels)
        path(msbooster_params)

    output:
        tuple val(sample_name), path(input_idxml), val(search_engine), path('run_*_edited.pin'), path(mapping)

    script:
    pin_files = input_pins instanceof List ? input_pins : [input_pins]
    """
    set -euo pipefail
    msbooster_jar=\$(find /fragpipe_bin -type f -name 'MSBooster-*.jar' -print -quit)
    batmass_jar=\$(find /fragpipe_bin -type f -name 'batmass-io-*.jar' -print -quit)
    diann_bin=\$(find /fragpipe_bin -type f -iname 'diann*' -perm -111 -print -quit)
    if [ -z "\$msbooster_jar" ] || [ -z "\$batmass_jar" ] || [ -z "\$diann_bin" ]; then
        echo 'MSBooster, batmass-io, or DIA-NN was not found in the configured FragPipe image' >&2
        exit 127
    fi
    for pin in ${pin_files.join(' ')}; do
        cp -L "\$pin" "\$pin.local"
        mv -f "\$pin.local" "\$pin"
    done
    cp -L ${msbooster_params} msbooster-run.params
    printf '\npinPepXMLDirectory=%s\nmzmlDirectory=%s\nDiaNN=%s\n' "\$PWD" "\$PWD" "\$diann_bin" >> msbooster-run.params
    java -Djava.awt.headless=true -cp "\$msbooster_jar:\$batmass_jar" mainsteps.MainClass --paramsList msbooster-run.params
    for pin in ${pin_files.join(' ')}; do
        test -s "\${pin%.pin}_edited.pin"
    done
    """
}
