process extract_PSMFeatures_MSBooster {
    container params.docker_image_fragpipe

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(input_pins), path(mapping), path(input_mzmls)
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
