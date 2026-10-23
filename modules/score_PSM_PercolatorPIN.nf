process score_PSM_PercolatorPIN {
    container params.docker_image_fragpipe

    input:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(input_pins), path(mapping)

    output:
        tuple val(sample_name), path(input_idxml), val(search_engine), path(mapping), path('target.pout'), path('decoy.pout')

    script:
    pin_files = input_pins instanceof List ? input_pins : [input_pins]
    tdc = params.percolator_post_processing_tdc ? '--post-processing-tdc' : ''
    decoy_pattern = params.percolator_decoy_pattern ?: params.decoy_string
    """
    set -euo pipefail
    percolator_bin=\$(find /fragpipe_bin -type f -name percolator -perm -111 -print -quit)
    if [ -z "\$percolator_bin" ]; then
        echo 'Percolator was not found in the configured FragPipe image' >&2
        exit 127
    fi
    first_pin=${pin_files[0]}
    head -n 1 "\$first_pin" > combined.pin
    for pin in ${pin_files.join(' ')}; do
        if [ "\$(head -n 1 "\$pin")" != "\$(head -n 1 "\$first_pin")" ]; then
            echo 'MSBooster PIN headers differ across fractions' >&2
            exit 1
        fi
        tail -n +2 "\$pin" >> combined.pin
    done
    "\$percolator_bin" --only-psms --no-terminate ${tdc} \\
        --num-threads ${task.cpus} \\
        --results-psms target.pout --decoy-results-psms decoy.pout \\
        --protein-decoy-pattern '${decoy_pattern}' combined.pin
    test -s target.pout
    test -s decoy.pout
    """
}
