process integrate_PSM_iProphet {
    container params.docker_image_fragpipe

    publishDir { "${params.final_output_dir}/${params.dataset_id}/ProteinProphet/output" },
        mode: 'copy', pattern: '*.pepXML'

    input:
        path(input_pepxmls)
        path(database_fasta)

    output:
        path("${params.dataset_id}-iProphet.pepXML"), emit: pepxml

    script:
    """
    set -euo pipefail
    philosopher_bin=''
    for candidate in /fragpipe_bin/fragpipe-${params.fragpipe_version}/fragpipe-${params.fragpipe_version}/tools/Philosopher/philosopher*; do
        if [ -f "\$candidate" ] && [ -x "\$candidate" ]; then
            philosopher_bin="\$candidate"
            break
        fi
    done
    if [ -z "\$philosopher_bin" ]; then
        echo 'Philosopher executable not found in the configured FragPipe image' >&2
        exit 127
    fi
    "\$philosopher_bin" workspace --init
    "\$philosopher_bin" database --annotate ${database_fasta} --prefix ${params.decoy_string}
    "\$philosopher_bin" iprophet \\
        --output ${params.dataset_id}-iProphet \\
        --threads ${task.cpus} \\
        ${input_pepxmls.join(' ')}
    mv ${params.dataset_id}-iProphet.pep.xml ${params.dataset_id}-iProphet.pepXML
    """
}
