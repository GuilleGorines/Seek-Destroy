process PREPARE_KRAKEN_REPORT {
    tag "$archive"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path("*.krona"), emit: untar
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def kraken_report_krona_rdy = "${meta.id}.krona"

    """
    cat ${kraken_report} | cut -f 2,3 > ${kraken_report_krona_rdy}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(echo \$(cut --version 2>&1) | head -n 1 | sed 's/^.*(GNU coreutils) //')
    END_VERSIONS
    """