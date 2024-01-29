process UNTAR {
    tag "$archive"
    label 'ALL'

    conda "conda-forge::sed=4.7 bioconda::grep=3.4 conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(name), path(archive)

    output:
    tuple val(name), path("${name}"), emit: untar
    path "versions.yml"             , emit: versions


    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    
    """
    mkdir ${name}

    ## Ensures --strip-components only applied when top level of tar contents is a directory
    ## If just files or multiple directories, place all in prefix
    if [[ \$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -C ${name} --strip-components 1 \\
            -xavf \\
            $args \\
            $archive \\
            $args2
    else
        tar \\
            -C ${name} \\
            -xavf \\
            $args \\
            $archive \\
            $args2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.toString().replaceFirst(/\.[^\.]+(.gz)?$/, ""))
    """
    mkdir ${name}
    touch ${prefix}/file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
