process BEDCLIP {
    label 'ALL'
    publishDir params.outdir
    container 'quay.io/biocontainers/ucsc-bedclip:377--h0b8a92a_2'

    input:
    tuple val(name), path(bedgraph)
    path  sizes

    output:
    tuple val(name), path("*.bedGraph"), emit: bedgraph

    """
    bedClip \\
        $bedgraph \\
        $sizes \\
        ${name}.bedGraph
    """
}