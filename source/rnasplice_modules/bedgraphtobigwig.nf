process BEDGRAPHTOBIGWIG {
    label 'process_single'
    publishDir params.outdir

    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1'

    input:
    tuple val(name), path(bedgraph)
    path  sizes

    output:
    tuple val(name), path("*.bigWig"), emit: bigwig

    script:
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${name}.bigWig
    """
}