process BEDTOOLS_GENOMECOV {
    label 'ALL'
    publishDir params.outdir

    container 'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

    input:
    tuple val(name), path(bam), val(condition)

    output:
    tuple val(name), path("*.forward.bedGraph"), emit: bedgraph_forward
    tuple val(name), path("*.reverse.bedGraph"), emit: bedgraph_reverse

    script:
    def prefix_forward = "${name}.forward"
    def prefix_reverse = "${name}.reverse"
    if (params.strand == 'reverse') {
        prefix_forward = "${name}.reverse"
        prefix_reverse = "${name}.forward"
    }
    """
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand + \\
        | bedtools sort > ${prefix_forward}.bedGraph
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand - \\
        | bedtools sort > ${prefix_reverse}.bedGraph
    """
}
