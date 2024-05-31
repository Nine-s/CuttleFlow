process GFFREAD_TX2GENE {
    tag "$gtf"
    publishDir params.outdir
    label 'ALL'

    container "zavolab/gffread:0.11.7"
    
    input:
    path gtf

    output:
    path "*.tx2gene.tsv" , emit: tx2gene

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--table transcript_id,gene_id,gene_name'
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    gffread $args $gtf | sort -u 1> ${prefix}.tx2gene.tsv
    """
}
