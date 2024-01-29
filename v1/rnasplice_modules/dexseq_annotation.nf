
process DEXSEQ_ANNOTATION {
    publishDir params.outdir
    tag "$gtf"
    label 'ALL'

    container 'quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0'

    input:
    path gtf         // path gtf file
    //val aggregation  // val params.aggregation

    output:
    path "*.gff"        , emit: gff

    script:

    """
    ${params.basedir}/bin/dexseq_prepare_annotation.py $gtf DEXSeq.gff
    """
}

