process FASTQC {
    label 'fastqc'
    publishDir params.outdir
    container "staphb/fastqc:0.11.9"

    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}", emit: zip

    script:
    """
    fastqc ${reads[0]} ${reads[1]} --thread ${params.threads}
    """
}
