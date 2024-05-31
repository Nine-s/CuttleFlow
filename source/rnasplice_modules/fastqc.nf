process FASTQC {
    label 'ALL'
    tag "fastqc $sample_id"
    publishDir params.outdir
    container "staphb/fastqc:0.11.9"

    input:
    tuple val(sample_id), path(reads), val(condition)

    output:
    path("*_fastqc.zip"), emit: zip


    script:
    """
    fastqc ${reads[0]} ${reads[1]} --thread ${params.threads}
    """
}