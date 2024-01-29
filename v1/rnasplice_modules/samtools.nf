process SAMTOOLS {
    label 'ALL'
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"

    input:
    tuple val(sample_name), path(sam_file), val(condition)
    
    output:
    tuple val(sample_name), path("${sam_file.baseName}*bam"), val(condition), emit: bam 
    path("alignment_summary.txt"), emit: log
    script:
    """
    samtools flagstat ${sam_file} > alignment_summary.txt
    samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sam_file.baseName}.sorted.bam -T tmp  -@ ${params.threads} 
    """
    
}