
process SAMTOOLS {
    label 'ALL'
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"
    
    input:
    tuple val(sample_name), path(sam_file)
    
    output:
    tuple val(sample_name), path("${sam_file}.sorted.bam"), emit: bam 
    path("alignment_summary.txt"), emit: log
    
    script:
    """
    samtools view -bS ${sam_file} | samtools sort -o ${sam_file}.sorted.bam -T tmp  
    samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sam_file.baseName}.sorted.bam -T tmp  -@ ${params.threads}
    """
    
}

// process SAMTOOLS {
//     label 'ALL'
//     publishDir params.outdir

//     container "biocontainers/samtools:v1.7.0_cv4"

//     input:
//     tuple val(sample_name), path(sam_file), val(condition)
    
//     output:
//     path("${sam_file.baseName}*bam"), emit: bam 
//     path("alignment_summary.txt"), emit: log
//     script:
//     """
//     samtools flagstat ${sam_file} > alignment_summary.txt
//     samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sam_file.baseName}.sorted.bam -T tmp  -@ ${params.threads} 
//     """
// }
