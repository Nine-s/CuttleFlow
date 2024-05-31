// process SAMTOOLS {
//     publishDir params.outdir

//     container "biocontainers/samtools:v1.7.0_cv4"
    
//     input:
//     tuple val(sample_name), path(sam_file)
    
//     output:
//     tuple val(sample_name), path("${sam_file}.sorted.bam"), emit: bam 
    
//     script:
//     """
//     samtools view -bS ${sam_file} | samtools sort -o ${sam_file}.sorted.bam -T tmp  
//     """
    
// }

process SAMTOOLS_MERGE {
    label 'ALL'
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"

    input:
    tuple val(sample_name), path(out_bam)
    
    output:
    tuple val(sample_name), path("${sample_name}.bam"), emit: merged
    
    script:
    """
    samtools merge ${sample_name}.bam ${out_bam}
    """
}