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
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"

    input:
    tuple val(sample_name), path(out_bam)
    
    output:
    tuple val("alignement_gathered.bam"), path("alignement_gathered.bam"), emit: merged
    
    script:
    """
    samtools merge alignement_gathered.bam ${out_bam}
    """
}