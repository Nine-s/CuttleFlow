
process SAMTOOLS {
    label 'ALL'
    publishDir params.outdir

    container "biocontainers/samtools:v1.7.0_cv4"
    
    input:
    tuple val(sample_name), path(sam_file)
    
    output:
    tuple val(sample_name), path("${sam_file.baseName}.sorted.bam"), emit: bam 
    path("alignment_summary.txt"), emit: log
    
    script:
    """
    samtools flagstat ${sam_file} > alignment_summary.txt
    samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sam_file.baseName}.sorted.bam -T tmp  -@ ${params.threads} 
    rm ${sam_file}
    """
    
}

process SAMTOOLS_MERGE {
    label 'samtools'
    publishDir params.outdir
    container "biocontainers/samtools:v1.7.0_cv4"

    input:
    tuple val(name), path(bam_files)
    
    output:
    tuple val(name), path("*.bam"), emit: merged
    
    script:
    """
    samtools merge ${bam_files[0].baseName}_merged.bam ${bam_files} -@ ${params.threads}
    """
}
