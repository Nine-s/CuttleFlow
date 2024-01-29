process SAMTOOLS {
    label 'samtools'
    publishDir params.outdir
    container "mgibio/samtools:1.9"
    
    input:
    tuple val(sample_name), path(sam_file), val(condition)
    
    output:
    path("${sample_name}.sorted.bam"), emit: sample_bam 
    
    script:
    """
    samtools view -bS ${sam_file} | samtools sort -o ${sam_file}.sorted.bam -T tmp  
    """
    
}

process SAMTOOLS_MERGE {
    label 'samtools'
    publishDir params.outdir
    container "mgibio/samtools:1.9"

    input:
    file out_bam
    
    output:
    tuple val("alignement_gathered.bam"), path("alignement_gathered.bam"), emit: merged
    
    script:
    """
    samtools merge alignement_gathered.bam ${out_bam}
    """
}

//samtools view -bS ${sam_file} -@ ${params.threads} | samtools sort -o ${sample_name}.sorted.bam -T tmp  -@ ${params.threads} 
