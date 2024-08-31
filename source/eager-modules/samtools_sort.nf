process SAMTOOLS_SORT {
    label 'bowtie2'
    publishDir params.outdir
    container "nfcore/eager:2.5.1"
    
    input:
    tuple val(sample_id), path(sam) 

    output:
    tuple val(sample_id), path("*.bam"), path("*.{bai,csi}"), emit: bam 

    script:
    //def size = params.large_ref ? '-c' : ''
    
    """
    samtools sort -@ ${params.threads} -O bam ${sam} > ${sam.baseName}.mapped.bam
    samtools index ${sam.baseName}.mapped.bam 
    """  
    //${size}
}
