process STAR_GENOMEGENERATE {
    label 'ALL'
    publishDir params.outdir
    container "mgibio/star:2.7.0f"
    
    input:
    path(reference)
    path(annotation)

    output:
    path("star/*"), emit: index

    script:
    """
    mkdir star
    STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles ${reference} \\
            --sjdbGTFfile ${annotation} \\
            --runThreadN ${params.threads} \\
	
    """
}