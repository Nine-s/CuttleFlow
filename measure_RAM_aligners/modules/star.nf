process STAR_INDEX_REFERENCE {
    label 'star'
    publishDir params.outdir
    
    input:
    path(reference)
    path(annotation)

    output:
    path("star/*")

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

process STAR_ALIGN {
    label 'star'
    publishDir params.outdir

    input:
    env STRANDNESS
    tuple val(sample_name), path(reads)
    path(index)
    path(annotation)

    output:
    tuple val(sample_name), path("${sample_name}*.sam"), emit: sample_sam 

    shell:
    '''
    if [[ ($STRANDNESS == "firststrand") || ($STRANDNESS == "secondstrand") ]]; then
		STAR \\
        		--genomeDir . \\
        		--readFilesIn !{reads[0]} !{reads[1]}  \\
        		--runThreadN !{params.threads} \\
        		--outFileNamePrefix !{sample_name}. \\
        		--sjdbGTFfile !{annotation} \\
			--alignSoftClipAtReferenceEnds No \\
			--outFilterIntronMotifs RemoveNoncanonical \\
			--outSAMattrIHstart 0

	elif [[ $STRANDNESS == "unstranded" ]]; then
		STAR \\
            --genomeDir . \\
            --readFilesIn !{reads[0]} !{reads[1]}  \\
			--alignSoftClipAtReferenceEnds No \\
			--outSAMstrandField intronMotif \\
			--outFilterIntronMotifs RemoveNoncanonical \\
        	--runThreadN !{params.threads} \\
        	--outFileNamePrefix !{sample_name}. \\
        	--sjdbGTFfile !{annotation} \\
			--outSAMattrIHstart 0

	else  
		echo $STRANDNESS > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi

   '''
}