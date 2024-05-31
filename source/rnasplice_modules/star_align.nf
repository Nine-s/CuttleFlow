
process STAR_ALIGN {
    label 'ALL'
    publishDir params.outdir
    container "mgibio/star:2.7.0f"

    input:
    tuple val(sample_name), path(reads), val(condition)
    path(index)
    path(annotation)

    output:
    tuple val(sample_name), path("${sample_name}*.sam"), val(condition), emit: sam 
    path('*Log.final.out')   , emit: log_final

    script:
    """
    if [[ (${params.strand} == "firststrand") || (${params.strand} == "secondstrand") ]]; then
		STAR \\
            --genomeDir . \\
            --readFilesIn ${reads[0]} ${reads[1]}  \\
            --runThreadN ${params.threads} \\
            --outFileNamePrefix ${sample_name}. \\
            --sjdbGTFfile ${annotation} \\
			--alignSoftClipAtReferenceEnds No \\
			--outFilterIntronMotifs RemoveNoncanonical \\
			--outSAMattrIHstart 0 \\
            --readFilesCommand zcat 

	elif [[ ${params.strand} == "unstranded" ]]; then
		STAR \\
            --genomeDir . \\
            --readFilesIn ${reads[0]} ${reads[1]}  \\
			--alignSoftClipAtReferenceEnds No \\
			--outSAMstrandField intronMotif \\
			--outFilterIntronMotifs RemoveNoncanonical \\
        	--runThreadN ${params.threads} \\
        	--outFileNamePrefix ${sample_name}. \\
        	--sjdbGTFfile ${annotation} \\
			--outSAMattrIHstart 0 \\
            --readFilesCommand zcat 
	else  
		echo ${params.strand} > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi
   """
}
