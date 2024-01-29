process CHECK_STRANDNESS {
	input:
		tuple val(sample_name), path(reads)
    		path(reference_cdna)
    		path(annotation)

	output: 
		env STRANDNESS

	shell:
	'''
	STRANDNESS="firststrand"
     '''
}

