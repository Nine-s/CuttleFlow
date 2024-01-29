process kallisto_index{

	label 'kallisto'
	publishDir params.outdir
	
	input: 
	path(fasta_input)
	
	output:
	path("${fasta_input.baseName}.index"), emit: index
	
	script:
	"""
	kallisto index -i ${fasta_input.baseName}.index ${fasta_input}
	"""
}

process kallisto_map{

	label 'kallisto'
	publishDir params.outdir
	
	input:
	env STRANDNESS
	tuple val(pair_id), path(reads)
	path(index)
	path(gtf)
	
	output:
	path("*.bam"), emit: bam
	path("*.tsv"), emit: tsv
	
	
 	shell:
 	'''
 	if [[ $STRANDNESS == "firststrand" ]]; then
		kallisto quant -i !{index} -o ./ --gtf !{gtf} --fr-stranded --genomebam --threads !{params.threads} !{reads[0]} !{reads[1]}
   	elif [[ $STRANDNESS == "secondstrand" ]]; then
		kallisto quant -i !{index} -o ./ --gtf !{gtf} --rf-stranded --genomebam --threads !{params.threads} !{reads[0]} !{reads[1]}
    	elif [[ $STRANDNESS == "unstranded" ]]; then
		kallisto quant -i !{index} -o ./ --gtf !{gtf} --genomebam --threads !{params.threads} !{reads[0]} !{reads[1]}
	else  
		echo $STRANDNESS > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
    	fi
 	'''

}