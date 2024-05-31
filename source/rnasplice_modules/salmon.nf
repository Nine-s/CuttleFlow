process SALMON_QUANT {
    label 'ALL'
    publishDir params.outdir
    container "zavolab/salmon:1.1.0"

    input:
    tuple val(sample_name), path(reads), val(condition)
    path(index)
    
    output:
    //tuple val(sample_name), path("transcripts_quant"), val(condition), emit: quantification 
    path(sample_name), emit: transcripts
    path("*info.json"), emit: json_info, optional: true

    script:
    """
    if [[ (${params.strand} == "firststrand") ]]; then 
	    salmon quant -i ${index} -l ISR -1 ${reads[0]} -2 ${reads[1]} --validateMappings -p ${params.threads} -o ${sample_name} 
    elif [[ (${params.strand} == "secondstrand") ]]; then 
        salmon quant -i ${index} -l ISF -1 ${reads[0]} -2 ${reads[1]} --validateMappings -p ${params.threads} -o ${sample_name}
	elif [[ ${params.strand} == "unstranded" ]]; then
		salmon quant -i ${index} -l IU -1 ${reads[0]} -2 ${reads[1]} --validateMappings -p ${params.threads} -o ${sample_name}
	else  
		echo ${params.strand} > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi
    cp ${sample_name}/aux_info/meta_info.json "${sample_name}_meta_info.json"
    """
    
}


