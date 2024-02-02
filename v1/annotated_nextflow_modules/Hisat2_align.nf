process HISAT2_ALIGN {
    label 'ALL'
    publishDir params.outdir
    container "nanozoo/hisat2:2.2.0--925b733"
    //container "biocontainers/hisat2:v2.1.0-2-deb_cv1"
 
    input:
    tuple val(sample_name), path(reads), val(condition)
    tuple path(reference), path(index)
    path(annotation)

    output:
    tuple val(sample_name), path("${reads[0].baseName}*.sam"), emit: sam
    path('*_summary.log')   , emit: log_final

    shell:
    '''
    if [[ (!{params.strand} == "firststrand") ]]; then
    
        hisat2 -x !{reference.baseName} -1 !{reads[0]} -2 !{reads[1]} --new-summary --summary-file !{reads[0].baseName}_summary.log --thread !{params.threads} --dta-cufflinks --rna-strandness FR -S !{reads[0].baseName}-!{condition}.sam

    elif [[ (!{params.strand} == "secondstrand") ]]; then
    
        hisat2 -x !{reference.baseName} -1 !{reads[0]} -2 !{reads[1]} --new-summary --summary-file !{reads[0].baseName}_summary.log --thread !{params.threads} --dta-cufflinks --rna-strandness RF -S !{reads[0].baseName}-!{condition}.sam

    elif [[ !{params.strand} == "unstranded" ]]; then
       
        hisat2 -x !{reference.baseName} -1 !{reads[0]} -2 !{reads[1]} --new-summary --summary-file !{reads[0].baseName}_summary.log --thread !{params.threads} --dta-cufflinks -S !{reads[0].baseName}-!{condition}.sam
    else  
		echo !{params.strand} > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi
    '''   
}
