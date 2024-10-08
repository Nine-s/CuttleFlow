workflow FASTQSPLIT{

    

    take:
    reads_channel 
    
    main:
    SPLIT_READS(reads_channel) \
	  	 | map { name, fastq, fastq1 -> tuple( groupKey( name, fastq.size()), fastq, fastq1 ) } \
       	 	 | transpose() \
       	 	 | map { name, fastq, fastq1 -> tuple(name, [fastq, fastq1] ) } \
       	 	 | view() \
       		 | set{ split_out }
    
    
    emit:
    split_reads = split_out
}
    
process SPLIT_READS{
    container "nfcore/rnaseq:1.4.2"
    label 'ALL'
    publishDir params.outdir

    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("*-*${name}*1*.f*q"), path("*-*${name}*2*.f*q"), emit: split_reads

    shell:
    '''
    zcat "!{fastq[0]}" > "!{name}_1.fastq"
    length=$(wc -l < "!{name}_1.fastq")
    length=$((length / 4))
    s=!{params.split}
    z=$((length / s))
    splitby=$((z + 1))
    /work/simon/bin/splitFastq -i "!{name}_1.fastq" -n "$splitby" -o "!{name}_1"
    rm !{name}_1.fastq
    
    zcat "!{fastq[1]}" > "!{name}_2.fastq"
    /work/simon/bin/splitFastq -i "!{name}_2.fastq" -n "$splitby" -o "!{name}_2"
    rm !{name}_2.fastq
    '''
}
// !{params.basedir}/../bin/splitFastq -i !{fastq[0]} -n ${splitby} -o !{name}_1
//     !{params.basedir}/../bin/splitFastq -i !{fastq[1]} -n ${splitby} -o !{name}_2
