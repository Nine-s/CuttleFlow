process CAT_FASTQ {
    label 'ALL'
    publishDir params.outdir

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    //input:
    //tuple val(sample), path(reads)

    // Define the input directory containing paired-end FastQ files
    input:
    path fastq_dir

    // Create channels for the two reads in the paired-end files
    reads1 = file(fastq_dir).listFiles { it.name =~ /_R1\.fastq/ }
    reads2 = file(fastq_dir).listFiles { it.name =~ /_R2\.fastq/ }

    // Define the output merged FastQ files
    output:
    path "merged_R1.fastq"
    path "merged_R2.fastq"

    // Merge the paired-end FastQ files using the 'cat' command
    script:
    """
    cat \${reads1.join(' ')} > merged_R1.fastq
    cat \${reads2.join(' ')} > merged_R2.fastq
    """
}