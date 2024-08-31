process BOWTIE2_INDEX {
    label 'bowtie2_index'
    publishDir path: params.outdir
    container "ummidock/bowtie2_samtools"
        

    input:
    path fasta

    output:
    path("${fasta.baseName}*"), emit: index

    script:
    """
    bowtie2-build -f ${fasta} --threads ${params.threads} ${fasta.baseName}
    """
}

