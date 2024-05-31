process CUSTOM_GETCHROMSIZES {
    label 'ALL'
    publishDir params.outdir
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
    path(reference_genome)

    output:
    path ("*.sizes"), emit: sizes
    path ("*.fai")  , emit: fai
    path ("*.gzi")  , emit: gzi, optional: true

    script:
    """
    samtools faidx ${reference_genome}
    cut -f 1,2 ${reference_genome}.fai > ${reference_genome}.sizes
    """

}

    // stub:
    // """
    // touch ${reference_genome}.fai
    // touch ${fasta}.sizes
    // """
