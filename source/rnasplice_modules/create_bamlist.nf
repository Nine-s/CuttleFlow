process CREATE_BAMLIST {
    label 'process_single'
    publishDir params.outdir

    container 'quay.io/biocontainers/sed:4.7.0'

    input:
    tuple val(contrast), val(cond1), val(cond2), path(bam1), path(bam2)

    output:
    tuple val(contrast), path("${cond1}_bamlist.txt"), path("${cond2}_bamlist.txt"), emit: bamlist

    script:
    """
    echo $bam1 | sed 's: :,:g' > ${cond1}_bamlist.txt
    echo $bam2 | sed 's: :,:g' > ${cond2}_bamlist.txt

    """
}