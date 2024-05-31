process DEXSEQ_COUNT {
    label 'ALL'
    publishDir params.outdir

    container "quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0"

    input:
    tuple val(name), path(bam), val(condition)
    path (gff)
    val alignment_quality                   // val params.alignment_quality

    output:
    path("*.clean.count.txt"), emit: dexseq_clean_txt

    script:

    def read_type = '-p yes'

    def alignment_quality = "-a ${alignment_quality}"

    def strand = ''
    if (params.strand == 'forward') {
        strand = '-s yes'
    } else if (params.strand == 'reverse') {
        strand = '-s reverse'
    } else if (params.strand == 'unstranded') {
        strand = '-s no'
    }

    script:
    """
    ${params.basedir}/bin/dexseq_count.py $gff $read_type -f bam $bam -r pos ${name}.clean.count.txt $alignment_quality $strand
    """
}