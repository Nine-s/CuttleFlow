process SALMON_GENOMEGENERATE {
    label 'ALL'
    publishDir params.outdir    
    container "zavolab/salmon:1.1.0"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon"      , emit: index

    script:
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    def gentrome      = "gentrome.fa"
    // if (genome_fasta.endsWith('.gz')) {
    //     get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    //     gentrome      = "gentrome.fa.gz"
    // }
    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    salmon index -t ${gentrome} -d decoys.txt -p ${params.threads} -i salmon
    """
//--decoys decoys.txt
}