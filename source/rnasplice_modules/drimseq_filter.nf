process DRIMSEQ_FILTER {
    label 'ALL'
    publishDir params.outdir
    container "zavolab/r_dge_dtu:3.5.1"

    input:
    path txi                    // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    path tximport_tx2gene       // path: tximport.tx2gene.tsv
    path samplesheet            // path: /path/to/samplesheet.csv
    val min_samps_gene_expr     // val params.min_samps_gene_expr
    val min_samps_feature_expr  // val params.min_samps_feature_expr
    val min_samps_feature_prop  // val params.min_samps_feature_prop
    val min_feature_expr        // val params.min_feature_expr
    val min_feature_prop        // val params.min_feature_prop
    val min_gene_expr           // val params.min_gene_expr

    output:
    path "dmDSdata.rds"  , emit: drimseq_dataset_rds
    path "samples.tsv"   , emit: drimseq_samples_tsv
    path "counts.tsv"    , emit: drimseq_counts_tsv

    script:

    """
    ${params.basedir}/bin/run_drimseq_filter.R $txi $tximport_tx2gene $samplesheet \\
        $min_samps_gene_expr \\
        $min_samps_feature_expr \\
        $min_samps_feature_prop \\
        $min_feature_expr \\
        $min_feature_prop \\
        $min_gene_expr

    """
}