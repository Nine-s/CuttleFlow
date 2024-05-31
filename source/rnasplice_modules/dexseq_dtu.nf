process DEXSEQ_DTU {
    label 'ALL'
    publishDir params.outdir

    //container "zavolab/r_dge_dtu:3.5.1"
    //container "biocontainers/htseq:v0.11.2-1-deb-py3_cv1"
    //container "biocontainers/htseq:v0.11.2-1-deb-py2_cv1"
    //container "machalen/dexseq:latest"
    container 'quay.io/biocontainers/bioconductor-dexseq:1.36.0--r40_0' //noit found
    //container "nfcore/rnaseq:1.4.2"
    //container "filipejesus/dexseq:3.8"

    input:
    path drimseq_sample_data
    path drimseq_d_counts
    path drimseq_contrast_data
    val ntop

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "DEXSeqResults.*.tsv"  , emit: dexseq_exon_results_tsv
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "perGeneQValue.*.tsv"  , emit: dexseq_gene_results_tsv

    script:
    """
    ${params.basedir}/bin/run_dexseq_dtu.R $drimseq_sample_data $drimseq_contrast_data $drimseq_d_counts $ntop
    """

}