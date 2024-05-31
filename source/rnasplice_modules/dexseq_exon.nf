process DEXSEQ_EXON {
    label 'ALL'
    publishDir params.outdir

    container 'quay.io/biocontainers/bioconductor-dexseq:1.36.0--r40_0' 

    input:
    path (dexseq_clean_counts)    // path: dexseq_clean_counts
    path gff                          // path: dexseq_gff
    path samplesheet                  // path: samplesheet
    path contrastsheet                // path: contrastsheet
    val ntop                          // val: n_dexseq_plot

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "DEXSeqResults.*.csv"  , emit: dexseq_exon_results_csv
    path "perGeneQValue.*.csv"  , emit: dexseq_gene_results_csv
    path "plotDEXSeq.*.pdf"     , emit: dexseq_plot_results_pdf

    script:
    //${params.baseDir}
    """

    ${params.basedir}/bin/run_dexseq_exon.R ${dexseq_clean_counts} $gff $samplesheet $contrastsheet $ntop

    """
}
