process MERGE_RESULTS_DEXSEQ {
    publishDir params.outdir
    container "zavolab/salmon:1.1.0"

    input:
    path(out_files)
    
    output:
    path("dexseq_clean_counts"), emit: clean_counts
    
    script:
    """
    mkdir dexseq_clean_counts
    mv  ${out_files} dexseq_clean_counts
    """
}