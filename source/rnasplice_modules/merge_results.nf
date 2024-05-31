
process MERGE_RESULTS_SALMON {
    publishDir params.outdir
    container "zavolab/salmon:1.1.0"

    input:
    path out_folders
    
    output:
    path("salmon"), emit:gathered_bam
    
    script:
    """
    mkdir salmon
    mv  ${out_folders} salmon
    """
}