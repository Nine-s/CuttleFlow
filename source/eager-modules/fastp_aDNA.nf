process FASTP {
    label 'fastp'
    publishDir params.outdir
    //container "nanozoo/fastp:0.20.0--093b790"
    container "nfcore/eager:2.0.5"
        
    input:
    tuple val(sample_id), path(reads), val(condition)

    output:
    tuple val(sample_id), path("*.fq.gz"), emit: reads
    path("*.json"), emit: json

    script:
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${reads[0].baseName}.pG.fq.gz" --out2 "${reads[1].baseName}.pG.fq.gz" -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${params.threads} --json "${sample_id}"_polyg_fastp.json 
    """
}
//    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${reads[0].baseName}.pG.fq.gz" --out2 "${reads[1].baseName}.pG.fq.gz" -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${params.threads} --json "${sample_id}"_polyg_fastp.json 

//-A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L 
