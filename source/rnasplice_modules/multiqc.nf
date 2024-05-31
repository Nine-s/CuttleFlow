process MULTIQC {
    label 'ALL'
    publishDir params.outdir
    //container "staphb/multiqc:1.8"
    container "nanozoo/multiqc:1.12--4f89fda"

    input:
    path(salmon)
    path(trim_g)                                             
    path(star)   
    path(fastqc)

    output:
    path "multiqc_report.html"

    script:
    """
    mkdir logs
    mv ${salmon} ${trim_g} ${star} ${fastqc} logs
    multiqc logs -o .
    """
}
    // path "*multiqc_report.html", emit: report
    // path "*_data"              , emit: data
    // path "*_plots"             , optional:true, emit: plots

    // mkdir -p multiqc_report
    // multiqc -o multiqc_report .


    //
        // echo ${salmon}
        // echo ${trim_galore}
        // echo ${star}
        // echo ${fastqc}
// process MULTIQC {
//     label 'ALL'
//     publishDir params.outdir
//     //container "staphb/multiqc:1.8"
//     container "nanozoo/multiqc:1.12--4f89fda"

//     input:
//     path(salmon_transcripts)
//     path(salmon_json_info)
//     tuple val(name_trim), path("*val*.f*q*.gz"), val(condition)
//     path(trim_log)                                             
//     tuple val(sample_name), path("${sample_name}*.sam"), val(condition)
//     path(star_log_final)   
//     path(star_log_out)   
//     path(fastqc)

//     output:
//     path "multiqc_report.html"

//     script:
//     """
//     mkdir logs
//     mv ${salmon_transcripts} ${salmon_json_info} ${trim_log} ${star_log_final} ${star_log_out} ${fastqc} logs
//     multiqc logs -o .
//     """
// }
//     // path "*multiqc_report.html", emit: report
//     // path "*_data"              , emit: data
//     // path "*_plots"             , optional:true, emit: plots

//     // mkdir -p multiqc_report
//     // multiqc -o multiqc_report .


//     //
//         // echo ${salmon}
//         // echo ${trim_galore}
//         // echo ${star}
//         // echo ${fastqc}