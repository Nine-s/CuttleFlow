process BOWTIE2 {
    label 'bowtie2'
    publishDir params.outdir
    container "ummidock/bowtie2_samtools"
    
    input:
    tuple val(sample_id), path(reads) 
    path(fasta)
    path(index) 

    output:
    tuple val(sample_id), path("*.sam"), emit: sam 

    script:
    def trim5 = params.bt2_trim5 != 0 ? "--trim5 ${params.bt2_trim5}" : ""
    def trim3 = params.bt2_trim3 != 0 ? "--trim3 ${params.bt2_trim3}" : ""
    def bt2n = params.bt2n != 0 ? "-N ${params.bt2n}" : ""
    def bt2l = params.bt2l != 0 ? "-L ${params.bt2l}" : ""

    if ( "${params.bt2_alignmode}" == "end-to-end"  ) {
      switch ( "${params.bt2_sensitivity}" ) {
        case "no-preset":
        sensitivity = ""; break
        case "very-fast":
        sensitivity = "--very-fast"; break
        case "fast":
        sensitivity = "--fast"; break
        case "sensitive":
        sensitivity = "--sensitive"; break
        case "very-sensitive":
        sensitivity = "--very-sensitive"; break
        default:
        sensitivity = ""; break
        }
      } else if ("${params.bt2_alignmode}" == "local") {
      switch ( "${params.bt2_sensitivity}" ) {
        case "no-preset":
        sensitivity = ""; break
        case "very-fast":
        sensitivity = "--very-fast-local"; break
        case "fast":
        sensitivity = "--fast-local"; break
        case "sensitive":
        sensitivity = "--sensitive-local"; break
        case "very-sensitive":
        sensitivity = "--very-sensitive-local"; break
        default:
        sensitivity = ""; break

        }
      } else {
        sensitivity = ""
      }

    """
    bowtie2 -x ${fasta.baseName} -1 ${reads[0]} -2 ${reads[1]} -p ${params.threads} ${sensitivity} ${bt2n} ${bt2l} ${trim5} ${trim3} --maxins ${params.bt2_maxins} -S ${reads[0].baseName}.sam
    """  
}
