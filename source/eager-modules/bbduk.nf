process BBDUK {
  label 'bbduk'
  publishDir params.outdir
  container "kincekara/bbduk:38.98"

  input:
  tuple val(sample_id), path(fastq)
  
  output:
  tuple val(sample_id), path("*_lowcomplexityremoved.fq.gz"), emit: fastq
  
  path("*_bbduk.stats") 

  script:
  """
  bbduk.sh in=${fastq} threads=${params.threads} entropymask=f entropy=${params.metagenomic_complexity_entropy} out=${fastq}_lowcomplexityremoved.fq.gz 2> ${fastq}_bbduk.stats
  """

}