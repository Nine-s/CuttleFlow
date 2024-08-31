process MALT {
  label 'malt'
  publishDir params.outdir

  input:
  path(fastqs) 
  path(db) 

  output:
  path("*.rma6"), emit: rma_malt_extract
  path("*.sam.gz"), optional: true
  path("malt.log"), emit: malt_log
  
  script:
  if ( "${params.malt_min_support_mode}" == "percent" ) {
    min_supp = "-supp ${params.malt_min_support_percent}" 
  } else if ( "${params.malt_min_support_mode}" == "reads" ) {
    min_supp = "-sup ${params.metagenomic_min_support_reads}"
  }
  def sam_out = params.malt_sam_output ? "-a . -f SAM" : ""
  """
  malt-run \
  -t ${params.threads} \
  -v \
  -o . \
  -d ${db} \
  ${sam_out} \
  -id ${params.percent_identity} \
  -m ${params.malt_mode} \
  -at ${params.malt_alignment_mode} \
  -top ${params.malt_top_percent} \
  ${min_supp} \
  -mq ${params.malt_max_queries} \
  --memoryMode ${params.malt_memory_mode} \
  -i ${fastqs.join(' ')} |&tee malt.log
  """
}


process MALTEXTRACT {
  label 'maltextract'
  publishDir params.outdir

  input:
  path(rma6)
  path(taxon_list)
  path(ncbifiles)
  
  output:
  path("results/*"), emit: results

  script:
  def destack = params.maltextract_destackingoff ? "--destackingOff" : ""
  def downsam = params.maltextract_downsamplingoff ? "--downSampOff" : ""
  def dupremo = params.maltextract_duplicateremovaloff ? "--dupRemOff" : ""
  def matches = params.maltextract_matches ? "--matches" : ""
  def megsum = params.maltextract_megansummary ? "--meganSummary" : ""
  def topaln = params.maltextract_topalignment ?  "--useTopAlignment" : ""
  def ss = params.single_stranded ? "--singleStranded" : ""
  """
  MaltExtract \
  -t ${taxon_list} \
  -i ${rma6.join(' ')} \
  -o results/ \
  -r ${ncbifiles} \
  -p ${params.threads} \
  -f ${params.maltextract_filter} \
  -a ${params.maltextract_toppercent} \
  --minPI ${params.maltextract_percentidentity} \
  ${destack} \
  ${downsam} \
  ${dupremo} \
  ${matches} \
  ${megsum} \
  ${topaln} \
  ${ss}

  postprocessing.AMPS.r -r results/ -m ${params.maltextract_filter} -t ${params.threads} -n ${taxon_list} -j
  """
}
