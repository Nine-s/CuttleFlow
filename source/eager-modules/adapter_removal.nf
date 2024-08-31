
process ADAPTER_REMOVAL {
    label 'adapter_removal'
    publishDir params.outdir
    container "jdidion/adapterremoval:2.2.2"

    input:
    tuple val(sample_id), path(reads)
    //path(adapterlist) 

    output:
    tuple val(sample_id), path("output/*pair*.truncated.gz"), emit: reads
    path("output/*.settings"), emit: settings

    script:
    def base = "${reads[0].baseName}"
    def adapters_to_remove = !params.clip_adapters_list ? "--adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor}" : "--adapter-list ${adapterlist}"
    //This checks whether we skip trimming and defines a variable respectively
    def preserve5p = params.preserve5p ? '--preserve5p' : '' // applies to any AR command - doesn't affect output file combination
    
    if ( !params.skip_collapse && !params.skip_trim  && !params.mergedonly && !params.preserve5p ) {
    """
    mkdir -p output

    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}

    cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
    
    mv *.settings output/
    
    """
    //PE mode, collapse and trim, outputting all reads, preserving 5p
    } else if (!params.skip_collapse && !params.skip_trim  && !params.mergedonly && params.preserve5p) {
    """
    mkdir -p output

    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}

    cat *.collapsed.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz

    mv *.settings output/

    """
    // PE mode, collapse and trim but only output collapsed reads
    } else if (  !params.skip_collapse && !params.skip_trim && params.mergedonly && !params.preserve5p ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe  --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}
    
    cat *.collapsed.gz *.collapsed.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
  
    mv *.settings output/
    """
    // PE mode, collapse and trim but only output collapsed reads, preserving 5p
    } else if ( !params.skip_collapse && !params.skip_trim && params.mergedonly && params.preserve5p ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe  --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}
    
    cat *.collapsed.gz > output/${base}.pe.combined.tmp.fq.gz
    

    mv *.settings output/
    """
    // PE mode, collapsing but skip trim, (output all reads). Note: seems to still generate `truncated` files for some reason, so merging for safety.
    // Will still do default AR length filtering I guess
    } else if ( !params.skip_collapse && params.skip_trim && !params.mergedonly ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --adapter1 "" --adapter2 ""
    
    cat *.collapsed.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
  
    mv *.settings output/
    """
    // PE mode, collapsing but skip trim, and only output collapsed reads. Note: seems to still generate `truncated` files for some reason, so merging for safety.
    // Will still do default AR length filtering I guess
    } else if ( !params.skip_collapse && params.skip_trim && params.mergedonly ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} --collapse ${preserve5p}  --adapter1 "" --adapter2 ""
    
    cat *.collapsed.gz > output/${base}.pe.combined.tmp.fq.gz

    mv *.settings output/
    """
    // PE mode, skip collapsing but trim (output all reads, as merging not possible) - activates paired-end mapping!
    } else if ( params.skip_collapse && !params.skip_trim ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base}.pe --gzip --threads ${params.threads} --qualitymax ${params.qualitymax} ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}
    
    mv ${base}.pe.pair*.truncated.gz *.settings output/
    """
    } 
       
}


//clip_adapters_list = '/workspace/ninon/description_prototype/modules_eager/adapter_list.txt'
	