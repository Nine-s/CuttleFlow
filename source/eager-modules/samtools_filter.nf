process SAMTOOLS_FILTER {
    label 'samtools_filter'
    publishDir params.outdir
    container "nfcore/eager:2.5.1"

    input: 
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), file("*filtered.bam"), file("*.{bai,csi}"), emit: filtered_bam
    tuple val(sample_id), file("*.unmapped.fastq.gz"), emit: unmapped_fastq
    tuple val(sample_id), file("*.unmapped.bam"), emit: unmapped_bam

    script:    
    def size = params.large_ref ? '-c' : ''
    
    // Unmapped/MAPQ Filtering WITHOUT min-length filtering
    if ( "${params.bam_unmapped_type}" == "keep"  && params.bam_filter_minreadlength == 0 ) {
        """
        samtools view -h ${bam} -@ ${params.threads} -q ${params.bam_mapping_quality_threshold} -b > ${sample_id}.filtered.bam
        samtools index ${sample_id}.filtered.bam ${size}
        """
    } else if ( "${params.bam_unmapped_type}" == "discard" && params.bam_filter_minreadlength == 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > ${sample_id}.filtered.bam
        samtools index ${sample_id}.filtered.bam ${size}
        """
    } else if ( "${params.bam_unmapped_type}" == "bam" && params.bam_filter_minreadlength == 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -f4 -b > ${sample_id}.unmapped.bam
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > ${sample_id}.filtered.bam
        samtools index ${sample_id}.filtered.bam ${size}
        """
    } else if ( "${params.bam_unmapped_type}" == "fastq" && params.bam_filter_minreadlength == 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -f4 -b > ${sample_id}.unmapped.bam
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > ${sample_id}.filtered.bam
        samtools index ${sample_id}.filtered.bam ${size}

        ## FASTQ
        samtools fastq -tN ${sample_id}.unmapped.bam | pigz -p ${params.threads} > ${sample_id}.unmapped.fastq.gz
        rm ${sample_id}.unmapped.bam
        """
    } else if ( "${params.bam_unmapped_type}" == "both" && params.bam_filter_minreadlength == 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -f4 -b > ${sample_id}.unmapped.bam
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > ${sample_id}.filtered.bam
        samtools index ${sample_id}.filtered.bam ${size}
        
        ## FASTQ
        samtools fastq -tN ${sample_id}.unmapped.bam | pigz -p ${params.threads} > ${sample_id}.unmapped.fastq.gz
        """
    // Unmapped/MAPQ Filtering WITH min-length filtering
    } else if ( "${params.bam_unmapped_type}" == "keep" && params.bam_filter_minreadlength != 0 ) {
        """
        samtools view -h ${bam} -@ ${params.threads} -q ${params.bam_mapping_quality_threshold} -b > tmp_mapped.bam
        filter_bam_fragment_length.py -a -l ${params.bam_filter_minreadlength} -o ${sample_id} tmp_mapped.bam
        samtools index ${sample_id}.filtered.bam ${size}
        """
    } else if ( "${params.bam_unmapped_type}" == "discard" && params.bam_filter_minreadlength != 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > tmp_mapped.bam
        filter_bam_fragment_length.py -a -l ${params.bam_filter_minreadlength} -o ${sample_id} tmp_mapped.bam
        samtools index ${sample_id}.filtered.bam ${size}
        """
    } else if ( "${params.bam_unmapped_type}" == "bam" && params.bam_filter_minreadlength != 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -f4 -b > ${sample_id}.unmapped.bam
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > tmp_mapped.bam
        filter_bam_fragment_length.py -a -l ${params.bam_filter_minreadlength} -o ${sample_id} tmp_mapped.bam
        samtools index ${sample_id}.filtered.bam ${size}
        """
    } else if ( "${params.bam_unmapped_type}" == "fastq" && params.bam_filter_minreadlength != 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -f4 -b > ${sample_id}.unmapped.bam
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > tmp_mapped.bam
        filter_bam_fragment_length.py -a -l ${params.bam_filter_minreadlength} -o ${sample_id} tmp_mapped.bam
        samtools index ${librsample_idaryid}.filtered.bam ${size}

        ## FASTQ
        samtools fastq -tN ${sample_id}.unmapped.bam | pigz -p ${params.threads} > ${sample_id}.unmapped.fastq.gz
        rm ${sample_id}.unmapped.bam
        """
    } else if ( "${params.bam_unmapped_type}" == "both" && params.bam_filter_minreadlength != 0 ){
        """
        samtools view -h ${bam} -@ ${params.threads} -f4 -b > ${sample_id}.unmapped.bam
        samtools view -h ${bam} -@ ${params.threads} -F4 -q ${params.bam_mapping_quality_threshold} -b > tmp_mapped.bam
        filter_bam_fragment_length.py -a -l ${params.bam_filter_minreadlength} -o ${sample_id} tmp_mapped.bam
        samtools index ${sample_id}.filtered.bam ${size}
        
        ## FASTQ
        samtools fastq -tN ${sample_id}.unmapped.bam | pigz -p ${params.threads} > ${sample_id}.unmapped.fastq.gz
        """
    }
}
