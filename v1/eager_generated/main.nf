nextflow.enable.dsl = 2


include { FASTQC as FASTQC_PRE_PREPROCESSING; FASTQC as FASTQC_POST_PREPROCESSING } from './modules_eager/fastqc.nf'
include { FASTP  } from './modules_eager/fastp_aDNA.nf'
include { ADAPTER_REMOVAL  } from './modules_eager/adapter_removal.nf'
include { BOWTIE2_INDEX  } from './modules_eager/bowtie2_index.nf'
include { BOWTIE2  } from './modules_eager/bowtie2.nf'
include { SAMTOOLS_FILTER  } from './modules_eager/samtools_filter.nf'
include { BBDUK  } from './modules_eager/bbduk.nf'
include { KRAKEN ; KRAKEN_PARSE ; KRAKEN_MERGE  } from './modules_eager/kraken.nf'

workflow{
        read_pairs_ch = Channel
            .fromPath( params.csv_input )
            .splitCsv(header: true, sep: ',')
            .map {row -> tuple(row.sample, [row.path_r1, row.path_r2], row.condition)}
            .view()
        
FASTQC_PRE_PREPROCESSING(read_pairs_ch)
FASTP(read_pairs_ch)
ADAPTER_REMOVAL(FASTP.out.reads, params.clip_adapters_list)
BOWTIE2_INDEX(params.genome)
FASTQC_POST_PREPROCESSING(ADAPTER_REMOVAL.out.reads)
BOWTIE2(ADAPTER_REMOVAL.out.reads, params.genome, BOWTIE2_INDEX.out.index)
SAMTOOLS_FILTER(BOWTIE2.out.bam)
BBDUK(SAMTOOLS_FILTER.out.unmapped_fastq)
KRAKEN(BBDUK.out.fastq)
KRAKEN_PARSE(KRAKEN.out.kraken_report)
KRAKEN_MERGE(KRAKEN_PARSE.out.kraken_parsed)

}