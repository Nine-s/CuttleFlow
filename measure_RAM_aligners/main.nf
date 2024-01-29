
nextflow.enable.dsl = 2

include { CHECK_STRANDNESS } from './modules/check_strandness.nf'
include { HISAT2_INDEX_REFERENCE ; HISAT2_INDEX_REFERENCE_MINIMAL ; HISAT2_ALIGN ; EXTRACT_SPLICE_SITES ; EXTRACT_EXONS } from './modules/hisat2.nf'
include { STAR_ALIGN ; STAR_INDEX_REFERENCE } from './modules/star.nf'
include { SALMON_ALIGN_QUANT ; SALMON_INDEX_REFERENCE ; GENERATE_DECOY_TRANSCIPTROME } from './modules/salmon.nf'
include {kallisto_index ; kallisto_map} from "./modules/kallisto"

log.info """\
         RNAseq analysis using NextFlow 
         =============================
         genome: ${params.reference_genome}
         annot : ${params.reference_annotation}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
params.outdir = 'results'

workflow {

    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
    CHECK_STRANDNESS( read_pairs_ch, params.reference_cdna, params.reference_annotation_ensembl )
    
    EXTRACT_EXONS( params.reference_annotation )
    EXTRACT_SPLICE_SITES( params.reference_annotation )
    HISAT2_INDEX_REFERENCE( params.reference_genome, EXTRACT_EXONS.out, EXTRACT_SPLICE_SITES.out )
    HISAT2_ALIGN( read_pairs_ch, HISAT2_INDEX_REFERENCE.out, CHECK_STRANDNESS.out.first() )
    
    STAR_INDEX_REFERENCE( params.reference_genome, params.reference_annotation )
    STAR_ALIGN( CHECK_STRANDNESS.out, read_pairs_ch, STAR_INDEX_REFERENCE.out, params.reference_annotation )
    
    GENERATE_DECOY_TRANSCIPTROME( params.reference_genome, params.reference_cdna )
    SALMON_INDEX_REFERENCE( GENERATE_DECOY_TRANSCIPTROME.out.decoy, GENERATE_DECOY_TRANSCIPTROME.out.gentrome )
    SALMON_ALIGN_QUANT( CHECK_STRANDNESS.out, read_pairs_ch, params.reference_annotation )

    kallisto_index(params.reference_cdna)
    kallisto_map(CHECK_STRANDNESS.out, read_pairs_ch, kallisto_index.out.index, params.gtf)
}