
process HISAT2_INDEX {
    label 'ALL'
    publishDir params.outdir
    container "nanozoo/hisat2:2.2.0--925b733"
    //container "biocontainers/hisat2:v2.1.0-2-deb_cv1"

    input:
    path(reference)
    path(annotation)

    output:
    tuple path(reference), path("${reference.baseName}*.ht2"), emit: index

    script:
    """
    hisat2-build ${reference} ${reference.baseName} -p ${params.threads} 
    """
    // #python /workspace/projects/Nine-s/description_prototype/bin/hisat2_extract_splice_sites.py ${annotation} > out_splice_sites.txt
    // #python /workspace/projects/Nine-s/description_prototype/bin/hisat2_extract_exons.py ${annotation}  > out_exons.txt
}