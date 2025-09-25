#!/usr/bin/env nextflow

process extract_mapped_ref {
    label "low"
    publishDir "${params.out_dir}/mapped_ref"
    
    input:
    tuple val(SampleName), path(amplicontxtfile)
    path (reference)
    
    output:
    val(SampleName)
    path("${SampleName}_mapped_ref.fasta"), emit:mapped_fasta
    
    script:
    """
    seqkit grep -f ${amplicontxtfile} ${reference} > ${SampleName}_mapped_ref.fasta
    """
}
