#!/usr/bin/env nextflow

process mapped_ref_bed {
    label "low"
    publishDir "${params.out_dir}/mapped_ref_bed"	
    
    input:
    path (reference)
    
    output:
    path("*.fai")
    path ("*.bed"),emit:bed
    
    script:
    """
    cp ${reference} reference.fasta
    samtools faidx reference.fasta

    awk 'BEGIN {FS=OFS="\t"}; {print \$1,0,\$2}' reference.fasta.fai > reference.bed
    """
}
