#!/usr/bin/env nextflow

process seq_length {
    label 'low'
    publishDir "${params.out_dir}/seqlengths/", mode: "copy"
    
    input:
    tuple val(SampleName), path(cons), path(orf)

    output:
    path("*.tsv")

    script:
    """
    cat ${cons} ${SampleName}_ORF.fasta | seqkit fx2tab -n -l > ${SampleName}_seqlengths.tsv
    """
}
