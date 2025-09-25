#!/usr/bin/env nextflow

process bedtools {
    label "high"
    publishDir "${params.out_dir}/bedtools/"
    
    input:
    tuple val(SampleName), path (bam)
    
    output:
    path ("${SampleName}*.bedgraph")
    
    script:
    """
    bedtools genomecov -ibam ${SampleName}.bam -bga > ${SampleName}.bedgraph
    """
}
