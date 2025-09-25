#!/usr/bin/env nextflow

process porechop {
    label "high"
    publishDir "${params.out_dir}/trimmed"
    
    input:
    tuple val(SampleName),path(SamplePath)
    
    output:
    tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
    
    script:
    """
    porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
    """
}
