#!/usr/bin/env nextflow

process minimap2 {
    publishDir "${params.out_dir}/minimap2/",mode:"copy"
    label "low"
    
    input:
    path (reference)
    tuple val(SampleName),path(SamplePath)
    
    output:
    tuple val(SampleName),path ("${SampleName}.sam")
    
    script:
    """
    minimap2 -ax map-ont ${reference} ${SamplePath} > ${SampleName}.sam
    """
}
