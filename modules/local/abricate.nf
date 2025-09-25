#!/usr/bin/env nextflow

process abricate {
    publishDir "${params.out_dir}/abricate/",mode:"copy"
    label "low"
    
    input:
    tuple val(SampleName),path(consensus)
    path(dbdir)
    
    output:
    path("${SampleName}_abricate.csv"),emit:abricate
    path("${SampleName}_withseq.csv"),emit:withseq
    
    script:
    """
    typing.sh ${SampleName} ${consensus} ${dbdir}
    """
}
