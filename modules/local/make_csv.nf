#!/usr/bin/env nextflow

process make_csv {
    publishDir "${params.out_dir}"
    
    input:
    path(fastq_input)
    
    output:
    path("samplelist.csv")
    
    script:
    """
    makecsv.sh ${fastq_input}
    """
}
