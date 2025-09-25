#!/usr/bin/env nextflow

process igvreports {
    label "high"
    publishDir "${params.out_dir}/igvreports/", mode: "copy"
    
    input:
    path(csv)
    path(reference)
    path(bed)
    path(bedgraph)
    
    output:
    path("*.html")
    
    script:
    """
    create_report ${bed} ${reference} --tracks *.bedgraph --output igv_coverage.html
    """
}
