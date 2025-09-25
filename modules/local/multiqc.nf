#!/usr/bin/env nextflow

process multiqc {
    publishDir "${params.out_dir}/multiqc/",mode:"copy"
    label "low"
    
    input:
    path '*'
    
    output:
    file ("multiqc_report.html")
    file ("multiqc_data")
    
    script:
    """
    multiqc .
    """
}
