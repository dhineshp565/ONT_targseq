#!/usr/bin/env nextflow

process krona_kraken {
    publishDir "${params.out_dir}/krona_kraken/",mode:"copy"
    label "low"
    
    input:
    path(raw)
    
    output:
    path ("rawreads_classified.html"),emit:raw
    
    script:
    """
    ktImportTaxonomy -t 5 -m 3 -o rawreads_classified.html ${raw}
    """
}
