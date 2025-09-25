#!/usr/bin/env nextflow

process kraken2 {
    publishDir "${params.out_dir}/kraken2/",mode:"copy"
    label "high"
    
    input:
    tuple val(SampleName),path (SamplePath)
    path(db_path)
    
    output:
    path ("${SampleName}_kraken.csv")
    path ("${SampleName}_kraken_report.csv"),emit:(kraken2_raw)
    
    script:
    """
    kraken2 --db $db_path --output ${SampleName}_kraken.csv --report ${SampleName}_kraken_report.csv --threads 1 ${SamplePath}
    """
}
