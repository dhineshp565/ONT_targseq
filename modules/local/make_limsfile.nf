#!/usr/bin/env nextflow

process make_limsfile {
    publishDir "${params.out_dir}/LIMS/",mode:"copy"
    label "low"
    
    input:
    path (withseq)
    path (software_version_file)
    
    output:
    path("LIMSfile_*.tsv")
    
    script:
    """
    date=\$(date '+%Y-%m-%d_%H-%M-%S')
    
    awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_withseq.csv > LIMSfile.tsv

    cat ${software_version_file} LIMSfile.tsv > LIMSfile_\${date}.tsv
    """
}
