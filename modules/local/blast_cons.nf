#!/usr/bin/env nextflow

process blast_cons {
    publishDir "${params.out_dir}/blast/", mode: "copy"
    containerOptions "-v ${params.blastdb_path}:${params.blastdb_path}"
    label "high"

    input:
    tuple val(SampleName), path(consensus)
    path(blastdb_path)
    val(blastdb_name)

    output:
    path("${SampleName}_report_blast.csv"), emit: blast_formatted
    path("${SampleName}_blast.csv")

    script:
    """
    # Create header for the report file
    echo -e "queryid\\tsubject_id\\talignment_length\\tquery_coverage\\t%identity\\tevalue\\tstaxids\\tsscinames\\tscomnames\\tstitle" > "${SampleName}_report_blast.csv"
    
    # Check if consensus file is empty or doesn't exist
    if [ ! -s "${consensus}" ]; then
        echo -e "none\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone" >> "${SampleName}_report_blast.csv"
        touch "${SampleName}_blast.csv"
        exit 0
    fi
    
    # Run BLAST
    blastn -db ${blastdb_path}/${blastdb_name} -query ${consensus} -out ${SampleName}_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue staxids ssciname scomnames stitle" -max_target_seqs 5
    
    # Process BLAST results
    if grep -q "# 0 hits found" "${SampleName}_blast.csv"; then
        # No hits found - add "none" row
        echo -e "none\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone\\tnone" >> "${SampleName}_report_blast.csv"
    else
        # Hits found - extract non-comment lines and add to report
        grep -v "^#" "${SampleName}_blast.csv" | sort | uniq >> "${SampleName}_report_blast.csv"
    fi
    """
}

