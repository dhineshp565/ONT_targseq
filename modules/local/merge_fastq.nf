#!/usr/bin/env nextflow

process merge_fastq {
    publishDir "${params.out_dir}/merged"
    label "low"
    
    input:
    tuple val(SampleName),path(SamplePath)
    val (qscore)
    
    output:
    tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}")
    
    shell:
    """
    count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
    
    
    if [[ "\${count}" != "0" ]]
    then
        cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
        nanoq -i ${SampleName}.fastq.gz -s -H > ${SampleName}_readstats.csv
        # nanoq -i ${SampleName}.fastq.gz -q ${qscore} -o ${SampleName}_filtered.fastq.gz
    
    else
        count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
        if [[ "\${count}" != "0" ]]
        then
            cat ${SamplePath}/*.fastq > ${SampleName}.fastq
            nanoq -i ${SampleName}.fastq -s -H > ${SampleName}_readstats.csv
            # nanoq -i ${SampleName}.fastq -q ${qscore}  -o ${SampleName}_filtered.fastq
        fi
    fi
    """
}
