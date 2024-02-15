#!/usr/bin/env nextflow
nextflow.enable.dsl=2



// make csv file with headers from the given input

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

//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.out_dir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}")
	
	shell:
	"""
	count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]]
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
			nanoq -i ${SampleName}.fastq.gz -s -H > ${SampleName}_readstats.csv
			nanoq -i ${SampleName}.fastq.gz -q 20 -o ${SampleName}_filtered.fastq.gz
		
		else
			count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
			if [[ "\${count}" != "0" ]]
			then
				cat ${SamplePath}/*.fastq > ${SampleName}.fastq
				nanoq -i ${SampleName}.fastq.gz -s -H > ${SampleName}_readstats.csv
				nanoq -i ${SampleName}.fastq -q 20 -o ${SampleName}_filtered.fastq
			fi
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "high"
	publishDir "${params.out_dir}/trimmed"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
	"""
}
// sequence alignment using minimap2
process minimap2 {
        publishDir "${params.out_dir}/minimap2/",mode:"copy"
		label "low"
        input:
        path (reference)
        tuple val(SampleName),path(SamplePath)
        output:
        val(SampleName)
		path ("${SampleName}.sam")
        script:
        """
        minimap2 -ax map-ont ${reference} ${SamplePath} > ${SampleName}.sam
        """
}

//convert minimap2 output sam to sorted bam and split bam files and create consensus
process splitbam {
	publishDir "${params.out_dir}/splitbam",mode:"copy"
	label "medium"
	input:
	val(SampleName)
	path(SamplePath)
	path (primerbed)
	output:
	val(SampleName),emit:SampleName
	path("${SampleName}_stats.txt"),emit:stats
	path("${SampleName}_mappedreads.txt"),emit:mapped
	path("${SampleName}_idxstats.txt"),emit:idxstats
	tuple val(SampleName),path("${SampleName}_consensus.fasta"),emit:consensus
	path("${SampleName}_consensus.fasta"),emit:cons_only
	path("${SampleName}_unfilt_stats.txt"),emit:unfilt_stats
	script:
	"""
	splitbam.sh ${SampleName} ${SamplePath} ${primerbed}

	"""
}


//multiqc generate mapped read statistics from samtools output

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

//kraken2 for classification
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
process kraken2_consensus {
	publishDir "${params.out_dir}/kraken2_cons/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
	
	output:
	path ("${SampleName}_cons_kraken.csv"),emit:(kraken2_cons)
	path ("${SampleName}_cons_kraken_report.csv")

	script:
	"""
	kraken2 --db $db_path --output ${SampleName}_cons_kraken.csv --report ${SampleName}_cons_kraken_report.csv --threads 3 ${SamplePath} --use-names --use-mpa-style
	"""
}

//krona plots
process krona_kraken {
	publishDir "${params.out_dir}/krona_kraken/",mode:"copy"
	label "low"
	input:
	path(raw)
	path(consensus)
	
	output:
	path ("rawreads_classified.html"),emit:raw
	path("consensus_classified.html"),emit:cons
	script:
	"""
	ktImportTaxonomy -t 5 -m 3 -o rawreads_classified.html ${raw}
	ktImportTaxonomy -t 5 -m 3 -o consensus_classified.html ${consensus}
	"""
}
process krona_centrifuge {
	publishDir "${params.out_dir}/krona_centrifuge/",mode:"copy"
	label "low"
	input:
	path(raw)
	path(consensus)
	
	output:
	path ("rawreads_classified.html"),emit:raw
	path("consensus_classified.html"),emit:cons
	script:
	"""
	ktImportTaxonomy -t 5 -m 3 -o rawreads_classified.html ${raw}
	ktImportTaxonomy -t 5 -m 3 -o consensus_classified.html ${consensus}
	"""
}
//make html report with rmarkdown
process make_report {
	publishDir "${params.out_dir}/",mode:"copy"
	
	label "low"
	input:
	path(csv)
	path(krona_reports_raw)
	path(mappedreads)
	path(orf)
	//path(kraken_cons)
	path(rmdfile)
	output:
	path("Ampliseq_results_report.html")
	script:
	"""
	
	cp ${csv} samples.csv
	cp ${krona_reports_raw} rawreads.html
	# handle empty mapped reads files
	for i in *mappedreads.txt
	do
	 	if [ \$(wc -l < "\${i}" ) -eq 0 ]
		 then
	 		echo "Amplicon_Name Size Reads" >> \${i}
			echo "NA NA NA" >> \${i}
	 	fi
	done
	


	cp ${rmdfile} report.Rmd
	

	Rscript -e 'rmarkdown::render(input="report.Rmd",params=list(csv="samples.csv",krona="rawreads.html"),output_file = "Ampliseq_results_report.html")'
	"""

}
// performs remote blast of the consensus sequences
process blast_cons {
	publishDir "${params.out_dir}/blast/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path(consensus)
	path (taxdb)
	output:
	path("${SampleName}_report_blast.csv")
	
	script:
	"""

	cp ${taxdb}/* ./
	blastn -db nt -query ${consensus} -out ${SampleName}_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue" -max_target_seqs 1 -remote
	grep -v "#" ${SampleName}_blast.csv|sort|uniq > ${SampleName}_report_blast.csv
	sed -i '1i queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tscinames' ${SampleName}_report_blast.csv
	
	"""

}



process orfipy {
	label "low"
	publishDir "${params.out_dir}/orf",mode:"copy"
	input:
	tuple val(SampleName),path("${SampleName}_consensus.fasta")
	output:
	path ("${SampleName}_ORF.fasta")
	script:
	"""
	orfipy ${SampleName}_consensus.fasta --dna ${SampleName}_ORF.fasta --min 700 --outdir ${SampleName}_ORF --start ATG
	mv ${SampleName}_ORF/${SampleName}_ORF.fasta ${SampleName}_ORF.fasta 
	if [ \$(wc -l < "${SampleName}_ORF.fasta") == "0" ]
		then 
			echo -e ">No_consensus/${SampleName}_ORF" > ${SampleName}_ORF.fasta
	else 
		sed -i '/>/ s/ORF.*/ORF/g' ${SampleName}_ORF.fasta
	fi

	"""

}

workflow {
	data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	reference=file("${baseDir}/reference.fasta")
	primerbed=file("${baseDir}/primer.bed")
	//trim barcodes and adapter sequences
	if (params.trim_barcodes){
		porechop(merge_fastq.out)
		minimap2(reference,porechop.out)
		 
	} else {
            minimap2(reference,merge_fastq.out)
		
        }
	// conditional for trim barcodes option
	if (params.trim_barcodes){
		if (params.kraken_db) {
			kraken=params.kraken_db
			kraken2(porechop.out,kraken)
		}
		         
			  
	 } else {
		if (params.kraken_db){
			kraken=params.kraken_db
			kraken2(merge_fastq.out,kraken)
		}
		
	}

	// create consensus
	splitbam(minimap2.out,primerbed)
	
	//condition for kraken2 classification
	if (params.kraken_db){
		kraken=params.kraken_db
		kraken2_consensus(splitbam.out.consensus,kraken)
		kraken_raw=kraken2.out.kraken2_raw
		kraken_cons=kraken2_consensus.out.kraken2_cons
		krona_kraken(kraken_raw.collect(),kraken_cons.collect())
		
	}
	
	
	// qc report using split bam out put
	stats=splitbam.out.unfilt_stats
	idxstats=splitbam.out.idxstats
	multiqc(stats.mix(idxstats).collect())
	
		
	//tax=("${baseDir}/taxdb")
	//blast_cons(splitbam.out.consensus,tax,db1)
	orfipy(splitbam.out.consensus)
	
	//generate report
	rmd_file=file("${baseDir}/Ampliseq.Rmd")
	if (params.kraken_db){
		make_report(make_csv.out,krona_kraken.out.raw,splitbam.out.mapped.collect(),orfipy.out.collect(),rmd_file)
	}
	
	
	
}

