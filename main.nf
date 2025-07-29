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
				nanoq -i ${SampleName}.fastq -s -H > ${SampleName}_readstats.csv
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
	path("${SampleName}_mappedreads.txt"),emit:mapped
	path("${SampleName}_idxstats.txt"),emit:idxstats
	tuple val(SampleName),path("${SampleName}_consensus.fasta"),emit:consensus
	path("${SampleName}_consensus.fasta"),emit:(cons_only)
	path("${SampleName}_unfilt_stats.txt"),emit:unfilt_stats
	path("${SampleName}_unfilt_idxstats.csv"),emit:unfilt_idx
	// tuple val(SampleName),path ("${SampleName}_amplicons.txt"),emit:amplicons
	tuple val(SampleName),path ("${SampleName}.bam"),emit:target_bam
	

	script:
	"""
	splitbam.sh ${SampleName} ${SamplePath} ${primerbed} ${params.read_count_threshold}

	"""
}

process medaka {
	publishDir "${params.out_dir}/medaka",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path(SamplePath)
	tuple val(SampleName),path(consensus)
	val (medaka_model)
	output:
	tuple val(SampleName),path("${SampleName}_medaka_consensus.fasta"),emit:consensus
	path("${SampleName}_medaka_consensus.fasta"),emit:cons_only
	script:
	"""
	
	if [ \$(wc -l < "${consensus}" ) -gt 1 ]
		then
		medaka_consensus -i ${SamplePath} -d ${consensus} -o ${SampleName}_medaka_consensus -m ${medaka_model}
		mv ${SampleName}_medaka_consensus/consensus.fasta ${SampleName}_medaka_consensus.fasta
	else 
		echo ">${SampleName} No consensus sequences to polish" > ${SampleName}_medaka_consensus.fasta
	fi
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

//krona plots
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

//make html report with rmarkdown
process make_report {
	publishDir "${params.out_dir}/",mode:"copy"
	
	label "medium"
	input:
	path(csv)
	path(krona)
	path(mappedreads)
	path(cons_only)
	path(abricate)
	path(blast_formatted)
	path(png)
	path(rmdfile)
	path(igv)
	path (orf)
	
	output:
	
	path("ONT_targseq*.html")
	script:
	"""
	
	cp ${rmdfile} rmdfile.Rmd
	

	Rscript -e "rmarkdown::render(input= 'rmdfile.Rmd', params=list(csv= '${csv}', krona= '${krona}', png='*.png', igv='${igv}'), output_file=paste0('ONT_targseq_results_report_',format(Sys.time(), '%Y-%m-%d_%H-%M-%S'),'.html'))"

	

    """

}


//performs blast of the consensus sequences
process blast_cons {
	publishDir "${params.out_dir}/blast/",mode:"copy"
	containerOptions "-v ${params.blastdb_path}:${params.blastdb_path}"

	label "high"
	input:
	tuple val(SampleName),path(consensus)
	//tuple val(SampleNAme),path(sanger_reads)
	path (blastdb_path)
	val(blastdb_name)

	output:
	path("${SampleName}_report_blast.csv"), emit:blast_formatted
	path("${SampleName}_blast.csv")
	
	script:
	"""
	

	# Run BLAST
	blastn -db ${blastdb_path}/${blastdb_name} -query ${consensus} -out ${SampleName}_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue staxids ssciname scomnames stitle" -max_target_seqs 5

	# Check for "# 0 hits found" and create an empty report if no hits found
	if grep -q "# 0 hits found" "${SampleName}_blast.csv"; then
		# Write header and 'none' row to report
		echo -e "queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tstaxids\tsscinames\tscomnames\tstitle" > "${SampleName}_report_blast.csv"
		echo -e "none\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone" >> "${SampleName}_report_blast.csv"
	else
		# remove lines starting with "#" and sort to get unique entries
		grep -v "#" "${SampleName}_blast.csv" | sort | uniq > "${SampleName}_report_blast.csv"
		# add header to the report
		sed -i '1i queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tstaxids\tsscinames\tscomnames\tstitle' "${SampleName}_report_blast.csv"
	fi
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
	orfy.sh ${SampleName} ${SampleName}_consensus.fasta

	"""

}

process abricate{
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

process make_LIMSfile {
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
process mafft {
	publishDir "${params.out_dir}/casefiles/",mode:"copy"
	label "low"
	input:
	path (csv)
	path (consensus)
    path (reference_sequences)
	output:
	path("*_msa.fasta")
	script:
	"""
	mafft.sh ${csv} ${reference_sequences}
	
	"""

}


// Process to run IQ-TREE for phylogenetic analysis
// This process takes multiple MSA files and runs IQ-TREE on each of them
process iqtree {
	publishDir "${params.out_dir}/iqtree/", mode: "copy"
	label "low"

	input:
	path (msa)

	output:
	path "*.treefile"

	script:
	"""
	#run iqtree2 for each msa file
    for file in ${msa};do

		# get the prefix of the file name to check if it is a no_msa file
		prefix=\$(basename "\${file}" .fasta)

		# check if the prefix is "no_msa" and skip iqtree analysis if it is
		if [ "\$prefix" = "no_msa" ]; then
			echo "no virus sequences found for this sample, skipping IQ-TREE analysis" > "\${prefix}.treefile"

		else
			# run iqtree2 with MFP model selection and 1000 bootstrap replicates
			iqtree2 -s "\${file}" -m MFP -bb 1000 -nt AUTO -pre "\${prefix}"
		fi
	done
    """

}

// this process generates a tree plot using ggtree
// It takes the tree files generated by iqtree and creates PNG files for visualization
process ggtree {
    publishDir "${params.out_dir}/ggtree/", mode: "copy"
    label "low"

    input:
    path (treefiles)

    output:
    path "*.png",emit:png
	

    script:
    """
	for treefile in ${treefiles}; do
    filename=\$(basename "\$treefile")
		if [[ "\$filename" != *no_msa* ]]; then
			plot_tree.R "\$treefile"
		else 
			echo "no virus sequences found for this sample, skipping ggtree analysis" > "\$filename.png"
		fi
	done
    """
}

process extract_mapped_ref {
	publishDir "${params.out_dir}/mapped_ref"
	input:
	tuple val(SampleName), path(amplicontxtfile)
	path (reference)
	output:
	val(SampleName)
	path("${SampleName}_mapped_ref.fasta"), emit:mapped_fasta
	script:
	"""
	seqkit grep -f ${amplicontxtfile} ${reference} > ${SampleName}_mapped_ref.fasta
	"""
}

//Generate fai and bed from mapped reference
process mapped_ref_bed {
	publishDir "${params.out_dir}/mapped_ref_bed"	
	input:
	path (reference)
	output:
	path ("*.bed")
	script:
	"""
	samtools faidx ${reference} > targseq_reference.fasta.fai
    awk 'BEGIN {FS=OFS="\t"}; {print \$1,0,\$2}' "targseq_reference.fasta.fai" > targseq_reference.bed
	"""
}
// generating bedgraph files from alignments
process bedtools {
	publishDir "${params.out_dir}/bedtools/"
	input:
	tuple val(SampleName), path (bam)
	output:
	path ("${SampleName}*.bedgraph")
	script:
	"""
	bedtools genomecov -ibam ${SampleName}.bam -bga > ${SampleName}.bedgraph
	"""

}
//igv reports from bedgraph
process igvreports {
	publishDir "${params.out_dir}/igvreports/", mode: "copy"
	input:
	path(csv)
	path(reference)
    path(bed)
	path(bedgraph)
	output:
	path{"*.html"}
	script:
	"""

	create_report ${bed} ${reference} --tracks *.bedgraph --output igv_coverage.html

	"""
}

		
workflow {
	data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	reference=file(params.reference)
	primerbed=file("${baseDir}/primer.bed")
	software_version_file=file("${baseDir}/software_version.tsv")
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

	// medaka polishing
	medaka(merge_fastq.out,splitbam.out.consensus,params.medaka_model)

	//condition for kraken2 classification
	if (params.kraken_db){
		kraken=params.kraken_db
		//kraken2_consensus(medaka.out.consensus,kraken)
		kraken_raw=kraken2.out.kraken2_raw
		//kraken_cons=kraken2_consensus.out.kraken2_cons
		krona_kraken(kraken_raw.collect())
		
	}
	
	
	// qc report using split bam out put
	stats=splitbam.out.unfilt_stats
	idxstats=splitbam.out.idxstats
	multiqc(stats.mix(idxstats).collect())
	dbdir=file("${baseDir}/targseq")
	
	abricate(splitbam.out.consensus,dbdir)
	make_LIMSfile(abricate.out.withseq.collect(),software_version_file)
	
	blast_cons(splitbam.out.consensus,params.blastdb_path,params.blastdb_name)

	refdir="${baseDir}/reference_sequences"
	mafft(make_csv.out,splitbam.out.cons_only.collect(),refdir)
	iqtree(mafft.out.collect())
	ggtree(iqtree.out.collect())
	orfipy(splitbam.out.consensus)
	
	//generate report


	
	bedtools(splitbam.out.target_bam)
	// extract_mapped_ref(splitbam.out.amplicons,reference)
	mapped_ref_bed(reference)
	igvreports(make_csv.out,reference,mapped_ref_bed.out,bedtools.out.collect())


	rmd_file=file("${baseDir}/targseq_rmdfile.Rmd")

	make_report(make_csv.out,krona_kraken.out.raw,splitbam.out.mapped.collect(),splitbam.out.cons_only.collect(),abricate.out.abricate.collect(),blast_cons.out.blast_formatted.collect(),ggtree.out.png,rmd_file,igvreports.out,orfipy.out.collect())
	

}

