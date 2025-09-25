#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { make_csv } from './modules/local/make_csv'
include { merge_fastq } from './modules/local/merge_fastq'
include { porechop } from './modules/local/porechop'
include { minimap2 } from './modules/local/minimap2'
include { splitbam } from './modules/local/splitbam'
include { medaka } from './modules/local/medaka'
include { multiqc } from './modules/local/multiqc'
include { kraken2 } from './modules/local/kraken2'
include { krona_kraken } from './modules/local/krona_kraken'
include { make_report } from './modules/local/make_report'
include { htmltopdf } from './modules/local/htmltopdf'
include { blast_cons } from './modules/local/blast_cons'
include { orfipy } from './modules/local/orfipy'
include { abricate } from './modules/local/abricate'
include { make_limsfile } from './modules/local/make_limsfile'
include { mafft } from './modules/local/mafft'
include { iqtree } from './modules/local/iqtree'
include { ggtree } from './modules/local/ggtree'
include { extract_mapped_ref } from './modules/local/extract_mapped_ref'
include { mapped_ref_bed } from './modules/local/mapped_ref_bed'
include { bedtools } from './modules/local/bedtools'
include { igvreports } from './modules/local/igvreports'
include { seq_length } from './modules/local/seq_length'

workflow {
	data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)}, params.qscore)
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
	splitbam(minimap2.out,primerbed,params.read_count_threshold,params.consensus_mode,params.qscore)

	// Pair fastq and consensus files by sample name
	paired_fastq_consensus = merge_fastq.out.join(splitbam.out.consensus)
    .map { sample, fastq, consensus -> tuple(sample, fastq, consensus) }

// medaka polishing
	medaka(paired_fastq_consensus)
	
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
	idxstats=splitbam.out.unfilt_idx
	multiqc(stats.mix(idxstats).collect())
	dbdir=file("${baseDir}/targseq")
	
	abricate(splitbam.out.consensus,dbdir)
	make_limsfile(abricate.out.withseq.collect(),software_version_file)
	
	blast_cons(splitbam.out.consensus,params.blastdb_path,params.blastdb_name)

	refdir="${baseDir}/reference_sequences"
	mafft(make_csv.out,splitbam.out.cons_only.collect(),refdir)
	iqtree(mafft.out.collect())
	ggtree(iqtree.out.collect())
	orfipy(medaka.out.consensus)

	// Pair consensus and ORF files by sample name
	paired_consensus_orf = splitbam.out.consensus.join(orfipy.out.orf)
    .map { sample, consensus, orf -> tuple(sample, consensus, orf) }
	seq_length(paired_consensus_orf)
	
	
	//generate report


	
	bedtools(splitbam.out.target_bam)
	// extract_mapped_ref(splitbam.out.amplicons,reference)
	mapped_ref_bed(reference)
	igvreports(make_csv.out,reference,mapped_ref_bed.out.bed,bedtools.out.collect())


	rmd_file=file("${baseDir}/targseq_rmdfile.Rmd")
	rmdfile_case=file("${baseDir}/targseq_rmdfile_case.Rmd")

	make_report(make_csv.out,krona_kraken.out.raw,splitbam.out.mapped.collect(),splitbam.out.cons_only.collect(),abricate.out.abricate.collect(),blast_cons.out.blast_formatted.collect(),ggtree.out.png,rmd_file,igvreports.out,orfipy.out.orf_only.collect(),seq_length.out.collect(),rmdfile_case)
	// htmltopdf(make_report.out.pdf)

}

