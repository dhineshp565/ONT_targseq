#!/bin/env bash

# This script takes the sam file from minimap2,
# filters q30 and primary alignmnets and creates consensus of each amplicon in a sample
# also generates read statistics on filtered and unfiltered reads read statistics
# ouput are used in almost all subsequet process in the workflow
# $1 = SampleName 
# $2 = Input path of sam file from minimap2, 
# $3 = primerbed file with primer coordinates for primer trimming $4 = threshold for minimum reads to be considered for consensus
# $4 = read count threshold
# $5 = samtools consensus mode "simple or bayesian"
# $6 = qscore

# generate stats prior to read filtering
samtools view -b -h $2|samtools sort > unfilt_$1.bam
samtools stats "unfilt_$1.bam" > $1_unfilt_stats.txt
samtools index "unfilt_$1.bam" > $1_unfilt.bai
samtools idxstats "unfilt_$1.bam" > $1_unfilt_idxstats.csv
#generate a  sorted bam file with primary alignments
samtools view -b -h -F 0x900 -q $6 $2|samtools sort > $1.bam	
samtools stats "$1.bam" > $1_stats.txt

#index sorted bam file and generate read counts for each amplicon
samtools index "$1.bam" > $1.bai
samtools idxstats "$1.bam" > $1_idxstats.txt
threshold=$4
awk -v var="$threshold" '{if ($3 >= var) print $1, $2, $3}' "${1}_idxstats.txt" > "${1}_mappedreads.txt"
# filter mapped reads to get only amplicon names and this is used for igvreport generation
#awk '{if ($3!=0) print $1}' ${1}_mappedreads.txt > ${1}_amplicons.txt

mode=$5
#conditional for empty mapped reads.txt file
if [ $(wc -l < "$1_mappedreads.txt") -ne 0 ]
then 
#using the list of mapped amplicons from text file, bam file is split based on amplicons and consensus is generated
	while read lines
	do 
			amp=$(echo $lines|cut -f1 -d' ')
			# split bam 
			samtools view -b "$1.bam" "${amp}" > $1_${amp}.bam
			#samtools ampliconclip --both-ends -b $3 "$1_${amp}.bam" > $1_trimmed_${amp}.bam
			# generate consensus for full length reads
			if [ "$mode" == "simple" ];then
				samtools consensus -f fasta -m simple -d "$threshold" -A -q "$1_${amp}.bam" > $1_${amp}.fasta
					# change fasta header with sample and amplicon names
				sed -i "s/>.*/>$1_${amp}_consensus/" $1_${amp}.fasta
			else
				samtools consensus -f fasta -m bayesian -A "$1_${amp}.bam" > $1_${amp}.fasta
					# change fasta header with sample and amplicon names
				sed -i "s/>.*/>$1_${amp}_consensus/" $1_${amp}.fasta
			fi
		
			
	done < "$1_mappedreads.txt"
	# merge consensus from all amplicons
	cat $1_*.fasta > $1_cons.fasta
	# convert multiline fasta to single line fasta
	awk '/^>/ {if (seq) print seq; print;seq =""} /^[^>]/{seq=seq$0} END {print seq}' $1_cons.fasta > $1_consensus.fasta
		
# handle empty consensus. when there are no mapped reads.add sequence header
else
		echo -e ">$1 No consensus" > $1_consensus.fasta
		echo -e "No reads found" >> $1_consensus.fasta
		
		echo "NA NA NA" >> "$1_mappedreads.txt"
		

fi
	# insert headers to mappedreads.txt
sed -i "1i Amplicon_Name Size $1" "$1_mappedreads.txt"

