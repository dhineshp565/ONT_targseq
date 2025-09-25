#!/usr/bin/env bash

# $1 = Sample name , $2 = Consensus sequence

# Run orfipy to predict ORFs from the consensus sequence ($2)
# Output DNA sequences to $1_ORF.fasta, minimum ORF length 600, output directory $1_ORF, start codon ATG

orfipy $2 --dna $1_ORF.fasta --min 600 --outdir $1_ORF --start ATG --include-stop

# Check if the output fasta file is empty
if [ "$(wc -l < "$1_ORF/$1_ORF.fasta")" -eq 0 ]
then 
	# If empty, write a placeholder fasta entry indicating no consensus ORF found
	echo -e ">No_600bp_orf_found_$1_consensus" > $1_ORF.fasta

else 
	# Otherwise, clean up the fasta header to standardize ORF naming
	sed -i '/>/ s/ORF.*/ORF/g' $1_ORF/$1_ORF.fasta
	# convert multiline fasta to single line fasta
	awk '/^>/ {if (seq) print seq; print;seq =""} /^[^>]/{seq=seq$0} END {print seq}' $1_ORF/$1_ORF.fasta > $1_ORF.fasta

fi
	
	