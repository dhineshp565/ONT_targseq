#!/usr/bin/env bash

# $1 = Sampple name , $2 = Consensus sequence

orfipy $2 --dna $1_ORF.fasta --min 700 --outdir $1_ORF --start ATG
mv $1_ORF/$1_ORF.fasta $1_ORF_.fasta 
if [ $(wc -l < "$1_ORF_.fasta") == "0" ]
then 
		echo -e ">No_consensus/$1_ORF" > $1_ORF.fasta
else 
	sed -i '/>/ s/ORF.*/ORF/g' $1_ORF_.fasta
	awk '/^>/ {if (NR>1) printf("\n"); printf("%s\n",$0); next;} {printf("%s",$0);} END {printf("\n");}' $1_ORF_.fasta | awk '{for(i=1; i<=length($0); i+=100) print substr($0, i, 100)}' > $1_ORF.fasta
fi
	