#!/usr/bin/env bash
# This script runs MAFFT on a given input file and outputs the result to a specified output file.
# list of viruses to process
# $1 = path to the directory containing reference sequences
virus_list=("FAdV" "IBDV" "IBV")

# concatenate all the consensus sequences into one file
cat *.fasta > all_consensus.fasta

# loop through each virus in the list
for virus in "${virus_list[@]}"; do
    # check if the virus is present in the concatenated file
    if grep -q "$virus" all_consensus.fasta; then
		# extract sequences for the virus and create a new file
		awk -v v="$virus" '/^>/ {keep = ($0 ~ v)} keep' all_consensus.fasta > "${virus}_only.fasta"
		# concatenate the extracted sequences with the reference sequences
        cat "${virus}_only.fasta" $1/${virus}_references.fasta > ${virus}_unaligned.fasta
		# run mafft on the concatenated file
        mafft --auto ${virus}_unaligned.fasta > ${virus}_msa.fasta
	# if the virus is not present, create an empty file
	else
		echo ">Novirusfromlist" > no_msa.fasta
    fi
done