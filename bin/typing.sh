#!/usr/bin/env bash

# This script performs typing analysis using the abricate tool and processes the results.
# Usage: typing.sh <SampleName> <consensus_sequence> <abricate_database_directory>
# Arguments:
#   $1 - SampleName: The name of the sample being analyzed.
#   $2 - consensus_sequence: The path to the consensus sequence file.
#   $3 - abricate_database_directory: The directory containing the abricate database.


# Run abricate with specified database and parameters, output results to a CSV file named <SampleName>_abricate.csv
abricate --datadir $3 --db targseq -minid 60  -mincov 60 --quiet $2 1> $1_abricate.csv

# Replace '_consensus' in the output CSV file with an empty string
sed -i -e "s/_consensus//g" -e "s/_medaka//g" "$1_abricate.csv"

# Check if the output CSV file has less than 2 lines, if so, append a line indicating no consensus and none for all other fields
if [[ $(wc -l < "$1_abricate.csv") -lt 2 ]]
then
	printf "$1.fasta\t$1_No_consensus\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNo_consensus" >> $1_abricate.csv	
fi

# Convert the consensus sequence file to a tab-separated format and save it as <SampleName>.tsv
awk '/^>/ {if (seq) {print seq}; printf "%s\t", substr($0, 2); seq=""} !/^>/ {seq = seq $0} END {if (seq) {print seq}}' $2 > $1.tsv	

# Add a header to the TSV file with sequence header and sequence columns
sed -i '1i SEQ HEADER\tSEQUENCE' $1.tsv

# Combine the abricate results CSV and the sequence TSV into a single CSV file named <SampleName>_withseq.csv
paste $1_abricate.csv $1.tsv > $1_withseq.csv