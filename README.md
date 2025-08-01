# ONT_targseq
Pipeline for analysing targeted amplicon sequencing using Oxford nanopore sequencing

Requires input directory with sub-directories with fastq file, reference sequence (fasta), kraken2 database and blast database
Outputs consensus sequences, kraken, krona, igv report and multiqc report



Usage:
```
nextflow run main.nf --input path_to_input --out_dir Results --kraken_db path_to_kraken_database --reference path_to_reference.fasta --blastdb_path /path/to/blastdb/nt --blastdb_name nt
```

```
Parameters:

--input        Path to input directory
--out_dir      Output directory
--reference    Path to fasta file with reference sequences
--kraken_db    Path to kraken database
--blastdv_path Path to blast database
--blastdb_name A string denoting the name of the blastdb

optional
--read_count_threshold  An interger denoting read depth. Default: 10
--trim_barcodes         barcode and adapter trimming using porechop
--medaka_model         Select basecalling model for polishing consensus sequence. Default mode:'r1041_e82_400bps_sup_g615'

```
## Dependencies
* nextflow
* docker
* wsl2
## Software Used
* nanoq (Steinig and Coin (2022). Nanoq: ultra-fast quality control for nanopore reads. Journal of Open Source Software, 7(69), 2991, https://doi.org/10.21105/joss.02991)
* porechop (https://github.com/rrwick/Porechop)
* minimap2 (Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191)
* samtools (https://github.com/samtools/samtools)
* kraken2 (https://ccb.jhu.edu/software/kraken2/)
* krona:2.7.1(Ondov, B.D., Bergman, N.H. & Phillippy, A.M. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12, 385 (2011). https://doi.org/10.1186/1471-2105-12-385)
* abricate (https://github.com/tseemann/abricate)
* rmarkdown (https://rmarkdown.rstudio.com/)
* bedtools (https://bedtools.readthedocs.io/en/latest/)
* igvreports (https://github.com/igvteam/igv-reports)
* Blast+ (https://doi.org/10.1186/1471-2105-10-421)
* mafft (https://mafft.cbrc.jp/alignment/server/index.html)
* iqtree (https://iqtree.github.io/)
* orfipy (https://github.com/urmi-21/orfipy)
