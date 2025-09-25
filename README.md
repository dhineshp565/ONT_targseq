# ONT_targseq

A comprehensive Nextflow pipeline for analyzing targeted amplicon sequencing data from Oxford Nanopore Technologies (ONT) sequencing platforms.

## Overview

This pipeline processes ONT targeted amplicon sequencing data to generate high-quality consensus sequences, taxonomic classifications, phylogenetic analyses, and comprehensive quality control reports. 


## Quick Start

### Prerequisites
- **Nextflow** (≥ 21.04.0), **Docker**, **WSL2** (Windows only)

### Basic Usage

```bash
nextflow run main.nf \
    --input /path/to/input \
    --out_dir Results \
    --reference /path/to/reference.fasta \
    --kraken_db /path/to/kraken2_db \
    --blastdb_path /path/to/blastdb/nt \
    --blastdb_name nt
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input directory containing sample subdirectories |
| `--out_dir` | Output directory for results |
| `--reference` | Path to reference FASTA file |
| `--kraken_db` | Path to Kraken2 database |
| `--blastdb_path` | Path to BLAST database directory |
| `--blastdb_name` | BLAST database name |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--read_count_threshold` | `10` | Minimum read depth for consensus generation |
| `--trim_barcodes` | `false` | Enable barcode/adapter trimming with Porechop |
| `--medaka_model` | `r1041_e82_400bps_sup_g615` | Medaka model for consensus polishing |
| `--qscore` | `10` | Minimum quality score threshold |
| `--consensus_mode`| `'simple' or 'bayesian'`| samtools option for generating consensus. Default: simple

## Input & Output

### Input Structure
```
input_directory/
├── Sample1/
│   └── *.fastq(.gz)
├── Sample2/
│   └── *.fastq(.gz)
└── ...
```

### Required Files
- Reference sequences (FASTA), `primer.bed` file, Kraken2 database, BLAST database

### Output Structure
Results organized by analysis type:
- `merged/` - Merged FASTQ files per sample
- `trimmed/` - Adapter-trimmed reads (if enabled)
- `splitbam/` - Split BAM files and initial consensus sequences
- `medaka/` - Polished consensus sequences
- `multiqc/` - Quality control summary reports
- `kraken2/` - Taxonomic classification results
- `blast/` - BLAST similarity search results
- `mafft/` - Multiple sequence alignments
- `iqtree/` - Phylogenetic tree files
- `ggtree/` - Tree visualization plots
- `igvreports/` - Interactive genome browser reports
- `*.html` - Final comprehensive HTML reports



## Software Dependencies

| Tool | Purpose | Description | Citation |
|------|---------|-------------|----------|
| **nanoq** | Quality control | Ultra-fast quality control tool for nanopore reads with comprehensive statistics and filtering capabilities | [Steinig & Coin (2022)](https://doi.org/10.21105/joss.02991) |
| **Porechop** | Adapter trimming | Tool for finding and removing adapters from Oxford Nanopore reads, with support for demultiplexing and chimeric read detection | [Wick et al.](https://github.com/rrwick/Porechop) |
| **minimap2** | Sequence alignment | Versatile sequence alignment program for aligning DNA or mRNA sequences against large reference databases | [Li (2018)](https://doi.org/10.1093/bioinformatics/bty191) |
| **SAMtools** | BAM processing | Suite of programs for interacting with high-throughput sequencing data in SAM/BAM format | [Danecek et al. (2021)](https://doi.org/10.1093/gigascience/giab008) |
| **Medaka** | Consensus polishing | Tool to create consensus sequences and call variants via neural networks from nanopore sequencing data | [Oxford Nanopore](https://github.com/nanoporetech/medaka) |
| **Kraken2** | Taxonomic classification | System for assigning taxonomic labels to DNA sequences using exact k-mer matches | [Wood et al. (2019)](https://doi.org/10.1186/s13059-019-1891-0) |
| **Krona** | Taxonomic visualization | Interactive metagenomic visualization tool that displays hierarchical data with multi-layered pie charts | [Ondov et al. (2011)](https://doi.org/10.1186/1471-2105-12-385) |
| **BLAST+** | Sequence similarity | Basic Local Alignment Search Tool for finding regions of local similarity between sequences | [Camacho et al. (2009)](https://doi.org/10.1186/1471-2105-10-421) |
| **ABRicate** | Sequence screening | Mass screening of contigs for antimicrobial resistance genes and virulence factors | [Seemann](https://github.com/tseemann/abricate) |
| **MAFFT** | Multiple alignment | Multiple sequence alignment program offering various algorithms including progressive, iterative, and structural methods | [Katoh & Standley (2013)](https://doi.org/10.1093/molbev/mst010) |
| **IQ-TREE** | Phylogenetic analysis | Efficient software for phylogenomic inference with model selection, bootstrap analysis, and tree topology tests | [Nguyen et al. (2015)](https://doi.org/10.1093/molbev/msu300) |
| **BEDTools** | Genomic intervals | Toolkit for genome arithmetic - intersecting, merging, counting, complementing, and shuffling genomic intervals | [Quinlan & Hall (2010)](https://doi.org/10.1093/bioinformatics/btq033) |
| **IGV Reports** | Genome visualization | Python application for generating self-contained HTML reports with embedded genome visualizations | [Robinson et al. (2023)](https://github.com/igvteam/igv-reports) |
| **R Markdown** | Report generation | Dynamic document format combining R code with narrative text to produce elegantly formatted output | [Xie et al. (2018)](https://bookdown.org/yihui/rmarkdown/) |
| **ORFiPy** | ORF prediction | Fast and flexible tool for extracting open reading frames (ORFs) from FASTA files with customizable parameters | [Singh & Wurtele (2021)](https://doi.org/10.1093/bioinformatics/btab090) |
| **MultiQC** | Quality control | Tool that aggregates results from multiple bioinformatics analyses across many samples into a single report | [Ewels et al. (2016)](https://doi.org/10.1093/bioinformatics/btw354) |
| **seqkit** | Sequence manipulation | Cross-platform and ultrafast toolkit for FASTA/Q file manipulation with rich functions | [Shen et al. (2016)](https://doi.org/10.1371/journal.pone.0163962) |
| **ggtree** | Tree visualization | R package for visualization and annotation of phylogenetic trees with associated data | [Yu et al. (2017)](https://doi.org/10.1111/2041-210X.12628) |




