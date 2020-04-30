CircMiner: Accurate and Rapid Detection of Circular RNA through Splice-Aware Pseudo-Alignment Scheme
===================
A sensitive and fast computational tool for detecting circular RNAs (circRNAs) from RNA-Seq data.

# Table of contents
1. [Installation](#installation)
2. [Commands Options](#commands-options)
3. [Test Run](#test-run)
4. [Example Commands](#example-commands)
5. [Output Files](#output-files)
6. [Output Format](#output-format)
7. [Contact & Support](#contact--support)

## Installation
### BIOCONDA

CircMiner can be installed using [conda](https://conda.io/) package manager via [bioconda](https://bioconda.github.io/) channel:
```
conda install -c bioconda circminer
```
### Installation from Source
> *Prerequisite.* You will need (i) GCC 4.9.4 or higher; and (ii) GNU sort to compile the source code.

To install, you need to first fetch the [git repository](https://github.com/vpc-ccg/circminer) or download the corresponding compressed files. 
```
git clone --recursive https://github.com/vpc-ccg/circminer.git
cd circminer
make 
```

Now you are ready to go!

## Commands Options

### Synopsis
	
	circminer --index -r FASTA_FILE [OPTIONS]
	circminer -r FASTA_FILE -g GTF_FILE -1 FASTQ_FILE_R1 -2 FASTQ_FILE_R2 [OPTIONS]

### OPTIONS
Run `circminer -h` to see available options.

#### Indexing options
	-i, --index:            Indicates the indexing stage.
	-m, --compact-index:	Use this option only while building the index to enable compact version of the index.
	-k, --kmer:		Kmer size [14..22] (default = 20).

#### General options:
	-r, --reference:	Reference file.
	-g, --gtf:	Gene model file.
	-s, --seq:	Single-end sequence file.
	-1, --seq1:	1st paired-end sequence file.
	-2, --seq2:	2nd paired-end sequence file.

#### Advanced options:
	-l, --rlen:		Max read length (default = 300).
	-e, --max-ed:		Max allowed edit distance on each mate (default = 4).
	-c, --max-sc:		Max allowed soft clipping on each mate (default = 7).
	-w, --band:		Band width for banded alignment (default = 3).
	-S, --seed-lim:		Skip seeds that have more than INT occurrences (default = 500).
	-T, --max-tlen:		Maximum template length of concordant mapping. Paired-end mode only (default = 500).
	-I, --max-intron:	Maximum length of an intron (default = 2000000).
	-C, --max-chain-list:	Maximum number of chained candidates to be processed (default = 30).
	-o, --output:		Prefix of output files (default = output).
	-t, --thread:		Number of threads (default = 1).
	-A, --sam:		Enables SAM output for aligned reads. Cannot be set along with --pam.
	-P, --pam:		Enables custom pam output for aligned reads. Cannot be set along with --sam.
	-d, --verbose:		Verbose mode: 0 to 1. Higher values output more information (default = 0).
	-a, --scan-lev:		Transcriptome/Genome scan level: 0 to 2. (default = 0)
				0: Report the first mapping.
				1: Continue processing the read unless it is perfectly mapped to cDNA.
				2: Report the best mapping.

#### Other options:
	-h, --help:	Shows help message.
	-v, --version:	Current version.

## Test Run
To make sure CircMiner is successfully installed and runable on your machine, you can download a small sample test [here](https://ndownloader.figshare.com/files/22423638). Unzip the downloaded package and use the following commands to run CircMiner for building the index and calling circRNAs.

### Example Commands
#### Indexing reference genome:

	$ ./circminer --index -r ref.fa -k 20 --thread 4

#### Mapping to reference genome and circRNA calling:
	
	$ ./circminer -r ref.fa -g ref.gtf -1 R1.fq -2 R2.fq -o output

After a successful run `output.circ_report` should contain the following two lines:

|   |        |        |   |     |       |       |      |                                                                                  |
|--:|:------:|:------:|:-:|:---:|:-----:|:-----:|:----:|:---------------------------------------------------------------------------------|
|1  |586821  |608056  |8  |STC  |CT-AG  |CT-AG  |Pass  |Circ1-12,Circ1-56,Circ1-18,Circ1-50,Circ1-110,Circ1-80,Circ1-2,Circ1-4            |
|1  |805799  |810170  |9  |STC  |TT-TA  |TT-TA  |Pass  |Circ2-32,Circ2-40,Circ2-36,Circ2-74,Circ2-76,Circ2-80,Circ2-60,Circ2-44,Circ2-2   |

### Annotate CircMiner's output file
To annotate `output.circ_report` with the corresponding gene, transcript, and exon numbers please run the following script:

    $ python scripts/annotate_transcript.py output.circ_report ref.gtf output.circ_report.annotated

### Converting GTF file
If you use UCSC GTF files please use the following script to convert them to Ensembl GTF style and use the converted GTF when running CircMiner.

    $ python scripts/convertGTF.py INPUT_GTF OUTPUT_GTF

## Output Files
When a successful run finishes, the structure of output directory (e.g. outdir) will be as follows:
```
outdir                                         # output directory
├── output.circ_report                         # detected circRNA report
├── output.candidates.pam                      # back-splice juntion read mappings
└── output.mapping.pam/output.mapping.sam      # pseudo-alignment mapping results (only available when corresponding argument is passed)
```

## Output Format
The information regarding the detected circRNAs is reported in `output.circ_report` file. It includes the exact breakpoint location, number of supporting back-splice junctions and their read names.

|Column|Type  |Description                                                              |
|-----:|:----:|:------------------------------------------------------------------------|
|1     |string|Chromosome name                                                          |
|2     |int   |Start genomic position of circRNA                                        |
|3     |int   |End genomic position of circRNA                                          |
|4     |int   |Number of supporting back-splice junction reads                          |
|5     |string|Type of circRNA                                                          |
|6     |string|Consensus of splice signal on supporting back-splice junction reads      |
|7     |string|Splice signal on reference                                               |
|8     |string|Pass/Fail (based on matching splice signal to reference)                 |
|9     |string|Supproting back-splice junction read names (comma-separated)             |

The back-splice juntion read mappings are stored in `output.candidates.pam`. If --pam/--sam is specified in the input arguments, the mapping results will be available in `output.mapping.pam/output.mapping.sam` file. 

Note: If the scan level parameter (`-a, --scan-lev`) is set to 2 while running the tool, the mapping with the smallest error and soft-clip values is reported.

PAM mapping format:

|Column|Type  |Description                                                |
|-----:|:----:|:----------------------------------------------------------|
|1     |string|Read name                                                  |
|2     |string|Chromosome name (R1)                                       |
|3     |int   |Start genomic position (R1)                                |
|4     |int   |End genomic position (R1)                                  |
|5     |int   |Number of aligned basepairs (R1)                           |
|6     |int   |Start of aligned position on read (R1)                     |
|7     |int   |End of aligned position on read (R1)                       |
|8     |char  |Relative strand: "+" or "-" (R1)                           |
|9     |int   |Edit distance (R1)                                         |
|10    |string|Chromosome name (R2)                                       |
|11    |int   |Start genomic position (R2)                                |
|12    |int   |End genomic position (R2)                                  |
|13    |int   |Number of aligned basepairs (R2)                           |
|14    |int   |Start of aligned position on read (R2)                     |
|15    |int   |End of aligned position on read (R2)                       |
|16    |char  |Relative strand: "+" or "-" (R2)                           |
|17    |int   |Edit distance (R2)                                         |
|18    |int   |Insert length                                              |
|19    |int   |Number of junctions happening between two mates            |
|20    |int   |Transcriptomic mapping: 1. Genomic mapping: 0              |


## Contact & Support
For any bug report, feature request, or questions please fill out an issue through CircMiner's [issue page](https://github.com/vpc-ccg/circminer/issues).

## Citation
If you use CircMiner please cite our paper:

Asghari H., Lin YY., Xu Y., Haghshenas E., Collins CC., Hach F. __[CircMiner: Accurate and Rapid Detection of Circular RNA through Splice-Aware Pseudo-Alignment Scheme](https://doi.org/10.1093/bioinformatics/btaa232)__. Bioinformatics (2020). btaa232

## Simulation Data
The simulated RNA-Seq reads used in this project can be accessed [here](https://figshare.com/projects/CircMiner/76488).

## Copyright and License
This software is released under GNU General Public License (v3.0).
