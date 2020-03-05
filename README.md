CircMiner: Accurate and Rapid Detection of Circular RNA through Splice-Aware Pseudo-Alignment Scheme
===================
## What is it?
It is a sensitive computational tool for detecting circular RNAs (circRNAs) from RNA-Seq data.

# Table of contents
1. [Installation](#installation)
2. [Commands Options](#commands-options)
3. [Example Commands](#example-commands)
4. [Output Format](#output-format)
5. [Contact & Support](#contact)

## Installation
The source code can be downloaded from [GitHub](https://github.com/vpc-ccg/circminer). Prerequisite are described below.

### Prerequisite
 - g++ 4.8.5 or higher
 - GNU sort

### Compilation and Configuration
To install, you need to first fetch the repository [git repository](https://github.com/vpc-ccg/circminer) or download the corresponding compressed files. 
```
git clone --recursive https://github.com/vpc-ccg/circminer.git
cd circminer
make -j
```

Now you are ready to go!

## Commands Options

### Synopsis
	
	circminer --index -r FASTA_FILE -k KMER_SIZE [OPTIONS]
	circminer -r FASTA_FILE -g GTF_FILE -f FASTQ_FILE -k KMER_SIZE [OPTIONS]

### OPTIONS
Run `circminer -h` to see available options.

### Example Commands
#### Indexing reference genome:

	$ ./circminer --index -r genome.fasta -k 20 --thread 4

#### Mapping to reference genome and circRNA calling:
	
	$ ./circminer -r genome.fasta -g ga.gtf -f reads.fastq -k 20 -o output [--pe] 

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
|9     |string|Supproting back-splice junction read names (comma-separated)|

The back-splice juntion read mappings are stored in output.candidates.pam. If --pam/--sam is specified in the input arguments, the mapping results will be available in output.mapping.pam/output.mapping.sam file. 
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
|8     |Char  |Relative strand: "+" or "-" (R1)                           |
|9     |int   |Edit distance (R1)                                         |
|10    |string|Chromosome name (R2)                                       |
|11    |int   |Start genomic position (R2)                                |
|12    |int   |End genomic position (R2)                                  |
|13    |int   |Number of aligned basepairs (R2)                           |
|14    |int   |Start of aligned position on read (R2)                     |
|15    |int   |End of aligned position on read (R2)                       |
|16    |Char  |Relative strand: "+" or "-" (R2)                           |
|17    |int   |Edit distance (R2)                                         |
|18    |int   |Insert length                                              |
|19    |int   |Number of junctions happening between two mates            |
|20    |int   |Transcriptomic mapping: 1. Genomic mapping: 0              |


## Simulation Data
The simulated RNA-Seq reads used in this project can be accessed [here](https://figshare.com/projects/CircMiner/76488).

## Contact & Support

Feel free to drop any inquiry to Hossein Asghari [hasghari at sfu dot ca].

## Copyright and License
This software is released under GNU General Public License (v3.0).
