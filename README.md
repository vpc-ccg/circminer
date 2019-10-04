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
git clone https://github.com/vpc-ccg/circminer.git
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

	$ ./circminer --index -r genome.fasta -k 20 --threads 4

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
|6     |string|Read names of the supproting back-splice junction reads (comma-separated)|

The mapping results are stored in .pam files. If the scan level parameter (`-s, --scan-lev`) is set to 2 while running the tool, the mapping with the smallest error and soft-clip values is reported.

## Contact & Support

Feel free to drop any inquiry to Hossein Asghari [hasghari at sfu dot ca].

## Copyright and License
This software is released under GNU General Public License (v3.0).
