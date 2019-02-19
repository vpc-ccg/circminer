Spliced Aligner for Transcriptome-Level Abberration Detection (circRNA, Fusion)
===================
## What is it?
It is a computational tool for detecting structural variations from RNA-Seq reads and gene models.

# Table of contents
1. [Installation](#installation)
2. [Commands Options](#commands-options)
3. [Example Commands](#example-commands)
4. [Output Format](#output-format)
5. [Contact & Support](#contact)

## Installation
The source code can be downloaded from [GitHub](https://github.com/vpc-ccg/circrna). Prerequisite are described below.

### Prerequisite
 - g++ 4.8.5 or higher
 - GNU sort

### Compilation and Configuration
To install, you need to first fetch the repository [git repository](https://github.com/vpc-ccg/circrna) or download the corresponding compressed files. 
```
git clone https://github.com/vpc-ccg/circrna.git
cd circrna
make -j
```

Now you are ready to go!

## Commands Options

### Synopsis
	
	circrna --index -r FASTA_FILE -k KMER_SIZE [OPTIONS]
	circrna -r FASTA_FILE -g GTF_FILE -f FASTQ_FILE -k KMER_SIZE [OPTIONS]

### OPTIONS
Run `circrna -h` to see available options.

### Example Commands
#### Indexing reference genome:

	$ ./circrna --index -r genome.fasta -k 20 --threads 4

#### Mapping to reference genome and circRNA calling:
	
	$ ./circrna -r genome.fasta -g ga.gtf -f reads.fastq -k 20 -o output [--pe] 

## Output Format
The mapping results are stored in .pam files if the coresponding option is used while running the tool.

## Contact & Support

Feel free to drop any inquiry to Hossein Asghari [hasghari at sfu dot ca].

## Copyright and License
This software is released under GNU General Public License (v3.0).
