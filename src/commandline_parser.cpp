#include <bitset>

#include "commandline_parser.h"

using namespace std;

uint32_t seedLim = FRAGLIM;
bool indexMode = false;
bool compactIndex = false;
bool pairedEnd = false;
bool finalCleaning = true;
bool internalSort = false;
int kmer = 19;
int maxReadLength = 120;
int verboseMode = 0;
int scanLevel = 0;
int maxEd = EDTH;
int maxSc = SOFTCLIPTH;
int bandWidth = INDELTH;
int maxTlen = MAXTLEN;
int maxIntronLen = MAXINTRON;
int maxChainLen = BESTCHAINLIM;
int threadCount = 1;
int stage = 2;
int maxCheckSumLen = sizeof(uint16_t) * 8 / 2;	// 2bit per bp
int reportMapping = DISCARDMAPREPORT;
int seqcnt = 0;

char gtfFilename[FILE_NAME_LENGTH];
char referenceFilename[FILE_NAME_LENGTH];
char fastqFilename[2][FILE_NAME_LENGTH];
char outputFilename[FILE_NAME_LENGTH]="output";
char outputDir[FILE_NAME_LENGTH] = "./";

char* contigName;
int contigNum;

uint32_t lookup_cnt;

vector <bitset <DEF_CONTIG_MAX_SIZE> > near_border_bs;
vector <bitset <DEF_CONTIG_MAX_SIZE> > intronic_bs;

/**********************************************/
// return:
// 0: normal execution
// 1: exit
int parse_command( int argc, char *argv[] )
{

	if (argc < 2) {
		printHELPShort();
		return 1;
	}

	int opt;
	int opt_index;

	static struct option long_opt[] = 
	{
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"index", no_argument, 0, 'i'},
		{"compact-index", no_argument, 0, 'm'},
		{"seq", required_argument, 0, 's'},
		{"seq1", required_argument, 0, '1'},
		{"seq2", required_argument, 0, '2'},
		{"reference", required_argument, 0, 'r'},
		{"gtf", required_argument, 0, 'g'},
		{"kmer", required_argument, 0, 'k'},
		{"rlen", required_argument, 0, 'l'},
		{"output", required_argument, 0, 'o'},
		{"verbose", required_argument, 0, 'd'},
		{"thread", required_argument, 0, 't'},
		{"scan-lev", required_argument, 0, 'a'},
		{"max-ed", required_argument, 0, 'e'},
		{"max-sc", required_argument, 0, 'c'},
		{"band", required_argument, 0, 'w'},
		{"seed-lim", required_argument, 0, 'S'},
		{"max-tlen", required_argument, 0, 'T'},
		{"max-intron", required_argument, 0, 'I'},
		{"max-chain-list", required_argument, 0, 'C'},
		{"stage", required_argument, 0, 'q'},
		{"keep-intermediate", no_argument, 0, 'z'},
		{"internal-sort", no_argument, 0, 'Z'},
		{"sam", no_argument, 0, 'A'},
		{"pam", no_argument, 0, 'P'},
		{0,0,0,0},
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "hvims:1:2:r:g:k:l:o:t:d:a:e:c:w:S:T:I:C:q:zZAP", long_opt, &opt_index )  ) ) 
	{
		switch(opt)
		{
			case 'h': {
				printHELP();
				return 1;
			}
			case 'v': {
				fprintf(stdout, "CircMiner %s.%s\n", versionNumberMajor, versionNumberMinor);
				return 1;
			}
			case 'i': {
				indexMode = true;
				break;
			}
			case 'm': {
				compactIndex = true;
				break;
			}
			case 's': {
				++seqcnt;
				pairedEnd = false;
				strncpy(fastqFilename[0], optarg, FILE_NAME_LENGTH );
				break;
			}
			case '1': {
				++seqcnt;
				pairedEnd = false;
				strncpy(fastqFilename[0], optarg, FILE_NAME_LENGTH );
				break;
			}
			case '2': {
				++seqcnt;
				pairedEnd = true;
				strncpy(fastqFilename[1], optarg, FILE_NAME_LENGTH );
				break;
			}
			case 'r': {
				strncpy(referenceFilename, optarg, FILE_NAME_LENGTH);
				break;
			}		
			case 'g': {
				strncpy(gtfFilename, optarg, FILE_NAME_LENGTH );
				//gtf_flag = 1;
				break;
			}
			case 'k': {
				kmer = atoi(optarg);
				break;
			}
			case 'l': {
				maxReadLength = atoi(optarg);
				break;
			}
			case 'o': {
				strncpy(outputFilename, optarg, FILE_NAME_LENGTH);
				break;
			}
			case 't': {
				threadCount = atoi(optarg);
				if (threadCount < 1 || threadCount > sysconf( _SC_NPROCESSORS_ONLN ))
					threadCount = sysconf( _SC_NPROCESSORS_ONLN );
				break;
			}
			case 'd': {
				verboseMode = atoi(optarg);
				break;
			}
			case 'a': {
				scanLevel = atoi(optarg);
				break;
			}
			case 'e': {
				maxEd = atoi(optarg);
				break;
			}
			case 'c': {
				maxSc = atoi(optarg);
				break;
			}
			case 'w': {
				bandWidth = atoi(optarg);
				break;
			}
			case 'S': {
				seedLim = atoi(optarg);
				break;
			}
			case 'T': {
				maxTlen = atoi(optarg);
				break;
			}
			case 'I': {
				maxIntronLen = atoi(optarg);
				break;
			}
			case 'C': {
				maxChainLen = atoi(optarg);
				break;
			}
			case 'q': {
				stage = atoi(optarg);
				if (stage > 2) {
					fprintf(stderr, "Invalid stage number: %d\nAborted\n", stage);
					return 1;
				}
				break;
			}
			case 'z': {
				finalCleaning = false;
				break;
			}
			case 'Z': {
				internalSort = true;
				break;
			}
			case 'A': {
				reportMapping = SAMFORMAT;
				break;
			}
			case 'P': {
				reportMapping = PAMFORMAT;
				break;
			}
			case '?': {
				//fprintf(stderr, "Unknown parameter: %s\n", long_opt[opt_index].name);
				exit(1);
				break;
			}
			default: {
				printHELPShort();
				return 1;
				}
		}
	}

	if (seqcnt == 0)
	{
		fprintf(stderr, "ERROR: no sequence file is provided.\n");
		return 1;
	}

	if (!pairedEnd and seqcnt != 1)
	{
		fprintf(stderr, "ERROR: %d single-end sequence files provided.\nOnly one is accepted.\n", seqcnt);
		return 1;
	}

	if (pairedEnd and seqcnt != 2)
	{
		fprintf(stderr, "ERROR: %d paired-end sequence file(s) provided.\nOnly two files are accepted.\n", seqcnt);
		return 1;
	}

	if (kmer > WINDOW_SIZE + maxCheckSumLen || kmer < WINDOW_SIZE)
	{
		fprintf(stderr, "ERROR: kmer size should be in [%d..%d]\n", WINDOW_SIZE, WINDOW_SIZE + maxCheckSumLen);
		return 1;
	}

	checkSumLength = (WINDOW_SIZE > kmer) ? 0 : kmer - WINDOW_SIZE;

	if (indexMode)
	{
		CONTIG_SIZE		= DEF_CONTIG_SIZE;
		CONTIG_MAX_SIZE	= DEF_CONTIG_MAX_SIZE;

		sprintf(fileName[0], "%s", referenceFilename);
		sprintf(fileName[1], "%s.index", fileName[0]);
	}

	initCommon();
	
	THREAD_COUNT = threadCount;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	return 0;
}

/**********************************************/
void printHELPShort(void) 
{
	fprintf(stdout, "usage: circminer --index ref.fa [options]\n");
	fprintf(stdout, "       circminer -r ref.fa -g gene_model.gtf -1 reads_1.fastq -2 reads_2.fastq [options]\n");
	fprintf(stderr, "For more details and command line options run \"circminer --help\"\n");
}
/**********************************************/
void printHELP(void)
{
	fprintf(stdout, "\nSYNOPSIS\n");
	fprintf(stdout, "\tcircminer --index FILE [options]\n");
	fprintf(stdout, "\tcircminer -r FILE -g FILE -s FILE [options]\n");
	fprintf(stdout, "\tcircminer -r FILE -g FILE -1 FILE -2 FILE [options]\n");
	fprintf(stdout, "\nGeneral options:\n");
	fprintf(stdout, "\t-r, --refernce:\tReference file.\n");
	fprintf(stdout, "\t-g, --gtf:\tGene model file.\n");
	fprintf(stdout, "\t-s, --seq:\tSingle-end sequence file.\n");
	fprintf(stdout, "\t-1, --seq1:\t1st paired-end sequence file.\n");
	fprintf(stdout, "\t-2, --seq2:\t2nd paired-end sequence file.\n");
	
	fprintf(stdout, "\nAdvanced options:\n");
	fprintf(stdout, "\t-m, --compact-index:\tUse this option only while building the index to enable compact version of the index.\n");
	fprintf(stdout, "\t-k, --kmer:\t\tKmer size [%d..%d] (default = 19).\n", WINDOW_SIZE, WINDOW_SIZE + maxCheckSumLen);
	fprintf(stdout, "\t-l, --rlen:\t\tMax read length (default = 120).\n");
	fprintf(stdout, "\t-e, --max-ed:\t\tMax allowed edit distance on each mate (default = %d).\n", EDTH);
	fprintf(stdout, "\t-c, --max-sc:\t\tMax allowed soft clipping on each mate (default = %d).\n", SOFTCLIPTH);
	fprintf(stdout, "\t-w, --band:\t\tBand width for banded alignment (default = %d).\n", INDELTH);
	fprintf(stdout, "\t-S, --seed-lim:\t\tSkip seeds that have more than INT occurrences (default = %d).\n", FRAGLIM);
	fprintf(stdout, "\t-T, --max-tlen:\t\tMaximum template length of concordant mapping. Paired-end mode only (default = %d).\n", MAXTLEN);
	fprintf(stdout, "\t-I, --max-intron:\tMaximum length of an intron (default = %d).\n", MAXINTRON);
	fprintf(stdout, "\t-C, --max-chain-list:\tMaximum number of chained candidates to be processed (default = %d).\n", BESTCHAINLIM);
	fprintf(stdout, "\t-o, --output:\t\tOutput file (default = output).\n");
	fprintf(stdout, "\t-t, --thread:\t\tNumber of threads (default = 1).\n");
	fprintf(stdout, "\t-A, --sam:\t\tEnables SAM output for aligned reads. Cannot be set along with --pam.\n");
	fprintf(stdout, "\t-P, --pam:\t\tEnables custom pam output for aligned reads. Cannot be set along with --sam.\n");
	fprintf(stdout, "\t-d, --verbose:\t\tVerbose mode: 0 to 1. Higher values output more information (default = 0).\n");
	fprintf(stdout, "\t-a, --scan-lev:\t\tTranscriptome/Genome scan level: 0 to 2. (default = 0)\n\t\t\t\t"
										"0: Report the first mapping.\n\t\t\t\t"
										"1: Continue processing the read unless it is perfectly mapped to cDNA.\n\t\t\t\t"
										"2: Report the best mapping.\n");
	
	fprintf(stdout, "\nOther options:\n");
	fprintf(stdout, "\t-h, --help:\tShows help message.\n");
	fprintf(stdout, "\t-v, --version:\tCurrent version.\n");

	fprintf(stdout, "\nExamples:\n");
	fprintf(stdout, "\tIndexing the reference genome:\n");
	fprintf(stdout, "\t$ ./circminer --index -r ref.fa -k 20\n");
	fprintf(stdout, "\tcircRNA detection of single-end RNA-Seq reads:\n");
	fprintf(stdout, "\t$ ./circminer -r ref.fa -g gene_model.gtf -s reads.fastq -k 20 -o output \n");
	fprintf(stdout, "\tcircRNA detection of paired-end RNA-Seq reads:\n");
	fprintf(stdout, "\t$ ./circminer -r ref.fa -g gene_model.gtf -1 reads_1.fastq -2 reads_2.fastq -k 20 -o output \n");

	fprintf(stdout, "\n");
}
/**********************************************/
