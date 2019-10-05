#include <bitset>

#include "commandline_parser.h"

using namespace std;

bool indexMode = false;
bool compactIndex = false;
bool pairedEnd = false;
bool finalCleaning = true;
int kmer = 19;
int maxReadLength = 120;
int verboseMode = 0;
int scanLevel = 0;
int maxEd = EDTH;
int maxSc = SOFTCLIPTH;
int bandWidth = INDELTH;
int seedLim = FRAGLIM;
int maxTlen = MAXTLEN;
int maxIntronLen = MAXINTRON;
int threadCount = 1;
int stage = 2;
int maxCheckSumLen = sizeof(uint16_t) * 8 / 2;	// 2bit per bp

char gtfFilename[FILE_NAME_LENGTH];
char referenceFilename[FILE_NAME_LENGTH];
char fastqFilename[FILE_NAME_LENGTH];
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
		{"fastq", required_argument, 0, 'f'},
		{"reference", required_argument, 0, 'r'},
		{"gtf", required_argument, 0, 'g'},
		{"pe", no_argument, 0, 'p'},
		{"kmer", required_argument, 0, 'k'},
		{"rlen", required_argument, 0, 'l'},
		{"output", required_argument, 0, 'o'},
		{"verbose", required_argument, 0, 'd'},
		{"thread", required_argument, 0, 't'},
		{"scan-lev", required_argument, 0, 's'},
		{"max-ed", required_argument, 0, 'e'},
		{"max-sc", required_argument, 0, 'c'},
		{"band", required_argument, 0, 'w'},
		{"seed-lim", required_argument, 0, 'S'},
		{"max-tlen", required_argument, 0, 'T'},
		{"max-intron", required_argument, 0, 'I'},
		{"stage", required_argument, 0, 'q'},
		{"keep-intermediate", required_argument, 0, 'z'},
		{0,0,0,0},
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "hvimf:r:g:pk:l:o:t:d:s:e:c:w:S:T:I:q:z", long_opt, &opt_index )  ) ) 
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
			case 'f': {
				strncpy(fastqFilename, optarg, FILE_NAME_LENGTH );
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
			case 'p': {
				pairedEnd = true;
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
			case 's': {
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
			case '?': {
				fprintf(stderr, "Unknown parameter: %s\n", long_opt[opt_index].name);
				exit(1);
				break;
			}
			default: {
				printHELPShort();
				return 1;
				}
		}
	}

	if (kmer > WINDOW_SIZE + maxCheckSumLen || kmer < WINDOW_SIZE)
	{
		fprintf(stdout, "ERROR: kmer size should be in [%d..%d]\n", WINDOW_SIZE, WINDOW_SIZE + maxCheckSumLen);
		return 1;
	}

	checkSumLength = (WINDOW_SIZE > kmer) ? 0 : kmer - WINDOW_SIZE;

	if (indexMode)
	{
		CONTIG_SIZE		= DEF_CONTIG_SIZE;
		CONTIG_MAX_SIZE	= DEF_CONTIG_MAX_SIZE;

		// if (referenceFilename == NULL)
		// {
		// 	fprintf(stdout, "ERROR: Reference(s) should be indicated for indexing\n");
		// 	return 0;
		// }
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
	fprintf(stdout, "       circminer -r ref.fa -g gene_model.gtf -f reads.fastq [options]\n");
	fprintf(stderr, "For more details and command line options run \"circminer --help\"\n");
}
/**********************************************/
void printHELP(void)
{
	fprintf(stdout, "\nSYNOPSIS\n");
	fprintf(stdout, "\tcircminer --index FILE [options]\n");
	fprintf(stdout, "\tcircminer -r FILE -g FILE -f FILE [options]\n");
	fprintf(stdout, "\nGeneral options:\n");
	fprintf(stdout, "\t-r, --refernce:\tReference file.\n");
	fprintf(stdout, "\t-g, --gtf:\tGene model file.\n");
	fprintf(stdout, "\t-f, --fastq:\tRead file (Only specify R1 in paired mode).\n");
	
	fprintf(stdout, "\nAdvanced options:\n");
	fprintf(stdout, "\t-m, --compact-index:\tUse this option only while building the index to enable compact version of the index.\n");
	fprintf(stdout, "\t-p, --pe:\t\tPaired end.\n");
	fprintf(stdout, "\t-k, --kmer:\t\tKmer size [%d..%d] (default = 19).\n", WINDOW_SIZE, WINDOW_SIZE + maxCheckSumLen);
	fprintf(stdout, "\t-l, --rlen:\t\tMax read length (default = 120).\n");
	fprintf(stdout, "\t-e, --max-ed:\t\tMax allowed edit distance on each mate (default = %d).\n", EDTH);
	fprintf(stdout, "\t-c, --max-sc:\t\tMax allowed soft clipping on each mate (default = %d).\n", SOFTCLIPTH);
	fprintf(stdout, "\t-w, --band:\t\tBand width for banded alignment (default = %d).\n", INDELTH);
	fprintf(stdout, "\t-S, --seed-lim:\t\tSkip seeds that have more than INT occurrences (default = %d).\n", FRAGLIM);
	fprintf(stdout, "\t-T, --max-tlen:\t\tMaximum template length of concordant mapping. Paired-end mode only (default = %d).\n", MAXTLEN);
	fprintf(stdout, "\t-I, --max-intron:\tMaximum length of an intron (default = %d).\n", MAXINTRON);
	fprintf(stdout, "\t-o, --output:\t\tOutput file (default = output).\n");
	fprintf(stdout, "\t-t, --thread:\t\tNumber of threads (default = 1).\n");
	fprintf(stdout, "\t-d, --verbose:\t\tVerbose mode: 0 to 1. Higher values output more information (default = 0).\n");
	fprintf(stdout, "\t-s, --scan-lev:\t\tTranscriptome/Genome scan level: 0 to 2. (default = 0)\n\t\t\t"
										"0: Report the first mapping.\n\t\t\t"
										"1: Continue processing the read unless it is perfectly mapped to cDNA.\n\t\t\t"
										"2: Report the best mapping.\n");
	
	fprintf(stdout, "\nOther options:\n");
	fprintf(stdout, "\t-h, --help:\tShows help message.\n");
	fprintf(stdout, "\t-v, --version:\tCurrent version.\n");

	fprintf(stdout, "\nExamples:\n");
	fprintf(stdout, "\tIndexing the reference genome:\n");
	fprintf(stdout, "\t$ ./circminer --index -r ref.fa -k 20\n");
	fprintf(stdout, "\tcircRNA detection of paired-end RNA-Seq reads:\n");
	fprintf(stdout, "\t$ ./circminer -r ref.fa -g gene_model.gtf -f reads_1.fastq -k 20 -o output --pe\n");

	fprintf(stdout, "\n");
}
/**********************************************/
