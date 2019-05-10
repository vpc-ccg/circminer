#include <bitset>

#include "commandline_parser.h"

using namespace std;

bool indexMode = false;
bool compactIndex = false;
bool pairedEnd = false;
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
int threads = 1;
int stage = 0;
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
	int opt;
	int opt_index;

	static struct option long_opt[] = 
	{
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"index", no_argument, 0, 'i'},
		{"compact_index", no_argument, 0, 'm'},
		{"fastq", required_argument, 0, 'f'},
		{"reference", required_argument, 0, 'r'},
		{"gtf", required_argument, 0, 'g'},
		{"pe", no_argument, 0, 'p'},
		{"kmer", required_argument, 0, 'k'},
		{"rlen", required_argument, 0, 'l'},
		{"output", required_argument, 0, 'o'},
		{"verbose", required_argument, 0, 'd'},
		{"thread", required_argument, 0, 't'},
		{"scan_lev", required_argument, 0, 's'},
		{"max_ed", required_argument, 0, 'e'},
		{"max_sc", required_argument, 0, 'c'},
		{"band", required_argument, 0, 'w'},
		{"seed_lim", required_argument, 0, 'S'},
		{"max_tlen", required_argument, 0, 'T'},
		{"max_intron", required_argument, 0, 'I'},
		{"stage", required_argument, 0, 'q'},
		{0,0,0,0},
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "hvimf:r:g:pk:l:o:t:d:s:e:c:w:S:T:I:q:", long_opt, &opt_index )  ) ) 
	{
		switch(opt)
		{
			case 'h': {
				printHELP();
				return 1;
			}
			case 'v': {
				fprintf(stdout, "%s.%s\n", versionNumberMajor, versionNumberMinor);
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
				threads = atoi(optarg);
				if (threads < 1 || threads > sysconf( _SC_NPROCESSORS_ONLN ))
					threads = sysconf( _SC_NPROCESSORS_ONLN );
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
			case '?': {
				fprintf(stderr, "Unknown parameter: %s\n", long_opt[opt_index].name);
				exit(1);
				break;
			}
			default:
				printHELP();
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
	
	THREAD_COUNT = threads;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	return 0;
}

/**********************************************/
void printHELP()
{
	fprintf(stdout, "\nGeneral Options:\n");
	fprintf(stdout, "-h|--help:\tShows help message.\n");
	fprintf(stdout, "-v|--version:\tCurrent version.\n");
	fprintf(stdout, "-f|--fastq:\tRead file (Only specify R1 in paired mode).\n");
	fprintf(stdout, "-r|--refernce:\tReference file.\n");
	fprintf(stdout, "-g|--gtf:\tGene model file.\n");
	
	fprintf(stdout, "\nAdvanced Options:\n");
	fprintf(stdout, "-p|--pe:\t\tPaired end.\n");
	fprintf(stdout, "-k|--kmer:\t\tKmer size [%d..%d] (default = 19).\n", WINDOW_SIZE, WINDOW_SIZE + maxCheckSumLen);
	fprintf(stdout, "-l|--rlen:\t\tMax read length (default = 120).\n");
	fprintf(stdout, "-e|--max_ed:\t\tMax allowed edit distance on each mate (default = %d).\n", EDTH);
	fprintf(stdout, "-c|--max_sc:\t\tMax allowed soft clipping on each mate (default = %d).\n", SOFTCLIPTH);
	fprintf(stdout, "-w|--band:\t\tBand width for banded alignment (default = %d).\n", INDELTH);
	fprintf(stdout, "-S|--seed_lim:\t\tSkip seeds that have more than INT occurrences (default = %d).\n", FRAGLIM);
	fprintf(stdout, "-T|--max_tlen:\t\tMaximum template length of concordant mapping. Paired-end mode only (default = %d).\n", MAXTLEN);
	fprintf(stdout, "-I|--max_intron:\tMaximum length of an intron (default = %d).\n", MAXINTRON);
	fprintf(stdout, "-o|--output:\t\tOutput file (default = output).\n");
	fprintf(stdout, "-t|--thread:\t\tNumber of threads (default = 1).\n");
	fprintf(stdout, "-d|--verbose:\t\tVerbose mode: 0 to 1. Higher values output more information (default = 0).\n");
	fprintf(stdout, "-s|--scan_lev:\t\tTranscriptome/Genome scan level: 0 to 2. (default = 0)\n\t\t\t"
										"0: Report the first mapping.\n\t\t\t"
										"1: Continue processing the read unless it is perfectly mapped to cDNA.\n\t\t\t"
										"2: Report the best mapping.\n");
	
	fprintf(stdout, "\nExample Command:\n");
	fprintf(stdout, "./circRNA -r hg19.fa -f reads_1.fastq -g gene_model.gtf -o output --pe\n");

	fprintf(stdout, "\n");
}
/**********************************************/
