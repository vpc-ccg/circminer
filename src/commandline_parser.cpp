#include "commandline_parser.h"

using namespace std;

bool pairedEnd = false;
int kmer = 19;
int maxReadLength = 76;
int verboseMode = 0;
char gtfFilename[FILE_NAME_LENGTH];
char referenceFilename[FILE_NAME_LENGTH];
char fastqFilename[FILE_NAME_LENGTH];
char outputFilename[FILE_NAME_LENGTH]="output";
char outputDir[FILE_NAME_LENGTH] = "./";

char* contigName;

uint32_t lookup_cnt;

uint8_t* near_border[3];

/**********************************************/
// return:
// 0: normal execution
// 1: exit
int parse_command( int argc, char *argv[] )
{
	int opt;
	int opt_index;
	int gtf_flag = 0;

	static struct option long_opt[] = 
	{
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"fastq", required_argument, 0, 'f'},
		{"reference", required_argument, 0, 'r'},
		{"gtf", required_argument, 0, 'g'},
		{"pe", no_argument, 0, 'p'},
		{"kmer", required_argument, 0, 'k'},
		{"rlen", required_argument, 0, 'l'},
		{"output", required_argument, 0, 'o'},
		{"verbose", required_argument, 0, 'd'},
		{0,0,0,0},
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "hvf:r:g:pk:l:o:d:", long_opt, &opt_index )  ) ) 
	{
		switch(opt)
		{
			case 'h':
				printHELP();
				return 1;
			case 'v':
				fprintf(stdout, "%s.%s\n", versionNumberMajor, versionNumberMinor);
				return 1;
			case 'f':
				strncpy(fastqFilename, optarg, FILE_NAME_LENGTH );
				break;
			case 'r':
				strncpy(referenceFilename, optarg, FILE_NAME_LENGTH);
				break;		
			case 'g':
				strncpy(gtfFilename, optarg, FILE_NAME_LENGTH );
				//gtf_flag = 1;
				break;
			case 'p':
				pairedEnd = true;
				break;
			case 'k':
				kmer = atoi(optarg);
				break;
			case 'l':
				maxReadLength = atoi(optarg);
				break;
			case 'o': {
				strncpy(outputFilename, optarg, FILE_NAME_LENGTH);
				break;
			}
			case 'd':
				verboseMode = atoi(optarg);
				break;
			case '?': 
				fprintf(stderr, "Unknown parameter: %s\n", long_opt[opt_index].name);
				exit(1);
				break;
			default:
				printHELP();
		}
	}
	
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
	fprintf(stdout, "-p|--pe:\tPaired end.\n");
	fprintf(stdout, "-k|--kmer:\tKmer size (default = 19).\n");
	fprintf(stdout, "-l|--rlen:\tMax read length (default = 76).\n");
	fprintf(stdout, "-o|--output:\tOutput file (default = output).\n");
	fprintf(stdout, "-d|--verbose:\tVerbose mode: 0 to 1. Higher values output more information (default = 0).\n");
	
	fprintf(stdout, "\nExample Commands:\n");
	fprintf(stdout, "./circRNA -r hg19.fa -f reads_1.fastq -g gene_model.gtf -o output --pe\n");

	fprintf(stdout, "\n");
}
/**********************************************/
