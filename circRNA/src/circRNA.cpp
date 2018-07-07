#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>

#include "common.h"
#include "commandline_parser.h"
#include "fastq_parser.h"
#include "filter.h"
#include "gene_annotation.h"
#include "align.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

char versionNumberMajor[10] = "0";
char versionNumberMinor[10] = "1";

using namespace std;

GTFParser gtf_parser;
Alignment alignment;

int main(int argc, char **argv) {
	time_t pre_time, curr_time;
	double diff_time;
	time(&pre_time);

	
	int exit_c = parse_command( argc, argv );
	if (exit_c == 1)
		return 0;

	char* ref_file = referenceFilename;
	char index_file [FILE_NAME_LENGTH];
	strcpy(index_file, ref_file);
	strcat(index_file, ".index");

	char* fq_file1 = fastqFilename;
	char fq_file2[FILE_NAME_LENGTH];
	bool is_pe = pairedEnd;
	if (pairedEnd) {
		get_mate_name(fq_file1, fq_file2);
	}
	
	/**********************/
	/**LOADING HASH TABLE**/
	int	flag;
	double tmpTime;
	checkSumLength = 5;

	initCommon();
	
	THREAD_COUNT = 8;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	if (!checkHashTable(index_file))
		return 1;

	if (!initLoadingHashTable(index_file))
		return 1;
	
	fprintf(stdout, "Started loading index...\n");
	flag = loadHashTable ( &tmpTime );  			// Reading a fragment

	time(&curr_time);
	diff_time = difftime(curr_time, pre_time);
	pre_time = curr_time;

	fprintf(stdout, "[P] Loaded genome index successfully in %.2f sec\n", diff_time);
	fprintf(stdout, "Winodw size: %d\nChecksum Len: %d\n", WINDOW_SIZE, checkSumLength);

	/***********************/
	
	alignment.init();

	gtf_parser.init(gtfFilename);
	if (! gtf_parser.load_gtf()) {
		fprintf(stderr, "Error in reading GTF file.\n");
		exit(1);
	}
	else 
		fprintf(stdout, "GTF file successfully loaded!\n");

	FASTQParser fq_parser1(fq_file1);
	Record* current_record1;

	FASTQParser fq_parser2;
	Record* current_record2;
	if (is_pe) {
		fq_parser2.init(fq_file2);
	}

	FilterRead filter_read(outputFilename, is_pe);

	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer)) - 1;	// considering both overlapping and non-overlapping kmers

	GIMatchedKmer* fl = (GIMatchedKmer*) malloc(max_seg_cnt * sizeof(GIMatchedKmer));
	GIMatchedKmer* bl = (GIMatchedKmer*) malloc(max_seg_cnt * sizeof(GIMatchedKmer));
	for (int i = 0; i < max_seg_cnt; i++) {
		fl[i].frag_count = 0;
		fl[i].frags = (GeneralIndex*) malloc(FRAGLIM * sizeof(GeneralIndex));
		fl[i].qpos = -1;
		bl[i].frag_count = 0;
		bl[i].frags = (GeneralIndex*) malloc(FRAGLIM * sizeof(GeneralIndex));
		bl[i].qpos = -1;
	}

	chain_list fbc_r1;
	chain_list bbc_r1;
	chain_list fbc_r2;
	chain_list bbc_r2;
	fbc_r1.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	bbc_r1.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	fbc_r2.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	bbc_r2.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	
	for (int i = 0; i < BESTCHAINLIM; i++) {
		fbc_r1.chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
		bbc_r1.chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
		fbc_r2.chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
		bbc_r2.chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
	}

	time(&curr_time);
	diff_time = difftime(curr_time, pre_time);
	fprintf(stdout, "[P] Loaded GTF in %.2f sec\n", diff_time);
	pre_time = curr_time;
	time_t fq_start_t = curr_time;

	fprintf(stdout, "Started reading FASTQ file...\n");
	int line = 0;

	while ( (current_record1 = fq_parser1.get_next()) != NULL ) { // go line by line on fastq file
		if (is_pe)
			current_record2 = fq_parser2.get_next();
		if (current_record1 == NULL)	// no new line
			break;

		line++;
		if (line % LINELOG == 0) {
			time(&curr_time);
			diff_time = difftime(curr_time, pre_time);
			pre_time = curr_time;
			fprintf(stdout, "[P] %d reads in %.2f sec\n", line, diff_time);
		}

		int state;
		if (is_pe) {
			state = filter_read.process_read(current_record1, current_record2, kmer, fl, bl, fbc_r1, bbc_r1, fbc_r2, bbc_r2);
			filter_read.write_read_category(current_record1, current_record2, state);
		}
		else {
			state = filter_read.process_read(current_record1, kmer, fl, bl, fbc_r1, bbc_r1);
			filter_read.write_read_category(current_record1, state);
		}
	}

	time(&curr_time);
	diff_time = difftime(curr_time, fq_start_t);
	fprintf(stdout, "[P] Mapping in %.2f sec\n", diff_time);

	free(fl);
	free(bl);

	return 0;
}
