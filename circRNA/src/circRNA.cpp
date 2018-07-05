#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include "common.h"
#include "commandline_parser.h"
#include "fastq_parser.h"
#include "filter.h"
#include "gene_annotation.h"
#include "align.h"

extern "C" {
#include "mrsfast/Common.h"
//#include "mrsfast/CommandLineParser.h"
//#include "mrsfast/Reads.h"
//#include "mrsfast/Output.h"
#include "mrsfast/HashTable.h"
//#include "mrsfast/MrsFAST.h"
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
	//sprintf(index_file, "%s.index", ref_file);
	strcpy(index_file, ref_file);
	strcat(index_file, ".index");

	char* fq_file1 = fastqFilename;
	char fq_file2[FILE_NAME_LENGTH];
	bool is_pe = pairedEnd;
	if (pairedEnd) {
		get_mate_name(fq_file1, fq_file2);
	}
	
	/*
	if (bwt_load(ref_file)) {
		fprintf(stdout, "Index not found!\n");
		//if (bwt_index(ref_file)) {
		//	fprintf(stderr, "Indexing failed!\n");
		//	return 1;
		//}
		//else
		//	fprintf(stderr, "Indexed successfully!\n");
	}
	else
		fprintf(stdout, "Index file successfully loaded!\n");
	
	*/

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
	
	gtf_parser.init(gtfFilename);
	alignment.init();
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

	//FragmentList fl(7);
	//FragmentList bl(7);
	GIMatchedKmer* fl = (GIMatchedKmer*) malloc(7 * sizeof(GIMatchedKmer));
	GIMatchedKmer* bl = (GIMatchedKmer*) malloc(7 * sizeof(GIMatchedKmer));
	for (int i = 0; i < 7; i++) {
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
	int max_frag_count = 7 + 1;
	for (int i = 0; i < BESTCHAINLIM; i++) {
		fbc_r1.chains[i].frags = (fragment_t*) malloc(max_frag_count * sizeof(fragment_t));
		bbc_r1.chains[i].frags = (fragment_t*) malloc(max_frag_count * sizeof(fragment_t));
		fbc_r2.chains[i].frags = (fragment_t*) malloc(max_frag_count * sizeof(fragment_t));
		bbc_r2.chains[i].frags = (fragment_t*) malloc(max_frag_count * sizeof(fragment_t));
	}

	time(&curr_time);
	diff_time = difftime(curr_time, pre_time);
	fprintf(stdout, "[P] Loaded GTF in %.2f sec\n", diff_time);
	pre_time = curr_time;
	time_t fq_start_t = curr_time;

	fprintf(stdout, "Started reading FASTQ file...\n");
	int line = 0;

	while ( fq_parser1.has_next() ) { // go line by line on fastq file
		current_record1 = fq_parser1.get_next();
		if (is_pe)
			current_record2 = fq_parser2.get_next();
		if (current_record1 == NULL)	// no new line
			break;

		line++;
		if (line % 100000 == 0) {
			time(&curr_time);
			diff_time = difftime(curr_time, pre_time);
			pre_time = curr_time;
			fprintf(stdout, "[P] %d reads in %.2f sec\n", line, diff_time);
		}

		//fprintf(stderr, "Seq: %s, len: %d\n", current_record1->seq, current_record1->seq_len);
		//fprintf(stderr, "RC Seq: %s, len: %d\n", current_record->rcseq, current_record->seq_len);
		//int occ = find_exact_positions(current_record->seq, current_record->seq_len, find_len);

		int is_chimeric;
		if (is_pe) {
			//unsigned long long occ1 = find_occ_sum(current_record1->seq, current_record1->seq_len, kmer);
			//unsigned long long occ2 = find_occ_sum(current_record2->seq, current_record2->seq_len, kmer);
			//continue;

			//is_chimeric = check_concordant_mates(current_record1, current_record2);
			//is_chimeric = filter_read.process_read(current_record1, current_record2, kmer);
			
			//is_chimeric = filter_read.process_read_chain(current_record1, current_record2, kmer);
			is_chimeric = filter_read.process_read_chain_hash(current_record1, current_record2, kmer, fl, bl, fbc_r1, bbc_r1, fbc_r2, bbc_r2);
			filter_read.write_read3(current_record1, current_record2, is_chimeric);
		}
		else {
			is_chimeric = filter_read.process_read(current_record1);
			filter_read.write_read(current_record1, is_chimeric);
		}
	}

	time(&curr_time);
	diff_time = difftime(curr_time, fq_start_t);
	fprintf(stdout, "[P] Mapping in %.2f sec\n", diff_time);

	free(fl);
	free(bl);

	return 0;
}
