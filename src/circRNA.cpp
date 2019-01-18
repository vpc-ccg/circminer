#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "common.h"
#include "commandline_parser.h"
#include "fastq_parser.h"
#include "filter.h"
#include "gene_annotation.h"
#include "align.h"
#include "process_circ.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

char versionNumberMajor[10] = "0";
char versionNumberMinor[10] = "1";

using namespace std;

GTFParser gtf_parser;
Alignment alignment;

int mapping(int& last_round_num);
void circ_detect(int last_round_num);

int main(int argc, char **argv) {
	int exit_c = parse_command( argc, argv );
	if (exit_c == 1)
		return 0;

	int last_round_num = 1;
	int map_ret = mapping(last_round_num);
	if (map_ret == 1)
		return 1;

	//circ_detect(last_round_num);

	return 0;
}

int mapping(int& last_round_num) {
	double cputime_start = get_cpu_time();
	double realtime_start = get_real_time();
	double cputime_curr;
	double realtime_curr;
	
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
	
	alignment.init();

	/*********************/
	/**Memory Allocation**/
	/*********************/

	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer)) - 1;	// considering both overlapping and non-overlapping kmers

	GIMatchedKmer* fl = (GIMatchedKmer*) malloc(max_seg_cnt * sizeof(GIMatchedKmer));
	GIMatchedKmer* bl = (GIMatchedKmer*) malloc(max_seg_cnt * sizeof(GIMatchedKmer));

	for (int i = 0; i < max_seg_cnt; i++) {
		fl[i].junc_dist = (JunctionDist*) malloc(FRAGLIM * sizeof(JunctionDist));
		bl[i].junc_dist = (JunctionDist*) malloc(FRAGLIM * sizeof(JunctionDist));
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

	/**********************/
	/**Loading Hash Table**/
	/**********************/

	int	flag;
	double tmpTime;

	checkSumLength = (WINDOW_SIZE > kmer) ? 0 : kmer - WINDOW_SIZE;

	initCommon();
	
	THREAD_COUNT = 8;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	if (!checkHashTable(index_file))
		return 1;

	ContigLen* orig_contig_len;
	int contig_cnt;
	if (!initLoadingHashTableMeta(index_file, &orig_contig_len, &contig_cnt))
		return 1;

	/*******************/
	/**GTF Parser Init**/
	/*******************/

	gtf_parser.init(gtfFilename, orig_contig_len, contig_cnt);
	if (! gtf_parser.load_gtf()) {
		fprintf(stdout, "Error in reading GTF file.\n");
		exit(1);
	}
	else 
		fprintf(stdout, "GTF file successfully loaded!\n");

	cputime_curr = get_cpu_time();
	realtime_curr = get_real_time();

	fprintf(stdout, "[P] Loaded GTF in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - cputime_start, realtime_curr - realtime_start);

	cputime_start = cputime_curr;
	realtime_start = realtime_curr;

	/*****************/
	/**Mapping Reads**/
	/*****************/

	//char* contig_name;
	int cat_count;
	bool is_first = true;
	bool is_last = false;
	Record* current_record1;
	Record* current_record2;

	do {
		fprintf(stdout, "Started loading index...\n");
	
		flag = loadHashTable ( &tmpTime );  			// Reading a fragment

		cputime_curr = get_cpu_time();
		realtime_curr = get_real_time();

		fprintf(stdout, "[P] Loaded genome index successfully in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - cputime_start, realtime_curr - realtime_start);
		fprintf(stdout, "Winodw size: %d\nChecksum Len: %d\n", WINDOW_SIZE, checkSumLength);

		cputime_start = cputime_curr;
		realtime_start = realtime_curr;

		double fq_cputime_start = cputime_curr;
		double fq_realtime_start = realtime_curr;

		fprintf(stdout, "Contig: %s\n", getRefGenomeName());
		contigName = getRefGenomeName();
		int contigNum = atoi(contigName);
		last_round_num = contigNum;

		FASTQParser fq_parser1(fq_file1);

		FASTQParser fq_parser2;
		if (is_pe) {
			fq_parser2.init(fq_file2);
		}

		is_last = !flag;
		FilterRead filter_read(outputFilename, is_pe, contigNum, is_first, is_last, fq_file1, fq_file2);
		is_first = false;

		int line = 0;
		lookup_cnt = 0;

		while ( (current_record1 = fq_parser1.get_next()) != NULL ) { // go line by line on fastq file
			if (is_pe)
				current_record2 = fq_parser2.get_next();

			if (current_record1 == NULL)	// no new line
				break;

			line++;
			if (line % LINELOG == 0) {
				cputime_curr = get_cpu_time();
				realtime_curr = get_real_time();

				fprintf(stdout, "[P] %d reads in %.2lf CPU sec (%.2lf real sec)\t Look ups: %u\n", line, cputime_curr - cputime_start, realtime_curr - realtime_start, lookup_cnt);
				fflush(stdout);
				
				cputime_start = cputime_curr;
				realtime_start = realtime_curr;

				lookup_cnt = 0;
			}

			int state;
			if (is_pe) {
				state = filter_read.process_read(current_record1, current_record2, kmer, fl, bl, fbc_r1, bbc_r1, fbc_r2, bbc_r2);
				if (current_record1->mr->type == CONCRD or is_last)
					filter_read.print_mapping(current_record1->rname, *(current_record1->mr));
				if ((!is_last and current_record1->mr->type != CONCRD) or (is_last and current_record1->mr->type == CHIBSJ))
					filter_read.write_read_category(current_record1, current_record2, *(current_record1->mr));
			}
			else {
				state = filter_read.process_read(current_record1, kmer, fl, bl, fbc_r1, bbc_r1);
				filter_read.write_read_category(current_record1, state);
			}
		}

		cputime_curr = get_cpu_time();
		realtime_curr = get_real_time();

		fprintf(stdout, "[P] Mapping in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

	} while (flag);

	/*************************/
	/**Free Allocated Memory**/
	/*************************/

	for (int i = 0; i < contig_cnt; i++) 
		freeMem(orig_contig_len[i].name, strlen(orig_contig_len[i].name));
	freeMem(orig_contig_len, contig_cnt * sizeof(ContigLen));

	finalizeLoadingHashTable();

	for (int i = 0; i < 3; i++)
		free(near_border[i]);

	for (int i = 0; i < max_seg_cnt; i++) {
		free(fl[i].junc_dist);
		free(bl[i].junc_dist);
	}
	free(fl);
	free(bl);

	for (int i = 0; i < BESTCHAINLIM; i++) {
		free(fbc_r1.chains[i].frags);
		free(bbc_r1.chains[i].frags);
		free(fbc_r2.chains[i].frags);
		free(bbc_r2.chains[i].frags);
	}
	free(fbc_r1.chains);
	free(bbc_r1.chains);
	free(fbc_r2.chains);
	free(bbc_r2.chains);
}

void circ_detect(int last_round_num) {
	int ws = 8;
	ProcessCirc process_circ(last_round_num, ws);
	process_circ.do_process();
}