#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <pthread.h>

#include "common.h"
#include "commandline_parser.h"
#include "fastq_parser.h"
#include "filter.h"
#include "gene_annotation.h"
#include "align.h"
#include "process_circ.h"
#include "hash_table.h"
#include "utils.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

char versionNumberMajor[10] = "0";
char versionNumberMinor[10] = "1";

using namespace std;

GTFParser gtf_parser;
FilterRead filter_read;
ScoreMatrix score_mat;
Alignment alignment;

int mapping(int& last_round_num);
void circ_detect(int last_round_num);

int main(int argc, char **argv) {
	int exit_c = parse_command( argc, argv );
	if (exit_c == 1)
		return 0;

	/****************************************************
	 * INDEXING
	 ***************************************************/
	if (indexMode)
	{
		if (compactIndex) {
			if (!generateHashTable(fileName[0], fileName[1]))
				return 1;
		}
		else {
			if (!generateHashTableOnDisk(fileName[0], fileName[1]))
				return 1;	
		}
	}

	/****************************************************
	 * SEARCHING
	 ***************************************************/
	else {
		score_mat.init();
		// 0 <= stage <= 2
		int last_round_num = 1;
		if (stage != 1) {
			int map_ret = mapping(last_round_num);
			if (map_ret == 1)
				return 1;
		}

		if (stage != 0) {
			// dirty turnaround --should be removed later
			if (stage == 1)
				last_round_num = 3;
			circ_detect(last_round_num);
		}
	}

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
	
	// alignment.init();

	/*********************/
	/**Memory Allocation**/
	/*********************/

	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer)) - 1;	// considering both overlapping and non-overlapping kmers

	vector <GIMatchedKmer*> fl(threadCount);
	vector <GIMatchedKmer*> bl(threadCount);

	// vector <chain_list> fbc_r1(threadCount);
	// vector <chain_list> bbc_r1(threadCount);
	// vector <chain_list> fbc_r2(threadCount);
	// vector <chain_list> bbc_r2(threadCount);
	
	vector <chain_list*> fbc_r1(threadCount);
	vector <chain_list*> bbc_r1(threadCount);
	vector <chain_list*> fbc_r2(threadCount);
	vector <chain_list*> bbc_r2(threadCount);

	for (int th = 0; th < threadCount; ++th) {
		fl[th] = (GIMatchedKmer*) malloc(max_seg_cnt * sizeof(GIMatchedKmer));
		bl[th] = (GIMatchedKmer*) malloc(max_seg_cnt * sizeof(GIMatchedKmer));

		fbc_r1[th] = (chain_list*) malloc(1 * sizeof(chain_list));
		bbc_r1[th] = (chain_list*) malloc(1 * sizeof(chain_list));
		fbc_r2[th] = (chain_list*) malloc(1 * sizeof(chain_list));
		bbc_r2[th] = (chain_list*) malloc(1 * sizeof(chain_list));

		fbc_r1[th]->chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
		bbc_r1[th]->chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
		fbc_r2[th]->chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
		bbc_r2[th]->chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
		
		for (int i = 0; i < BESTCHAINLIM; i++) {
			fbc_r1[th]->chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
			bbc_r1[th]->chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
			fbc_r2[th]->chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
			bbc_r2[th]->chains[i].frags = (fragment_t*) malloc(max_seg_cnt * sizeof(fragment_t));
		}
	}

	pthread_t *cm_threads = (pthread_t*) malloc(threadCount * sizeof(pthread_t));

	FilterArgs* filter_args[threadCount];
	for (int th = 0; th < threadCount; ++th)
		filter_args[th] = new FilterArgs(kmer);

	/**********************/
	/**Loading Hash Table**/
	/**********************/

	int	flag;
	double tmpTime;

	if (!checkHashTable(index_file))
		return 1;

	fprintf(stdout, "Load full table? %d\n", loadFullHashTable);

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
	// Record* current_record1;
	// Record* current_record2;

	Record** current_records1;
	Record** current_records2;

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
		contigNum = atoi(contigName) - 1;
		last_round_num = contigNum + 1;

		FASTQParser fq_parser1(fq_file1);

		FASTQParser fq_parser2;
		if (is_pe) {
			fq_parser2.init(fq_file2);
		}

		is_last = !flag;

		if (!is_first) {
			filter_read.finalize();
		}
		filter_read.init(outputFilename, is_pe, contigNum + 1, is_first, is_last, fq_file1, fq_file2);
		is_first = false;

		for (int th = 0; th < threadCount; ++th)
			filter_args[th]->set(fl[th], bl[th], fbc_r1[th], bbc_r1[th], fbc_r2[th], bbc_r2[th]);

		int line = 0;
		int block_size = 0;
		lookup_cnt = 0;

		while ( (current_records1 = fq_parser1.get_next_block()) != NULL ) { // go line by line on fastq file
			if (is_pe)
				current_records2 = fq_parser2.get_next_block();

			if (current_records1 == NULL)	// no new line
				break;

			block_size = fq_parser1.get_block_size();
			if (line % LINELOG == 0) {
				cputime_curr = get_cpu_time();
				realtime_curr = get_real_time();

				fprintf(stdout, "[P] %d reads in %.2lf CPU sec (%.2lf real sec)\t Look ups: %u\n", 
								line, cputime_curr - cputime_start, realtime_curr - realtime_start, lookup_cnt);
				fflush(stdout);
				
				cputime_start = cputime_curr;
				realtime_start = realtime_curr;

				lookup_cnt = 0;
			}

			for (int th = 0; th < threadCount; ++th) {
				filter_args[th]->current_records1 = current_records1;
				filter_args[th]->current_records2 = current_records2;
				filter_args[th]->id = th;
				filter_args[th]->block_size = block_size;
				// pthread_create(cm_threads + th, NULL, process_block, filter_args[th]);
				process_block(filter_args[th]);
			}
			// for (int th = 0; th < threadCount; ++th)
			// 	pthread_join(cm_threads[th], NULL);

			bool skip;
			if (is_pe) {
				for (int i = 0; i < block_size; ++i) {
					skip = (scanLevel == 0 and current_records1[i]->mr->type == CONCRD) or 
							(scanLevel == 1 and current_records1[i]->mr->type == CONCRD and current_records1[i]->mr->gm_compatible and
							(current_records1[i]->mr->ed_r1 + current_records1[i]->mr->ed_r2 == 0) and 
							(current_records1[i]->mr->mlen_r1 + current_records1[i]->mr->mlen_r2 == current_records1[i]->seq_len + current_records2[i]->seq_len));

					if (skip or is_last)
						filter_read.print_mapping(current_records1[i]->rname, *(current_records1[i]->mr));
					if ((!is_last and !skip) or (is_last and (current_records1[i]->mr->type == CHIBSJ or current_records1[i]->mr->type == CHI2BSJ)))
						filter_read.write_read_category(current_records1[i], current_records2[i], *(current_records1[i]->mr));
				}
			}
			else {
				for (int i = 0; i < block_size; ++i) {
					filter_read.write_read_category(current_records1[i], current_records1[i]->mr->type);
				}
			}

			line += block_size;
		}

		cputime_curr = get_cpu_time();
		realtime_curr = get_real_time();

		fprintf(stdout, "[P] Mapping in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

	} while (flag);

	/*************************/
	/**Free Allocated Memory**/
	/*************************/

	filter_read.finalize();

	for (int i = 0; i < contig_cnt; i++) 
		freeMem(orig_contig_len[i].name, strlen(orig_contig_len[i].name));
	freeMem(orig_contig_len, contig_cnt * sizeof(ContigLen));

	finalizeLoadingHashTable();

	for (int th = 0; th < threadCount; ++th) {
		free(fl[th]);
		free(bl[th]);

		for (int i = 0; i < BESTCHAINLIM; i++) {
			free(fbc_r1[th]->chains[i].frags);
			free(bbc_r1[th]->chains[i].frags);
			free(fbc_r2[th]->chains[i].frags);
			free(bbc_r2[th]->chains[i].frags);
		}
		free(fbc_r1[th]->chains);
		free(bbc_r1[th]->chains);
		free(fbc_r2[th]->chains);
		free(bbc_r2[th]->chains);

		free(fbc_r1[th]);
		free(bbc_r1[th]);
		free(fbc_r2[th]);
		free(bbc_r2[th]);		

		delete filter_args[th];
	}

	free(cm_threads);

	return 0;
}

void circ_detect(int last_round_num) {
	int ws = 8;
	ProcessCirc process_circ(last_round_num, ws);
	process_circ.do_process();
}