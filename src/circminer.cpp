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

using namespace std;

char versionNumberMajor[10] = "0";
char versionNumberMinor[10] = "1";

pthread_mutex_t write_lock;
pthread_mutex_t pmap_lock;
pthread_mutex_t read_lock;

GTFParser gtf_parser;
FilterRead filter_read;
ScoreMatrix score_mat;
FASTQParser fq_parser1;
FASTQParser fq_parser2;

int mapping(int& last_round_num);
void circ_detect(int last_round_num);
void* map_reads (void* args);

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
		if (pairedEnd)
			fq_parser1.set_mate(&fq_parser2);

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

	char* fq_file1 = fastqFilename[0];
	char* fq_file2;
	bool is_pe = pairedEnd;
	if (is_pe) {
		fq_file2 = fastqFilename[1];
	}
	
	/*********************/
	/**Memory Allocation**/
	/*********************/

	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer)) - 1;	// considering both overlapping and non-overlapping kmers

	vector <GIMatchedKmer*> fl(threadCount);
	vector <GIMatchedKmer*> bl(threadCount);
	
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

		fbc_r1[th]->chains = (chain_t*) malloc(maxChainLen * sizeof(chain_t));
		bbc_r1[th]->chains = (chain_t*) malloc(maxChainLen * sizeof(chain_t));
		fbc_r2[th]->chains = (chain_t*) malloc(maxChainLen * sizeof(chain_t));
		bbc_r2[th]->chains = (chain_t*) malloc(maxChainLen * sizeof(chain_t));
		
		for (int i = 0; i < maxChainLen; i++) {
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

	fprintf(stdout, "%s mode\n", pairedEnd ? "paired-end" : "single-end");
	fprintf(stdout, "%s\n", loadFullHashTable ? "Load full hash table from index" : "Create hash table on the fly");

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

	int cat_count;
	bool is_first = true;
	bool is_last = false;

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

		fq_parser1.reset(fq_file1);
		if (is_pe) {
			fq_parser2.reset(fq_file2);
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

		for (int th = 0; th < threadCount; ++th) {
			filter_args[th]->id = th;
			pthread_create(cm_threads + th, NULL, map_reads, filter_args[th]);
		}

		for (int th = 0; th < threadCount; ++th) {
			pthread_join(cm_threads[th], NULL);
		}

		cputime_curr = get_cpu_time();
		realtime_curr = get_real_time();

		fprintf(stdout, "[P] Mapping in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

		cputime_start = cputime_curr;
		realtime_start = realtime_curr;

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

		for (int i = 0; i < maxChainLen; i++) {
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

void* map_reads (void* args) {
	FilterArgs* fa = (struct FilterArgs*) args;
	printf("--- thread #%d\n", fa->id);

	filter_print_func print_func_ptr;
	if (reportMapping == SAMFORMAT)
		print_func_ptr = &FilterRead::print_sam;
	else if (reportMapping == PAMFORMAT)
		print_func_ptr = &FilterRead::print_pam;
	else
		print_func_ptr = &FilterRead::print_nothing;
	

	int rid;
	Record* current_record1;
	Record* current_record2;
	int state;
	int is_last = filter_read.get_last_round();
	while ( (rid = fq_parser1.get_next_rec_id()) >= 0 ) { // go line by line on fastq file
		current_record1 = fq_parser1.get_next(rid);
		if (pairedEnd) {
			current_record2 = fq_parser2.get_next(rid);
			if (current_record1 == NULL or current_record2 == NULL)	// no new line
				break;
		
			state = filter_read.process_read(fa->id, current_record1, current_record2, 
					fa->kmer_size, fa->fl, fa->bl, *(fa->fbc_r1), *(fa->bbc_r1), *(fa->fbc_r2), *(fa->bbc_r2));
			bool skip = (scanLevel == 0 and state == CONCRD) or 
						(scanLevel == 1 and state == CONCRD and current_record1->mr->gm_compatible and
						(current_record1->mr->ed_r1 + current_record1->mr->ed_r2 == 0) and 
						(current_record1->mr->mlen_r1 + current_record1->mr->mlen_r2 == current_record1->seq_len + current_record2->seq_len));

			if (skip or is_last) {
				(filter_read.*(print_func_ptr))(current_record1, current_record2);
			}
			if ((!is_last and !skip) or (is_last and (current_record1->mr->type == CHIBSJ or current_record1->mr->type == CHI2BSJ)))
				filter_read.write_read_category(current_record1, current_record2, *(current_record1->mr));
		}
		else {
			state = filter_read.process_read(fa->id, current_record1, fa->kmer_size, fa->fl, fa->bl, 
											 *(fa->fbc_r1), *(fa->bbc_r1));
			filter_read.write_read_category(current_record1, state);
		}
	}

}
