#ifndef __MATCHREAD_H__
#define __MATCHREAD_H__

#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <inttypes.h>

#include <string>
#include <vector>
#include <algorithm>

#include "bwt.h"
#include "fastq_parser.h"
#include "common.h"
#include "fragment_list.h"

using namespace std;

//#define GENETHRESH 160000
#define RANGELIM 1000
#define REGIONSIZELIM 1000
#define MRLSIZELIM 20
#define MAXMISSKMER 1

typedef struct {
	bool is_concord;
	int type;
	char* chr;
	bwtint_t start_pos;
	bwtint_t end_pos;
	int matched_len;
	int dir;
	string gene_id;
} MatchedRead;

typedef struct {
	bwtint_t sp;
	bwtint_t ep;
	int matched_len;
	int end_ind;
} ExtendMatch;

typedef struct ExactMatchRes {
	bwtint_t sp;
	bwtint_t ep;
	int q_ind;
	int matched_len;
	int occ;

	bool operator < (const ExactMatchRes& mr) const {
		return occ < mr.occ;
	}
} ExactMatchRes;

typedef struct {
	uint32_t* locs;
	int q_ind;
	int matched_len;
	int occ;
} ExactMatchHash;

//typedef struct {
//	MatchedRead* match_list;
//	int size;
//} MatchedReadList;

void sort_positions(const bwtint_t& sp, const bwtint_t& ep, const int& len, vector<bwtint_t>& forward_list, bwtint_t& flist_size, vector<bwtint_t>& backward_list, bwtint_t& blist_size);
bwtint_t binary_search(const vector<bwtint_t>& list, bwtint_t beg, bwtint_t end, const bwtint_t& target);
void get_mate_name(char* fq1, char* fq2);

int split_match(const char* rseq, int rseq_len, int kmer_size, vector<fragment_t>& forward_fragments, int& forward_fragment_count, vector<fragment_t>& backward_fragments, int& backward_fragment_count);
int split_match_ll(const char* rseq, int rseq_len, int kmer_size, FragmentList& forward_fragments, FragmentList&  backward_fragments);
int split_match_hash(char* rseq, int rseq_len, int kmer_size, FragmentList& forward_fragments);
//int split_match_hash(char* rseq, int rseq_len, int kmer_size, MatchedKmer* start);
int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* start);

void get_reference_chunk_left(uint32_t pos, int len, char* res_str);
void get_reference_chunk_right(uint32_t pos, int len, char* res_str);
void get_reference_chunk(uint32_t pos, int len, char* res_str);
void get_reference_chunk2(uint32_t pos, int len, char* res_str);

void print_location_list(int verbosity, const bwtint_t& sp, const bwtint_t& ep, const int& len);
//void vafprintf(int verbosity, FILE *stream, const char *format, ...);

#endif // __MATCHREAD_H__
