#ifndef __MATCHREAD_H__
#define __MATCHREAD_H__

#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "fastq_parser.h"
#include "common.h"

using namespace std;

typedef struct {
	bool is_concord;
	int type;
	char* chr;
	uint64_t start_pos;
	uint64_t end_pos;
	int matched_len;
	int dir;
	string gene_id;
} MatchedRead;

typedef struct {
	uint32_t* locs;
	int q_ind;
	int matched_len;
	int occ;
} ExactMatchHash;

void get_mate_name(char* fq1, char* fq2);

int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* start);

void get_reference_chunk_left(uint32_t pos, int len, char* res_str);
void get_reference_chunk_right(uint32_t pos, int len, char* res_str);

#endif // __MATCHREAD_H__
