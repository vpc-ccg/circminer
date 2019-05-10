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

int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* start);

bool pac2char(uint32_t start, int len, char* str);

void get_reference_chunk_left(uint32_t pos, int len, char* res_str);
void get_reference_chunk_right(uint32_t pos, int len, char* res_str);

#endif // __MATCHREAD_H__
