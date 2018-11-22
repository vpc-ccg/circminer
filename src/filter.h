#ifndef __READFILTER_H__
#define __READFILTER_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "fastq_parser.h"
#include "match_read.h"
#include "chain.h"

using namespace std;

class FilterRead {
private:
	FILE* cat_file_pam[CATNUM];
	FILE* temp_fq_r1;
	FILE* temp_fq_r2;

	bool is_pe;
	int cat_count;

	bool first_round;
	bool last_round;

	char comment[400];

public:
	FilterRead (char* save_fname, bool pe, int round, bool first_round, bool last_round, char* fq_file1, char* fq_file2);
	~FilterRead (void);

	int process_read (	Record* current_record, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
						chain_list& forward_best_chain, chain_list& backward_best_chain);

	int process_read (	Record* current_record1, Record* current_record2, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
						chain_list& forward_best_chain_r1, chain_list& backward_best_chain_r1, 
						chain_list& forward_best_chain_r2, chain_list& backward_best_chain_r2);

	void write_read_category (Record* current_record, int is_chimeric);
	void write_read_category (Record* current_record1, Record* current_record2, const MatchedRead& mr);

	void print_mapping (char* rname, const MatchedRead& mr);

};

#endif //__READFILTER_H__
