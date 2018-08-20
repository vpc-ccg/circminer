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

#define CATNUM 10	// number of output categories

class FilterRead {
private:
	FILE* cat_file_r1[CATNUM];
	FILE* cat_file_r2[CATNUM];

	bool is_pe;
	int cat_count;

public:
	FilterRead (char* save_fname, bool pe, char* filter_temp_name, int num_files, char* fq_file1, char* fq_file2);
	~FilterRead (void);

	int process_read (	Record* current_record, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
						chain_list& forward_best_chain, chain_list& backward_best_chain);

	int process_read (	Record* current_record1, Record* current_record2, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
						chain_list& forward_best_chain_r1, chain_list& backward_best_chain_r1, 
						chain_list& forward_best_chain_r2, chain_list& backward_best_chain_r2);

	void write_read_category (Record* current_record, int is_chimeric);
	void write_read_category (Record* current_record1, Record* current_record2, int state);

};

#endif //__READFILTER_H__
