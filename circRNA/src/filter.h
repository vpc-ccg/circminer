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

#define CATNUM 4	// number of output categories

// the order matters:
#define CONCRD 0
#define CANDID 1
#define OEANCH 2
#define ORPHAN 3
#define CHIBSJ 4
#define CHIFUS 5

class FilterRead {
private:
	FILE* cat_file_r1[CATNUM];
	FILE* cat_file_r2[CATNUM];

	bool is_pe;

public:
	FilterRead (char* save_fname, bool pe);
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
