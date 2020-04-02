#ifndef __READFILTER_H__
#define __READFILTER_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "fastq_parser.h"
#include "match_read.h"
#include "chain.h"
#include "extend.h"
#include "output.h"

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

	TransExtension* extension;
	SAMOutput sam_output;

public:
	FilterRead (void);
	FilterRead (char* save_fname, bool pe, int round, bool first_round, bool last_round, char* fq_file1, char* fq_file2, const vector <ContigLen>& chr_info);
	~FilterRead (void);

	void init (char* save_fname, bool pe, int round, bool first_round, bool last_round, char* fq_file1, char* fq_file2, const vector <ContigLen>& chr_info);
	void finalize (void);

	int process_read (int thid, Record* current_record, int kmer_size, GIMatchedKmer* fl, GIMatchedKmer* bl, 
						chain_list& forward_best_chain, chain_list& backward_best_chain);

	int process_read (int thid, Record* current_record1, Record* current_record2, int kmer_size, GIMatchedKmer* fl, GIMatchedKmer* bl, 
						chain_list& forward_best_chain_r1, chain_list& backward_best_chain_r1, 
						chain_list& forward_best_chain_r2, chain_list& backward_best_chain_r2);

	int process_mates(int thid, const chain_list& forward_chain, const Record* forward_rec, const chain_list& backward_chain, const Record* backward_rec, 
						MatchedRead& mr, bool r1_forward);

	void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer* frag_l, int& high_hits);
	void pair_chains(const chain_list& forward_chain, const chain_list& reverse_chain, vector <MatePair>& mate_pairs, bool* forward_paired, bool* reverse_paired, int saved_type);

	void write_read_category (Record* current_record, int is_chimeric);
	void write_read_category (Record* current_record1, Record* current_record2, const MatchedRead& mr);

	void print_sam (Record* rec);
	void print_pam (Record* rec);
	inline void print_nothing (Record* rec) {}

	void print_sam (Record* rec1, Record* rec2);
	void print_pam (Record* rec1, Record* rec2);
	inline void print_nothing (Record* rec1, Record* rec2) {}

	int get_last_round (void) { return last_round; };

};

typedef void (FilterRead::*filter_print_func)(Record*, Record*);

void* process_block (void* args);

extern FilterRead filter_read;

#endif //__READFILTER_H__
