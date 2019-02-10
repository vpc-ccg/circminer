#ifndef __PROCESSCIRC_H__
#define __PROCESSCIRC_H__

#include <cstdio>
#include "common.h"
#include "hash_table.h"
#include "fastq_parser.h"

class ProcessCirc {
private:
	char fq_file1[FILE_NAME_LENGTH];
	char fq_file2[FILE_NAME_LENGTH];

	int window_size;
	int step;
	
	MatchedRead mr;
	chain_list bc;

public:
	ProcessCirc (int last_round_num, int ws);
	~ProcessCirc (void);

	void sort_fq(char* fqname);

	void do_process (void);

	void call_circ(Record* current_record1, Record* current_record2);

	void binning(uint32_t qspos, uint32_t qepos, const RegionalHashTable& regional_ht, char* remain_seq, uint32_t gene_len);
	void chaining(uint32_t qspos, uint32_t qepos, const RegionalHashTable& regional_ht, char* remain_seq, uint32_t gene_len, uint32_t shift, uint32_t& rspos, uint32_t& repos);

	bool find_exact_coord(MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, 
							int dir, uint32_t qspos, char* rseq, int rlen, int whole_len);

	int get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos);


};

#endif
