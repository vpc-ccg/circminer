#ifndef __PROCESSCIRC_H__
#define __PROCESSCIRC_H__

#include <cstdio>
#include "common.h"
#include "hash_table.h"

class ProcessCirc {
private:
	char fq_file1[FILE_NAME_LENGTH];
	char fq_file2[FILE_NAME_LENGTH];

	int window_size;
	int step;
	
	MatchedRead mr;

public:
	ProcessCirc (int last_round_num, int ws);
	~ProcessCirc (void);

	void do_process (void);

	void binning(uint32_t qspos, uint32_t qepos, const RegionalHashTable& regional_ht, char* remain_seq, uint32_t gene_len);
	void chaining(uint32_t qspos, uint32_t qepos, const RegionalHashTable& regional_ht, char* remain_seq, uint32_t gene_len, uint32_t shift);

	int get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos);
};

#endif
