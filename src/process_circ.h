#ifndef __PROCESSCIRC_H__
#define __PROCESSCIRC_H__

#include <cstdio>
#include "common.h"

class ProcessCirc {
private:
	char fq_file1[FILE_NAME_LENGTH];
	char fq_file2[FILE_NAME_LENGTH];

	int window_size;
	MatchedRead mr;

public:
	ProcessCirc (int last_round_num, int ws);
	~ProcessCirc (void);

	void do_process (void);

	int get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos);
};

#endif
