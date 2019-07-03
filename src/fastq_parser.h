#ifndef __FASTQPARSER_H__
#define __FASTQPARSER_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>

#include "common.h"

#define BLOCKSIZE 1000
#define FQCOMMENTCNT 22

class FASTQParser {
private:
	FILE* input;
	char comp[ASCISIZE];

	// Record* current_record[BLOCKSIZE];
	Record** current_record;
	size_t max_line_size;
	int curr_read;
	int filled_size;
	int block_size;

	char tokens[FQCOMMENTCNT][100];

	bool has_next (void);
	bool read_block (void);

	void set_comp (void);
	void set_reverse_comp (int r_ind);

	int extract_map_info(char* str, int r_ind);
	void fill_map_info(int cnt, int r_ind);

public:
	FASTQParser (void);
	FASTQParser (char* filename);
	~FASTQParser (void);

	void init (char* filename);
	Record* get_next (void);
	Record** get_next_block (void);
	int get_block_size (void);
};

#endif	//__FASTQPARSER_H__
