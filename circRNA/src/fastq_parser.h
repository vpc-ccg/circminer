#ifndef __FASTQPARSER_H__
#define __FASTQPARSER_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>

#include "common.h"

typedef struct {
	char *rname;
	char *seq;
	char *rcseq;
	char *comment;
	char *qual;

	int seq_len;
	int state;
} Record;


class FASTQParser {
private:
	FILE* input;
	char comp[ASCISIZE];

	Record* current_record;
	size_t max_line_size;
	int size;
	bool read_state;

	bool has_next (void);

public:
	FASTQParser (bool read_state);
	FASTQParser (char* filename, bool read_state);
	~FASTQParser (void);

	void init (char* filename);
	Record* get_next (void);
	bool read_next (void);

	void set_reverse_comp (void);

	void set_comp (void);
};

#endif	//__FASTQPARSER_H__
