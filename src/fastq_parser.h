#ifndef __FASTQPARSER_H__
#define __FASTQPARSER_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>

#include "common.h"

#define FQCOMMENTCNT 22

struct Record {
	char* rname;
	char* seq;
	char* rcseq;
	char* comment;
	char* qual;

	int seq_len;

	MatchedRead* mr;

	Record(void) {
		mr = new MatchedRead;
	}

	~Record(void) {
		delete mr;
	}
};


class FASTQParser {
private:
	FILE* input;
	char comp[ASCISIZE];

	Record* current_record;
	size_t max_line_size;
	int size;

	char tokens[FQCOMMENTCNT][100];

	bool has_next (void);

public:
	FASTQParser (void);
	FASTQParser (char* filename);
	~FASTQParser (void);

	void init (char* filename);
	Record* get_next (void);
	bool read_next (void);

	void set_reverse_comp (void);

	void set_comp (void);

	int extract_map_info(char* str);
	void fill_map_info(int cnt);
};

#endif	//__FASTQPARSER_H__
