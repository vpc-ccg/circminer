#ifndef __FASTQPARSER_H__
#define __FASTQPARSER_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>

typedef struct {
	char *rname;
	char *seq;
	char *rcseq;
	char *comment;
	char *qual;

	int seq_len;
} Record;


class FASTQParser {
private:
	FILE* input;
	size_t file_size;

	Record* current_record;
	size_t max_line_size;
	int size;

	short comp[30];

public:
	FASTQParser ();
	FASTQParser (char* filename);
	~FASTQParser (void);

	void init (char* filename);
	Record* get_next (void);
	bool has_next (void);
	bool read_next (void);

	void set_reverse_comp (void);

	void set_comp (void);
	char get_comp (char nt);
};

#endif	//__FASTQPARSER_H__
