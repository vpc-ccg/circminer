#ifndef __FASTQParser_H__
#define __FASTQParser_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>

typedef struct {
	char *rname;
	char *seq;
	char *comment;
	char *qual;

	int seq_len;
} Record;


class FASTQParser {
private:
	FILE *input;
	size_t file_size;

	Record* current_record;
	size_t max_line_size;
	int size;

public:
	FASTQParser (char* filename);
	~FASTQParser (void);

	Record* get_next (void);
	bool has_next (void);
	bool read_next (void);
};

#endif
