#ifndef __PROCESSCIRC_H__
#define __PROCESSCIRC_H__

#include <cstdio>
#include "common.h"

class ProcessCirc {
private:
	FILE* circ_file;

	int window_size;
	size_t max_line_size;
	
	char* line;
	char tokens[12][100];

	MatchedRead mr;

public:
	ProcessCirc (char* fname, int ws);
	~ProcessCirc (void);

	bool read_next (void);
	void tokenize (void);
	void do_process (void);
};

#endif
