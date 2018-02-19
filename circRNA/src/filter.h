#ifndef __READFILTER_H__
#define __READFILTER_H__

#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <inttypes.h>

#include <string>
#include <vector>
#include <algorithm>

#include "bwt.h"
#include "fastq_parser.h"

using namespace std;

#define GENETHRESH 100000
#define RANGELIM 1000
#define REGIONSIZELIM 2e5

typedef struct {
	bool is_concord;
	bwtint_t start_pos;
	int matched_len;
	int dir;
} MatchedRead;

class FilterRead {
private:
	FILE* ignore_r1;
	FILE* keep_r1;

	FILE* ignore_r2;
	FILE* keep_r2;

	bool is_pe;

public:
	FilterRead (char* save_fname, bool pe);
	~FilterRead (void);

	void write_read (Record* current_record, bool is_chimeric);
	void write_read (Record* current_record1, Record* current_record2, bool is_chimeric);

};

#endif //__READFILTER_H__
