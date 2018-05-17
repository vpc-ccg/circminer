#ifndef __READFILTER_H__
#define __READFILTER_H__

#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <inttypes.h>

#include "bwt.h"
#include "fastq_parser.h"
#include "match_read.h"
#include "chain.h"

using namespace std;

//#define GENETHRESH 100000
//#define RANGELIM 1000
//#define REGIONSIZELIM 2e5

class FilterRead {
private:
	FILE* ignore_r1;
	FILE* keep_r1;
	FILE* chimeric_bsj_r1;
	FILE* chimeric_fusion_r1;
	FILE* partly_unmappable_r1;
	FILE* unmappable_r1;

	FILE* ignore_r2;
	FILE* keep_r2;
	FILE* chimeric_bsj_r2;
	FILE* chimeric_fusion_r2;
	FILE* partly_unmappable_r2;
	FILE* unmappable_r2;

	bool is_pe;

public:
	FilterRead (char* save_fname, bool pe);
	~FilterRead (void);

	int process_read (Record* current_record);
	int process_read (Record* current_record1, Record* current_record2, int kmer_size);
	int process_read_chain (Record* current_record1, Record* current_record2, int kmer_size);

	void write_read (Record* current_record, int is_chimeric);
	void write_read (Record* current_record1, Record* current_record2, int is_chimeric);
	void write_read2 (Record* current_record1, Record* current_record2, int is_chimeric);
	void write_read3 (Record* current_record1, Record* current_record2, int state);

};

#endif //__READFILTER_H__
