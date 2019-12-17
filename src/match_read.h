#ifndef __MATCHREAD_H__
#define __MATCHREAD_H__

#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "fastq_parser.h"
#include "common.h"

using namespace std;

#define DONOTLOADCONTIG 0
#define LOADCONTIGSTRINMEM 1

class GenomeSeeder {
private:
	int mode;
	char lookup[8];
	char* genome_str;

	int frag_binary_search(const vector<fragment_t>& list, int beg, int end, uint32_t target);
	int get_exact_locs_hash(char* seq, int32_t qpos, uint32_t len, GIMatchedKmer* mk);
	bool reduce_hits_behind(GIMatchedKmer* sl, GIMatchedKmer* ll);
	bool reduce_hits_ahead(GIMatchedKmer* sl, GIMatchedKmer* ll);

	int kmer_match_skip_hash(char* rseq, int rseq_len, int kmer_size, int shift, int skip, int ll_step, GIMatchedKmer* mk_res, int& em_count);

	void print_hits(GIMatchedKmer*, int);

public:
	GenomeSeeder(void);
	~GenomeSeeder(void);
	
	void init(int mode);
	void finalize(void);

	int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* start);
	
	bool pac2char(uint32_t start, int len, char* str);
	bool pac2char_otf(uint32_t start, int len, char* str);
	bool pac2char_whole_contig(void);
	
	void get_reference_chunk_left(uint32_t pos, int len, char* res_str);
	void get_reference_chunk_right(uint32_t pos, int len, char* res_str);
};

extern GenomeSeeder genome_seeder;

#endif // __MATCHREAD_H__
