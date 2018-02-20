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

void sort_positions(const bwtint_t& sp, const bwtint_t& ep, const int& len, vector<bwtint_t>& forward_list, bwtint_t& flist_size, vector<bwtint_t>& backward_list, bwtint_t& blist_size);
bwtint_t binary_search(vector<bwtint_t>& list, bwtint_t beg, bwtint_t end, bwtint_t& target);
bool is_concordant_sorted(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh);
bool is_concordant_sorted2(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh, MatchedRead& mr);
bool is_concordant(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh);
int find_expanded_positions(const char* rseq, const char* rcseq, const int& rseq_len);
int find_expanded_positions2(const char* rseq, const char* rcseq, const int& rseq_len, MatchedRead& mr);
int check_concordant_mates(const Record* m1, const Record* m2);
int find_exact_positions(const char* rseq, int rseq_len, int window_size);
void get_mate_name(char* fq1, char* fq2);

#endif //__READFILTER_H__
