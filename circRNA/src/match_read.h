#ifndef __MATCHREAD_H__
#define __MATCHREAD_H__

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
#include "common.h"
#include "fragment_list.h"

using namespace std;

//#define GENETHRESH 160000
#define RANGELIM 1000
#define REGIONSIZELIM 1000
#define MRLSIZELIM 20
#define MAXMISSKMER 1

typedef struct {
	bool is_concord;
	int type;
	char* chr;
	bwtint_t start_pos;
	bwtint_t end_pos;
	int matched_len;
	int dir;
	string gene_id;
} MatchedRead;

typedef struct {
	bwtint_t sp;
	bwtint_t ep;
	int matched_len;
	int end_ind;
} ExtendMatch;

typedef struct ExactMatchRes{
	bwtint_t sp;
	bwtint_t ep;
	int q_ind;
	int matched_len;
	int occ;

	bool operator < (const ExactMatchRes& mr) const {
		return occ < mr.occ;
	}
} ExactMatchRes;

//typedef struct {
//	MatchedRead* match_list;
//	int size;
//} MatchedReadList;

void sort_positions(const bwtint_t& sp, const bwtint_t& ep, const int& len, vector<bwtint_t>& forward_list, bwtint_t& flist_size, vector<bwtint_t>& backward_list, bwtint_t& blist_size);
bwtint_t binary_search(const vector<bwtint_t>& list, bwtint_t beg, bwtint_t end, const bwtint_t& target);
bool is_concordant_sorted(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh);
bool is_concordant_sorted2(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh, MatchedRead& mr);
bool is_concordant(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh);
int find_expanded_positions(const char* rseq, const char* rcseq, const int& rseq_len);
int find_expanded_positions2(const char* rseq, const char* rcseq, const int& rseq_len, MatchedRead& mr);
int check_concordant_mates(const Record* m1, const Record* m2);
int find_exact_positions(const char* rseq, int rseq_len, int window_size);
void get_mate_name(char* fq1, char* fq2);

int check_concordant_mates_noexpand(const Record* m1, const Record* m2);
int find_exact_positions_slide(const char* rseq, int rseq_len, const int& window_size, const int& shift_step, MatchedRead& mr);

bool is_chimeric_intersect(const vector<bwtint_t>& forwardlist_f, const bwtint_t& flist_size_f, const vector<bwtint_t>& backwardlist_f, const bwtint_t& blist_size_f, const int& len_f,
							const vector<bwtint_t>& forwardlist_b, const bwtint_t& flist_size_b, const vector<bwtint_t>& backwardlist_b, const bwtint_t& blist_size_b, const int& len_b);
int intersect(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, vector <MatchedRead>& mrl, int& mrl_size, bool same_strand);
int find_expanded_sliding_positions(const char* rseq, const char* rcseq, const int& rseq_len, const int& window_size, const int& step, const int& junction_detect_size_lim, vector <MatchedRead>& mrl, int& mrl_size);
int find_expanded_sliding_positions2(const char* rseq, const char* rcseq, const int& rseq_len, const int& window_size, const int& step, const int& junction_detect_size_lim, vector <MatchedRead>& mrl, int& mrl_size, bwtint_t& sp_b, bwtint_t& ep_b, int& exp_len_back);
int check_concordant_mates_expand(const Record* m1, const Record* m2, int kmer_size);

int chop_read_match(const char* rseq, int rseq_len, int kmer_size, int shift, bool recursive, vector<fragment_t>& forward_fragments, int& forward_fragment_count, vector<fragment_t>& backward_fragments, int& backward_fragment_count);
int split_match(const char* rseq, int rseq_len, int kmer_size, vector<fragment_t>& forward_fragments, int& forward_fragment_count, vector<fragment_t>& backward_fragments, int& backward_fragment_count);
int split_match_ll(const char* rseq, int rseq_len, int kmer_size, FragmentList& forward_fragments, FragmentList&  backward_fragments);

void get_reference_chunk(uint32_t pos, int len, char* res_str);
void get_reference_chunk2(uint32_t pos, int len, char* res_str);

void print_location_list(int verbosity, const bwtint_t& sp, const bwtint_t& ep, const int& len);
//void vafprintf(int verbosity, FILE *stream, const char *format, ...);

// kmer analysis
unsigned long long find_occ_sum(const char* rseq, int rseq_len, const int& kmer_size);

#endif // __MATCHREAD_H__
