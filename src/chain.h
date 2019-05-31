#ifndef __CHAIN_H__
#define __CHAIN_H__

#include "common.h"

using namespace std;

typedef struct {
	double score;
	int prev_list;
	int prev_ind;
} chain_cell;

typedef struct {
	chain_cell chain_list[BESTCHAINLIM];
	uint32_t count;
} chain_cell_list;

class Chaining {
private:
	int thid;
	// chain_cell dp[kmer_cnt][max_frag_cnt + 1];
	chain_cell** dp;

	bool compare_frag(fragment_t a, fragment_t b);
	bool check_junction(uint32_t s1, uint32_t s2, const IntervalInfo<UniqSeg>* ol_exons, int kmer, int read_dist, int& trans_dist);
public:
	Chaining(void);
	Chaining(int kmer_cnt, int max_frag_cnt);
	~Chaining(void);

	void init(int kmer_cnt, int max_frag_cnt);

	void chain_seeds_sorted_kbest(int seq_len, GIMatchedKmer* fragment_list, chain_list& best_chain);
	void chain_seeds_sorted_kbest2(int seq_len, GIMatchedKmer* fragment_list, chain_list& best_chain, 
									int kmer, int kmer_cnt, int shift);
};

#endif	//__CHAIN_H__
