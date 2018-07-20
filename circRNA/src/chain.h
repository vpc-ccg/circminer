#ifndef __CHAIN_H__
#define __CHAIN_H__

#include "common.h"

using namespace std;

typedef struct {
	double score;
	int prev_list;
	int prev_ind;
} chain_cell;

void chain_seeds_sorted_kbest(int seq_len, GIMatchedKmer*& fragment_list, chain_list& best_chain);
void chain_seeds_sorted_kbest_old(GIMatchedKmer*& fragment_list, chain_list& best_chain);

#endif	//__CHAIN_H__
