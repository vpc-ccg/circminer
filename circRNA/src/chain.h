#ifndef __CHAIN_H__
#define __CHAIN_H__

#include <vector>
#include "common.h"
#include "fragment_list.h"

using namespace std;

typedef struct {
	double score;
	int prev_list;
	int prev_ind;
} chain_cell;

void chain_seeds_n2(vector<fragment_t>& fragment_list, uint32_t fragment_count, chain_t &best_chain);
void chain_seeds_n2_kbest(vector<fragment_t>& fragment_list, uint32_t fragment_count, vector <chain_t>& best_chain);

void chain_seeds_n2_kbest(FragmentList& fragment_list, vector <chain_t>& best_chain);
//void chain_seeds_sorted_kbest(FragmentList& fragment_list, vector <chain_t>& best_chain);
//void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, vector <chain_t>& best_chain);
void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, chain_list& best_chain);

#endif
