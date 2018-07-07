#ifndef __CHAIN_H__
#define __CHAIN_H__

#include <vector>
#include "common.h"

using namespace std;

typedef struct {
	double score;
	int prev_list;
	int prev_ind;
} chain_cell;

void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, chain_list& best_chain);

#endif	//__CHAIN_H__
