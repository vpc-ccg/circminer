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

void chain_seeds_sorted_kbest(int seq_len, GIMatchedKmer *fragment_list, chain_list &best_chain);
void chain_seeds_sorted_kbest2(int seq_len, GIMatchedKmer *fragment_list, chain_list &best_chain,
                               int kmer, int kmer_cnt, int shift);

#endif    //__CHAIN_H__
