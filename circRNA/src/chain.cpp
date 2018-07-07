#include <cstdio>

#include "chain.h"

#define INF 100000

inline double score_beta(int distr, int distt, int frag_len) {
	int maxd = distr < distt ? distt : distr;
	int mind = distr < distt ? distr : distt;

	return 0.15 * (maxd - mind);
}

inline double score_alpha(int distr, int distl, int frag_len) {
	return 7.0 * frag_len;
}

bool compare_frag(fragment_t a, fragment_t b) {
	return a.qpos < b.qpos;
}

// Assumption: fragment list sorted by rpos
//
// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, chain_list& best_chain) {
	int i, j;

	int distr, distt;
	double a_score, b_score;
	int kmer_cnt = 7;
	int max_frag_cnt = FRAGLIM;
	chain_cell dp[kmer_cnt][max_frag_cnt + 1];

	double best_score = -1;
	int max_best = BESTCHAINLIM;
	chain_cell best_indices[max_best];
	int best_count = 0;
	double temp_score;

	GIMatchedKmer* cur_mk;
	GIMatchedKmer* pc_mk;	// previously calculated matched kmer
	int ii;
	int jj;
	int max_rpos_lim;
	int pre_start_ind;
	for (ii = kmer_cnt - 1; ii >= 0; ii--) {
		cur_mk = fragment_list + ii;
		for (i = 0; i < cur_mk->frag_count; i++) {
			dp[ii][i].score = kmer;
			dp[ii][i].prev_list = -1;
			dp[ii][i].prev_ind = -1;
		}

		for (jj = ii+1; jj < kmer_cnt; jj++) {
			pc_mk = fragment_list + jj;
			if (pc_mk->frag_count <= 0)  
				continue;

			pre_start_ind = 0;
			for (i = 0; i < cur_mk->frag_count; i++) {
				max_rpos_lim = cur_mk->frags[i].info + GENETHRESH;
				if (max_rpos_lim < pc_mk->frags[pre_start_ind].info)
					continue;

				while ((pre_start_ind < pc_mk->frag_count) and (pc_mk->frags[pre_start_ind].info <= cur_mk->frags[i].info))	// skip out of range
					pre_start_ind++;

				if (pre_start_ind >= pc_mk->frag_count)	// no more fragment from pre are in range for the remaining of the t
					break;

				j = pre_start_ind;

				while ((j < pc_mk->frag_count) and (pc_mk->frags[j].info <= max_rpos_lim)) {
					distr = pc_mk->qpos - cur_mk->qpos;

					distt = pc_mk->frags[j].info - cur_mk->frags[i].info;

					//a_score = score_alpha(distr, distt, t->frags[i].len);
					//b_score = score_beta (distr, distt, t->frags[i].len);

					temp_score = dp[jj][j].score + score_alpha(distr, distt, kmer) - score_beta(distr, distt, kmer);
					if (temp_score > dp[ii][i].score) {
						dp[ii][i].score = temp_score;
						dp[ii][i].prev_list = jj;
						dp[ii][i].prev_ind = j;
					}
					j++;
				}
			}
		}

		for (i = 0; i < cur_mk->frag_count; i++) {
			if (dp[ii][i].score > best_score) {
				best_score = dp[ii][i].score;
				best_count = 1;
				best_indices[0].prev_list = ii;
				best_indices[0].prev_ind = i;
			}
			else if (dp[ii][i].score == best_score and best_count < max_best) {
				best_indices[best_count].prev_list = ii;
				best_indices[best_count].prev_ind = i;
				best_count++;
			}
		}
	}

	// back-tracking
	chain_cell best_index;	
	best_chain.best_chain_count = best_count;
	int tmp_list_ind;
	for (j = 0; j < best_count; j++) {
		best_index = best_indices[j];
		i = 0;
		while (best_index.prev_list != -1) {
			best_chain.chains[j].frags[i].rpos = (fragment_list + best_index.prev_list)->frags[best_index.prev_ind].info;
			best_chain.chains[j].frags[i].qpos = (fragment_list + best_index.prev_list)->qpos;
			best_chain.chains[j].frags[i].len = kmer;
			
			tmp_list_ind = best_index.prev_list;
			best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
			best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;
			
			i++;
		}

		best_chain.chains[j].score = best_score;
		best_chain.chains[j].chain_len = i;

	}
}
