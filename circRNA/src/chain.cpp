#include <cstdio>
#include <cmath>
#include <cstring>

#include "chain.h"
#include "gene_annotation.h"

#define REWARD_COEF		2e4
#define PENALTY_COEF	0.1

#define MAX_INTRON	2000000

inline double score_beta(int distr, int distt, int frag_len) {
	int maxd = distr < distt ? distt : distr;
	int mind = distr < distt ? distr : distt;

	return PENALTY_COEF * (maxd - mind);
}

inline double score_alpha(int distr, int distl, int frag_len) {
	return REWARD_COEF * frag_len;
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
void chain_seeds_sorted_kbest(int seq_len, GIMatchedKmer*& fragment_list, chain_list& best_chain) {
	best_chain.best_chain_count = 0;

	int kmer_cnt = 2 * ceil(1.0 * seq_len / kmer) - 1;
	int max_frag_cnt = FRAGLIM;
	chain_cell dp[kmer_cnt][max_frag_cnt + 1];

	int max_best = BESTCHAINLIM;
	int best_count = 0;
	double best_score = -1;
	chain_cell best_indices[max_best];

	int distr, distt;
	int genome_dist;
	int trans_dist;
	int read_dist;

	double a_score, b_score;
	double temp_score;

	int skipped = 0;
	int i, j;
	int ii, jj;
	int lb_ind[kmer_cnt];
	uint32_t max_lpos_lim;
	uint32_t read_remain;
	uint32_t seg_start;
	uint32_t seg_end;
	uint32_t max_exon_end;
	GIMatchedKmer* cur_mk;
	GIMatchedKmer* pc_mk;	// previously calculated matched kmer

	// Ignore empty fragment list at the back
	while ((kmer_cnt >= 1) and (fragment_list + kmer_cnt - 1)->frag_count <= 0)
		kmer_cnt--;

	if (kmer_cnt <= 0)
		return;

	// Initialize dp array
	for (ii = kmer_cnt - 1; ii >= 0; ii--) {
		cur_mk = fragment_list + ii;
		for (i = 0; i < cur_mk->frag_count; i++) {
			dp[ii][i].score = kmer;
			dp[ii][i].prev_list = -1;
			dp[ii][i].prev_ind = -1;
		}
	}
		
	// Updating stage
	for (ii = kmer_cnt - 2; ii >= 0; ii--) {
		cur_mk = fragment_list + ii;
		read_remain = seq_len - cur_mk->qpos - kmer;
		//memset(lb_ind, 0, kmer_cnt * sizeof(int));
		for (int k = 0; k < kmer_cnt; k++)
			lb_ind[k] = 0;
		
		for (i = 0; i < cur_mk->frag_count; i++) {
			seg_start = cur_mk->frags[i].info;
			seg_end = cur_mk->frags[i].info + kmer - 1;
					
			//max_lpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end);	// = 0 means not found
			max_lpos_lim = -1;

			for (jj = ii+1; jj < kmer_cnt; jj++) {
				pc_mk = fragment_list + jj;
				if (pc_mk->frag_count <= 0 or lb_ind[jj] >= pc_mk->frag_count)	// no fragment left in pc_mk to chain
					continue;

				//if (max_lpos_lim < pc_mk->frags[lb_ind[jj]].info)	// will bot chain to any fragment in thie list
				if (cur_mk->frags[i].info + MAX_INTRON < pc_mk->frags[lb_ind[jj]].info)	// will not chain to any fragment in thie list
					continue;

				while ((lb_ind[jj] < pc_mk->frag_count) and (pc_mk->frags[lb_ind[jj]].info <= cur_mk->frags[i].info))	// skip fragments starting before target
				{
					lb_ind[jj]++;
				}

				if (lb_ind[jj] >= pc_mk->frag_count)	// no more fragment from pre are in range for the remaining of the t
					continue;

				j = lb_ind[jj];

				if (max_lpos_lim == -1)
				{
					gtf_parser.get_upper_bound_alu(seg_start, kmer, read_remain, cur_mk->junc_dist[i]);
					max_lpos_lim = cur_mk->junc_dist[i].range;
					max_exon_end = cur_mk->junc_dist[i].max_end;
					//max_lpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end);	// = 0 means not found
					//vafprintf(2, stderr, "[%d-%d] -> upper bound: %u\n", seg_start, seg_end, max_lpos_lim);
				}

				distr = pc_mk->qpos - cur_mk->qpos - kmer;
				read_dist = distr;
				while ((j < pc_mk->frag_count) and (pc_mk->frags[j].info <= max_lpos_lim)) {
					if (max_exon_end != 0 and (pc_mk->frags[j].info <= max_exon_end) and (pc_mk->frags[j].info + kmer - 1) > max_exon_end) {	// fragment on exon boundary
						j++;
						continue;
					}
					//distt = pc_mk->frags[j].info - cur_mk->frags[i].info;
					gtf_parser.get_upper_bound_alu(pc_mk->frags[j].info, kmer, seq_len - pc_mk->qpos - kmer, pc_mk->junc_dist[j]);
					genome_dist = pc_mk->frags[j].info - seg_end - 1;
					trans_dist = cur_mk->junc_dist[i].dr + pc_mk->junc_dist[j].dl;
					if (abs(genome_dist - read_dist) <= 10) {
						distt = genome_dist;
					}
					else if (abs(trans_dist - read_dist) <= 10) {
						distt = trans_dist;
					}
					else {
						j++;
						continue;
					}

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
	}

		// Finding best score
	for (ii = kmer_cnt - 1; ii >= 0; ii--) {
		cur_mk = fragment_list + ii;
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

	// Back-tracking
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

void chain_seeds_sorted_kbest_old(GIMatchedKmer*& fragment_list, chain_list& best_chain) {
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
	uint32_t max_rpos_lim;
	int pre_start_ind;
	uint32_t read_remain;
	uint32_t seg_start;
	uint32_t seg_end;
	uint32_t ub;
	uint32_t max_exon_end;

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
				read_remain = 76 - cur_mk->qpos - kmer;
				seg_start = cur_mk->frags[i].info;
				seg_end = cur_mk->frags[i].info + kmer - 1;
				
				max_rpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end);	// ub = 0 means not found
				//max_rpos_lim = (ub == 0) ? seg_end + read_remain : ub + read_remain;
				//vafprintf(2, stderr, "[%d-%d] -> upper bound: %u\n", seg_start, seg_end, max_rpos_lim);
				if (max_rpos_lim < pc_mk->frags[pre_start_ind].info)
					continue;

				while ((pre_start_ind < pc_mk->frag_count) and (pc_mk->frags[pre_start_ind].info <= cur_mk->frags[i].info))	// skip out of range
					pre_start_ind++;

				if (pre_start_ind >= pc_mk->frag_count)	// no more fragment from pre are in range for the remaining of the t
					break;

				j = pre_start_ind;

				while ((j < pc_mk->frag_count) and (pc_mk->frags[j].info <= max_rpos_lim)) {
					if (max_exon_end != 0 and (pc_mk->frags[j].info <= max_exon_end) and (pc_mk->frags[j].info + kmer - 1) > max_exon_end) {
						j++;
						continue;
					}
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
