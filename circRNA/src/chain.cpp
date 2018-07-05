#include <algorithm>
#include <deque>
#include <cstdio>

#include "chain.h"
#include "fragment_list.h"

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

// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
void chain_seeds_n2(vector<fragment_t>& fragment_list, uint32_t fragment_count, chain_t &best_chain) {
	int i, j;

	int distr, distt;
	double a_score, b_score;
	double dp[fragment_count + 1];
	int prev [fragment_count + 1];

	double best_score = -1;
	int best_index = -1;

	for (i = 0; i < fragment_count; i++) {
		//fprintf(stderr, "qpos: %d\trpos: %zu\tlen: %zu\n", fragment_list[i].qpos, fragment_list[i].rpos, fragment_list[i].len);
	}
	sort(fragment_list.begin(), fragment_list.begin() + fragment_count, compare_frag);

	for (i = 0; i < fragment_count; i++) {
		dp[i] = fragment_list[i].len;
		//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
		//fprintf(stderr, "Frag Len: %f\n", dp[i]);
		prev[i] = -1;
		
		for (j = i - 1; j >= 0; j--) {
			distr = fragment_list[i].qpos - (fragment_list[j].qpos + fragment_list[j].len - 1);
			if (distr <= 0)
				continue;

			distt = fragment_list[i].rpos - (fragment_list[j].rpos + fragment_list[j].len - 1);
			if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
				continue;

			a_score = score_alpha(distr, distt, fragment_list[i].len);
			b_score = score_beta (distr, distt, fragment_list[i].len);

			if (dp[j] + a_score - b_score > dp[i]) {
				dp[i] = dp[j] + a_score - b_score;
				prev[i] = j;
			}
		}

		if (dp[i] > best_score) {
			//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[i], i);
			best_score = dp[i];
			best_index = i;
		}
	}

	// back-tracking
	//
	deque<int> chain_index;

	while (best_index != -1) {
		chain_index.push_front(best_index);
		best_index = prev[best_index];
	}

	best_chain.score = best_score;
	best_chain.chain_len = chain_index.size();

	for (i = 0; i < chain_index.size(); i++) {
		best_chain.frags[i] = fragment_list[chain_index[i]];
	}
}

// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
void chain_seeds_n2_kbest(vector<fragment_t>& fragment_list, uint32_t fragment_count, vector <chain_t>& best_chain) {
	int i, j;

	int distr, distt;
	double a_score, b_score;
	double dp[fragment_count + 1];
	int prev [fragment_count + 1];

	double best_score = -1;
	int max_best = best_chain.size();
	vector <int> best_indices(max_best);
	int best_count = 0;

	//sort(fragment_list.begin(), fragment_list.begin() + fragment_count, compare_frag);

	for (i = 0; i < fragment_count; i++) {
		dp[i] = fragment_list[i].len;
		//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
		//fprintf(stderr, "Frag Len: %f\n", dp[i]);
		prev[i] = -1;
		
		for (j = i - 1; j >= 0; j--) {
			distr = fragment_list[i].qpos - (fragment_list[j].qpos + fragment_list[j].len - 1);
			if (distr <= 0)
				continue;

			distt = fragment_list[i].rpos - (fragment_list[j].rpos + fragment_list[j].len - 1);
			if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
				continue;

			a_score = score_alpha(distr, distt, fragment_list[i].len);
			b_score = score_beta (distr, distt, fragment_list[i].len);

			if (dp[j] + a_score - b_score > dp[i]) {
				dp[i] = dp[j] + a_score - b_score;
				prev[i] = j;
			}
		}

		if (dp[i] > best_score) {
			//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[i], i);
			best_score = dp[i];
			best_count = 1;
			best_indices[0] = i;
		}
		else if (dp[i] == best_score and best_count < max_best) {
			best_indices[best_count] = i;
			best_count++;
		}
	}

	// back-tracking
	//
	deque<int> chain_index;
	int best_index;	
	best_chain.resize(best_count);
	for (int j = 0; j < best_count; j++) {
		chain_index.clear();
		best_index = best_indices[j];
		while (best_index != -1) {
			chain_index.push_front(best_index);
			best_index = prev[best_index];
		}

		best_chain[j].score = best_score;
		best_chain[j].chain_len = chain_index.size();

		for (i = 0; i < chain_index.size(); i++) {
			best_chain[j].frags[i] = fragment_list[chain_index[i]];
		}
	}
}

// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
void chain_seeds_n2_kbest(FragmentList& fragment_list, vector <chain_t>& best_chain) {
	int i, j;

	int distr, distt;
	double a_score, b_score;
	chain_cell dp[fragment_list.get_size()][fragment_list.get_max_frag_size() + 1];

	double best_score = -1;
	int max_best = best_chain.size();
	vector <chain_cell> best_indices(max_best);
	int best_count = 0;

	//sort(fragment_list.begin(), fragment_list.begin() + fragment_count, compare_frag);

	MatchedKmer* t;
	MatchedKmer* pre;
	MatchedKmer* frag_tmp[fragment_list.get_size()];
	int ii = 0;
	int jj = 0;
	for (t = fragment_list.get_head(); t != NULL; t = t->next) {
		frag_tmp[ii] = t;
		for (i = 0; i < t->frag_count; i++) {
			dp[ii][i].score = t->frags[i].len;
			//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
			//fprintf(stderr, "Frag Len: %f\n", dp[i]);
			dp[ii][i].prev_list = -1;
			dp[ii][i].prev_ind = -1;

			jj = ii-1;
			for (pre = t->prev; pre != NULL; pre = pre->prev) {
				for (j = pre->frag_count-1; j >= 0; j--) {
					distr = t->frags[i].qpos - pre->frags[j].qpos;								// also consider overlapping
					//distr = t->frags[i].qpos - (pre->frags[j].qpos + pre->frags[j].len - 1);	// for non-overlapping
					if (distr <= 0)
						continue;

					distt = t->frags[i].rpos - pre->frags[j].rpos;								// also consider overlapping
					//distt = t->frags[i].rpos - (pre->frags[j].rpos + pre->frags[j].len - 1);	// for non-overlapping
					if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
						continue;

					a_score = score_alpha(distr, distt, t->frags[i].len);
					b_score = score_beta (distr, distt, t->frags[i].len);

					if (dp[jj][j].score + a_score - b_score > dp[ii][i].score) {
						dp[ii][i].score = dp[jj][j].score + a_score - b_score;
						dp[ii][i].prev_list = jj;
						dp[ii][i].prev_ind = j;
					}
				}
				jj--;
			}

			if (dp[ii][i].score > best_score) {
				//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[i], i);
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
		ii++;
	}

	// back-tracking
	deque<chain_cell> chain_index;
	chain_cell best_index;	
	best_chain.resize(best_count);
	int tmp_list_ind;
	for (int j = 0; j < best_count; j++) {
		chain_index.clear();
		best_index = best_indices[j];
		while (best_index.prev_list != -1) {
			chain_index.push_front(best_index);
			tmp_list_ind = best_index.prev_list;
			best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
			best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;
		}

		best_chain[j].score = best_score;
		best_chain[j].chain_len = chain_index.size();

		for (i = 0; i < chain_index.size(); i++) {
			best_chain[j].frags[i] = frag_tmp[chain_index[i].prev_list]->frags[chain_index[i].prev_ind];
		}
	}
}

// Assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
int frag_binary_search(const fragment_t* list, int beg, int end, uint32_t target) {
	if (end - beg <= 1)
		return end;
	//fprintf(stderr, "In Binary Search: (%d, %d) looking for %lu\n", beg, end, target);
	int mid = (beg + end) / 2;
	if (target < list[mid].rpos)
		return frag_binary_search(list, beg, mid, target);
	else 
		return frag_binary_search(list, mid, end, target);
}

// Assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
int frag_binary_search(const GeneralIndex* list, int beg, int end, uint32_t target) {
	if (end - beg <= 1)
		return end;
	//fprintf(stderr, "In Binary Search: (%d, %d) looking for %lu\n", beg, end, target);
	int mid = (beg + end) / 2;
	if (target < list[mid].info)
		return frag_binary_search(list, beg, mid, target);
	else 
		return frag_binary_search(list, mid, end, target);
}

// Assumption: fragment list sorted by rpos
//
// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
void chain_seeds_sorted_kbest(FragmentList& fragment_list, vector <chain_t>& best_chain) {
	int i, j;

	int distr, distt;
	double a_score, b_score;
	chain_cell dp[fragment_list.get_size()][fragment_list.get_max_frag_size() + 1];

	double best_score = -1;
	int max_best = best_chain.size();
	vector <chain_cell> best_indices(max_best);
	int best_count = 0;
	double temp_score;

	MatchedKmer* t;
	MatchedKmer* pre;
	MatchedKmer* frag_tmp[fragment_list.get_size()];
	int ii = 0;
	int jj = 0;
	int max_rpos_lim;
	for (t = fragment_list.get_head(); t != NULL; t = t->next) {
		frag_tmp[ii] = t;
		for (i = 0; i < t->frag_count; i++) {
			dp[ii][i].score = t->frags[i].len;
			//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
			//fprintf(stderr, "Frag Len: %f\n", dp[i]);
			dp[ii][i].prev_list = -1;
			dp[ii][i].prev_ind = -1;
		}

		jj = ii;
		for (pre = t->prev; pre != NULL; pre = pre->prev) {
			jj--;
			if (pre->frag_count <= 0)  
				continue;

			for (i = 0; i < t->frag_count; i++) {
				if ((t->frags[i].rpos - 1) < pre->frags[0].rpos)
					continue;

				max_rpos_lim = t->frags[i].rpos - GENETHRESH;
				//j = frag_binary_search(pre->frags, 0, pre->frag_count, t->frags[i].rpos - GENETHRESH - 1);
				j = frag_binary_search(pre->frags, 0, pre->frag_count, t->frags[i].rpos - 1) - 1;
				//fprintf(stderr, "found index on flist[%d]: %d\n", jj, j);
				while ((j >= 0) and (pre->frags[j].rpos >= max_rpos_lim)) {
					//fprintf(stderr, "working on t[%d], pre[%d]\n", i, j);
					distr = t->frags[i].qpos - pre->frags[j].qpos;
					//if (distr <= 0)
					//	continue;

					distt = t->frags[i].rpos - pre->frags[j].rpos;
					//if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
					//	continue;

					//a_score = score_alpha(distr, distt, t->frags[i].len);
					//b_score = score_beta (distr, distt, t->frags[i].len);

					temp_score = dp[jj][j].score + score_alpha(distr, distt, t->frags[i].len) - score_beta(distr, distt, t->frags[i].len);
					if (temp_score > dp[ii][i].score) {
						dp[ii][i].score = temp_score;
						dp[ii][i].prev_list = jj;
						dp[ii][i].prev_ind = j;
					}
					j--;
				}
			}
		}

		for (i = 0; i < t->frag_count; i++) {
			if (dp[ii][i].score > best_score) {
				//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[i], i);
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
		ii++;
	}

	// back-tracking
	deque<chain_cell> chain_index;
	chain_cell best_index;	
	best_chain.resize(best_count);
	int tmp_list_ind;
	for (int j = 0; j < best_count; j++) {
		chain_index.clear();
		best_index = best_indices[j];
		while (best_index.prev_list != -1) {
			chain_index.push_front(best_index);
			tmp_list_ind = best_index.prev_list;
			best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
			best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;
		}

		best_chain[j].score = best_score;
		best_chain[j].chain_len = chain_index.size();

		for (i = 0; i < chain_index.size(); i++) {
			best_chain[j].frags[i] = frag_tmp[chain_index[i].prev_list]->frags[chain_index[i].prev_ind];
		}
	}
}

// Assumption: fragment list sorted by rpos
//
// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
//void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, chain_list& best_chain) {
//	int i, j;
//
//	int distr, distt;
//	double a_score, b_score;
//	int kmer_cnt = 7;
//	int max_frag_cnt = FRAGLIM;
//	chain_cell dp[kmer_cnt][max_frag_cnt + 1];
//
//	double best_score = -1;
//	int max_best = BESTCHAINLIM;
//	vector <chain_cell> best_indices(max_best);
//	int best_count = 0;
//	double temp_score;
//
//	GIMatchedKmer* t;
//	GIMatchedKmer* pre;
//	int ii = 0;
//	int jj = 0;
//	int max_rpos_lim;
//	int pre_start_ind;
//	for (ii = 0; ii < kmer_cnt; ii++) {
//		t = fragment_list + ii;
//		for (i = 0; i < t->frag_count; i++) {
//			dp[ii][i].score = kmer;
//			//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
//			//fprintf(stderr, "Frag Len: %f\n", dp[i]);
//			dp[ii][i].prev_list = -1;
//			dp[ii][i].prev_ind = -1;
//		}
//
//		for (jj = ii-1; jj >= 0; jj--) {
//			pre = fragment_list + jj;
//			if (pre->frag_count <= 0)  
//				continue;
//
//			pre_start_ind = 0;
//			for (i = 0; i < t->frag_count; i++) {
//				if ((t->frags[i].info - 1) < pre->frags[pre_start_ind].info)
//					continue;
//
//				max_rpos_lim = t->frags[i].info - GENETHRESH;
//				//j = frag_binary_search(pre->frags, 0, pre->frag_count, t->frags[i].rpos - GENETHRESH - 1);
//				j = frag_binary_search(pre->frags, pre_start_ind, pre->frag_count, t->frags[i].info - 1) - 1;
//				//fprintf(stderr, "found index on flist[%d]: %d\n", jj, j);
//				while ((j >= pre_start_ind) and (pre->frags[j].info >= max_rpos_lim)) {
//					//fprintf(stderr, "working on t[%d], pre[%d]\n", i, j);
//					distr = t->qpos - pre->qpos;
//					//if (distr <= 0)
//					//	continue;
//
//					distt = t->frags[i].info - pre->frags[j].info;
//					//if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
//					//	continue;
//
//					//a_score = score_alpha(distr, distt, t->frags[i].len);
//					//b_score = score_beta (distr, distt, t->frags[i].len);
//
//					temp_score = dp[jj][j].score + score_alpha(distr, distt, kmer) - score_beta(distr, distt, kmer);
//					if (temp_score > dp[ii][i].score) {
//						dp[ii][i].score = temp_score;
//						dp[ii][i].prev_list = jj;
//						dp[ii][i].prev_ind = j;
//					}
//					j--;
//				}
//				pre_start_ind = j + 1;
//			}
//		}
//
//		for (i = 0; i < t->frag_count; i++) {
//			if (dp[ii][i].score > best_score) {
//				//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[ii][i].score, i);
//				best_score = dp[ii][i].score;
//				best_count = 1;
//				best_indices[0].prev_list = ii;
//				best_indices[0].prev_ind = i;
//			}
//			else if (dp[ii][i].score == best_score and best_count < max_best) {
//				best_indices[best_count].prev_list = ii;
//				best_indices[best_count].prev_ind = i;
//				best_count++;
//			}
//		}
//	}
//
//	// back-tracking
//	deque<chain_cell> chain_index;
//	chain_cell best_index;	
//	best_chain.best_chain_count = best_count;
//	int tmp_list_ind;
//	for (int j = 0; j < best_count; j++) {
//		chain_index.clear();
//		best_index = best_indices[j];
//		while (best_index.prev_list != -1) {
//			chain_index.push_front(best_index);
//			tmp_list_ind = best_index.prev_list;
//			best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
//			best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;
//		}
//
//		best_chain.chains[j].score = best_score;
//		best_chain.chains[j].chain_len = chain_index.size();
//
//		for (i = 0; i < chain_index.size(); i++) {
//			//best_chain[j].frags[i] = frag_tmp[chain_index[i].prev_list]->frags[chain_index[i].prev_ind];
//			best_chain.chains[j].frags[i].rpos = (fragment_list + chain_index[i].prev_list)->frags[chain_index[i].prev_ind].info;
//			best_chain.chains[j].frags[i].qpos = (fragment_list + chain_index[i].prev_list)->qpos;
//			best_chain.chains[j].frags[i].len = kmer;
//		}
//	}
//}

// merging instead of binary search
// Assumption: fragment list sorted by rpos
//
// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
//void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, chain_list& best_chain) {
//	int i, j;
//
//	int distr, distt;
//	double a_score, b_score;
//	int kmer_cnt = 7;
//	int max_frag_cnt = FRAGLIM;
//	chain_cell dp[kmer_cnt][max_frag_cnt + 1];
//
//	double best_score = -1;
//	int max_best = BESTCHAINLIM;
//	vector <chain_cell> best_indices(max_best);
//	int best_count = 0;
//	double temp_score;
//
//	GIMatchedKmer* t;
//	GIMatchedKmer* pre;
//	int ii = 0;
//	int jj = 0;
//	int max_rpos_lim;
//	int pre_start_ind;
//	for (ii = 0; ii < kmer_cnt; ii++) {
//		t = fragment_list + ii;
//		for (i = 0; i < t->frag_count; i++) {
//			dp[ii][i].score = kmer;
//			//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
//			//fprintf(stderr, "Frag Len: %f\n", dp[i]);
//			dp[ii][i].prev_list = -1;
//			dp[ii][i].prev_ind = -1;
//		}
//
//		for (jj = ii-1; jj >= 0; jj--) {
//			pre = fragment_list + jj;
//			if (pre->frag_count <= 0)  
//				continue;
//
//			pre_start_ind = 0;
//			for (i = 0; i < t->frag_count; i++) {
//				if ((t->frags[i].info - 1) < pre->frags[pre_start_ind].info)
//					continue;
//
//				max_rpos_lim = t->frags[i].info - GENETHRESH;
//
//				while ((pre_start_ind < pre->frag_count) and (pre->frags[pre_start_ind].info < max_rpos_lim))	// skip out of range
//					pre_start_ind++;
//
//				if (pre_start_ind >= pre->frag_count)	// no more fragment from pre are in range for the remaining of the t
//					break;
//
//				//j = frag_binary_search(pre->frags, 0, pre->frag_count, t->frags[i].rpos - GENETHRESH - 1);
//				//j = frag_binary_search(pre->frags, pre_start_ind, pre->frag_count, t->frags[i].info - 1) - 1;
//				j = pre_start_ind;
//
//				//fprintf(stderr, "found index on flist[%d]: %d\n", jj, j);
//				//while ((j >= pre_start_ind) and (pre->frags[j].info >= max_rpos_lim)) {
//				while ((j < pre->frag_count) and (pre->frags[j].info < t->frags[i].info)) {
//					//fprintf(stderr, "working on t[%d], pre[%d]\n", i, j);
//					distr = t->qpos - pre->qpos;
//					//if (distr <= 0)
//					//	continue;
//
//					distt = t->frags[i].info - pre->frags[j].info;
//					//if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
//					//	continue;
//
//					//a_score = score_alpha(distr, distt, t->frags[i].len);
//					//b_score = score_beta (distr, distt, t->frags[i].len);
//
//					temp_score = dp[jj][j].score + score_alpha(distr, distt, kmer) - score_beta(distr, distt, kmer);
//					if (temp_score > dp[ii][i].score) {
//						dp[ii][i].score = temp_score;
//						dp[ii][i].prev_list = jj;
//						dp[ii][i].prev_ind = j;
//					}
//					//j--;
//					j++;
//				}
//				//pre_start_ind = j + 1;
//			}
//		}
//
//		for (i = 0; i < t->frag_count; i++) {
//			if (dp[ii][i].score > best_score) {
//				//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[ii][i].score, i);
//				best_score = dp[ii][i].score;
//				best_count = 1;
//				best_indices[0].prev_list = ii;
//				best_indices[0].prev_ind = i;
//			}
//			else if (dp[ii][i].score == best_score and best_count < max_best) {
//				best_indices[best_count].prev_list = ii;
//				best_indices[best_count].prev_ind = i;
//				best_count++;
//			}
//		}
//	}
//
//	// back-tracking
//	deque<chain_cell> chain_index;
//	chain_cell best_index;	
//	best_chain.best_chain_count = best_count;
//	int tmp_list_ind;
//	for (int j = 0; j < best_count; j++) {
//		chain_index.clear();
//		best_index = best_indices[j];
//		while (best_index.prev_list != -1) {
//			chain_index.push_front(best_index);
//			tmp_list_ind = best_index.prev_list;
//			best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
//			best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;
//		}
//
//		best_chain.chains[j].score = best_score;
//		best_chain.chains[j].chain_len = chain_index.size();
//
//		for (i = 0; i < chain_index.size(); i++) {
//			//best_chain[j].frags[i] = frag_tmp[chain_index[i].prev_list]->frags[chain_index[i].prev_ind];
//			best_chain.chains[j].frags[i].rpos = (fragment_list + chain_index[i].prev_list)->frags[chain_index[i].prev_ind].info;
//			best_chain.chains[j].frags[i].qpos = (fragment_list + chain_index[i].prev_list)->qpos;
//			best_chain.chains[j].frags[i].len = kmer;
//		}
//	}
//}

// removed queue
// Assumption: fragment list sorted by rpos
//
// f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
//void chain_seeds_sorted_kbest(GIMatchedKmer*& fragment_list, vector <chain_t>& best_chain) {
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
			//fprintf(stderr, "Frag Len: %zu\n", fragment_list[i].len);
			//fprintf(stderr, "Frag Len: %f\n", dp[i]);
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

				//j = frag_binary_search(pre->frags, 0, pre->frag_count, t->frags[i].rpos - GENETHRESH - 1);
				//j = frag_binary_search(pre->frags, pre_start_ind, pre->frag_count, t->frags[i].info - 1) - 1;
				j = pre_start_ind;

				//fprintf(stderr, "found index on flist[%d]: %d\n", jj, j);
				//while ((j >= pre_start_ind) and (pre->frags[j].info >= max_rpos_lim)) {
				while ((j < pc_mk->frag_count) and (pc_mk->frags[j].info <= max_rpos_lim)) {
					//fprintf(stderr, "working on t[%d], pre[%d]\n", i, j);
					distr = pc_mk->qpos - cur_mk->qpos;
					//if (distr <= 0)
					//	continue;

					distt = pc_mk->frags[j].info - cur_mk->frags[i].info;
					//if (distt <= 0 or distt > GENETHRESH)	// should be in gene size range
					//	continue;

					//a_score = score_alpha(distr, distt, t->frags[i].len);
					//b_score = score_beta (distr, distt, t->frags[i].len);

					temp_score = dp[jj][j].score + score_alpha(distr, distt, kmer) - score_beta(distr, distt, kmer);
					if (temp_score > dp[ii][i].score) {
						dp[ii][i].score = temp_score;
						dp[ii][i].prev_list = jj;
						dp[ii][i].prev_ind = j;
					}
					//j--;
					j++;
				}
				//pre_start_ind = j + 1;
			}
		}

		for (i = 0; i < cur_mk->frag_count; i++) {
			if (dp[ii][i].score > best_score) {
				//fprintf(stderr, "Best score updated to %f\ti = %d\n", dp[ii][i].score, i);
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
	//deque<chain_cell> chain_index;
	chain_cell best_index;	
	best_chain.best_chain_count = best_count;
	int tmp_list_ind;
	for (j = 0; j < best_count; j++) {
		//chain_index.clear();
		best_index = best_indices[j];
		i = 0;
		while (best_index.prev_list != -1) {
			//chain_index.push_front(best_index);
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
