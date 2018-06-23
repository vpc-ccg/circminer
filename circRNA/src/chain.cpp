#include <algorithm>
#include <deque>
#include <cstdio>

#include "chain.h"
#include "fragment_list.h"

#define INF 100000

inline double score_beta(int distr, int distt, int frag_len) {
	int maxd = distr < distt ? distt : distr;
	int mind = distr < distt ? distr : distt;

	return 0.15 * (maxd - mind) + 0 * mind;
}

inline double score_alpha(int distr, int distl, int frag_len) {
	return 7 * frag_len;
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
					distr = t->frags[i].qpos - (pre->frags[j].qpos + pre->frags[j].len - 1);
					if (distr <= 0)
						continue;

					distt = t->frags[i].rpos - (pre->frags[j].rpos + pre->frags[j].len - 1);
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
