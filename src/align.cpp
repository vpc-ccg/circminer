#include <cstdio>
#include <cstring>
#include <string>
#include <set>

#include "align.h"
#include "common.h"
#include "utils.h"

using namespace std;

#define DPTINF 10000000

#define compBP(A,B) ((A ^ B) ? (1) : (0))

Alignment::Alignment(void) {
	init();
}

Alignment::~Alignment(void) {
}

void Alignment::init(void) {
	diff_ch = edit_mat.diff_ch;
	score_ch = score_mat.diff_ch;
}

// dir > 0
// return edit distance, soft clipped length
int Alignment::hamming_distance(char* s, char* t, int n) {
	int ed = 0;

	for (int i = 0; i < n; i++) {
		ed += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t>(t[i]) ];
		if (ed > maxEd)
			return ed;
	}

	return ed;
}

// dir > 0
bool Alignment::hamming_match_right(char* s, int n, char* t, int m) {
	int mm_th = maxSc / 2;
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist = 0;

	int mid_range = min_len - maxSc;
	for (int i = 0; i < mid_range; i++) {
		mid_dist += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];
		if (mid_dist > maxEd)
			return false;
	}
	for (int i = mid_range; i < min_len; i++)
		side_dist += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];

	return (side_dist >= mm_th) ? true : ((mid_dist + side_dist) <= maxEd);
}

// dir > 0
// return edit distance, soft clipped length
int Alignment::hamming_distance_right(char* s, int n, char* t, int m, int& sclen) {
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist[maxSc+1];
	memset(side_dist, 0, (maxSc+1) * sizeof(int));

	int mid_range = min_len - maxSc;
	for (int i = 0; i < mid_range; i++) {
		mid_dist += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];
		if (mid_dist > maxEd)
			return mid_dist;
	}

	for (int i = min_len - 1; i >= mid_range; i--)
		side_dist[min_len - i] = side_dist[min_len - i - 1] + 
					diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];

	int i = maxSc;
	while (i > 0 and (side_dist[i] == side_dist[i-1] or (2 * side_dist[i]) < i))
		i--;

	sclen = i;
	// mismatch exactly after the soft clip ends =>  not extendable
	if (sclen == maxSc and 
		diff_ch[ static_cast<uint8_t> (s[mid_range-1]) ][ static_cast<uint8_t> (t[mid_range-1]) ] == 1)
		return maxEd+1;

	return mid_dist + side_dist[maxSc] - side_dist[i];
}

// dir < 0
bool Alignment::hamming_match_left(char* s, int n, char* t, int m) {
	int mm_th = maxSc / 2;
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist = 0;

	for (int i = maxSc; i < min_len; i++) {
		mid_dist += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];
		if (mid_dist > maxEd)
			return false;
	}

	for (int i = 0; i < maxSc; i++)
		side_dist += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];

	return (side_dist >= mm_th) ? true : ((mid_dist + side_dist) <= maxEd);
}

// dir < 0
int Alignment::hamming_distance_left(char* s, int n, char* t, int m, int& sclen) {
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist[maxSc+1];
	memset(side_dist, 0, (maxSc+1) * sizeof(int));

	for (int i = maxSc; i < min_len; i++) {
		mid_dist += diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];
		if (mid_dist > maxEd)
			return mid_dist;
	}

	for (int i = 0; i < maxSc; i++)
		side_dist[i + 1] = side_dist[i] + diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[i]) ];

	int i = maxSc;
	while (i > 0 and (side_dist[i] == side_dist[i-1] or (2 * side_dist[i]) < i))
		i--;

	sclen = i;
	// mismatch exactly after the soft clip ends =>  not extendable
	if (sclen == maxSc and 
		diff_ch[ static_cast<uint8_t> (s[maxSc]) ][ static_cast<uint8_t> (t[maxSc]) ] == 1)
		return maxEd+1;

	return mid_dist + side_dist[maxSc] - side_dist[i];
}

int Alignment::global_alignment(char* s, int n, char* t, int m, int gap_pen, int mm_pen) {
	for (int i = 0; i <= n; i++)
		dp[i][0] = i * gap_pen;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j * gap_pen;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			int penalty = mm_pen * diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1])];
			dp[i][j] = dp[i-1][j-1] + penalty;
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + gap_pen), dp[i][j-1] + gap_pen);
		}
	}

	//fprintf(stderr, "DP\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

	return dp[n][m];
}

int Alignment::global_alignment(char* s, int n, char* t, int m) {
	for (int i = 0; i <= n; i++)
		dp[i][0] = i;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			dp[i][j] = dp[i-1][j-1] + diff_ch[static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ];
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + 1), dp[i][j-1] + 1);
		}
	}

	//fprintf(stderr, "DP\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

	return dp[n][m];
}

int Alignment::global_alignment_reverse(char* s, int n, char* t, int m) {
	for (int i = 0; i <= n; i++)
		dp[i][0] = i;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			dp[i][j] = dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[n-i]) ][ static_cast<uint8_t> (t[m-j]) ];
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + 1), dp[i][j-1] + 1);
		}
	}

	//fprintf(stderr, "DP\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

	return dp[n][m];
}

// w is the band size
// m = n + w
// w >= 0 : m >= n
// one sided means we are not allowing complex insertion and deletions
// returns edit distance
int Alignment::global_one_side_banded_alignment(char* s, int n, char* t, int m, int w) {
	if (w < 0 or n <= w) {
		global_alignment(s, n, t, m);
		return dp[n][m];
	}

	int i, j;
	
	// memset(dp, 127, MAXSTRSIZE * MAXSTRSIZE * sizeof(dp[0][0]));
	j = 0;
	for (i = 1; i <= n; i++) {
		dp[i][j++] = DPTINF;
	}

	i = 0;
	for (j = w+1; j <= m; j++) {
		dp[i++][j] = DPTINF;
	}

	// uint32_t dp[n+1][m+1];
	// memset(dp, 127, (n+1) * (m+1) * sizeof(uint32_t));

	for (j = 0; j <= w; j++) {
		dp[0][j] = j;
	}

	for (i = 1; i <= n; i++) {
		for (j = i; j <= i + w; j++) {
			dp[i][j] = dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t>(t[j-1]) ];
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + 1), dp[i][j-1] + 1);
		}
	}
	return dp[n][m];
}

void Alignment::global_banded_alignment_drop(char* s, int n, char* t, int m, int w, int& on_s, int& on_t) {
	int32_t pre_optimum = 0;
	int32_t cur_optimum = 0;

	int i, j, k;

	//for (int i = 0; i <= n; ++i)
	//	for (int j = 0; j <= m; ++j)
	//		dpx[i][j] = 100;

	j = 0;
	for (i = w+1; i <= n; i++) {
		dpx[i][j++] = -DPTINF;
	}

	i = 0;
	for (j = w+1; j <= m; j++)
		dpx[i++][j] = -DPTINF;

	for (i = 0; i <= w; i++) {
		dpx[i][0] = i * score_mat.ind_sc;
		dpx[0][i] = i * score_mat.ind_sc;
	}

	on_s = 0;
	on_t = 0;

	if (m <= 0 or n <= 0)
		return;

	// lower and upper bound 
	// lb[0], ub[0] -> current anti diagonal
	// lb[1], ub[1] -> past anti diagonal
	// lb[2], ub[2] -> the anti diagonal before last
	int lb[3] = {1, 1, 1}, ub[3] = {1, DPTINF, DPTINF};

	// int new_lb;
	int new_ub, pre_ub = 0;
	int best_i = 0, best_j = 0;
	for (k = 2; k <= m+n; ++k) {
		//new_lb = -1;
		new_ub = -1;
		for (i = lb[0]; i <= ub[0]; ++i) {
			j = k - i;

			dpx[i][j] = maxM3(dpx[i-1][j-1] + score_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ], 
							 dpx[i-1][j] + score_mat.ind_sc, 
							 dpx[i][j-1] + score_mat.ind_sc);

// 			fprintf(stderr, "i: %d, j: %d\n", i, j);
			cur_optimum = maxM(cur_optimum, dpx[i][j]);
			if (dpx[i][j] >= cur_optimum) {
				cur_optimum = dpx[i][j];
				best_i = i;
				best_j = j;
			}

			if (dpx[i][j] + score_mat.xd < pre_optimum)
				dpx[i][j] = -DPTINF;

			// get first and last non-minus-infinity values in current anti-diagonal
			//if (dpx[i][j] > -DPTINF) {
			//	new_ub = i;
			//	if (new_lb == -1)
			//		new_lb = i;
			//}

			if (dpx[i][j] > -DPTINF) {
				new_ub = i;
			}
		}
		
		// update lower-bound and upper-bound based on band
		//lb[2] = lb[1];
		//ub[2] = ub[1];
		//lb[1] = new_lb;
		//ub[1] = new_ub;

		int lb_t = k - lb[0];
		if ( lb_t == m or (k > w and ((k - w) % 2 == 0)) )
			++lb[0];

		if ( ub[0] < n and (k <= w or (k > w and ((k - w) % 2 == 1))) )
			++ub[0];

		//lb[0] = maxM(lb[0], minM(lb[1], lb[2] + 1));
		//ub[0] = minM(ub[0], maxM(lb[1], lb[2]));

		if ((pre_ub == -1 and new_ub == -1) or lb[0] > ub[0]) {
			//fprintf(stderr, "k: %d\n", k);
			break;
		}

		pre_ub = new_ub;

		//pre_optimum = cur_optimum;
		pre_optimum = maxM(pre_optimum, cur_optimum);
	}

	on_s = best_i;
	on_t = best_j;

	// fprintf(stderr, "DP\n");
	// fprintf(stderr, "s:%s\nt:%s\n", s, t);
	// fprintf(stderr, "Best score: dpx[%d][%d] = %d\n", best_i, best_j, dpx[best_i][best_j]);

	// for (int i = 0; i <= n; i++) {
	// 	if (i == 0)
	// 		fprintf(stderr, "          ");
	// 	else
	// 		fprintf(stderr, "%4d ", i);
	// }
	// fprintf(stderr, "\n");
	// for (int i = 0; i <= n; i++) {
	// 	if (i == 0)
	// 		fprintf(stderr, "          ");
	// 	else
	// 		fprintf(stderr, "%4c ", s[i-1]);
	// }
	// fprintf(stderr, "\n");

	// for (int j = 0; j <= m; j++) {
	// 	if (j == 0)
	// 		fprintf(stderr, "     ");
	// 	else
	// 		fprintf(stderr, "%2d %c ", j, t[j-1]);

	// 	for (int i = 0; i <= n; i++) {
	// 		if (dpx[i][j] > -DPTINF)
	// 			fprintf(stderr, "%4d ", dpx[i][j]);
	// 		else
	// 			fprintf(stderr, "---- ");
	// 	}
	// 	fprintf(stderr, "\n");
	// }
}

// w is the band size
// w < 0 disables band
// n >= m
void Alignment::global_banded_alignment(char* s, int n, char* t, int m, int w) {
	if (w < 0 or n <= 2 * w or m <= w) {
		global_alignment(s, n, t, m);
		return;
	}

	int i, j;

	j = 0;
	for (i = w+1; i <= n; i++) {
		dp[i][j++] = DPTINF;
	}

	i = 0;
	for (j = w+1; j <= m; j++)
		dp[i++][j] = DPTINF;

	for (i = 0; i <= w; i++) {
		dp[i][0] = i;
		dp[0][i] = i;
	}

	for (j = 1; j <= w; j++) {
		for (i = 1; i <= j + w; i++) {
			dp[i][j] = minM3(dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ], 
							 dp[i-1][j] + 1, 
							 dp[i][j-1] + 1);
		}
	}

	for (j = w + 1; j <= n - w; j++) {
		for (i = j - w; i <= j + w; i++) {
			dp[i][j] = minM3(dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ], 
							 dp[i-1][j] + 1, 
							 dp[i][j-1] + 1);
		}	
	}
	for (j = n - w + 1; j <= m; j++) {
		for (i = j - w; i <= n; i++) {
			dp[i][j] = minM3(dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ], 
							 dp[i-1][j] + 1, 
							 dp[i][j-1] + 1);
		}	
	}
	
	// fprintf(stderr, "DP\n");
	// for (int j = 0; j <= m; j++) {
	// 	for (int i = 0; i <= n; i++)
	// 		fprintf(stderr, "%2d ", dp[i][j]);
	// 	fprintf(stderr, "\n");
	// }

}

// n >= m
void Alignment::global_banded_alignment_reverse(char* s, int n, char* t, int m, int w) {
	if (w < 0 or n <= 2 * w or m <= w) {
		global_alignment_reverse(s, n, t, m);
		return;
	}

	int i, j;

	j = 0;
	for (i = w+1; i <= n; i++) {
		dp[i][j++] = DPTINF;
	}

	i = 0;
	for (j = w+1; j <= m; j++)
		dp[i++][j] = DPTINF;

	for (i = 0; i <= w; i++) {
		dp[i][0] = i;
		dp[0][i] = i;
	}

	for (j = 1; j <= w; j++) {
		for (i = 1; i <= j + w; i++) {
			dp[i][j] = minM3(dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[n-i]) ][ static_cast<uint8_t> (t[m-j]) ],
							 dp[i-1][j] + 1,
							 dp[i][j-1] + 1);
		}
	}

	for (j = w + 1; j <= n - w; j++) {
		for (i = j - w; i <= j + w; i++) {
			dp[i][j] = minM3(dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[n-i]) ][ static_cast<uint8_t> (t[m-j]) ],
							 dp[i-1][j] + 1,
							 dp[i][j-1] + 1);
		}	
	}
	for (j = n - w + 1; j <= m; j++) {
		for (i = j - w; i <= n; i++) {
			dp[i][j] = minM3(dp[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[n-i]) ][ static_cast<uint8_t> (t[m-j]) ],
							 dp[i-1][j] + 1,
							 dp[i][j-1] + 1);
		}	
	}
	

	// fprintf(stderr, "DP\n");
	// for (int j = 0; j <= m; j++) {
	// 	for (int i = 0; i <= n; i++)
	// 		fprintf(stderr, "%2u ", dp[i][j]);
	// 	fprintf(stderr, "\n");
	// }

}


// TO DO: fix n < m + w
void Alignment::hamming_distance_bottom(char* s, int n, char* t, int m, int max_sclen) {
	for (int i = maxM(0, n - 2*bandWidth); i <= n; i++)
		hamm[i][m] = 0;

	for (int j = m - 1; j >= max_sclen; j--) {
		for (int i = maxM(0, j - bandWidth); i <= minM(j + bandWidth, n); i++) {
			hamm[i][j] = hamm[i+1][j+1] + diff_ch[ static_cast<uint8_t> (s[i]) ][ static_cast<uint8_t> (t[j]) ];
		}
	}
}

void Alignment::hamming_distance_top(char* s, int n, char* t, int m, int max_sclen) {
	for (int i = 0; i <= minM(2*bandWidth, n); i++)
		hamm[i][0] = 0;

	for (int j = 1; j <= max_sclen; j++) {
		for (int i = j; i <= minM(j + 2*bandWidth, n); i++) {
			hamm[i][j] = hamm[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ];
		}
	}
}

void Alignment::hamming_distance(char* s, int n, char* t, int m) {
	for (int i = 0; i <= bandWidth; i++)
		hamm[i][0] = 0;

	for (int j = 0; j <= bandWidth; j++)
		hamm[0][j] = 0;

	for (int j = 1; j <= m; j++) {
		for (int i = maxM(1, j - bandWidth); i <= j + bandWidth; i++) {
			hamm[i][j] = hamm[i-1][j-1] + diff_ch[ static_cast<uint8_t> (s[i-1]) ][ static_cast<uint8_t> (t[j-1]) ];
		}
	}

	//fprintf(stderr, "Hamm\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", hamm[i][j]);
	//	fprintf(stderr, "\n");
	//}
}

int Alignment::local_alignment_right(char* s, int n, char* t, int m, int& indel) {
	int max_sclen = maxSc;
	int max_indel = bandWidth;
	uint32_t max_edit  = maxEd;

	global_banded_alignment(s, n, t, m, bandWidth);

	AlignCandid best(max_edit + 1, max_sclen + 1, max_indel + 1);

	for (int i = maxM(0, m - max_indel); i <= minM(m + max_indel, n); i++) {
		if (dp[i][m] <= max_edit)
			//best = minM(best, AlignCandid(dp[i][m], 0, m - i));
			best.update(AlignCandid(dp[i][m], 0, m - i));
	}

	//fprintf(stderr, "Right: Best: ed = %d\tindel = %d\n", best.ed, best.indel);

	align_score = -1 * best.ed;
	indel = best.indel;
	return best.ed;
}

int Alignment::local_alignment_left(char* s, int n, char* t, int m, int& indel) {
	int max_sclen = maxSc;
	int max_indel = bandWidth;
	uint32_t max_edit  = maxEd;

	global_banded_alignment_reverse(s, n, t, m, bandWidth);
	
	AlignCandid best(max_edit + 1, max_sclen + 1, max_indel + 1);

	for (int i = maxM(0, m - max_indel); i <= minM(m + max_indel, n); i++) {
		if (dp[i][m] <= max_edit) {
			// fprintf(stderr, "[%d][%d] -> (%d, %d, %d)\n", i, m, dp[i][m], 0, m-i);
			//best = minM(best, AlignCandid(dp[i][m], 0, m - i));
			best.update(AlignCandid(dp[i][m], 0, m - i));
		}
	}

	//fprintf(stderr, "Left: Best: ed = %d\tindel = %d\n", best.ed, best.indel);

	align_score = -1 * best.ed;
	indel = best.indel;
	return best.ed;
}

int EditDistAlignment::local_alignment_right_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_sclen = minM(maxSc, m);
	int max_indel = bandWidth;
	uint32_t max_edit = maxEd;

	global_banded_alignment(s, n, t, m, bandWidth);
	//hamming_distance_bottom(s, n, t, m, max_sclen);

	AlignCandid best(max_edit + 1, maxSc + 1, max_indel + 1);

	for (int j = m; j >= m - max_sclen; j--) {
		for (int i = maxM(0, j - max_indel); i <= minM(j + max_indel, n); i++) {
			//if (dp[i][j] <= max_edit and 3*dp[i][j] < static_cast<uint32_t> (j)) {
			if (dp[i][j] <= max_edit) {
				//best = minM(best, AlignCandid(dp[i][j], m - j, j - i));
				best.update(AlignCandid(dp[i][j], m - j, j - i));
			}
		}
	}

// 	fprintf(stderr, "EDIT ALIGNMENT Right: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	// very small seq with mismatch that is not going to go over edit dist threshold
	if (m <= maxEd) {
		best.update(AlignCandid(m, 0, 0));
	}

	align_score = m - best.sclen - 2 * best.ed;
	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
}

int EditDistAlignment::local_alignment_left_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_sclen = minM(maxSc, m);
	int max_indel = bandWidth;
	uint32_t max_edit = maxEd;

	global_banded_alignment_reverse(s, n, t, m, bandWidth);
	//hamming_distance_top(s, n, t, m, max_sclen);
	
	AlignCandid best(max_edit + 1, maxSc + 1, max_indel + 1);

	for (int j = m; j >= m - max_sclen; j--) {
		for (int i = maxM(0, j - max_indel); i <= minM(j + max_indel, n); i++) {
			//if (dp[i][j] <= max_edit and 3*dp[i][j] < static_cast<uint32_t> (j)) {
			if (dp[i][j] <= max_edit) {
				// fprintf(stderr, "[%d][%d] -> (%d, %d, %d)\n", i, j, dp[i][j], m-j, j-i);
				//best = minM(best, AlignCandid(dp[i][j], m - j, j - i));
				best.update(AlignCandid(dp[i][j], m - j, j - i));
			}
		}
	}

// 	fprintf(stderr, "EDIT ALIGNMENT Left: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	// very small seq with mismatch that is not going to go over edit dist threshold
	if (m <= maxEd) {
		best.update(AlignCandid(m, 0, 0));
	}

	align_score = m - best.sclen - 2 * best.ed;
	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
}

int DropAlignment::local_alignment_right_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_indel = bandWidth;
	uint32_t max_edit  = maxEd;

	int on_s, on_t;
	global_banded_alignment_drop(s, n, t, m, bandWidth, on_s, on_t);
	int32_t score = dpx[on_s][on_t];
	uint32_t ed = (score_mat.mat_sc * maxM(on_s, on_t) - score) / (score_mat.mat_sc - score_mat.mis_sc);
	int indel_cnt = on_t - on_s;
	int clip = m - on_t;

	AlignCandid best(max_edit + 1, maxM(maxSc, m) + 1, max_indel + 1, 0);

	if (ed <= max_edit) {
		best.update(AlignCandid(ed, clip, indel_cnt, score));
	}

// 	fprintf(stderr, "DROP ALIGNMENT\nRight: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	align_score = score;
	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
}

int DropAlignment::local_alignment_left_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_indel = bandWidth;
	uint32_t max_edit  = maxEd;

	int on_s, on_t;
	char revs[1000];
	char revt[1000];

	reverse_str(s, n, revs);
	reverse_str(t, m, revt);

	global_banded_alignment_drop(revs, n, revt, m, bandWidth, on_s, on_t);
	int32_t score = dpx[on_s][on_t];
	uint32_t ed = (score_mat.mat_sc * maxM(on_s, on_t) - score) / (score_mat.mat_sc - score_mat.mis_sc);
	int indel_cnt = on_t - on_s;
	int clip = m - on_t;

	AlignCandid best(max_edit + 1, maxM(m, maxSc) + 1, max_indel + 1, 0);

	if (ed <= max_edit) {
		best.set(AlignCandid(ed, clip, indel_cnt, score));
	}

	//fprintf(stderr, "Left: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	align_score = score;
	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
}

ScoreMatrix::ScoreMatrix(void) {
	diff_ch = (int**) malloc(ASCISIZE * sizeof(int*));
	diff_ch[0] = (int*) malloc(ASCISIZE * ASCISIZE * sizeof(int));
	for (int i = 1; i < ASCISIZE; ++i)
		diff_ch[i] = diff_ch[0] + i * ASCISIZE;

	init(0, -1, -1, 0);
}

ScoreMatrix::~ScoreMatrix(void) {
	free(diff_ch[0]);
	free(diff_ch);
}

void ScoreMatrix::init(int mat, int mis, int ind, int xval) {
	mat_sc = mat;
	mis_sc = mis;
	ind_sc = ind;
	xd = xval;

	string nucs = "ACGTacgt";
	int nlen = nucs.length();
	int i;

	for (i = 0; i < ASCISIZE; i++)
		fill(diff_ch[i], diff_ch[i] + ASCISIZE, mis);

	for (i = 0; i < nlen; i++)
		diff_ch[ static_cast<uint8_t> (nucs[i]) ][ static_cast<uint8_t> (nucs[i]) ] = mat;

	int half_nlen = nlen / 2;
	for (i = 0; i < half_nlen; i++) {
		diff_ch[ static_cast<uint8_t> (nucs[i]) ][ static_cast<uint8_t> (nucs[i+half_nlen]) ] = mat;
		diff_ch[ static_cast<uint8_t> (nucs[i+half_nlen]) ][ static_cast<uint8_t>(nucs[i]) ] = mat;
	}
}
