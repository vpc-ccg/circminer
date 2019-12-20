#include <cstdio>
#include <cstring>
#include <string>
#include <set>

#include "align.h"
#include "common.h"

using namespace std;

#define DPTINF 1000000

#define compBP(A,B) ((A ^ B) ? (1) : (0))

Alignment::Alignment(void) {
	init();
}

Alignment::~Alignment(void) {
}

void Alignment::init(void) {
	diff_ch = score_mat.diff_ch;
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

int Alignment::local_alignment_right_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_sclen = minM(maxSc, m);
	int max_indel = bandWidth;
	uint32_t max_edit  = maxEd;

	global_banded_alignment(s, n, t, m, bandWidth);
	//hamming_distance_bottom(s, n, t, m, max_sclen);

	AlignCandid best(max_edit + 1, maxSc + 1, max_indel + 1);

	for (int j = m; j >= m - max_sclen; j--) {
		for (int i = maxM(0, j - max_indel); i <= minM(j + max_indel, n); i++) {
			if (dp[i][j] <= max_edit)
				//best = minM(best, AlignCandid(dp[i][j], m - j, j - i));
				best.update(AlignCandid(dp[i][j], m - j, j - i));
		}
	}

	//fprintf(stderr, "Right: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
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

	indel = best.indel;
	return best.ed;
}

int Alignment::local_alignment_left_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_sclen = minM(maxSc, m);
	int max_indel = bandWidth;
	uint32_t max_edit  = maxEd;

	global_banded_alignment_reverse(s, n, t, m, bandWidth);
	//hamming_distance_top(s, n, t, m, max_sclen);
	
	AlignCandid best(max_edit + 1, maxSc + 1, max_indel + 1);

	for (int j = m; j >= m - max_sclen; j--) {
		for (int i = maxM(0, j - max_indel); i <= minM(j + max_indel, n); i++) {
			if (dp[i][j] <= max_edit) {
				// fprintf(stderr, "[%d][%d] -> (%d, %d, %d)\n", i, j, dp[i][j], m-j, j-i);
				//best = minM(best, AlignCandid(dp[i][j], m - j, j - i));
				best.update(AlignCandid(dp[i][j], m - j, j - i));
			}
		}
	}

	//fprintf(stderr, "Left: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	sc_len = best.sclen;
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

	indel = best.indel;
	return best.ed;
}

ScoreMatrix::ScoreMatrix(void) {
	diff_ch = (int**) malloc(ASCISIZE * sizeof(int*));
	diff_ch[0] = (int*) malloc(ASCISIZE * ASCISIZE * sizeof(int));
	for (int i = 1; i < ASCISIZE; ++i)
		diff_ch[i] = diff_ch[0] + i * ASCISIZE;

	init();
}

ScoreMatrix::~ScoreMatrix(void) {
	free(diff_ch[0]);
	free(diff_ch);
}

void ScoreMatrix::init(void) {
	string nucs = "ACGTacgt";
	int nlen = nucs.length();
	int i;

	for (i = 0; i < ASCISIZE; i++)
		fill(diff_ch[i], diff_ch[i] + ASCISIZE, 1);

	for (i = 0; i < nlen; i++)
		diff_ch[ static_cast<uint8_t> (nucs[i]) ][ static_cast<uint8_t> (nucs[i]) ] = 0;

	int half_nlen = nlen / 2;
	for (i = 0; i < half_nlen; i++) {
		diff_ch[ static_cast<uint8_t> (nucs[i]) ][ static_cast<uint8_t> (nucs[i+half_nlen]) ] = 0;
		diff_ch[ static_cast<uint8_t> (nucs[i+half_nlen]) ][ static_cast<uint8_t>(nucs[i]) ] = 0;
	}
}
