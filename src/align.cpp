#include <cstdio>
#include <cstring>
#include <string>
#include <set>

#include "align.h"
#include "common.h"

using namespace std;

Alignment::Alignment(void) {
	init();
}

Alignment::~Alignment(void) {
}

void Alignment::init(void) {
	string nucs = "ACGTacgt";
	int nlen = nucs.length();
	int i;

	for (i = 0; i < ASCISIZE; i++)
		fill(diff_ch[i], diff_ch[i] + ASCISIZE, 1);

	for (i = 0; i < nlen; i++)
		diff_ch[nucs[i]][nucs[i]] = 0;

	int half_nlen = nlen / 2;
	for (i = 0; i < half_nlen; i++) {
		diff_ch[nucs[i]][nucs[i+half_nlen]] = 0;
		diff_ch[nucs[i+half_nlen]][nucs[i]] = 0;
	}
}

// dir > 0
bool Alignment::hamming_match_right(char* s, int n, char* t, int m) {
	int mm_th = SOFTCLIPTH / 2;
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist = 0;

	int mid_range = min_len - SOFTCLIPTH;
	for (int i = 0; i < mid_range; i++) {
		mid_dist += diff_ch[s[i]][t[i]];
		if (mid_dist > EDTH)
			return false;
	}
	for (int i = mid_range; i < min_len; i++)
		side_dist += diff_ch[s[i]][t[i]];

	return (side_dist >= mm_th) ? true : ((mid_dist + side_dist) <= EDTH);
}

// dir > 0
// return edit distance, soft clipped length
int Alignment::hamming_distance_right(char* s, int n, char* t, int m, int& sclen) {
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist[SOFTCLIPTH+1] = {0};

	int mid_range = min_len - SOFTCLIPTH;
	for (int i = 0; i < mid_range; i++) {
		mid_dist += diff_ch[s[i]][t[i]];
		if (mid_dist > EDTH)
			return mid_dist;
	}

	for (int i = min_len - 1; i >= mid_range; i--)
		side_dist[min_len - i] = side_dist[min_len - i - 1] + diff_ch[s[i]][t[i]];

	int i = SOFTCLIPTH;
	while (i > 0 and (side_dist[i] == side_dist[i-1] or (2 * side_dist[i]) < i))
		i--;

	sclen = i;
	if (sclen == SOFTCLIPTH and diff_ch[s[mid_range-1]][t[mid_range-1]] == 1)	// mismatch exactly after the soft clip ends =>  not extendable
		return EDTH+1;

	return mid_dist + side_dist[SOFTCLIPTH] - side_dist[i];
}

// dir < 0
bool Alignment::hamming_match_left(char* s, int n, char* t, int m) {
	int mm_th = SOFTCLIPTH / 2;
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist = 0;

	for (int i = SOFTCLIPTH; i < min_len; i++) {
		mid_dist += diff_ch[s[i]][t[i]];
		if (mid_dist > EDTH)
			return false;
	}

	for (int i = 0; i < SOFTCLIPTH; i++)
		side_dist += diff_ch[s[i]][t[i]];

	return (side_dist >= mm_th) ? true : ((mid_dist + side_dist) <= EDTH);
}

// dir < 0
int Alignment::hamming_distance_left(char* s, int n, char* t, int m, int& sclen) {
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist[SOFTCLIPTH+1] = {0};

	for (int i = SOFTCLIPTH; i < min_len; i++) {
		mid_dist += diff_ch[s[i]][t[i]];
		if (mid_dist > EDTH)
			return mid_dist;
	}

	for (int i = 0; i < SOFTCLIPTH; i++)
		side_dist[i + 1] = side_dist[i] + diff_ch[s[i]][t[i]];

	int i = SOFTCLIPTH;
	while (i > 0 and (side_dist[i] == side_dist[i-1] or (2 * side_dist[i]) < i))
		i--;

	sclen = i;
	if (sclen == SOFTCLIPTH and diff_ch[s[SOFTCLIPTH]][t[SOFTCLIPTH]] == 1)	// mismatch exactly after the soft clip ends =>  not extendable
		return EDTH+1;

	return mid_dist + side_dist[SOFTCLIPTH] - side_dist[i];
}

int Alignment::global_alignment(char* s, int n, char* t, int m, int gap_pen, int mm_pen) {
	for (int i = 0; i <= n; i++)
		dp[i][0] = i * gap_pen;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j * gap_pen;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			int penalty = mm_pen * diff_ch[s[i-1]][t[j-1]];
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

void Alignment::global_banded_alignment(char* s, int n, char* t, int m) {
	int i, j;
	int penalty;

	for (i = 0; i <= INDELTH; i++)
		dp[i][0] = i;

	for (j = 0; j <= INDELTH; j++)
		dp[0][j] = j;

	for (j = INDELTH + 1; j <= m; j++) {
		i = j - INDELTH;
		dp[i][j] = dp[i-1][j-1] + diff_ch[s[i-1]][t[j-1]];
	}
	
	for (j = 1; j <= m; j++) {
		i = j + INDELTH;
		dp[i][j] = dp[i-1][j-1] + diff_ch[s[i-1]][t[j-1]];
	}

	for (j = 1; j < INDELTH; j++) {
		for (i = 1; i < j + INDELTH; i++) {
			dp[i][j] = dp[i-1][j-1] + diff_ch[s[i-1]][t[j-1]];
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + 1), dp[i][j-1] + 1);
		}
	}

	for (j = INDELTH; j <= m; j++) {
		for (i = j - INDELTH + 1; i < j + INDELTH; i++) {
			dp[i][j] = dp[i-1][j-1] + diff_ch[s[i-1]][t[j-1]];
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + 1), dp[i][j-1] + 1);
		}
	}

	//fprintf(stderr, "DP\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

}

void Alignment::global_banded_alignment_reverse(char* s, int n, char* t, int m) {
	int i, j;
	int penalty;

	for (i = n - INDELTH; i <= n; i++)
		dp[i][m] = n - i;

	for (j = m - INDELTH; j <= m; j++)
		dp[n][j] = m - j;

	int _2indelth = 2 * INDELTH;
	for (j = m - INDELTH - 1; j >= 0; j--) {
		i = j + _2indelth;
		dp[i][j] = dp[i+1][j+1] + diff_ch[s[i]][t[j]];
	}
	
	for (j = m - 1; j >= 0; j--) {
		i = j;
		dp[i][j] = dp[i+1][j+1] + diff_ch[s[i]][t[j]];
	}

	for (j = m - 1; j > m - INDELTH; j--) {
		for (i = n - 1; i > j; i--) {
			dp[i][j] = dp[i+1][j+1] + diff_ch[s[i]][t[j]];
			dp[i][j] = minM(minM(dp[i][j], dp[i+1][j] + 1), dp[i][j+1] + 1);
		}
	}

	for (j = m - INDELTH; j >= 0; j--) {
		for (i = j + _2indelth - 1; i > j; i--) {
			dp[i][j] = dp[i+1][j+1] + diff_ch[s[i]][t[j]];
			dp[i][j] = minM(minM(dp[i][j], dp[i+1][j] + 1), dp[i][j+1] + 1);
		}
	}

	//fprintf(stderr, "DP\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

}

void Alignment::hamming_distance_bottom(char* s, int n, char* t, int m) {
	for (int i = n - 2*INDELTH; i <= n; i++)
		hamm[i][m] = 0;

	for (int j = m - 1; j >= m - SOFTCLIPTH; j--) {
		for (int i = j - INDELTH; i <= j + INDELTH; i++) {
			hamm[i][j] = hamm[i+1][j+1] + diff_ch[s[i]][t[j]];
		}
	}
}

void Alignment::hamming_distance_top(char* s, int n, char* t, int m) {
	for (int i = 0; i <= 2*INDELTH; i++)
		hamm[i][0] = 0;

	for (int j = 1; j <= SOFTCLIPTH; j++) {
		for (int i = j; i <= j + 2*INDELTH; i++) {
			hamm[i][j] = hamm[i-1][j-1] + diff_ch[s[i-1]][t[j-1]];
		}
	}
}

void Alignment::hamming_distance(char* s, int n, char* t, int m) {
	for (int i = 0; i <= INDELTH; i++)
		hamm[i][0] = 0;

	for (int j = 0; j <= INDELTH; j++)
		hamm[0][j] = 0;

	for (int j = 1; j <= m; j++) {
		for (int i = maxM(1, j - INDELTH); i <= j + INDELTH; i++) {
			hamm[i][j] = hamm[i-1][j-1] + diff_ch[s[i-1]][t[j-1]];
		}
	}

	//fprintf(stderr, "Hamm\n");
	//for (int j = 0; j <= m; j++) {
	//	for (int i = 0; i <= n; i++)
	//		fprintf(stderr, "%2d ", hamm[i][j]);
	//	fprintf(stderr, "\n");
	//}
}

int Alignment::local_alignment_right(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_sclen = SOFTCLIPTH;
	int max_indel = INDELTH;
	int max_edit  = EDTH;

	global_banded_alignment(s, n, t, m);
	hamming_distance_bottom(s, n, t, m);

	AlignCandid best(max_edit + 1, max_sclen + 1, max_indel + 1);

	for (int j = m; j >= m - max_sclen; j--) {
		for (int i = j - max_indel; i <= j + max_indel; i++) {
			if (dp[i][j] <= max_edit and (2 * hamm[i][j]) >= (m - j) and (!diff_ch[s[i-1]][t[j-1]]))
				best = minM(best, AlignCandid(dp[i][j], m - j, j - i));
		}
	}

	//fprintf(stderr, "Right: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
}

int Alignment::local_alignment_left(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	int max_sclen = SOFTCLIPTH;
	int max_indel = INDELTH;
	int max_edit  = EDTH;

	global_banded_alignment_reverse(s, n, t, m);
	hamming_distance_top(s, n, t, m);

	AlignCandid best(max_edit + 1, max_sclen + 1, max_indel + 1);

	for (int j = 0; j <= max_sclen; j++) {
		for (int i = j; i <= j + 2*max_indel; i++) {
			if (dp[i][j] <= max_edit and (2 * hamm[i][j]) >= j and (!diff_ch[s[i]][t[j]]))
				best = minM(best, AlignCandid(dp[i][j], j, i - j - max_indel));
		}
	}

	//fprintf(stderr, "Left: Best: ed = %d\tsclen = %d\t indel = %d\n", best.ed, best.sclen, best.indel);

	sc_len = best.sclen;
	indel = best.indel;
	return best.ed;
}

int Alignment::local_alignment_left2(char* s, int n, char* t, int m, int& sc_len, int& indel) {
	char sp[n];
	char tp[m];
	for (int i = 0; i < n; i++)
		sp[i] = s[n-1-i];
	
	for (int j = 0; j < m; j++)
		tp[j] = t[m-1-j];

	return local_alignment_right(sp, n, tp, m, sc_len, indel);
}

