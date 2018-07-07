#include <cstdio>
#include <cstring>
#include <string>

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

	return (side_dist > mm_th) ? true : ((mid_dist + side_dist) <= EDTH);
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

	return (side_dist > mm_th) ? true : ((mid_dist + side_dist) <= EDTH);
}

// dir > 0
int Alignment::hamming_distance_right(char* s, int n, char* t, int m, int max_sc) {
	int mm_th = max_sc / 2;
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist = 0;

	int mid_range = min_len - max_sc;
	for (int i = 0; i < mid_range; i++)
		mid_dist += diff_ch[s[i]][t[i]];

	for (int i = mid_range; i < min_len; i++)
		side_dist += diff_ch[s[i]][t[i]];

	return (side_dist > mm_th) ? (mid_dist) : (mid_dist + side_dist);
}

// dir < 0
int Alignment::hamming_distance_left(char* s, int n, char* t, int m, int max_sc) {
	int mm_th = max_sc / 2;
	int min_len = minM(n, m);
	int mid_dist = 0;
	int side_dist = 0;

	for (int i = 0; i < max_sc; i++)
		side_dist += diff_ch[s[i]][t[i]];

	for (int i = max_sc; i < min_len; i++)
		mid_dist += diff_ch[s[i]][t[i]];

	return (side_dist > mm_th) ? (mid_dist) : (mid_dist + side_dist);
}

int Alignment::alignment(char* s, int n, char* t, int m, int gap_pen, int mm_pen) {
	for (int i = 0; i <= n; i++)
		dp[i][0] = i * gap_pen;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j * gap_pen;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			int penalty = (s[i-1] == t[j-1]) ? 0 : mm_pen;
			dp[i][j] = dp[i-1][j-1] + penalty;
			dp[i][j] = minM(minM(dp[i][j], dp[i-1][j] + gap_pen), dp[i][j-1] + gap_pen);
		}
	}

	return dp[n][m];
}

loc_align Alignment::local_alignment(char* ref, int n, char* query, int m, int match_score, int gap_pen, int mm_pen) {
	loc_align max_score;
	max_score.q_matched = 0;
	max_score.r_matched = 0;
	max_score.score = 0;

	for (int i = 0; i <= n; i++)
		dp[i][0] = i * gap_pen;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j * gap_pen;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			int score = (ref[i-1] == query[j-1]) ? match_score : mm_pen;
			dp[i][j] = dp[i-1][j-1] + score;
			dp[i][j] = maxM(maxM(dp[i][j], dp[i-1][j] + gap_pen), dp[i][j-1] + gap_pen);
			if (dp[i][j] >= max_score.score) {	// in case of same max score, greater index is prefered
				max_score.r_matched = i;
				max_score.q_matched = j;
				max_score.score = dp[i][j];
			}
		}
	}

	return max_score;
}

loc_align Alignment::local_alignment_reverse(char* ref, int n, char* query, int m, int match_score, int gap_pen, int mm_pen) {
	loc_align max_score;
	max_score.q_matched = 0;
	max_score.r_matched = 0;
	max_score.score = 0;

	for (int i = 0; i <= n; i++)
		dp[i][m] = (n - i) * gap_pen;

	for (int j = 0; j <= m; j++)
		dp[n][j] = (m - j) * gap_pen;

	for (int i = n-1; i >= 0; i--) {
		for (int j = m-1; j >= 0; j--) {
			int score = (ref[i] == query[j]) ? match_score : mm_pen;
			dp[i][j] = dp[i+1][j+1] + score;
			dp[i][j] = maxM(maxM(dp[i][j], dp[i+1][j] + gap_pen), dp[i][j+1] + gap_pen);
			if (dp[i][j] >= max_score.score) {	// in case of same max score, greater index is prefered
				max_score.r_matched = n - i;
				max_score.q_matched = m - j;
				max_score.score = dp[i][j];
			}
		}
	}

	return max_score;
}

