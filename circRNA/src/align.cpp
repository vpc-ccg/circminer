#include <vector>
#include <cstdio>
#include "align.h"

using namespace std;

template <typename T>
inline T const& minimum(T const& a, T const& b, T const& c) {
	return (a < b) ? ((a < c) ? (a) : (c)) : ((b < c) ? (b) : (c));
}

template <typename T>
inline T const& maximum(T const& a, T const& b, T const& c) {
	return (a > b) ? ((a > c) ? (a) : (c)) : ((b > c) ? (b) : (c));
}

int alignment(char* s, int n, char* t, int m, int gap_pen, int mm_pen) {
	vector <vector <int> > dp(n+1);
	for (int i = 0; i <= n; i++)
		dp[i].resize(m+1);

	for (int i = 0; i <= n; i++)
		dp[i][0] = i * gap_pen;

	for (int j = 0; j <= m; j++)
		dp[0][j] = j * gap_pen;

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			int penalty = (s[i-1] == t[j-1]) ? 0 : mm_pen;
			dp[i][j] = dp[i-1][j-1] + penalty;
			dp[i][j] = minimum(dp[i][j], dp[i-1][j] + gap_pen, dp[i][j-1] + gap_pen);
		}
	}

	//for (int i = 0; i <= n; i++) {
	//	for (int j = 0; j <= m; j++) 
	//		fprintf(stderr, "%d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

	return dp[n][m];
}

loc_align local_alignment(char* ref, int n, char* query, int m, int match_score, int gap_pen, int mm_pen) {
	vector <vector <int> > dp(n+1);
	for (int i = 0; i <= n; i++)
		dp[i].resize(m+1);

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
			dp[i][j] = maximum(dp[i][j], dp[i-1][j] + gap_pen, dp[i][j-1] + gap_pen);
			if (dp[i][j] >= max_score.score) {	// in case of same max score, greater index is prefered
				max_score.r_matched = i;
				max_score.q_matched = j;
				max_score.score = dp[i][j];
			}
		}
	}

	//for (int i = 0; i <= n; i++) {
	//	for (int j = 0; j <= m; j++) 
	//		fprintf(stderr, "%d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

	return max_score;
}

loc_align local_alignment_reverse(char* ref, int n, char* query, int m, int match_score, int gap_pen, int mm_pen) {
	vector <vector <int> > dp(n+1);
	for (int i = 0; i <= n; i++)
		dp[i].resize(m+1);

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
			dp[i][j] = maximum(dp[i][j], dp[i+1][j] + gap_pen, dp[i][j+1] + gap_pen);
			if (dp[i][j] >= max_score.score) {	// in case of same max score, greater index is prefered
				max_score.r_matched = n - i;
				max_score.q_matched = m - j;
				max_score.score = dp[i][j];
			}
		}
	}

	//for (int i = 0; i <= n; i++) {
	//	for (int j = 0; j <= m; j++) 
	//		fprintf(stderr, "%d ", dp[i][j]);
	//	fprintf(stderr, "\n");
	//}

	return max_score;
}

