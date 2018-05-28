#include <vector>
#include <cstdio>
#include "align.h"

using namespace std;

template <typename T>
inline T const& minimum(T const& a, T const& b, T const& c) {
	return (a < b) ? ((a < c) ? (a) : (c)) : ((b < c) ? (b) : (c));
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
