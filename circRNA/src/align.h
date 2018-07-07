#ifndef __ALIGN_H__
#define __ALIGN_H__

#define maxM(A,B) (((A) > (B)) ? (A) : (B))
#define minM(A,B) (((A) < (B)) ? (A) : (B))

#define MAXSTRSIZE 100

#include "common.h"

typedef struct {
	int q_matched;
	int r_matched;
	int score;
	int ed_dist;
} loc_align;

class Alignment {
public:
	Alignment(void);
	~Alignment(void);
	void init(void);

	bool hamming_match_right(char* s, int n, char* t, int m);
	bool hamming_match_left (char* s, int n, char* t, int m);

	int  hamming_distance_right(char* s, int n, char* t, int m, int max_sc);
	int  hamming_distance_left (char* s, int n, char* t, int m, int max_sc);
	
	int alignment(char* s, int n, char* t, int m, int gap_pen, int mm_pen);
	
	loc_align local_alignment(char* s, int n, char* t, int m, int match_score, int gap_pen, int mm_pen);
	loc_align local_alignment_reverse(char* ref, int n, char* query, int m, int match_score, int gap_pen, int mm_pen);

private:
	int dp[MAXSTRSIZE][MAXSTRSIZE];
	int diff_ch[ASCISIZE][ASCISIZE];
};

extern Alignment alignment;

#endif	//__ALIGN_H__
