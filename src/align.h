#ifndef __ALIGN_H__
#define __ALIGN_H__

#define MAXSTRSIZE 600

#define DROP_ALIGNMENT 0
#define EDIT_ALIGNMENT 1

#include <cmath>
#include "common.h"

struct AlignRes {
	uint32_t pos;
	int ed;
	int sclen;
	int indel;
	int qcovlen;	// on read
	int rcovlen;	// on ref
	int score;

	AlignRes (uint32_t p, int e, int s, int i, int qc, int scr) : pos(p), ed(e), sclen(s), indel(i), qcovlen(qc), rcovlen(qc - i), score(scr) {}
	AlignRes (uint32_t p) : pos(p), ed(0), sclen(0), indel(0), qcovlen(0), rcovlen(0), score(-INF) {}

	void set (uint32_t p, int e, int s, int i, int qc, int scr) {
		this->pos = p;
		this->ed = e;
		this->sclen = s;
		this->indel = i;
		this->qcovlen = qc;
		this->rcovlen = qc - i;
		this->score = scr;
	}

	void update (int edit_dist, int sclength, uint32_t newpos, int indel, int qcovlen, int scr) {
		this->pos = newpos;
		this->ed += edit_dist;
		this->sclen = sclength;
		this->indel += indel;
		this->qcovlen += qcovlen;
		this->rcovlen += qcovlen - indel;
		this->score = scr;
	}

	bool update_by_score (AlignRes& r) {
// 		print();
// 		r.print();
		
		if (score < r.score) {
// 			fprintf(stderr, "Updated\n");
			this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			return true;
		}

		return false;
	}
		
	void update_right (AlignRes& r) {
		// fprintf(stderr, "Right UPDATE--->\n");
		// print();
		// r.print();
		
		if (r.qcovlen > qcovlen) {
			// int pre_ed = (qcovlen == 0) ? 0 : ed;
			int pre_ed = ed;
			if (r.ed <= maxEd and r.sclen <= maxSc and 2 * (r.ed - pre_ed) < (r.qcovlen - qcovlen)) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			}
		}
		else if (r.qcovlen < qcovlen) {
			if (r.ed <= maxEd and r.sclen <= maxSc and 2 * (ed - r.ed) >= (qcovlen - r.qcovlen)) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			}
		}
		else if (r.qcovlen == qcovlen) {
			if ((r.ed < ed) or (r.ed == ed and r.sclen < sclen) or (r.ed == ed and r.sclen == sclen and r.pos < pos)) {	
				// less edit distance, then less soft clip, then smaller junction distance
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			}
		}
	}

	void update_left (AlignRes& r) {
		// fprintf(stderr, "Left UPDATE--->\n");
		// print();
		// r.print();
		if (r.qcovlen > qcovlen) {
			// int pre_ed = (qcovlen == 0) ? 0 : ed;
			int pre_ed = ed;
			if (r.ed <= maxEd and r.sclen <= maxSc and 2 * (r.ed - pre_ed) < (r.qcovlen - qcovlen)) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			}
		}
		else if (r.qcovlen < qcovlen) {
			if (r.ed <= maxEd and r.sclen <= maxSc and 2 * (ed - r.ed) >= (qcovlen - r.qcovlen)) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			}
		}
		else if (r.qcovlen == qcovlen) {
			if ((r.ed < ed) or (r.ed == ed and r.sclen < sclen) or (r.ed == ed and r.sclen == sclen and r.pos > pos)) {	
				// less edit distance, then less soft clip, then smaller junction distance
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen, r.score);
			}
		}
	}

	void print() {
		fprintf(stderr, "Pos: %d\nEdit: %d\nSClen: %d\nIndel: %d\nRef Covlen: %d\nRead Covlen: %d\nScore: %d\n", pos, ed, sclen, indel, rcovlen, qcovlen, score);
	}
};

struct AlignCandid{
	int ed;
	int sclen;
	int indel;	// > 0: (ins) -> extra consumption on read, < 0: (del) -> fewer consumption on read
	int score;

	AlignCandid (int e, int s, int i) : ed(e), sclen(s), indel(i), score(-1*s - 2*e) {}
	AlignCandid (int e, int s, int i, int scr) : ed(e), sclen(s), indel(i), score(scr) {}

	bool operator < (const AlignCandid& r) const {
		// score = MatchedLen - 2 * EditDist ~ -1 * SCLen - 2 * EditDist
		if (score != r.score) 
			return score > r.score;
		if (ed != r.ed)
			return ed < r.ed;
		return abs(indel) < abs(r.indel);
	}

	void update(const AlignCandid& r) {
		if (r < *this) {
			set(r);
		}
	}

	void set(const AlignCandid& r) {
		ed = r.ed;
		sclen = r.sclen;
		indel = r.indel;
		score = r.score;
	}
};

class Alignment {
public:
	Alignment(void);
	virtual ~Alignment(void);
	void init(void);
	
	void hamming_distance(char* s, int n, char* t, int m);
	void hamming_distance_bottom(char* s, int n, char* t, int m, int max_sclen);
	void hamming_distance_top(char* s, int n, char* t, int m, int max_sclen);
	
	int global_one_side_banded_alignment(char* s, int n, char* t, int m, int w);
	
	// in the following are prefix to global alignment 
	// prefix on s
	// global on t
	virtual int local_alignment_right_sc(char* s, int n, char* t, int m, int& sc_len, int& indel) = 0;
	virtual int local_alignment_left_sc (char* s, int n, char* t, int m, int& sc_len, int& indel) = 0;

	int local_alignment_right(char* s, int n, char* t, int m, int& indel);
	int local_alignment_left (char* s, int n, char* t, int m, int& indel);

	void global_banded_alignment_drop(char* s, int n, char* t, int m, int w, int& on_s, int& on_t);
	void global_banded_alignment_reverse_drop(char* s, int n, char* t, int m, int w);

	int32_t get_score() { return align_score; }

protected:
	int hamming_distance(char* s, char* t, int n);

	bool hamming_match_right(char* s, int n, char* t, int m);
	bool hamming_match_left (char* s, int n, char* t, int m);

	int  hamming_distance_right(char* s, int n, char* t, int m, int& sc_len);
	int  hamming_distance_left (char* s, int n, char* t, int m, int& sc_len);

	int global_alignment(char* s, int n, char* t, int m, int gap_pen, int mm_pen);
	int global_alignment(char* s, int n, char* t, int m);
	int global_alignment_reverse(char* s, int n, char* t, int m);

	void global_banded_alignment(char* s, int n, char* t, int m, int w);
	void global_banded_alignment_reverse(char* s, int n, char* t, int m, int w);

protected:
	uint32_t dp[MAXSTRSIZE][MAXSTRSIZE];
	int32_t dpx[MAXSTRSIZE][MAXSTRSIZE];
	int hamm[MAXSTRSIZE][MAXSTRSIZE];
	int** diff_ch;
	int** score_ch;
	int32_t align_score;
};

class EditDistAlignment : public Alignment {
public:
	EditDistAlignment(void) {}
	~EditDistAlignment(void) {}

	// in the following are prefix to global alignment 
	// prefix on s
	// global on t
	virtual int local_alignment_right_sc(char* s, int n, char* t, int m, int& sc_len, int& indel);
	virtual int local_alignment_left_sc (char* s, int n, char* t, int m, int& sc_len, int& indel);

private:

};

class DropAlignment : public Alignment {
public:
	DropAlignment(void) {}
	~DropAlignment(void) {}

	// in the following are prefix to global alignment 
	// prefix on s
	// global on t
	virtual int local_alignment_right_sc(char* s, int n, char* t, int m, int& sc_len, int& indel);
	virtual int local_alignment_left_sc (char* s, int n, char* t, int m, int& sc_len, int& indel);

private:

};

class ScoreMatrix {
public:
	ScoreMatrix(void);
	~ScoreMatrix(void);

	void init(int mat, int mis, int ind, int xval);

	int mat_sc;
	int mis_sc;
	int ind_sc;
	int xd;

	int** diff_ch;
};

extern ScoreMatrix score_mat;
extern ScoreMatrix edit_mat;

#endif	//__ALIGN_H__
