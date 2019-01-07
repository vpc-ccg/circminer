#ifndef __EXTEND_H__
#define __EXTEND_H__

#include <cstdio>
#include "common.h"

struct AlignRes {
	uint32_t pos;
	int ed;
	int sclen;
	int indel;
	int covlen;	// on ref

	AlignRes (uint32_t p, int e, int s, int i, int c) : pos(p), ed(e), sclen(s), indel(i), covlen(c) {}
	AlignRes (uint32_t p) : pos(p), ed(0), sclen(0), indel(0), covlen(0) {}

	void set (uint32_t p, int e, int s, int i, int c) {
		this->pos = p;
		this->ed = e;
		this->sclen = s;
		this->indel = i;
		this->covlen = c;
	}

	void update (int edit_dist, int sclength, uint32_t newpos, int indel, int covlen) {
		this->pos = newpos;
		this->ed += edit_dist;
		this->sclen = sclength;
		this->indel += indel;
		this->covlen += covlen;
	}

	void update_right (AlignRes& r) {
		if (r.covlen > covlen) {
			if (r.ed <= EDTH and r.sclen <= SOFTCLIPTH) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.covlen);
			}
		}
		else if (r.covlen == covlen)
			if ((r.ed < ed) or (r.ed == ed and r.sclen < sclen) or (r.ed == ed and r.sclen == sclen and r.pos < pos)) {	
				// less hamming distance, then less soft clip, then smaller junction distance
				this->set(r.pos, r.ed, r.sclen, r.indel, r.covlen);
			}
	}

	void update_left (AlignRes& r) {
		if (r.covlen > covlen) {
			if (r.ed <= EDTH and r.sclen <= SOFTCLIPTH) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.covlen);
			}
		}
		else if (r.covlen == covlen) {
			if ((r.ed < ed) or (r.ed == ed and r.sclen < sclen) or (r.ed == ed and r.sclen == sclen and r.pos > pos)) {	
				// less hamming distance, then less soft clip, then smaller junction distance
				this->set(r.pos, r.ed, r.sclen, r.indel, r.covlen);
			}
		}
	}

	void print() {
		fprintf(stderr, "Pos: %d\nEdit: %d\nSClen: %d\nIndel: %d\nCovlen: %d\n", pos, ed, sclen, indel, covlen);
	}
};

bool extend_right(char* seq, uint32_t& pos, int len, uint32_t ub, AlignRes& best_alignment);
bool extend_left(char* seq, uint32_t& pos, int len, uint32_t lb, AlignRes& best_alignment);

#endif