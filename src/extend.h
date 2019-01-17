#ifndef __EXTEND_H__
#define __EXTEND_H__

#include <cstdio>
#include "common.h"

struct AlignRes {
	uint32_t pos;
	int ed;
	int sclen;
	int indel;
	int rcovlen;	// on ref
	int qcovlen;	// on read

	AlignRes (uint32_t p, int e, int s, int i, int qc) : pos(p), ed(e), sclen(s), indel(i), qcovlen(qc), rcovlen(qc - i) {}
	AlignRes (uint32_t p) : pos(p), ed(0), sclen(0), indel(0), rcovlen(0), qcovlen(0) {}

	void set (uint32_t p, int e, int s, int i, int qc) {
		this->pos = p;
		this->ed = e;
		this->sclen = s;
		this->indel = i;
		this->qcovlen = qc;
		this->rcovlen = qc - i;
	}

	void update (int edit_dist, int sclength, uint32_t newpos, int indel, int qcovlen) {
		this->pos = newpos;
		this->ed += edit_dist;
		this->sclen = sclength;
		this->indel += indel;
		this->qcovlen += qcovlen;
		this->rcovlen += qcovlen - indel;
	}

	void update_right (AlignRes& r) {
		if (r.qcovlen > qcovlen) {
			if (r.ed <= EDTH and r.sclen <= SOFTCLIPTH) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen);
			}
		}
		else if (r.qcovlen == qcovlen)
			if ((r.ed < ed) or (r.ed == ed and r.sclen < sclen) or (r.ed == ed and r.sclen == sclen and r.pos < pos)) {	
				// less hamming distance, then less soft clip, then smaller junction distance
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen);
			}
	}

	void update_left (AlignRes& r) {
		if (r.qcovlen > qcovlen) {
			if (r.ed <= EDTH and r.sclen <= SOFTCLIPTH) {
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen);
			}
		}
		else if (r.qcovlen == qcovlen) {
			if ((r.ed < ed) or (r.ed == ed and r.sclen < sclen) or (r.ed == ed and r.sclen == sclen and r.pos > pos)) {	
				// less hamming distance, then less soft clip, then smaller junction distance
				this->set(r.pos, r.ed, r.sclen, r.indel, r.qcovlen);
			}
		}
	}

	void print() {
		fprintf(stderr, "Pos: %d\nEdit: %d\nSClen: %d\nIndel: %d\nRef Covlen: %d\nRead Covlen: %d\n", pos, ed, sclen, indel, rcovlen, qcovlen);
	}
};

bool extend_right(char* seq, uint32_t& pos, int len, uint32_t ub, AlignRes& best_alignment);
bool extend_left(char* seq, uint32_t& pos, int len, uint32_t lb, AlignRes& best_alignment);

#endif