#ifndef __EXTEND_H__
#define __EXTEND_H__

#include <map>
#include <vector>
#include <cstdio>
#include "common.h"
#include "align.h"

class TransExtension {
private:
	int thid;
	Alignment alignment;

	bool extend_right(const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, int ed_th, uint32_t ub, AlignRes& best_alignment);
	bool extend_left (const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, int ed_th, uint32_t lb, AlignRes& best_alignment);

	void extend_right_trans(uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
							int ed_th, uint32_t ub,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res);
	void extend_left_trans (uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
							int ed_th, uint32_t lb,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res);

	bool extend_right_middle(uint32_t pos, char* ref_seq, uint32_t exon_len, char* qseq, int qseq_len, 
							int ed_th, AlignRes& best, AlignRes& curr, AlignRes& exon_res);

	void extend_right_end(uint32_t pos, char* ref_seq, uint32_t ref_len, char* qseq, int qseq_len, 
						int ed_th, AlignRes& best, AlignRes& curr, AlignRes& exon_res);

	bool extend_left_middle(uint32_t pos, char* ref_seq, uint32_t exon_len, char* qseq, int qseq_len, 
							int ed_th, AlignRes& best, AlignRes& curr, AlignRes& exon_res);

	void extend_left_end(uint32_t pos, char* ref_seq, uint32_t ref_len, char* qseq, int qseq_len, 
						int ed_th, AlignRes& best, AlignRes& curr, AlignRes& exon_res);

public:
	TransExtension(void);
	TransExtension(int id);
	~TransExtension(void);
	void init(int id);

	bool extend_both_mates(const chain_t& lch, const chain_t& rch, const vector<uint32_t>& common_tid, char* lseq, char* rseq, 
							int lqspos, int rqspos, int lseq_len, int rseq_len, MatchedMate& lmm, MatchedMate& rmm);

	int extend_chain_both_sides(const chain_t& ch, char* seq, int seq_len, MatchedMate& mr, int dir);

	bool extend_chain_right(const vector <uint32_t>& common_tid, const chain_t& ch, char* seq, int seq_len, int ub, MatchedMate& mr, int& err);
	bool extend_chain_left (const vector <uint32_t>& common_tid, const chain_t& ch, char* seq, int32_t qspos, int lb, MatchedMate& mr, int& err);

	int calc_middle_ed(const chain_t& ch, int edth, char* qseq, int qseq_len);

};

#endif