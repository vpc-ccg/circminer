#ifndef __EXTEND_H__
#define __EXTEND_H__

#include <map>
#include <vector>
#include <cstdio>
#include "common.h"
#include "align.h"

bool extend_chain_right(const vector <uint32_t>& common_tid, const chain_t& ch, char* seq, int seq_len, int ub, MatchedMate& mr, int& err);
bool extend_chain_left (const vector <uint32_t>& common_tid, const chain_t& ch, char* seq, int32_t qspos, int lb, MatchedMate& mr, int& err);

bool extend_right(const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, int ed_th, uint32_t ub, AlignRes& best_alignment);
bool extend_left (const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, int ed_th, uint32_t lb, AlignRes& best_alignment);

void extend_right_trans(uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
						int ed_th, uint32_t ub,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res);
void extend_left_trans (uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
						int ed_th, uint32_t lb,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res);



#endif