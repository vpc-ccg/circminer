#ifndef __EXTEND_H__
#define __EXTEND_H__

#include <cstdio>
#include "common.h"
#include "align.h"

bool extend_right(char* seq, uint32_t& pos, int len, uint32_t ub, AlignRes& best_alignment);
bool extend_right(const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, uint32_t ub, AlignRes& best_alignment);

bool extend_left(char* seq, uint32_t& pos, int len, uint32_t lb, AlignRes& best_alignment);
bool extend_left(const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, uint32_t lb, AlignRes& best_alignment);

#endif