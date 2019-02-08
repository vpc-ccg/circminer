#ifndef __UTILS_H__
#define __UTILS_H__

#include "common.h"

void update_match_mate_info(bool lok, bool rok, int err, MatchedMate& mm);

int estimate_middle_error(const chain_t& ch);
int calc_middle_ed(const chain_t& ch, int edth, char* qseq, int qseq_len);
bool check_middle_ed(const chain_t& ch, int edth, char* qseq, int qseq_len);

int calc_tlen(const MatchedMate& sm, const MatchedMate& lm, int& intron_num);

bool is_concord(const chain_t& a, int seq_len, MatchedMate& mr);
bool concordant_explanation(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm);
bool check_chimeric(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm);
bool check_bsj(MatchedMate& sm, MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm);

bool same_transcript(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r, MatePair& mp);
bool same_gene(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r);
bool same_gene(const IntervalInfo<UniqSeg>* mate, uint32_t s, uint32_t e);
bool same_gene(const MatchedMate& mm, const MatchedMate& other);
bool same_gene(uint32_t sme, const IntervalInfo<GeneInfo>* smg, uint32_t lms, const IntervalInfo<GeneInfo>* lmg);

void overlap_to_epos(MatchedMate& mr);
void overlap_to_spos(MatchedMate& mr);
void gene_overlap(MatchedMate& mr);

#endif	// __UTILS_H__