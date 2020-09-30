#ifndef __UTILS_H__
#define __UTILS_H__

#include "common.h"


void get_mate_name(char *fq1, char *fq2);

void update_match_mate_info(bool lok, bool rok, int err, MatchedMate &mm);

int estimate_middle_error(const chain_t &ch);

int calc_tlen(const MatchedMate &sm, const MatchedMate &lm, int &intron_num);

void get_junctions(MatchedMate &mm);

bool is_concord(const chain_t &a, uint32_t seq_len, MatchedMate &mr);
bool is_concord2(const chain_t &a, uint32_t seq_len, MatchedMate &mr);

bool concordant_explanation(const MatchedMate &sm, const MatchedMate &lm, MatchedRead &mr, const string &chr,
                            uint32_t shift, bool r1_sm, int pair_type);

bool check_chimeric(const MatchedMate &sm, const MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift,
                    bool r1_sm);

bool check_bsj(MatchedMate &sm, MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift, bool r1_sm);
bool check_2bsj(MatchedMate &sm, MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift, bool r1_sm);

void intersect_trans(const vector <uint32_t> &tid1, const vector <uint32_t> &tid2, vector <uint32_t> &common_tid);

bool same_transcript(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r, vector <uint32_t> &common_tid);

bool same_transcript(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r, const IntervalInfo<UniqSeg> *q,
                     vector <uint32_t> &common_tid);

bool same_transcript(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r,
                     const IntervalInfo<UniqSeg> *q, const IntervalInfo<UniqSeg> *p, vector <uint32_t> &common_tid);

bool same_transcript(vector <MatchedMate> &segments, vector <uint32_t> &common_tid);
bool same_transcript(vector <MatchedMate> &segments, int size, vector <uint32_t> &common_tid);

bool same_gene(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r);
bool same_gene(const IntervalInfo<UniqSeg> *mate, uint32_t s, uint32_t e);
bool same_gene(const MatchedMate &mm, const MatchedMate &other);
bool same_gene(uint32_t sme, const IntervalInfo<GeneInfo> *smg, uint32_t lms, const IntervalInfo<GeneInfo> *lmg);

void overlap_to_epos(MatchedMate &mr);
void overlap_to_spos(MatchedMate &mr);

const IntervalInfo<UniqSeg> *overlap_to_mpos(MatchedMate &mr);

void gene_overlap(MatchedMate &mr);

string get_consensus(const string &s1, const string &s2);
string get_consensus(const vector <string> &vseq);

void reverse_str(char *s, int n, char *revs);

bool is_left_chain(chain_t a, chain_t b, int read_length);

void remove_side_introns(MatchedMate &mm, int rlen);

#endif    // __UTILS_H__
