#include "utils.h"
#include "common.h"
#include "align.h"
#include "gene_annotation.h"
#include "match_read.h"


void get_mate_name(char *fq1, char *fq2) {
    strcpy(fq2, fq1);
    int i = strlen(fq1) - 1;
    while (fq1[i--] != '.' and i >= 1);
    if (fq1[i] == '1')
        fq2[i] = '2';
    else if (fq1[i] == '2')
        fq2[i] = '1';
    else {
        fprintf(stderr, "Error: PE FASTQ names are not in the correct format\n");
        exit(1);
    }
}

void update_match_mate_info(bool lok, bool rok, int err, MatchedMate &mm) {
    mm.left_ok = lok and (mm.sclen_left <= maxSc);
    mm.right_ok = rok and (mm.sclen_right <= maxSc);
    if (lok and rok and (err <= maxEd) and (mm.sclen_right <= maxSc) and (mm.sclen_left <= maxSc)) {
        mm.is_concord = true;
        mm.type = CONCRD;
    } else if (lok or rok) {
        mm.type = CANDID;
    } else
        mm.type = ORPHAN;
}


int estimate_middle_error(const chain_t &ch) {
    int mid_err = 0;
    for (uint32_t i = 0; i < ch.chain_len - 1; i++) {
        if (ch.frags[i + 1].qpos > int32_t(ch.frags[i].qpos + ch.frags[i].len)) {
            int diff = (ch.frags[i + 1].rpos - ch.frags[i].rpos) - (ch.frags[i + 1].qpos - ch.frags[i].qpos);
            if (diff == 0)
                mid_err++;
            else if (diff > 0 and diff <= bandWidth)
                mid_err += diff;
            else if (diff < 0 and diff >= (-1 * bandWidth))
                mid_err -= diff;
        }
    }
    return mid_err;
}


// calculate tlen including sm.epos and lm.spos
int calc_tlen(const MatchedMate &sm, const MatchedMate &lm, int &intron_num) {
    const IntervalInfo<UniqSeg> *this_region;
    uint32_t tid;
    int start_ind;
    uint32_t start_table_ind;
    uint32_t end_table_ind;
    int tlen;
    int min_tlen = INF;
    int this_it_ind;
    int in;

    for (unsigned int i = 0; i < sm.exons_epos->seg_list.size(); i++) {
        for (unsigned int j = 0; j < sm.exons_epos->seg_list[i].trans_id.size(); j++) {
            tid = sm.exons_epos->seg_list[i].trans_id[j];
            start_ind = gtf_parser.get_trans_start_ind(contigNum, tid);
            start_table_ind = sm.exon_ind_epos - start_ind;
            if (start_table_ind < 0)    // assert
                continue;

            end_table_ind = lm.exon_ind_spos - start_ind;
            // transcript does not contain lm exon
            if (lm.exon_ind_spos < start_ind or
                end_table_ind >= gtf_parser.trans2seg[contigNum][tid].size() or
                gtf_parser.trans2seg[contigNum][tid][end_table_ind] == 0)

                continue;

            if (start_table_ind == end_table_ind) {
                in = 0;
                tlen = lm.spos - sm.epos + 1;
            } else {
                bool pre_zero = false;
                in = 0;
                tlen = sm.exons_epos->epos - sm.epos + 1;
                this_it_ind = sm.exon_ind_epos;
                for (uint32_t k = start_table_ind + 1; k < end_table_ind; k++) {
                    this_it_ind++;
                    if (gtf_parser.trans2seg[contigNum][tid][k] != 0) {
                        this_region = gtf_parser.get_interval(this_it_ind);
                        tlen += this_region->epos - this_region->spos + 1;
                        pre_zero = false;
                    } else {
                        if (!pre_zero)
                            in++;
                        pre_zero = true;
                    }
                }
                tlen += lm.spos - lm.exons_spos->spos + 1;
            }

            //fprintf(stdout, "tr[%d]: %s\ttlen: %d\tintrons: %d\n", tid, gtf_parser.transcript_ids[contigNum][tid].c_str(), tlen, in);

            if (tlen < min_tlen) {
                intron_num = in;
                min_tlen = tlen;
            }
        }
    }

    return (min_tlen == INF) ? -1 : min_tlen + sm.matched_len - 1 + lm.matched_len - 1;
}


bool is_concord(const chain_t &a, uint32_t seq_len, MatchedMate &mr) {
    if (a.chain_len < 2) {
        mr.is_concord = false;
    } else if ((a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len - a.frags[0].qpos) >= seq_len) {
        mr.is_concord = true;
        mr.type = CONCRD;

        mr.spos = a.frags[0].rpos;
        mr.epos = a.frags[a.chain_len - 1].rpos + a.frags[a.chain_len - 1].len - 1;
        mr.matched_len = a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len - a.frags[0].qpos;
        mr.qspos = a.frags[0].qpos;
        mr.qepos = a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len - 1;
    } else {
        mr.is_concord = false;
    }
    return mr.is_concord;
}

bool is_concord2(const chain_t &a, uint32_t seq_len, MatchedMate &mr) {
    if (a.chain_len < 2) {
        mr.is_concord = false;
    } else if ((a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len - a.frags[0].qpos) >= seq_len) {
        mr.is_concord = true;
        mr.type = CONCRD;

        mr.spos = a.frags[0].rpos;
        mr.epos = a.frags[a.chain_len - 1].rpos + a.frags[a.chain_len - 1].len - 1;
        mr.matched_len = a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len - a.frags[0].qpos;
        mr.qspos = a.frags[0].qpos;
        mr.qepos = a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len - 1;
    } else {
        mr.is_concord = false;
        if (a.frags[0].qpos == 0 or a.frags[a.chain_len - 1].qpos + a.frags[a.chain_len - 1].len == seq_len) {
            mr.type = CANDID;
        }
    }
    return mr.is_concord;
}


// sm should start before lm
bool
concordant_explanation(const MatchedMate &sm, const MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift,
                       bool r1_sm, int pair_type) {
    if (sm.spos > lm.spos)
        return false;

    int32_t tlen;
    bool on_cdna;
    on_cdna =
            (sm.exons_spos != NULL) and (sm.exons_epos != NULL) and (lm.exons_spos != NULL) and (lm.exons_epos != NULL);
    if (sm.exons_spos == NULL or lm.exons_spos == NULL) {
        tlen = lm.spos - sm.epos - 1 + lm.matched_len + sm.matched_len;
        if (tlen <= maxTlen)
            mr.update(sm, lm, chr, shift, tlen, 0, false, CONGNM, r1_sm);
        else if (tlen <= MAXDISCRDTLEN)
            mr.update(sm, lm, chr, shift, tlen, 0, false, CONGNM, r1_sm);
    } else {
        //fprintf(stderr, "Left Mate [%u-%u] dir=%d, type=%d, Right Mate[%u-%u] dir=%d, type=%d\n", sm.spos, sm.epos, sm.dir, sm.type, lm.spos, lm.epos, lm.dir, lm.type);
        // starts on same exon
        for (unsigned int i = 0; i < sm.exons_spos->seg_list.size(); i++)
            for (unsigned int j = 0; j < lm.exons_spos->seg_list.size(); j++)
                if (sm.exons_spos->seg_list[i].same_exon(lm.exons_spos->seg_list[j])) {
                    // => assume genomic locations
                    tlen = lm.spos + lm.matched_len - sm.spos;
                    if (tlen <= maxTlen)
                        mr.update(sm, lm, chr, shift, tlen, 0, on_cdna, ((pair_type == 0) ? CONCRD : CONGEN), r1_sm);
                    else
                        mr.update(sm, lm, chr, shift, tlen, 0, on_cdna, DISCRD, r1_sm);
                }
    }

    if (sm.exons_epos == NULL or lm.exons_spos == NULL) {
        tlen = lm.spos - sm.epos - 1 + sm.matched_len + lm.matched_len;
        if (tlen <= maxTlen)
            mr.update(sm, lm, chr, shift, tlen, 0, false, CONGNM, r1_sm);
        else if (tlen <= MAXDISCRDTLEN)
            mr.update(sm, lm, chr, shift, tlen, 0, false, CONGNM, r1_sm);
    } else {
        int intron_num;
        tlen = calc_tlen(sm, lm, intron_num);
        //fprintf(stdout, "tlen: %d\n", tlen);
        if (tlen >= 0 and tlen <= maxTlen) {
            mr.update(sm, lm, chr, shift, tlen, intron_num, on_cdna, ((pair_type == 0) ? CONCRD : CONGEN), r1_sm);
        } else {
            if (tlen < 0) {
                tlen = lm.spos - sm.epos - 1 + sm.matched_len + lm.matched_len;
                intron_num = 0;
            }
            mr.update(sm, lm, chr, shift, tlen, intron_num, on_cdna, DISCRD, r1_sm);
        }
    }

    if (mr.type == CONCRD)
        return true;
    else
        return false;
}

bool check_chimeric(const MatchedMate &sm, const MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift,
                    bool r1_sm) {
    if (mr.type == CONCRD)
        return false;

    if (sm.exons_spos == NULL or lm.exons_spos == NULL)
        return false;

    //fprintf(stderr, "Left Mate [%u-%u] dir=%d, type=%d, Right Mate[%u-%u] dir=%d, type=%d\n", sm.spos, sm.epos, sm.dir, sm.type, lm.spos, lm.epos, lm.dir, lm.type);
    for (unsigned int i = 0; i < sm.exons_spos->seg_list.size(); i++)
        for (unsigned int j = 0; j < lm.exons_spos->seg_list.size(); j++)
            if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j]) and sm.spos < lm.spos) {
                mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIORF, r1_sm);
                return true;
            }
    return false;
}

// check for one valid split mate
// -> one fully mapped mate and one partially mapped
bool check_bsj(MatchedMate &sm, MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift, bool r1_sm) {
    if (mr.type == CONCRD or mr.type == DISCRD)
        return false;

    if ((!sm.right_ok) or (!lm.left_ok))
        return false;

    if (sm.exons_spos == NULL or lm.exons_spos == NULL) {
        if ((sm.exons_spos != NULL and same_gene(sm, lm)) or (lm.exons_spos != NULL and same_gene(lm, sm))) {
            mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ, r1_sm);
            return true;
        }

        // checking for ciRNA
        //fprintf(stderr, "R1 start ind: %d\tR2 end ind: %d\n To beg of intron: %d\n", sm.exon_ind_spos, lm.exon_ind_epos, sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos));
        if ((intronic_bs[contigNum][sm.spos]) and ((intronic_bs[contigNum][lm.spos])) and
            (sm.exon_ind_spos >= 0) and (lm.exon_ind_epos >= 0) and (sm.exon_ind_spos == lm.exon_ind_epos) and
            (sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos) <= LARIAT2BEGTH)) {
            mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ, r1_sm);
            return true;
        }
        return false;
    }

    for (unsigned int i = 0; i < sm.exons_spos->seg_list.size(); i++)
        for (unsigned int j = 0; j < lm.exons_spos->seg_list.size(); j++)
            if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j])) {
                mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ, r1_sm);
                return true;
            }
    return false;
}

// check for two valid split mates
// -> two partially mapped mates
bool check_2bsj(MatchedMate &sm, MatchedMate &lm, MatchedRead &mr, const string &chr, uint32_t shift, bool r1_sm) {
    if (mr.type < CHI2BSJ)
        return false;

    if (sm.spos > lm.spos)
        return false;

    // ...<----
    // ...-->
    if (sm.right_ok and lm.right_ok and (sm.spos != lm.spos))
        return false;

    //  <--...
    // --->...
    if (sm.left_ok and lm.left_ok and (sm.epos != lm.epos))
        return false;

    // ...<---   --->...
    // OR
    // ...-->     <--...
    //if sm.right_ok and lm.left_ok -> continue

    // o.w.
    if (sm.left_ok and lm.right_ok)
        return false;

    if (sm.exons_spos == NULL or lm.exons_spos == NULL) {
        if ((sm.exons_spos != NULL and same_gene(sm, lm)) or (lm.exons_spos != NULL and same_gene(lm, sm))) {
            mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHI2BSJ, r1_sm);
            return true;
        }

        // checking for ciRNA
        //fprintf(stderr, "R1 start ind: %d\tR2 end ind: %d\n To beg of intron: %d\n", sm.exon_ind_spos, lm.exon_ind_epos, sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos));
        if ((intronic_bs[contigNum][sm.spos]) and ((intronic_bs[contigNum][lm.spos])) and
            (sm.exon_ind_spos >= 0) and (lm.exon_ind_epos >= 0) and (sm.exon_ind_spos == lm.exon_ind_epos) and
            (sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos) <= LARIAT2BEGTH)) {
            mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHI2BSJ, r1_sm);
            return true;
        }
        return false;
    }

    for (unsigned int i = 0; i < sm.exons_spos->seg_list.size(); i++)
        for (unsigned int j = 0; j < lm.exons_spos->seg_list.size(); j++)
            if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j])) {
                mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHI2BSJ, r1_sm);
                return true;
            }
    return false;
}

void intersect_trans(const vector <uint32_t> &tid_l1, const vector <uint32_t> &tid_l2, vector <uint32_t> &common_tid) {
    for (unsigned int i = 0; i < tid_l1.size(); i++) {
        for (unsigned int j = 0; j < tid_l2.size(); j++) {
            uint32_t tid1 = tid_l1[i];
            uint32_t tid2 = tid_l2[j];
            if (tid1 == tid2) {
                common_tid.push_back(tid1);
                break;
            }
        }
    }
}

bool same_transcript(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r, vector <uint32_t> &common_tid) {
    common_tid.clear();
    if (s == NULL or r == NULL)
        return false;

    vector <uint32_t> seg1_tid;
    vector <uint32_t> seg2_tid;

    for (unsigned int i = 0; i < s->seg_list.size(); i++)
        for (unsigned int k = 0; k < s->seg_list[i].trans_id.size(); k++)
            seg1_tid.push_back(s->seg_list[i].trans_id[k]);

    for (unsigned int j = 0; j < r->seg_list.size(); j++)
        for (unsigned int l = 0; l < r->seg_list[j].trans_id.size(); l++)
            seg2_tid.push_back(r->seg_list[j].trans_id[l]);

    intersect_trans(seg1_tid, seg2_tid, common_tid);

    return common_tid.size() != 0;
}

bool same_transcript(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r, const IntervalInfo<UniqSeg> *q,
                     vector <uint32_t> &common_tid) {
    common_tid.clear();
    if (s == NULL or r == NULL or q == NULL)
        return false;

    vector <uint32_t> sr_common_tid;
    bool sr_intersect = same_transcript(s, r, sr_common_tid);
    if (!sr_intersect)
        return false;

    vector <uint32_t> seg_tid;

    for (unsigned int i = 0; i < s->seg_list.size(); i++)
        for (unsigned int k = 0; k < s->seg_list[i].trans_id.size(); k++)
            seg_tid.push_back(s->seg_list[i].trans_id[k]);

    intersect_trans(sr_common_tid, seg_tid, common_tid);

    return common_tid.size() != 0;
}

bool same_transcript(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r,
                     const IntervalInfo<UniqSeg> *q, const IntervalInfo<UniqSeg> *p, vector <uint32_t> &common_tid) {

    common_tid.clear();
    if (s == NULL or r == NULL or q == NULL or p == NULL)
        return false;

    vector <uint32_t> sr_common_tid;
    bool sr_intersect = same_transcript(s, r, sr_common_tid);
    if (!sr_intersect)
        return false;

    vector <uint32_t> qp_common_tid;
    bool qp_intersect = same_transcript(q, p, qp_common_tid);
    if (!qp_intersect)
        return false;

    intersect_trans(sr_common_tid, qp_common_tid, common_tid);

    return common_tid.size() != 0;
}

bool same_transcript(vector <MatchedMate> &segments, vector <uint32_t> &common_tid) {
    common_tid.clear();

    vector<const IntervalInfo<UniqSeg> *> mpos_ol; // mid_pos_overlap
    for (unsigned int i = 0; i < segments.size(); ++i) {
        const IntervalInfo<UniqSeg> *s = overlap_to_mpos(segments[i]);
        mpos_ol.push_back(s);
    }

    if (segments.size() == 4)
        return same_transcript(mpos_ol[0], mpos_ol[1], mpos_ol[2], mpos_ol[3], common_tid);
    else if (segments.size() == 3)
        return same_transcript(mpos_ol[0], mpos_ol[1], mpos_ol[2], common_tid);
    else if (segments.size() == 2)
        return same_transcript(mpos_ol[0], mpos_ol[1], common_tid);
    return false;

}

bool same_transcript(vector <MatchedMate> &segments, int size, vector <uint32_t> &common_tid) {
    bool success;
    if (size == 2) {
        overlap_to_spos(segments[0]);
        overlap_to_spos(segments[1]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[1]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[0]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, common_tid);
        if (success)
            return success;
    }

    if (size == 3) {
        overlap_to_spos(segments[0]);
        overlap_to_spos(segments[1]);
        overlap_to_spos(segments[2]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, segments[2].exons_spos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[2]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, segments[2].exons_epos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[1]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, segments[2].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, segments[2].exons_epos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[0]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, segments[2].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, segments[2].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, segments[2].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, segments[2].exons_epos, common_tid);
        if (success)
            return success;

    }

    if (size == 4) {
        overlap_to_spos(segments[0]);
        overlap_to_spos(segments[1]);
        overlap_to_spos(segments[2]);
        overlap_to_spos(segments[3]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, segments[2].exons_spos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[2]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, segments[2].exons_epos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[1]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, segments[2].exons_spos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, segments[2].exons_epos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[0]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, segments[2].exons_spos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, segments[2].exons_epos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, segments[2].exons_spos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, segments[2].exons_epos,
                                  segments[3].exons_spos, common_tid);
        if (success)
            return success;

        overlap_to_epos(segments[3]);
        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, segments[2].exons_spos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_spos, segments[2].exons_epos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, segments[2].exons_spos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_spos, segments[1].exons_epos, segments[2].exons_epos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, segments[2].exons_spos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_spos, segments[2].exons_epos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, segments[2].exons_spos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

        common_tid.clear();
        success = same_transcript(segments[0].exons_epos, segments[1].exons_epos, segments[2].exons_epos,
                                  segments[3].exons_epos, common_tid);
        if (success)
            return success;

    }
    return false;
}

bool same_gene(const IntervalInfo<UniqSeg> *s, const IntervalInfo<UniqSeg> *r) {
    if (s == NULL or r == NULL)
        return false;

    for (unsigned int i = 0; i < s->seg_list.size(); i++)
        for (unsigned int j = 0; j < r->seg_list.size(); j++)
            if (s->seg_list[i].gene_id == r->seg_list[j].gene_id)
                return true;

    return false;
}

bool same_gene(const IntervalInfo<UniqSeg> *mate, uint32_t s, uint32_t e) {
    GeneInfo *ginfo;
    for (unsigned int i = 0; i < mate->seg_list.size(); i++) {
        ginfo = gtf_parser.get_gene_info(mate->seg_list[i].gene_id);
        if (ginfo->start <= s and e <= ginfo->end)
            return true;
    }


    return false;
}

bool same_gene(const MatchedMate &mm, const MatchedMate &other) {
    GeneInfo *ginfo;
    for (unsigned int i = 0; i < mm.exons_spos->seg_list.size(); i++) {
        ginfo = gtf_parser.get_gene_info(mm.exons_spos->seg_list[i].gene_id);
        //fprintf(stderr, "Gene[%d][%s]: [%d - %d], [%d - %d]\n", i, mm.exons_spos->seg_list[i].gene_id.c_str(), ginfo->start, ginfo->end, other.spos, other.epos);
        if (ginfo->start <= other.spos and other.epos <= ginfo->end)
            return true;
    }

    return false;
}

bool same_gene(uint32_t sme, const IntervalInfo<GeneInfo> *smg, uint32_t lms, const IntervalInfo<GeneInfo> *lmg) {
    if (smg == NULL or lmg == NULL)
        return false;

    if (smg->seg_list.size() == 0 or lmg->seg_list.size() == 0)
        return false;

    bool same_intron;
    int step = 10;
    for (unsigned int i = 0; i < smg->seg_list.size(); i++)
        for (unsigned int j = 0; j < lmg->seg_list.size(); j++)
            if (smg->seg_list[i].start == lmg->seg_list[j].start and smg->seg_list[i].end == lmg->seg_list[j].end) {
                same_intron = true;
                for (uint32_t k = sme; k <= lms; k += step)
                    if (!(intronic_bs[contigNum][k])) {
                        same_intron = false;
                        break;
                    }
                if (same_intron)
                    return true;
            }

    return false;
}


void overlap_to_epos(MatchedMate &mr) {
    if (mr.looked_up_epos or mr.exons_epos != NULL)
        return;

    mr.exons_epos = gtf_parser.get_location_overlap_ind(mr.epos, false, mr.exon_ind_epos);
    mr.looked_up_epos = true;
    //fprintf(stdout, "End Seg list size: %d\n", mr.exons_epos->seg_list.size());
}

void overlap_to_spos(MatchedMate &mr) {
    if (mr.looked_up_spos or mr.exons_spos != NULL)
        return;

    mr.exons_spos = gtf_parser.get_location_overlap_ind(mr.spos, false, mr.exon_ind_spos);
    mr.looked_up_spos = true;
    //fprintf(stdout, "Start Seg list size: %d\n", mr.exons_spos->seg_list.size());
}

const IntervalInfo<UniqSeg> *overlap_to_mpos(MatchedMate &mr) {
    int ind;
    return gtf_parser.get_location_overlap_ind((mr.spos + mr.epos) / 2, false, ind);
}

void gene_overlap(MatchedMate &mr) {
    if (mr.looked_up_gene or mr.gene_info != NULL)
        return;
    mr.gene_info = gtf_parser.get_gene_overlap(mr.spos, false);
    mr.looked_up_gene = true;
}

void get_junctions(MatchedMate &mm) {
    overlap_to_spos(mm);
    overlap_to_epos(mm);

    mm.junc_info.clear();

    if (mm.exons_spos == NULL or mm.exons_epos == NULL)
        return;

    for (unsigned int i = 0; i < mm.exons_spos->seg_list.size(); ++i) {
        for (unsigned int j = 0; j < mm.exons_spos->seg_list[i].trans_id.size(); ++j) {
            uint32_t covered = 0;
            uint32_t tid = mm.exons_spos->seg_list[i].trans_id[j];
            int start_ind = gtf_parser.get_trans_start_ind(contigNum, tid);
            uint32_t start_table_ind = mm.exon_ind_spos - start_ind;
            if (start_table_ind < 0)    // assert
                continue;

            uint32_t end_table_ind = mm.exon_ind_epos - start_ind;

            // transcript does not contain lm exon
            if (mm.exon_ind_epos < start_ind or
                end_table_ind >= gtf_parser.trans2seg[contigNum][tid].size() or
                gtf_parser.trans2seg[contigNum][tid][end_table_ind] == 0)

                continue;

            // no junction
            if (start_table_ind == end_table_ind) {
                return;
            }

            uint32_t junc_start = mm.exons_spos->epos;
            const IntervalInfo<UniqSeg> *this_region;
            covered = mm.exons_spos->epos - mm.spos + 1;
            int this_it_ind = mm.exon_ind_spos;
            for (uint32_t k = start_table_ind + 1; k < end_table_ind; k++) {
                this_it_ind++;
                if (gtf_parser.trans2seg[contigNum][tid][k] != 0) {
                    this_region = gtf_parser.get_interval(this_it_ind);

                    //fprintf(stderr, "[%u-%u]\n", junc_start, this_region->spos);

                    mm.junc_info.push_back(junc_start, this_region->spos, covered);
                    covered += this_region->epos - this_region->spos + 1;
                    junc_start = this_region->epos;
                }
            }
            mm.junc_info.push_back(junc_start, mm.exons_epos->spos, covered);
            covered += mm.epos - mm.exons_epos->spos + 1;

            //fprintf(stderr, "While building: Covered = %d\n", covered);
            //fprintf(stderr, "on read: [%d-%d] matched_len: %d\n", mm.qspos, mm.qepos, mm.matched_len);
            //mm.junc_info.print();
            if (abs(static_cast<int32_t> (covered - mm.matched_len)) <= INDELTH)
                return;
            else
                mm.junc_info.clear();
        }
    }
}

string get_consensus(const string &s1, const string &s2) {
    string res = "";
    if (s1.length() != s2.length())
        return res;

    for (unsigned int i = 0; i < s1.length(); ++i) {
        res += (s1[i] == s2[i]) ? s1[i] : 'N';
    }

    return res;
}

string get_consensus(const vector <string> &vseq) {
    string res = "";

    if (vseq.size() == 0)
        return res;

    for (unsigned int i = 1; i < vseq.size(); ++i) {
        if (vseq[i].length() != vseq[i - 1].length())
            return res;
    }

    unsigned int counts[ASCISIZE];
    char nuc[4] = {'A', 'C', 'G', 'T'};
    for (unsigned int i = 0; i < vseq[0].length(); ++i) {
        memset(counts, 0, sizeof(int) * ASCISIZE);
        for (unsigned int j = 0; j < vseq.size(); ++j)
            ++counts[static_cast<uint8_t> (vseq[j][i])];

        counts[static_cast<uint8_t> ('A')] += counts[static_cast<uint8_t> ('a')];
        counts[static_cast<uint8_t> ('C')] += counts[static_cast<uint8_t> ('c')];
        counts[static_cast<uint8_t> ('G')] += counts[static_cast<uint8_t> ('g')];
        counts[static_cast<uint8_t> ('T')] += counts[static_cast<uint8_t> ('t')];
        counts[static_cast<uint8_t> ('N')] += counts[static_cast<uint8_t> ('n')];

        counts[static_cast<uint8_t> ('a')] = 0;
        counts[static_cast<uint8_t> ('c')] = 0;
        counts[static_cast<uint8_t> ('g')] = 0;
        counts[static_cast<uint8_t> ('t')] = 0;
        counts[static_cast<uint8_t> ('n')] = 0;

        unsigned int max_cnt = 0;
        char ch = 'N';
        for (int k = 0; k < 4; ++k) {
            if (counts[static_cast<uint8_t> (nuc[k])] > max_cnt) {
                max_cnt = counts[static_cast<uint8_t> (nuc[k])];
                ch = nuc[k];
            }
        }

        if (max_cnt >= (vseq.size() / 2))
            res += ch;
        else
            res += 'N';
    }

    return res;
}

void reverse_str(char *s, int n, char *revs) {
    for (int i = 0; i < n; ++i)
        revs[i] = s[n - i - 1];

    revs[n] = '\0';
}

// is a on the left side?
bool is_left_chain(chain_t a, chain_t b, int read_length) {
    uint32_t a_beg = a.frags[0].rpos;
    uint32_t b_beg = b.frags[0].rpos;
    uint32_t a_end = a.frags[a.chain_len - 1].rpos + a.frags[a.chain_len - 1].len - 1;
    uint32_t b_end = b.frags[b.chain_len - 1].rpos + b.frags[b.chain_len - 1].len - 1;

    // check whether they are overlapping
    bool non_overlaping = (b_beg > a_end) or (a_beg > b_end);

    if (non_overlaping) {
// 		fprintf(stderr, "NOV\n");
        return a_beg < b_beg;
    } else {
        uint32_t i = 0, j = 0;
        int best_distance = INF;
        int best_i = -1;
        int best_j = -1;
        while (i < a.chain_len and j < b.chain_len) {
            uint32_t bj_beg = b.frags[j].rpos;
            uint32_t ai_end = a.frags[i].rpos + a.frags[i].len - 1;
            if (ai_end < bj_beg) {
                int distance = bj_beg - ai_end;
                if (distance < best_distance) {
                    best_distance = distance;
                    best_i = i;
                    best_j = j;
                }
                ++i;
                continue;
            }
            uint32_t ai_beg = a.frags[i].rpos;
            uint32_t bj_end = b.frags[j].rpos + b.frags[j].len - 1;
            if (bj_end < ai_beg) {
                int distance = ai_beg - bj_end;
                if (distance < best_distance) {
                    best_distance = distance;
                    best_i = i;
                    best_j = j;
                }
                ++j;
                continue;
            }
            best_i = i;
            best_j = j;
            break;
        }


        uint32_t common_bp = maxM(a.frags[best_i].rpos, b.frags[best_j].rpos);
        int32_t a_ov_qpos = a.frags[best_i].qpos + (common_bp - a.frags[best_i].rpos);
        int32_t b_ov_qpos = b.frags[best_j].qpos + (common_bp - b.frags[best_j].rpos);

// 			fprintf(stderr, "OV -> Decision made\n");
        if (a_ov_qpos < read_length and b_ov_qpos < read_length)
            return a_ov_qpos >= b_ov_qpos;

// 		fprintf(stderr, "OV -> Ambiguous\n");

        return a_beg < b_beg;
    }
}

void remove_side_introns(MatchedMate &mm, int rlen) {
    overlap_to_spos(mm);
    if (mm.exons_spos == NULL) {
        const IntervalInfo<UniqSeg> *it_seg = gtf_parser.get_interval(mm.exon_ind_spos + 1);
        if (it_seg == NULL) {
// 			fprintf(stderr, "Failed to remove intron left\n");
            return;
        }

        int diff = it_seg->spos - mm.spos;
// 		fprintf(stderr, "left intron retention: %d bp\n", diff);
        if (diff > 0 and (static_cast<uint32_t> (diff)) < mm.matched_len) {
            // update mm
            mm.spos = it_seg->spos;
            mm.qspos += diff;
            mm.matched_len -= diff;

            if ((mm.qspos - 1) < (rlen - mm.qepos)) // left-side matched
                mm.sclen_left += diff;

            mm.looked_up_spos = false;
            mm.exons_spos = NULL;
// 			fprintf(stderr, "Removed left intron retention: %d bp\n", diff);
        }
    }

    overlap_to_epos(mm);
    if (mm.exons_epos == NULL) {
        const IntervalInfo<UniqSeg> *it_seg = gtf_parser.get_interval(mm.exon_ind_epos);
        if (it_seg == NULL) {
// 			fprintf(stderr, "Failed to remove intron right\n");
            return;
        }

        int diff = mm.epos - it_seg->epos;
// 		fprintf(stderr, "right intron retention: %d bp\n", diff);
        if (diff > 0 and (static_cast<uint32_t> (diff)) < mm.matched_len) {
            // update mm
            mm.epos = it_seg->epos;
            mm.qepos -= diff;
            mm.matched_len -= diff;

            if ((mm.qspos - 1) > (rlen - mm.qepos)) // right-side matched
                mm.sclen_right += diff;

            mm.looked_up_epos = false;
            mm.exons_epos = NULL;
// 			fprintf(stderr, "Removed right intron retention: %d bp\n", diff);
        }
    }
}

