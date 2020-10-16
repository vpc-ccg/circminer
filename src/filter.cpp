#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <set>

#include "filter.h"
#include "align.h"
#include "common.h"
#include "gene_annotation.h"
#include "extend.h"
#include "match_read.h"
#include "utils.h"

extern "C" {
#include "mrsfast/Common.h"
}


FilterRead::FilterRead(void) {
}

FilterRead::FilterRead(char *save_fname, bool pe, int round, bool first_round, bool last_round,
                       char *fq_file1, char *fq_file2, const vector <ContigLen> &chr_info) {

    init(save_fname, pe, round, first_round, last_round, fq_file1, fq_file2, chr_info);
}

FilterRead::~FilterRead(void) {
    // finalize();
}

void FilterRead::init(char *save_fname, bool pe, int round, bool first_round, bool last_round,
                      char *fq_file1, char *fq_file2, const vector <ContigLen> &chr_info) {

    char mode[2];
    sprintf(mode, "%s", (first_round) ? "w" : "a");

    sam_output.init(save_fname, mode, chr_info);
    extension = new TransExtension[threadCount];

    for (int i = 0; i < threadCount; ++i)
        extension[i].init(i, DROP_ALIGNMENT);

    this->is_pe = pe;
    this->first_round = first_round;
    this->last_round = last_round;

    // temp fastq file(s) to be read in next round
    char *temp_fname = (char *) malloc(FILE_NAME_MAX_LEN);
    char *wmode = (char *) malloc(FILE_NAME_MAX_LEN);
    sprintf(wmode, "%c", 'w');
    if (is_pe) {
        sprintf(temp_fname, "%s_%d_remain_R1.fastq", save_fname, round);
        temp_fq_r1 = open_file(temp_fname, wmode);

        sprintf(temp_fname, "%s_%d_remain_R2.fastq", save_fname, round);
        temp_fq_r2 = open_file(temp_fname, wmode);
    } else {
        sprintf(temp_fname, "%s_%d_remain.fastq", save_fname, round);
        temp_fq_r1 = open_file(temp_fname, wmode);
    }
    free(temp_fname);
    free(wmode);

    // updating fq file to be read in next round
    if (is_pe) {
        sprintf(fq_file1, "%s_%d_remain_R1.fastq", save_fname, round);
        sprintf(fq_file2, "%s_%d_remain_R2.fastq", save_fname, round);
    } else {
        sprintf(fq_file1, "%s_%d_remain.fastq", save_fname, round);
    }

}

void FilterRead::finalize(void) {
    sam_output.finalize();
    delete[] extension;

    close_file(temp_fq_r1);
    close_file(temp_fq_r2);
}

// SE mode
int FilterRead::process_read(int thid, Record *current_record, int kmer_size, GIMatchedKmer *fl, GIMatchedKmer *bl,
                             chain_list &forward_best_chain, chain_list &backward_best_chain) {

    vafprintf(1, stderr, "%s\n", current_record->rname);

    MatchedMate mr;
    int ex_ret, min_ret = ORPHAN;
    int fhh, bhh;

    get_best_chains(current_record->seq, current_record->seq_len, kmer_size, forward_best_chain, fl, fhh);
    for (int i = 0; i < forward_best_chain.best_chain_count; i++) {
        ex_ret = extension[thid].extend_chain_both_sides(forward_best_chain.chains[i], current_record->seq,
                                                         current_record->seq_len, mr, 1);

        if (ex_ret == CONCRD)
            return CONCRD;

        if (ex_ret < min_ret)
            min_ret = ex_ret;
    }

    get_best_chains(current_record->rcseq, current_record->seq_len, kmer_size, backward_best_chain, bl, bhh);
    for (int i = 0; i < backward_best_chain.best_chain_count; i++) {
        ex_ret = extension[thid].extend_chain_both_sides(backward_best_chain.chains[i], current_record->rcseq,
                                                         current_record->seq_len, mr, -1);

        if (ex_ret == CONCRD)
            return CONCRD;

        if (ex_ret < min_ret)
            min_ret = ex_ret;
    }

    return min_ret;

}

// PE mode
int
FilterRead::process_read(int thid, Record *current_record1, Record *current_record2, int kmer_size,
                         GIMatchedKmer *fl, GIMatchedKmer *bl,
                         chain_list &forward_best_chain_r1, chain_list &backward_best_chain_r1,
                         chain_list &forward_best_chain_r2, chain_list &backward_best_chain_r2) {

    int fhh_r1, bhh_r1;
    int fhh_r2, bhh_r2;

    // R1
    vafprintf(1, stderr, "R1/%s\n", current_record1->rname);

    get_best_chains(current_record1->seq, current_record1->seq_len, kmer_size, forward_best_chain_r1, fl, fhh_r1);
    get_best_chains(current_record1->rcseq, current_record1->seq_len, kmer_size, backward_best_chain_r1, bl, bhh_r1);

    vafprintf(1, stderr, "R1/%s\n", current_record1->rname);
    vafprintf(1, stderr, "R1 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r1.chains[0].score,
              (unsigned long) forward_best_chain_r1.best_chain_count);
    for (int j = 0; j < forward_best_chain_r1.best_chain_count; j++)
        for (unsigned int i = 0; i < forward_best_chain_r1.chains[j].chain_len; i++) {
            vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, forward_best_chain_r1.chains[j].frags[i].rpos,
                      forward_best_chain_r1.chains[j].frags[i].qpos, forward_best_chain_r1.chains[j].frags[i].len);
        }

    vafprintf(1, stderr, "R1 Reverse score:%.4f,\t len: %lu\n", backward_best_chain_r1.chains[0].score,
              (unsigned long) backward_best_chain_r1.best_chain_count);
    for (int j = 0; j < backward_best_chain_r1.best_chain_count; j++)
        for (unsigned int i = 0; i < backward_best_chain_r1.chains[j].chain_len; i++) {
            vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, backward_best_chain_r1.chains[j].frags[i].rpos,
                      backward_best_chain_r1.chains[j].frags[i].qpos, backward_best_chain_r1.chains[j].frags[i].len);
        }

    // R2
    vafprintf(1, stderr, "R2/%s\n", current_record2->rname);

    get_best_chains(current_record2->seq, current_record2->seq_len, kmer_size, forward_best_chain_r2, fl, fhh_r2);
    get_best_chains(current_record2->rcseq, current_record2->seq_len, kmer_size, backward_best_chain_r2, bl, bhh_r2);

    vafprintf(1, stderr, "R2/%s\n", current_record2->rname);
    vafprintf(1, stderr, "R2 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r2.chains[0].score,
              (unsigned long) forward_best_chain_r2.best_chain_count);
    for (int j = 0; j < forward_best_chain_r2.best_chain_count; j++)
        for (unsigned int i = 0; i < forward_best_chain_r2.chains[j].chain_len; i++) {
            vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, forward_best_chain_r2.chains[j].frags[i].rpos,
                      forward_best_chain_r2.chains[j].frags[i].qpos, forward_best_chain_r2.chains[j].frags[i].len);
        }

    vafprintf(1, stderr, "R2 Reverse score:%.4f,\t len: %lu\n", backward_best_chain_r2.chains[0].score,
              (unsigned long) backward_best_chain_r2.best_chain_count);
    for (int j = 0; j < backward_best_chain_r2.best_chain_count; j++)
        for (unsigned int i = 0; i < backward_best_chain_r2.chains[j].chain_len; i++) {
            vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, backward_best_chain_r2.chains[j].frags[i].rpos,
                      backward_best_chain_r2.chains[j].frags[i].qpos, backward_best_chain_r2.chains[j].frags[i].len);
        }

    // Orphan / OEA
    if (forward_best_chain_r1.best_chain_count + backward_best_chain_r1.best_chain_count +
        forward_best_chain_r2.best_chain_count + backward_best_chain_r2.best_chain_count <= 0) {
        if ((fhh_r1 + bhh_r1 > 0) and (fhh_r2 + bhh_r2 > 0)) {
            current_record1->mr->update_type(NOPROC_MANYHIT);
            return NOPROC_MANYHIT;
        } else {
            current_record1->mr->update_type(NOPROC_NOMATCH);
            return NOPROC_NOMATCH;
        }
    }
    if ((forward_best_chain_r1.best_chain_count + backward_best_chain_r1.best_chain_count <= 0) or
        (forward_best_chain_r2.best_chain_count + backward_best_chain_r2.best_chain_count <= 0)) {
        current_record1->mr->update_type(OEANCH);
        return OEANCH;
    }

    // checking read concordancy
    // any concordancy evidence/explanation ?
    float fc_score_r1 = (forward_best_chain_r1.best_chain_count > 0) ? forward_best_chain_r1.chains[0].score : 0;
    float bc_score_r1 = (backward_best_chain_r1.best_chain_count > 0) ? backward_best_chain_r1.chains[0].score : 0;
    float fc_score_r2 = (forward_best_chain_r2.best_chain_count > 0) ? forward_best_chain_r2.chains[0].score : 0;
    float bc_score_r2 = (backward_best_chain_r2.best_chain_count > 0) ? backward_best_chain_r2.chains[0].score : 0;
    vafprintf(2, stderr, "Scores: fc1=%f, bc1=%f, fc2=%f, bc2=%f\n", fc_score_r1, bc_score_r1, fc_score_r2,
              bc_score_r2);

    int attempt1, attempt2;
    if (fc_score_r1 + bc_score_r2 >= fc_score_r2 + bc_score_r1) {
        vafprintf(1, stderr, "Forward R1 / Backward R2\n");
        attempt1 = process_mates(thid, forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2,
                                 *(current_record1->mr), true);
        if (scanLevel == 0 and attempt1 == CONCRD) {
            return CONCRD;
        }

        vafprintf(1, stderr, "Backward R1 / Forward R2\n");
        attempt2 = process_mates(thid, forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1,
                                 *(current_record1->mr), false);
        if (scanLevel == 0 and attempt2 == CONCRD) {
            return CONCRD;
        }

        // return (attempt1 < attempt2) ? attempt1 : attempt2;
        return current_record1->mr->type;
    } else {
        vafprintf(1, stderr, "Backward R1 / Forward R2\n");
        attempt1 = process_mates(thid, forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1,
                                 *(current_record1->mr), false);
        if (scanLevel == 0 and attempt1 == CONCRD) {
            return CONCRD;
        }

        vafprintf(1, stderr, "Forward R1 / Backward R2\n");
        attempt2 = process_mates(thid, forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2,
                                 *(current_record1->mr), true);
        if (scanLevel == 0 and attempt2 == CONCRD) {
            return CONCRD;
        }

        // return (attempt1 < attempt2) ? attempt1 : attempt2;
        return current_record1->mr->type;
    }
}

// r1 is the forward mate and r2 is backward
int FilterRead::process_mates(int thid, const chain_list &forward_chain, const Record *forward_rec,
                              const chain_list &backward_chain, const Record *backward_rec, MatchedRead &mr,
                              bool r1_forward) {
    vector <MatePair> mate_pairs;

    bool forward_paired[maxChainLen];
    bool backward_paired[maxChainLen];
    pair_chains(forward_chain, backward_chain, mate_pairs, forward_paired, backward_paired, mr.type);

    int min_ret1 = ORPHAN;
    int min_ret2 = ORPHAN;

    bool r1_genic = false;
    bool r2_genic = false;

    // concordant?
    vafprintf(1, stderr, "#pairs = %d\n", mate_pairs.size());
    for (unsigned int i = 0; i < mate_pairs.size(); i++) {
        vafprintf(2, stderr, "Mate[%d]: %d, %d\n", i, mate_pairs[i].forward.frags[0].rpos,
                  mate_pairs[i].reverse.frags[0].rpos);
        // fprintf(stderr, "# common_tid: %d\n", mate_pairs[i].common_tid.size());
        // for (int k = 0; k < mate_pairs[i].common_tid.size(); k++)
        // 	fprintf(stderr, "%d, ", mate_pairs[i].common_tid[k]);
        // fprintf(stderr, "---\n");
        MatchedMate r1_mm;
        MatchedMate r2_mm;

        r1_mm.dir = 1;
        r2_mm.dir = -1;

        bool success;
// 		uint32_t forward_start = mate_pairs[i].forward.frags[0].rpos;
// 		uint32_t reverse_start = mate_pairs[i].reverse.frags[0].rpos;
// 		uint32_t reverse_end   = mate_pairs[i].reverse.frags[mate_pairs[i].reverse.chain_len-1].rpos + mate_pairs[i].reverse.frags[mate_pairs[i].reverse.chain_len-1].len - 1;

        bool is_forward_left = is_left_chain(mate_pairs[i].forward, mate_pairs[i].reverse, forward_rec->seq_len);

// 		if (forward_start <= reverse_end) {
        if (is_forward_left) {
            success = extension[thid].extend_both_mates(mate_pairs[i].forward, mate_pairs[i].reverse,
                                                        mate_pairs[i].common_tid, forward_rec->seq,
                                                        backward_rec->rcseq, 1, 1, forward_rec->seq_len,
                                                        backward_rec->seq_len, r1_mm, r2_mm);

            if (success) {
                ConShift con_shift = gtf_parser.get_shift(contigNum, r1_mm.spos);

                overlap_to_epos(r1_mm);
                overlap_to_spos(r1_mm);

                overlap_to_epos(r2_mm);
                overlap_to_spos(r2_mm);

                if (r1_mm.type == CONCRD and r2_mm.type == CONCRD) {
                    if (concordant_explanation(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward,
                                               mate_pairs[i].type) and scanLevel == 0) {
                        return CONCRD;
                    }
                }

                // potentially back splice junction?
                else if ((r1_mm.type == CANDID and r2_mm.type == CONCRD) or
                         (r1_mm.type == CONCRD and r2_mm.type == CANDID)) {
                    check_bsj(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward);
                } else if (r1_mm.type == CANDID and r2_mm.type == CANDID) {
                    check_2bsj(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward);
                }
            }
        }

// 		if (forward_start > reverse_start) {
        else {
            success = extension[thid].extend_both_mates(mate_pairs[i].reverse, mate_pairs[i].forward,
                                                        mate_pairs[i].common_tid, backward_rec->rcseq,
                                                        forward_rec->seq, 1, 1, backward_rec->seq_len,
                                                        forward_rec->seq_len, r2_mm, r1_mm);

            if (success) {
                ConShift con_shift = gtf_parser.get_shift(contigNum, r2_mm.spos);

                overlap_to_epos(r1_mm);
                overlap_to_spos(r1_mm);

                overlap_to_epos(r2_mm);
                overlap_to_spos(r2_mm);

                if (r1_mm.type == CONCRD and r2_mm.type == CONCRD) {
                    check_chimeric(r2_mm, r1_mm, mr, con_shift.contig, con_shift.shift, !r1_forward);
                }

                // potentially back splice junction?
                else if ((r1_mm.type == CANDID and r2_mm.type == CONCRD) or
                         (r1_mm.type == CONCRD and r2_mm.type == CANDID)) {
                    check_bsj(r2_mm, r1_mm, mr, con_shift.contig, con_shift.shift, !r1_forward);
                } else if (r1_mm.type == CANDID and r2_mm.type == CANDID) {
                    check_2bsj(r2_mm, r1_mm, mr, con_shift.contig, con_shift.shift, !r1_forward);
                }
            }
        }

        min_ret1 = minM(r1_mm.type, min_ret1);
        min_ret2 = minM(r2_mm.type, min_ret2);

        r1_genic = (r1_mm.exons_spos != NULL) or (r1_mm.exons_epos != NULL);
        r2_genic = (r2_mm.exons_spos != NULL) or (r2_mm.exons_epos != NULL);

    }

    if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ or mr.type == CHI2BSJ) {
        return mr.type;
    }

    MatchedMate mm1;
    int ex_ret;
    if (min_ret1 != CONCRD)
        for (int i = 0; i < forward_chain.best_chain_count; i++) {
            if (!forward_paired[i]) {
                ex_ret = extension[thid].extend_chain_both_sides(forward_chain.chains[i], forward_rec->seq,
                                                                 forward_rec->seq_len, mm1, 1);
                min_ret1 = minM(ex_ret, min_ret1);

                overlap_to_spos(mm1);
                overlap_to_epos(mm1);

                r1_genic = (mm1.exons_spos != NULL) or (mm1.exons_epos != NULL);
            }
        }

    MatchedMate mm2;
    if (min_ret2 != CONCRD)
        for (int i = 0; i < backward_chain.best_chain_count; i++) {
            if (!backward_paired[i]) {
                ex_ret = extension[thid].extend_chain_both_sides(backward_chain.chains[i], backward_rec->rcseq,
                                                                 backward_rec->seq_len, mm2, -1);
                min_ret2 = minM(ex_ret, min_ret2);

                overlap_to_spos(mm2);
                overlap_to_epos(mm2);

                r2_genic = (mm2.exons_spos != NULL) or (mm2.exons_epos != NULL);
            }
        }

    int new_type = (((min_ret1 == ORPHAN) and (min_ret2 == CONCRD)) or ((min_ret1 == CONCRD) and (min_ret2 == ORPHAN)))
                   ? OEANCH : ((min_ret1 == ORPHAN) or (min_ret2 == ORPHAN))
                   ? ORPHAN : ((min_ret1 == CONCRD) and (min_ret2 == CONCRD) and (r1_genic and r2_genic))
                   ? CHIFUS: ((min_ret1 == CONCRD) and (min_ret2 == CONCRD))
                   ? OEA2 : CANDID;

    mr.update_type(new_type);
    return mr.type;
}


// write reads SE mode
void FilterRead::write_read_category(Record *current_record, int state) {
    //state = minM(state, current_record->state);
    //int cat = (state >= cat_count) ? DISCRD : state;
    if (!last_round and state != CONCRD) {
        mutex_lock(&write_lock);

        fprintf(temp_fq_r1, "%s\n%s\n%s%d\n%s\n", current_record->rname, current_record->seq, current_record->comment,
                state, current_record->qual);

        mutex_unlock(&write_lock);
    }
}

// write reads PE mode
void FilterRead::write_read_category(Record *current_record1, Record *current_record2, const MatchedRead &mr) {
    char r1_dir = (mr.r1_forward) ? '+' : '-';
    char r2_dir = (mr.r2_forward) ? '+' : '-';

    mutex_lock(&write_lock);

    fprintf(temp_fq_r1, "@%s", current_record1->rname);
    fprintf(temp_fq_r2, "@%s", current_record2->rname);

    if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or
        mr.type == CHIBSJ or mr.type == CHI2BSJ or mr.type == CONGNM or mr.type == CONGEN) {

        string con = mr.chr_r1;
        uint32_t con_spos = mr.spos_r1;
        uint32_t con_epos = mr.epos_r1;
        gtf_parser.chrloc2conloc(con, con_spos, con_epos);

        uint64_t gspos = uint64_t(mr.contig_num) * DEF_CONTIG_SIZE + con_spos;

        fprintf(temp_fq_r1, " %" PRId64 " %d %s %u %u %d %u %u %c %d %s %u %u %d %u %u %c %d %d %d %d %d",
                gspos, mr.type,
                mr.chr_r1.c_str(), mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, r1_dir, mr.ed_r1,
                mr.chr_r2.c_str(), mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, r2_dir, mr.ed_r2,
                mr.tlen, mr.junc_num, mr.gm_compatible, mr.contig_num);
        fprintf(temp_fq_r2, " %" PRId64 " %d %s %u %u %d %u %u %c %d %s %u %u %d %u %u %c %d %d %d %d %d",
                gspos, mr.type,
                mr.chr_r1.c_str(), mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, r1_dir, mr.ed_r1,
                mr.chr_r2.c_str(), mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, r2_dir, mr.ed_r2,
                mr.tlen, mr.junc_num, mr.gm_compatible, mr.contig_num);
    } else {
        fprintf(temp_fq_r1, " * %d * * * * * * * * * * * * * * * * * * * *", mr.type);
        fprintf(temp_fq_r2, " * %d * * * * * * * * * * * * * * * * * * * *", mr.type);
    }

    char sep = '\n';
    fprintf(temp_fq_r1, "%c%s%c%s%c%s\n", sep, current_record1->seq, sep,
            current_record1->comment, sep, current_record1->qual);
    fprintf(temp_fq_r2, "%c%s%c%s%c%s\n", sep, current_record2->seq, sep,
            current_record2->comment, sep, current_record2->qual);

    mutex_unlock(&write_lock);

}

void FilterRead::print_sam(Record *rec) {
    sam_output.write_sam_rec_se(rec);
}

void FilterRead::print_sam(Record *rec1, Record *rec2) {
    sam_output.write_sam_rec_pe(rec1, rec2);
}

void FilterRead::print_pam(Record *rec1, Record *rec2) {
    sam_output.write_pam_rec_pe(rec1, rec2);
}

void
FilterRead::get_best_chains(char *read_seq, int seq_len, int kmer_size, chain_list &best_chain, GIMatchedKmer *frag_l,
                            int &high_hits) {
    // considering both overlapping and non-overlapping kmers
    int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer_size)) - 1;

    genome_seeder.split_match_hash(read_seq, seq_len, kmer_size, frag_l);
    chain_seeds_sorted_kbest(seq_len, frag_l, best_chain);

    high_hits = 0;
    for (int i = 0; i < max_seg_cnt; i += 2)
        if (((frag_l + i)->frags != NULL) and ((frag_l + i)->frag_count == 0))
            high_hits++;
}

void
FilterRead::pair_chains(const chain_list &forward_chain, const chain_list &reverse_chain, vector <MatePair> &mate_pairs,
                        bool *forward_paired, bool *reverse_paired, int saved_type) {
    vector<const IntervalInfo<UniqSeg> *> forward_exon_list(forward_chain.best_chain_count);
    vector<const IntervalInfo<UniqSeg> *> reverse_exon_list(reverse_chain.best_chain_count);

    uint32_t pos;

    for (int i = 0; i < forward_chain.best_chain_count; i++) {
        pos = forward_chain.chains[i].frags[0].rpos;
        forward_exon_list[i] = gtf_parser.get_location_overlap(pos, false);
    }

    for (int j = 0; j < reverse_chain.best_chain_count; j++) {
        pos = reverse_chain.chains[j].frags[0].rpos;
        reverse_exon_list[j] = gtf_parser.get_location_overlap(pos, false);
    }

    mate_pairs.clear();

    memset(forward_paired, 0, maxChainLen * sizeof(bool));
    memset(reverse_paired, 0, maxChainLen * sizeof(bool));

    uint32_t fs, fe;
    uint32_t rs, re;
    int tlen;

    for (int i = 0; i < forward_chain.best_chain_count; i++) {
        for (int j = 0; j < reverse_chain.best_chain_count; j++) {

            fs = forward_chain.chains[i].frags[0].rpos;
            rs = reverse_chain.chains[j].frags[0].rpos;

            fe = forward_chain.chains[i].frags[forward_chain.chains[i].chain_len - 1].rpos +
                 forward_chain.chains[i].frags[forward_chain.chains[i].chain_len - 1].len;
            re = reverse_chain.chains[j].frags[reverse_chain.chains[j].chain_len - 1].rpos +
                 reverse_chain.chains[j].frags[reverse_chain.chains[j].chain_len - 1].len;

            tlen = (fs < rs) ? (re - fs) : (fe - rs);

            MatePair temp;
            bool same_tr = false, same_gen = false;
            if (forward_exon_list[i] != NULL and reverse_exon_list[j] != NULL)
                same_tr = same_transcript(forward_exon_list[i], reverse_exon_list[j], temp.common_tid);
            if (!same_tr and forward_exon_list[i] != NULL and
                ((scanLevel == 0 and saved_type > CONGEN) or (scanLevel > 0 and saved_type >= CONGEN)))
                same_gen = same_gene(forward_exon_list[i], rs, re);
            if (!same_gen and reverse_exon_list[j] != NULL and (saved_type >= CONGEN))
                same_gen = same_gene(reverse_exon_list[j], fs, fe);

            if (same_tr or same_gen
                or ((tlen <= MAXDISCRDTLEN) and (saved_type >= CONGNM))) {
                //or (forward_exon_list[i] == NULL and reverse_exon_list[j] == NULL and tlen <= maxTlen)) {
                temp.forward = forward_chain.chains[i];
                temp.reverse = reverse_chain.chains[j];
                temp.score = forward_chain.chains[i].score + reverse_chain.chains[j].score;
                temp.type = (same_tr) ? 0 : ((same_gen) ? 1 : 2);
                mate_pairs.push_back(temp);

                forward_paired[i] = true;
                reverse_paired[j] = true;
            }
        }
    }

    //sort(mate_pairs.begin(), mate_pairs.end());

}

