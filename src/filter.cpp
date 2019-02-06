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

extern "C" {
#include "mrsfast/Common.h"
}

#define MINLB 0
#define MAXUB 4294967295	//2^32 - 1

//set <GenRegion> trans_extensions;

void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer* frag_l, int& high_hits);
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedMate& mr, int dir);
int process_mates(const chain_list& forward_chain, const Record* forward_rec, const chain_list& backward_chain, const Record* backward_rec, MatchedRead& mr, bool r1_forward);

void overlap_to_epos(MatchedMate& mr);
void overlap_to_spos(MatchedMate& mr);
void gene_overlap(MatchedMate& mr);


FilterRead::FilterRead (char* save_fname, bool pe, int round, bool first_round, bool last_round, char* fq_file1, char* fq_file2) {
	is_pe = pe;
	this->first_round = first_round;
	this->last_round = last_round;
	
	// temp fastq file(s) to be read in next round
	//if (! last_round) {
	char temp_fname [FILE_NAME_LENGTH];
	if (is_pe) {
		sprintf(temp_fname, "%s_%d_remain_R1.fastq", save_fname, round);
		temp_fq_r1 = open_file(temp_fname, "w");

		sprintf(temp_fname, "%s_%d_remain_R2.fastq", save_fname, round);
		temp_fq_r2 = open_file(temp_fname, "w");
	}
	else {
		sprintf(temp_fname, "%s_%d_remain.fastq", save_fname, round);
		temp_fq_r1 = open_file(temp_fname, "w");
	}
	//}

	// updating fq file to be read in next round
	if (is_pe) {
		sprintf(fq_file1, "%s_%d_remain_R1.fastq", save_fname, round);
		sprintf(fq_file2, "%s_%d_remain_R2.fastq", save_fname, round);
	}
	else {
		sprintf(fq_file1, "%s_%d_remain.fastq", save_fname, round);
	}

	// openning pam files
	char* output_names[11] = { "concordant", "discordant", "circ_RF", "circ_bsj", "fusion", "OEA2", "keep", "OEA", "orphan", "many_hits", "no_hit" };
	char cat_fname [FILE_NAME_LENGTH];

	char mode[2];
	sprintf(mode, "%s", (first_round) ? "w" : "a");

	sprintf(cat_fname, "%s.%s.pam", save_fname, output_names[0]);
	cat_file_pam[0] = open_file(cat_fname, mode);

	if (! last_round)
		return;

	for (int i = 1; i < CATNUM; i++) {
		sprintf(cat_fname, "%s.%s.pam", save_fname, output_names[i]);
		cat_file_pam[i] = open_file(cat_fname, "w");
	}

}

FilterRead::~FilterRead (void) {
	close_file(cat_file_pam[0]);

	//if (! last_round) {
	close_file(temp_fq_r1);
	close_file(temp_fq_r2);
	//}
	if (last_round) {
		for (int i = 1; i < CATNUM; i++) {
			close_file(cat_file_pam[i]);
		}
	}
}

// SE mode
int FilterRead::process_read (	Record* current_record, int kmer_size, GIMatchedKmer* fl, GIMatchedKmer* bl, 
								chain_list& forward_best_chain, chain_list& backward_best_chain) {
	
	vafprintf(1, stderr, "%s\n", current_record->rname);

	MatchedMate mr;
	int ex_ret, min_ret = ORPHAN;
	int fhh, bhh;

	get_best_chains(current_record->seq, current_record->seq_len, kmer_size, forward_best_chain, fl, fhh);
	for (int i = 0; i < forward_best_chain.best_chain_count; i++) {
		ex_ret = extend_chain(forward_best_chain.chains[i], current_record->seq, current_record->seq_len, mr, 1); 
		
		if (ex_ret == CONCRD)
			return CONCRD;

		if (ex_ret < min_ret)
			min_ret = ex_ret;
	}

	get_best_chains(current_record->rcseq, current_record->seq_len, kmer_size, backward_best_chain, bl, bhh);
	for (int i = 0; i < backward_best_chain.best_chain_count; i++) {
		ex_ret = extend_chain(backward_best_chain.chains[i], current_record->rcseq, current_record->seq_len, mr, -1); 
		
		if (ex_ret == CONCRD)
			return CONCRD;

		if (ex_ret < min_ret)
			min_ret = ex_ret;
	}

	return min_ret;

}

// PE mode
int FilterRead::process_read (	Record* current_record1, Record* current_record2, int kmer_size, GIMatchedKmer* fl, GIMatchedKmer* bl, 
								chain_list& forward_best_chain_r1, chain_list& backward_best_chain_r1, 
								chain_list& forward_best_chain_r2, chain_list& backward_best_chain_r2) {

	int max_frag_count = current_record1->seq_len / kmer_size + 1;
	int kmer_count = 2 * (current_record1->seq_len / kmer_size);
	int fhh_r1, bhh_r1;
	int fhh_r2, bhh_r2;

	// R1
	vafprintf(1, stderr, "R1/%s\n", current_record1->rname);

	get_best_chains(current_record1->seq, current_record1->seq_len, kmer_size, forward_best_chain_r1, fl, fhh_r1);
	get_best_chains(current_record1->rcseq, current_record1->seq_len, kmer_size, backward_best_chain_r1, bl, bhh_r1);

	vafprintf(1, stderr, "R1/%s\n", current_record1->rname);
	vafprintf(1, stderr, "R1 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r1.chains[0].score, (unsigned long)forward_best_chain_r1.best_chain_count);
	for (int j = 0; j < forward_best_chain_r1.best_chain_count; j++)
		for (int i = 0; i < forward_best_chain_r1.chains[j].chain_len; i++) {
			vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, forward_best_chain_r1.chains[j].frags[i].rpos, forward_best_chain_r1.chains[j].frags[i].qpos, forward_best_chain_r1.chains[j].frags[i].len);
		}

	vafprintf(1, stderr, "R1 Reverse score:%.4f,\t len: %lu\n", backward_best_chain_r1.chains[0].score, (unsigned long)backward_best_chain_r1.best_chain_count);
	for (int j = 0; j < backward_best_chain_r1.best_chain_count; j++)
		for (int i = 0; i < backward_best_chain_r1.chains[j].chain_len; i++) {
			vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, backward_best_chain_r1.chains[j].frags[i].rpos, backward_best_chain_r1.chains[j].frags[i].qpos, backward_best_chain_r1.chains[j].frags[i].len);
		}

	// R2
	vafprintf(1, stderr, "R2/%s\n", current_record2->rname);

	get_best_chains(current_record2->seq, current_record2->seq_len, kmer_size, forward_best_chain_r2, fl, fhh_r2);
	get_best_chains(current_record2->rcseq, current_record2->seq_len, kmer_size, backward_best_chain_r2, bl, bhh_r2);

	vafprintf(1, stderr, "R2/%s\n", current_record2->rname);
	vafprintf(1, stderr, "R2 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r2.chains[0].score, (unsigned long)forward_best_chain_r2.best_chain_count);
	for (int j = 0; j < forward_best_chain_r2.best_chain_count; j++)
		for (int i = 0; i < forward_best_chain_r2.chains[j].chain_len; i++) {
			vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, forward_best_chain_r2.chains[j].frags[i].rpos, forward_best_chain_r2.chains[j].frags[i].qpos, forward_best_chain_r2.chains[j].frags[i].len);
		}

	vafprintf(1, stderr, "R2 Reverse score:%.4f,\t len: %lu\n", backward_best_chain_r2.chains[0].score, (unsigned long)backward_best_chain_r2.best_chain_count);
	for (int j = 0; j < backward_best_chain_r2.best_chain_count; j++)
		for (int i = 0; i < backward_best_chain_r2.chains[j].chain_len; i++) {
			vafprintf(2, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, backward_best_chain_r2.chains[j].frags[i].rpos, backward_best_chain_r2.chains[j].frags[i].qpos, backward_best_chain_r2.chains[j].frags[i].len);
		}

	// Orphan / OEA
	if (forward_best_chain_r1.best_chain_count + backward_best_chain_r1.best_chain_count + forward_best_chain_r2.best_chain_count + backward_best_chain_r2.best_chain_count <= 0) {
		if ((fhh_r1 + bhh_r1 > 0) and (fhh_r2 + bhh_r2 > 0)) {
			current_record1->mr->update_type(NOPROC_MANYHIT);
			return NOPROC_MANYHIT;
		}
		else {
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
	float fc_score_r1 = (forward_best_chain_r1.best_chain_count  > 0) ? forward_best_chain_r1.chains[0].score  : 0;
	float bc_score_r1 = (backward_best_chain_r1.best_chain_count > 0) ? backward_best_chain_r1.chains[0].score : 0;
	float fc_score_r2 = (forward_best_chain_r2.best_chain_count  > 0) ? forward_best_chain_r2.chains[0].score  : 0;
	float bc_score_r2 = (backward_best_chain_r2.best_chain_count > 0) ? backward_best_chain_r2.chains[0].score : 0;
	vafprintf(2, stderr, "Scores: fc1=%f, bc1=%f, fc2=%f, bc2=%f\n", fc_score_r1, bc_score_r1, fc_score_r2, bc_score_r2);
	
	int attempt1, attempt2;
	if (fc_score_r1 + bc_score_r2 >= fc_score_r2 + bc_score_r1) {
		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt1 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2, *(current_record1->mr), true);
		if (scanLevel == 0 and attempt1 == CONCRD) {
			return CONCRD;
		}

		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt2 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1, *(current_record1->mr), false);
		if (scanLevel == 0 and attempt2 == CONCRD) {
			return CONCRD;
		}

		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
	else {
		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt1 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1, *(current_record1->mr), false);
		if (scanLevel == 0 and attempt1 == CONCRD) {
			return CONCRD;
		}

		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt2 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2, *(current_record1->mr), true);
		if (scanLevel == 0 and attempt2 == CONCRD) {
			return CONCRD;
		}

		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
}

// write reads SE mode
void FilterRead::write_read_category (Record* current_record, int state) {
	//state = minM(state, current_record->state);
	//int cat = (state >= cat_count) ? DISCRD : state;
	if (!last_round and state != CONCRD) {
		fprintf(temp_fq_r1, "%s\n%s%s%d\n%s", current_record->rname, current_record->seq, current_record->comment, state, current_record->qual);
	}
}

// write reads PE mode
void FilterRead::write_read_category (Record* current_record1, Record* current_record2, const MatchedRead& mr) {
	char r1_dir = (mr.r1_forward) ? '+' : '-';
	char r2_dir = (mr.r2_forward) ? '+' : '-';
	if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ) {
		sprintf(comment, " %d %s %u %u %d %u %u %c %d %s %u %u %d %u %u %c %d %d %d %d", 
						mr.type, 
						mr.chr_r1.c_str(), mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, r1_dir, mr.ed_r1,
						mr.chr_r2.c_str(), mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, r2_dir, mr.ed_r2,
						mr.tlen, mr.junc_num, mr.gm_compatible);
	}
	else {
		sprintf(comment, " %d * * * * * * * * * * * * * * * * * * *", mr.type);
	}

	fprintf(temp_fq_r1, "@%s%s\n%s%s%s", current_record1->rname, comment, current_record1->seq, current_record1->comment, current_record1->qual);
	fprintf(temp_fq_r2, "@%s%s\n%s%s%s", current_record2->rname, comment, current_record2->seq, current_record2->comment, current_record2->qual);

}

void FilterRead::print_mapping (char* rname, const MatchedRead& mr) {
	char r1_dir = (mr.r1_forward) ? '+' : '-';
	char r2_dir = (mr.r2_forward) ? '+' : '-';
	if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ) {
		fprintf(cat_file_pam[mr.type], "%s\t%s\t%u\t%u\t%d\t%u\t%u\t%c\t%d\t%s\t%u\t%u\t%d\t%u\t%u\t%c\t%d\t%d\t%d\t%d\t%d\n", 
										rname, 
										mr.chr_r1.c_str(), mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, r1_dir, mr.ed_r1, 
										mr.chr_r2.c_str(), mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, r2_dir, mr.ed_r2,
										mr.tlen, mr.junc_num, mr.gm_compatible, mr.type);
	}

	else {
		fprintf(cat_file_pam[mr.type], "%s\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\n", rname);
	}
}

bool is_concord(const chain_t& a, int seq_len, MatchedMate& mr) {
	if (a.chain_len < 2) {
		mr.is_concord = false;
	}
	else if ((a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos) >= seq_len) {
		mr.is_concord = true;
		mr.type = CONCRD;
		
		mr.spos = a.frags[0].rpos;
		mr.epos = a.frags[a.chain_len-1].rpos + a.frags[a.chain_len-1].len - 1;
		mr.matched_len = a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos;
		mr.qspos = a.frags[0].qpos;
		mr.qepos = a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - 1;
	}
	else {
		mr.is_concord = false;
	}
	return mr.is_concord;
}

void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer* frag_l, int& high_hits) {
	int kmer_count = ceil(seq_len / kmer_size);
	int forward_fragment_count, backward_fragment_count;
	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer_size)) - 1;	// considering both overlapping and non-overlapping kmers

	split_match_hash(read_seq, seq_len, kmer_size, frag_l);
	chain_seeds_sorted_kbest(seq_len, frag_l, best_chain);

	high_hits = 0;
	for (int i = 0; i < max_seg_cnt; i+=2)
		if (((frag_l+i)->frags != NULL) and ((frag_l+i)->frag_count == 0))
			high_hits++;
}

bool extend_chain_left(const vector <uint32_t>& common_tid, const chain_t& ch, char* seq, int seq_len, int lb, MatchedMate& mr, int& err) {
	bool left_ok = true;

	uint32_t lm_pos = ch.frags[0].rpos;
	int remain_beg = ch.frags[0].qpos;

	left_ok = (remain_beg <= 0);
	AlignRes best_alignment(lb);
	
	char remain_str_beg[remain_beg+5];
	if (remain_beg > 0) {
		left_ok = extend_left(common_tid, seq, lm_pos, remain_beg, maxEd - err, lb, best_alignment);
	}
	
	int sclen_left = best_alignment.sclen;
	int err_left = best_alignment.ed;
	remain_beg -= best_alignment.rcovlen;

	mr.spos = lm_pos;
	mr.matched_len -= (left_ok) ? sclen_left : remain_beg;
	mr.qspos += (left_ok) ? sclen_left : remain_beg;
	mr.sclen_left = sclen_left;
	mr.left_ed = best_alignment.ed;

	err += err_left;

	return left_ok;
}

bool extend_chain_right(const vector <uint32_t>& common_tid, const chain_t& ch, char* seq, int seq_len, int ub, MatchedMate& mr, int& err) {
	bool right_ok = true;

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len);

	right_ok = (remain_end <= 0);
	AlignRes best_alignment(ub);

	char remain_str_end[remain_end+5];
	if (remain_end > 0) {
		right_ok = extend_right(common_tid, seq + seq_len - remain_end, rm_pos, remain_end, maxEd - err, ub, best_alignment);
	}

	int sclen_right = best_alignment.sclen;
	int err_right = best_alignment.ed;
	remain_end -= best_alignment.rcovlen;

	mr.epos = rm_pos;
	mr.matched_len -= (right_ok)? sclen_right : remain_end;
	mr.qepos -= (right_ok)? sclen_right : remain_end;
	mr.sclen_right = sclen_right;
	mr.right_ed = best_alignment.ed;

	err += err_right;

	return right_ok;
}

void update_match_mate_info(bool lok, bool rok, int err, MatchedMate& mm) {
	mm.left_ok = lok;
	mm.right_ok = rok;
	if (lok and rok and (err <= maxEd)) {
		mm.is_concord = true;
		mm.type = CONCRD;
	}
	else if (lok or rok) {
		mm.type = CANDID;
	}
	else
		mm.type = ORPHAN;
}

int estimate_middle_error(const chain_t& ch) {
	int mid_err = 0;
	for (int i = 0; i < ch.chain_len - 1; i++) {
		if (ch.frags[i+1].qpos > ch.frags[i].qpos + ch.frags[i].len) {
			int diff = (ch.frags[i+1].rpos - ch.frags[i].rpos) - (ch.frags[i+1].qpos - ch.frags[i].qpos);
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

int calc_middle_ed(const chain_t& ch, int edth, char* qseq, int qseq_len) {
	char rseq[qseq_len + 4*bandWidth];
	int mid_err = 0;
	int32_t qspos;
	uint32_t rspos;
	int qlen;
	int rlen;
	for (int i = 0; i < ch.chain_len - 1; i++) {
		if (ch.frags[i+1].qpos > ch.frags[i].qpos + ch.frags[i].len) {
			int diff = (ch.frags[i+1].rpos - ch.frags[i].rpos) - (ch.frags[i+1].qpos - ch.frags[i].qpos);

			qspos = ch.frags[i].qpos + ch.frags[i].len;
			qlen = ch.frags[i+1].qpos - qspos;
			rspos = ch.frags[i].rpos + ch.frags[i].len;
			rlen = qlen + diff;

			// fprintf(stderr, "Diff: %d\n", diff);
			if (diff >= 0 and diff <= bandWidth) {
				pac2char(rspos, rlen, rseq);
				// fprintf(stderr, "qlen: %d\nrlen: %d\nQ str: %s\nT str: %s\n", qlen, rlen, qseq+qspos, rseq);
				mid_err += alignment.global_one_side_banded_alignment(qseq + qspos, qlen, rseq, rlen, diff);
			}
			else if (diff < 0 and diff >= (-1 * bandWidth)) {
				pac2char(rspos, rlen, rseq);
				// fprintf(stderr, "qlen: %d\nrlen: %d\nQ str: %s\nT str: %s\n", qlen, rlen, qseq+qspos, rseq);
				mid_err += alignment.global_one_side_banded_alignment(rseq, rlen, qseq + qspos, qlen, -1*diff);
			}
			if (mid_err > edth)
				return edth+1;
		}
	}
	return mid_err;
}

bool check_middle_ed(const chain_t& ch, int edth, char* qseq, int qseq_len) {
	int mid_err = calc_middle_ed(ch, edth, qseq, qseq_len);
	return (mid_err <= edth);
}

bool extend_both_mates(const chain_t& lch, const chain_t& rch, const vector<uint32_t>& common_tid, char* lseq, char* rseq, 
						int lseq_len, int rseq_len, MatchedMate& lmm, MatchedMate& rmm) {
	
	// lmm.middle_ed = estimate_middle_error(lch);
	// rmm.middle_ed = estimate_middle_error(rch);

	lmm.middle_ed = calc_middle_ed(lch, maxEd, lseq, lseq_len);
	// fprintf(stderr, "Left middle ed: %d\n", lmm.middle_ed);
	if (lmm.middle_ed > maxEd)
		return false;

	rmm.middle_ed = calc_middle_ed(rch, maxEd, rseq, lseq_len);
	// fprintf(stderr, "Right middle ed: %d\n", rmm.middle_ed);
	if (rmm.middle_ed > maxEd)
		return false;

	bool l_extend = true;
	lmm.is_concord = false;
	if (lch.chain_len <= 0) {
		lmm.type = ORPHAN;
		lmm.matched_len = 0;
		l_extend = false;
	}

	bool r_extend = true;
	rmm.is_concord = false;
	if (rch.chain_len <= 0) {
		rmm.type = ORPHAN;
		rmm.matched_len = 0;
		r_extend = false;
	}

	bool llok;
	bool lrok;
	bool rlok;
	bool rrok;

	int lerr = lmm.middle_ed;
	int rerr = rmm.middle_ed;

	if (l_extend) {
		lmm.matched_len = lseq_len;
		lmm.qspos = 1;
		lmm.qepos = lseq_len;
		llok = extend_chain_left(common_tid, lch, lseq, lseq_len, MINLB, lmm, lerr);
	}

	if (r_extend) {
		rmm.matched_len = rseq_len;
		rmm.qspos = 1;
		rmm.qepos = rseq_len;
		rlok = extend_chain_left(common_tid, rch, rseq, rseq_len, (l_extend) ? lmm.spos : MINLB, rmm, rerr);
	}
	
	if (r_extend) {
		rrok = extend_chain_right(common_tid, rch, rseq, rseq_len, MAXUB, rmm, rerr);
	}
	
	if (l_extend) {
		lrok = extend_chain_right(common_tid, lch, lseq, lseq_len, (r_extend) ? rmm.epos : MAXUB, lmm, lerr);
	}
	
	if (l_extend) {
		update_match_mate_info(llok, lrok, lerr, lmm);
	}

	if (r_extend) {
		update_match_mate_info(rlok, rrok, rerr, rmm);
	}

	// check accurate edit distance for middle part
	return true;

	// int ledth = maxEd - (lmm.left_ed + lmm.right_ed);
	// lmm.middle_ed = calc_middle_ed(lch, ledth, lseq, lseq_len);
	// if (lmm.middle_ed + lmm.right_ed + lmm.left_ed > maxEd)
	// 	return false;
	
	// int redth = maxEd - (rmm.left_ed + rmm.right_ed);
	// rmm.middle_ed = calc_middle_ed(rch, redth, rseq, lseq_len);
	// return (rmm.middle_ed + rmm.right_ed + rmm.left_ed <= maxEd);
}

// return:
// CONCRD:0 if mapped to both sides
// CANDID:1 if at least one side is mapped
// ORPHAN:3 if not mappable (not reaching either side)
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedMate& mr, int dir) {
	mr.is_concord = false;
	if (ch.chain_len <= 0) {
		mr.type = ORPHAN;
		return mr.type;
	}

	mr.middle_ed = estimate_middle_error(ch);

	if (is_concord(ch, seq_len, mr)) {
		mr.dir = dir;
		return mr.type;
	}

	bool left_ok = true;
	bool right_ok = true;
	int sclen_left = 0;
	int sclen_right = 0;
	
	int err_left = 0;
	int err_right = 0;

	uint32_t lm_pos = ch.frags[0].rpos;
	int remain_beg = ch.frags[0].qpos;

	left_ok = (remain_beg <= 0);
	AlignRes best_alignment_left(MINLB);
	vector <uint32_t> empty;

	char remain_str_beg[remain_beg+5];
	if (remain_beg > 0) {
		left_ok = extend_left(empty, seq, lm_pos, remain_beg, maxEd - mr.middle_ed, MINLB, best_alignment_left);
	}

	err_left = best_alignment_left.ed;
	sclen_left = best_alignment_left.sclen;
	remain_beg -= best_alignment_left.rcovlen;

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len);

	right_ok = (remain_end <= 0);
	AlignRes best_alignment(MAXUB);

	char remain_str_end[remain_end+5];
	if (remain_end > 0) {
		right_ok = extend_right(empty, seq + seq_len - remain_end, rm_pos, remain_end, maxEd - mr.middle_ed - err_left, MAXUB, best_alignment);
	}

	err_right = best_alignment.ed;
	sclen_right = best_alignment.sclen;
	remain_end -= best_alignment.rcovlen;

	mr.spos = lm_pos;
	mr.epos = rm_pos;
	mr.matched_len = seq_len;
	mr.matched_len -= (left_ok) ? sclen_left : remain_beg;
	mr.matched_len -= (right_ok) ? sclen_right: remain_end;

	mr.qspos = 1 + (left_ok) ? sclen_left : remain_beg;
	mr.qepos = seq_len - (right_ok) ? sclen_right: remain_end;

	mr.right_ed = best_alignment.ed;
	mr.left_ed  = best_alignment_left.ed;

	mr.dir = dir;
	
	//fprintf(stderr, "############################# left_ok: %d\tright_ok: %d\terr_left: %d\terr_right: %d\n", left_ok, right_ok, err_left, err_right);
	if (left_ok and right_ok and (err_left + err_right <= maxEd)) {
		mr.is_concord = true;
		mr.type = CONCRD;
	}
	else if (left_ok or right_ok) {
		mr.type = CANDID;
	}
	else
		mr.type = ORPHAN;

	return mr.type;
}

bool same_gene(const MatchedMate& mm, const MatchedMate& other) {
	GeneInfo* ginfo;
	for (int i = 0; i < mm.exons_spos->seg_list.size(); i++) {
		ginfo = gtf_parser.get_gene_info(mm.exons_spos->seg_list[i].gene_id);
		//fprintf(stderr, "Gene[%d][%s]: [%d - %d], [%d - %d]\n", i, mm.exons_spos->seg_list[i].gene_id.c_str(), ginfo->start, ginfo->end, other.spos, other.epos);
		if (ginfo->start <= other.spos and other.epos <= ginfo->end)
			return true;
	}

	return false;
}

bool same_gene(uint32_t sme, const IntervalInfo<GeneInfo>* smg, uint32_t lms, const IntervalInfo<GeneInfo>* lmg) {
	if (smg == NULL or lmg == NULL)
		return false;

	if (smg->seg_list.size() == 0 or lmg->seg_list.size() == 0)
		return false;

	bool same_intron;
	int step = 10;
	for (int i = 0; i < smg->seg_list.size(); i++)
		for (int j = 0; j < lmg->seg_list.size(); j++)
			if (smg->seg_list[i].start == lmg->seg_list[j].start and smg->seg_list[i].end == lmg->seg_list[j].end) {
				same_intron = true;
				for (int k = sme; k <= lms; k += step)
					if (!(intronic_bs[contigNum][k])) {
						same_intron = false;
						break;
					}
				if (same_intron)
					return true;
			}
	
	return false;
}

// calculate tlen including sm.epos and lm.spos
int calc_tlen(const MatchedMate& sm, const MatchedMate& lm, int& intron_num) {
	const IntervalInfo<UniqSeg>* this_region;
	uint32_t tid;
	int start_ind;
	int start_table_ind;
	int end_table_ind;
	int tlen;
	int min_tlen = INF;
	int this_it_ind;
	int in;

	for (int i = 0; i < sm.exons_epos->seg_list.size(); i++) {
		for (int j = 0; j < sm.exons_epos->seg_list[i].trans_id.size(); j++) {
			tid = sm.exons_epos->seg_list[i].trans_id[j];
			start_ind = gtf_parser.get_trans_start_ind(contigNum, tid);
			start_table_ind = sm.exon_ind_epos - start_ind;
			if (start_table_ind < 0)	// assert
				continue;
			
			end_table_ind = lm.exon_ind_spos - start_ind;
			if (end_table_ind >= gtf_parser.trans2seg[contigNum][tid].size() or gtf_parser.trans2seg[contigNum][tid][end_table_ind] == 0)	// transcript does not contain lm exon
				continue;

			if (start_table_ind == end_table_ind) {
				in = 0;
				tlen = lm.spos - sm.epos + 1;
			}
			else {
				bool pre_zero = false;
				in = 0;
				tlen = sm.exons_epos->epos - sm.epos + 1;
				this_it_ind = sm.exon_ind_epos;
				for (int k = start_table_ind + 1; k < end_table_ind; k++) {
					this_it_ind++;
					if (gtf_parser.trans2seg[contigNum][tid][k] != 0) {
						this_region = gtf_parser.get_interval(this_it_ind);
						tlen += this_region->epos - this_region->spos + 1;
						pre_zero = false;
					}
					else {
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

// sm should start before lm
bool concordant_explanation(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm) {
	if (sm.spos > lm.spos)
		return false;

	int32_t tlen;
	int32_t intron_gap;
	bool on_cdna;
	on_cdna = (sm.exons_spos != NULL) and (sm.exons_epos != NULL) and (lm.exons_spos != NULL) and (lm.exons_epos != NULL);
	if (sm.exons_spos == NULL or lm.exons_spos == NULL) {
		tlen = lm.spos - sm.epos - 1 + lm.matched_len + sm.matched_len;
		if (tlen <= maxTlen)
			mr.update(sm, lm, chr, shift, tlen, 0, false, CONCRD, r1_sm);
		else if (tlen <= MAXDISCRDTLEN)
			mr.update(sm, lm, chr, shift, tlen, 0, false, DISCRD, r1_sm);
	}
	else {
		//fprintf(stderr, "Left Mate [%u-%u] dir=%d, type=%d, Right Mate[%u-%u] dir=%d, type=%d\n", sm.spos, sm.epos, sm.dir, sm.type, lm.spos, lm.epos, lm.dir, lm.type);
	// starts on same exon
		for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
			for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
				if (sm.exons_spos->seg_list[i].same_exon(lm.exons_spos->seg_list[j])) {
					// => assume genomic locations
					tlen = lm.spos + lm.matched_len - sm.spos;
					if (tlen <= maxTlen)
						mr.update(sm, lm, chr, shift, tlen, 0, on_cdna, CONCRD, r1_sm);
					else
						mr.update(sm, lm, chr, shift, tlen, 0, on_cdna, DISCRD, r1_sm);
				}
	}

	if (sm.exons_epos == NULL or lm.exons_spos == NULL) {
		tlen = lm.spos - sm.epos - 1 + sm.matched_len + lm.matched_len;
		if (tlen <= maxTlen)
			mr.update(sm, lm, chr, shift, tlen, 0, false, CONCRD, r1_sm);
		else if (tlen <= MAXDISCRDTLEN)
			mr.update(sm, lm, chr, shift, tlen, 0, false, DISCRD, r1_sm);
	}
	else {
		int intron_num;
		tlen = calc_tlen(sm, lm, intron_num);
		//fprintf(stdout, "tlen: %d\n", tlen);
		if (tlen >= 0 and tlen <= maxTlen) {
			mr.update(sm, lm, chr, shift, tlen, intron_num, on_cdna, CONCRD, r1_sm);
		}
		else {
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

bool check_chimeric(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm) {
	if (mr.type == CONCRD)
		return false;

	if (sm.exons_spos == NULL or lm.exons_spos == NULL)
		return false;

	//fprintf(stderr, "Left Mate [%u-%u] dir=%d, type=%d, Right Mate[%u-%u] dir=%d, type=%d\n", sm.spos, sm.epos, sm.dir, sm.type, lm.spos, lm.epos, lm.dir, lm.type);
	for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
		for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
			if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j]) and sm.spos < lm.spos) {
				mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIORF, r1_sm);
				return true;
			}
	return false;
}

bool check_bsj(MatchedMate& sm, MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm) {
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
		//overlap_to_spos(sm);
		//overlap_to_epos(lm);
		//fprintf(stderr, "R1 start ind: %d\tR2 end ind: %d\n To beg of intron: %d\n", sm.exon_ind_spos, lm.exon_ind_epos, sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos));
		if ((intronic_bs[contigNum][sm.spos]) and ((intronic_bs[contigNum][lm.spos])) and 
			(sm.exon_ind_spos >= 0) and (lm.exon_ind_epos >= 0) and (sm.exon_ind_spos == lm.exon_ind_epos) and 
			(sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos) <= LARIAT2BEGTH)) {
			mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ, r1_sm);
			return true;
		}
		
		//gene_overlap(sm);
		//gene_overlap(lm);
		//if (same_gene(sm.epos, sm.gene_info, lm.spos, lm.gene_info)) {
		//	mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ);
		//	return true;
		//}

		return false;
	}

	for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
		for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
			if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j])) {
				mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ, r1_sm);
				return true;
			}
	return false;
}

bool same_transcript(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r, MatePair& mp) {
	mp.common_tid.clear();
	if (s == NULL or r == NULL)
		return false;

	// fprintf(stderr, "In same_transcript\nseg size1 = %d\tseg size2 = %d\n", s->seg_list.size(), r->seg_list.size());
	vector<uint32_t> seg1_tid;
	vector<uint32_t> seg2_tid;

	for (int i = 0; i < s->seg_list.size(); i++)
		for (int k = 0; k < s->seg_list[i].trans_id.size(); k++)
			seg1_tid.push_back(s->seg_list[i].trans_id[k]);

	for (int j = 0; j < r->seg_list.size(); j++)
		for (int l = 0; l < r->seg_list[j].trans_id.size(); l++)
			seg2_tid.push_back(r->seg_list[j].trans_id[l]);

	for (int i = 0; i < seg1_tid.size(); i++)
		for (int j = 0; j < seg2_tid.size(); j++) {
			uint32_t tid1 = seg1_tid[i];
			uint32_t tid2 = seg2_tid[j];
			// fprintf(stderr, "tr[%d]: %s\ttr[%d]: %s", tid1, gtf_parser.transcript_ids[contigNum][tid1].c_str(), tid2, gtf_parser.transcript_ids[contigNum][tid2].c_str());
			if (tid1 == tid2) {
				//fprintf(stderr, "\ttid: %d\n", tid1);
				mp.common_tid.push_back(tid1);
				break;
			}
			// fprintf(stderr, "\n" );
		}

	// for (int i = 0; i < s->seg_list.size(); i++)
	// 	for (int j = 0; j < r->seg_list.size(); j++) {
	// 		fprintf(stderr, "trans size[%d] = %d\ttrans size[%d] = %d\n", i, s->seg_list[i].trans_id.size(), j, r->seg_list[j].trans_id.size());
	// 		for (int k = 0; k < s->seg_list[i].trans_id.size(); k++) {
	// 			for (int l = 0; l < r->seg_list[j].trans_id.size(); l++) {
	// 				int tid1 = s->seg_list[i].trans_id[k];
	// 				int tid2 = r->seg_list[j].trans_id[l];
	// 				fprintf(stderr, "tr[%d][%d]: %s\ttr[%d][%d]: %s", i, tid1, gtf_parser.transcript_ids[contigNum][tid1].c_str(), j, tid2, gtf_parser.transcript_ids[contigNum][tid2].c_str());
	// 				if (tid1 == tid2) {
	// 					fprintf(stderr, "\ttid: %d\n", tid1);
	// 					mp.common_tid.push_back(tid1);
	// 					break;
	// 				}
	// 				fprintf(stderr, "\n" );
	// 			}
	// 		}
	// 	}

	return mp.common_tid.size() != 0;
}

bool same_gene(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r) {
	if (s == NULL or r == NULL)
		return false;

	for (int i = 0; i < s->seg_list.size(); i++)
		for (int j = 0; j < r->seg_list.size(); j++)
			if (s->seg_list[i].gene_id == r->seg_list[j].gene_id)
				return true;

	return false;
}

bool same_gene(const IntervalInfo<UniqSeg>* mate, uint32_t s, uint32_t e) {
	GeneInfo* ginfo;
	for (int i = 0; i < mate->seg_list.size(); i++) {
		ginfo = gtf_parser.get_gene_info(mate->seg_list[i].gene_id);
		if (ginfo->start <= s and e <= ginfo->end)
			return true;
	}
		

	return false;
}

void pair_chains(const chain_list& forward_chain, const chain_list& reverse_chain, vector <MatePair>& mate_pairs, bool* forward_paired, bool* reverse_paired) {
	vector <const IntervalInfo<UniqSeg>*> forward_exon_list(forward_chain.best_chain_count);
	vector <const IntervalInfo<UniqSeg>*> reverse_exon_list(reverse_chain.best_chain_count);

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

	memset(forward_paired, 0, BESTCHAINLIM * sizeof(bool));
	memset(reverse_paired, 0, BESTCHAINLIM * sizeof(bool));

	uint32_t fs, fe;
	uint32_t rs, re;
	int tlen;

	for (int i = 0; i < forward_chain.best_chain_count; i++) {
		for (int j = 0; j < reverse_chain.best_chain_count; j++) {
			
			fs = forward_chain.chains[i].frags[0].rpos;
			rs = reverse_chain.chains[j].frags[0].rpos;

			fe = forward_chain.chains[i].frags[forward_chain.chains[i].chain_len-1].rpos + forward_chain.chains[i].frags[forward_chain.chains[i].chain_len-1].len;
			re = reverse_chain.chains[j].frags[reverse_chain.chains[j].chain_len-1].rpos + reverse_chain.chains[j].frags[reverse_chain.chains[j].chain_len-1].len;

			tlen = (fs < rs) ? (re - fs) : (fe - rs);

			MatePair temp;
			if ((forward_exon_list[i] != NULL and reverse_exon_list[j] != NULL and same_transcript(forward_exon_list[i], reverse_exon_list[j], temp))
				or (forward_exon_list[i] != NULL and same_gene(forward_exon_list[i], rs, re))
				or (reverse_exon_list[j] != NULL and same_gene(reverse_exon_list[j], fs, fe))
				or (tlen <= MAXDISCRDTLEN)) {
				//or (forward_exon_list[i] == NULL and reverse_exon_list[j] == NULL and tlen <= maxTlen)) {
				temp.forward = forward_chain.chains[i];
				temp.reverse = reverse_chain.chains[j];
				temp.score = forward_chain.chains[i].score + reverse_chain.chains[j].score;
				mate_pairs.push_back(temp);

				forward_paired[i] = true;
				reverse_paired[j] = true;
			}
		}
	}

	//sort(mate_pairs.begin(), mate_pairs.end());
	
}

// r1 is the forward mate and r2 is backward
int process_mates(const chain_list& forward_chain, const Record* forward_rec, const chain_list& backward_chain, const Record* backward_rec, MatchedRead& mr, bool r1_forward) {
	vector <MatePair> mate_pairs;

	bool forward_paired[BESTCHAINLIM];
	bool backward_paired[BESTCHAINLIM];
	pair_chains(forward_chain, backward_chain, mate_pairs, forward_paired, backward_paired);

	int min_ret1 = ORPHAN;
	int min_ret2 = ORPHAN;

	bool r1_genic = false;
	bool r2_genic = false;

	// concordant?
	vafprintf(1, stderr, "#pairs = %d\n", mate_pairs.size());
	for (int i = 0; i < mate_pairs.size(); i++) {
		vafprintf(2, stderr, "Mate[%d]: %d, %d\n", i, mate_pairs[i].forward.frags[0].rpos, mate_pairs[i].reverse.frags[0].rpos);
		// fprintf(stderr, "# common_tid: %d\n", mate_pairs[i].common_tid.size());
		// for (int k = 0; k < mate_pairs[i].common_tid.size(); k++)
		// 	fprintf(stderr, "%d, ", mate_pairs[i].common_tid[k]);
		// fprintf(stderr, "---\n");
		MatchedMate r1_mm;
		MatchedMate r2_mm;

		r1_mm.dir = 1;
		r2_mm.dir = -1;

		bool success;
		uint32_t forward_start = mate_pairs[i].forward.frags[0].rpos;
		uint32_t reverse_start = mate_pairs[i].reverse.frags[0].rpos;
		uint32_t reverse_end   = mate_pairs[i].reverse.frags[mate_pairs[i].reverse.chain_len-1].rpos + mate_pairs[i].reverse.frags[mate_pairs[i].reverse.chain_len-1].len - 1;
		
		if (forward_start <= reverse_end) {
			//extend_both_mates(mate_pairs[i].forward, mate_pairs[i].reverse, forward_rec->seq, backward_rec->rcseq, forward_rec->seq_len, backward_rec->seq_len, r1_mm, r2_mm);
			success = extend_both_mates(mate_pairs[i].forward, mate_pairs[i].reverse, mate_pairs[i].common_tid, forward_rec->seq, 
								backward_rec->rcseq, forward_rec->seq_len, backward_rec->seq_len, r1_mm, r2_mm);
			
			if (success and r1_mm.type == CONCRD and r2_mm.type == CONCRD) {
				ConShift con_shift = gtf_parser.get_shift(contigNum, r1_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				if (concordant_explanation(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward) and scanLevel == 0) {
					return CONCRD;
				}
			}
			
			// potentially back splice junction?
			else if (success and ((r1_mm.type == CANDID and r2_mm.type == CONCRD) or (r1_mm.type == CONCRD and r2_mm.type == CANDID))) {
				ConShift con_shift = gtf_parser.get_shift(contigNum, r1_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				check_bsj(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward);
			}
		}

		if (forward_start > reverse_start) {
			//extend_both_mates(mate_pairs[i].reverse, mate_pairs[i].forward, backward_rec->rcseq, forward_rec->seq, backward_rec->seq_len, forward_rec->seq_len, r2_mm, r1_mm);
			success = extend_both_mates(mate_pairs[i].reverse, mate_pairs[i].forward, mate_pairs[i].common_tid, backward_rec->rcseq, forward_rec->seq, backward_rec->seq_len, forward_rec->seq_len, r2_mm, r1_mm);
			
			if (success and r1_mm.type == CONCRD and r2_mm.type == CONCRD) {
				ConShift con_shift = gtf_parser.get_shift(contigNum, r2_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				check_chimeric(r2_mm, r1_mm, mr, con_shift.contig, con_shift.shift, !r1_forward);
			}
			
			// potentially back splice junction?
			else if (success and ((r1_mm.type == CANDID and r2_mm.type == CONCRD) or (r1_mm.type == CONCRD and r2_mm.type == CANDID))) {
				ConShift con_shift = gtf_parser.get_shift(contigNum, r2_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				check_bsj(r2_mm, r1_mm, mr, con_shift.contig, con_shift.shift, !r1_forward);
			}
		}

		min_ret1 = minM(r1_mm.type, min_ret1);
		min_ret2 = minM(r2_mm.type, min_ret2);

		r1_genic = (r1_mm.exons_spos != NULL) or (r1_mm.exons_epos != NULL);
		r2_genic = (r2_mm.exons_spos != NULL) or (r2_mm.exons_epos != NULL);

	}

	if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ) {
		return mr.type;
	}

	MatchedMate mm1;
	int ex_ret;
	if (min_ret1 != CONCRD)
		for (int i = 0; i < forward_chain.best_chain_count; i++) {
			if (!forward_paired[i]) {
				ex_ret = extend_chain(forward_chain.chains[i], forward_rec->seq, forward_rec->seq_len, mm1, 1);
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
				ex_ret = extend_chain(backward_chain.chains[i], backward_rec->rcseq, backward_rec->seq_len, mm2, -1);
				min_ret2 = minM(ex_ret, min_ret2);
				
				overlap_to_spos(mm2);
				overlap_to_epos(mm2);

				r2_genic = (mm2.exons_spos != NULL) or (mm2.exons_epos != NULL);
			}
		}

	int new_type = (((min_ret1 == ORPHAN) and (min_ret2 == CONCRD)) or ((min_ret1 == CONCRD) and (min_ret2 == ORPHAN))) ? OEANCH 
				: ((min_ret1 == ORPHAN) or (min_ret2 == ORPHAN)) ? ORPHAN 
				: ((min_ret1 == CONCRD) and (min_ret2 == CONCRD) and (r1_genic and r2_genic)) ? CHIFUS
				: ((min_ret1 == CONCRD) and (min_ret2 == CONCRD)) ? OEA2
				: CANDID;

	mr.update_type(new_type);
	return mr.type;
}

void overlap_to_epos(MatchedMate& mr) {
	if (mr.looked_up_epos or mr.exons_epos != NULL)
		return;

	mr.exons_epos = gtf_parser.get_location_overlap_ind(mr.epos, false, mr.exon_ind_epos);
	mr.looked_up_epos = true;
	//fprintf(stdout, "End Seg list size: %d\n", mr.exons_epos->seg_list.size());
}

void overlap_to_spos(MatchedMate& mr) {
	if (mr.looked_up_spos or mr.exons_spos != NULL)
		return;

	mr.exons_spos = gtf_parser.get_location_overlap_ind(mr.spos, false, mr.exon_ind_spos);
	mr.looked_up_spos = true;
	//fprintf(stdout, "Start Seg list size: %d\n", mr.exons_spos->seg_list.size());
}

void gene_overlap(MatchedMate& mr) {
	if (mr.looked_up_gene or mr.gene_info != NULL)
		return;
	mr.gene_info = gtf_parser.get_gene_overlap(mr.spos, false);
	mr.looked_up_gene = true;
}
