#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <set>

#include "filter.h"
#include "align.h"
#include "common.h"
#include "gene_annotation.h"

extern "C" {
#include "mrsfast/Common.h"
}

#define MINLB 0
#define MAXUB 4294967295	//2^32 - 1

//set <GenRegion> trans_extensions;

void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer*& frag_l, int& high_hits);
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedMate& mr, int dir);
int process_mates(const chain_list& forward_chain, const Record* forward_rec, const chain_list& backward_chain, const Record* backward_rec, MatchedRead& mr, bool r1_forward);

bool extend_right(char* seq, uint32_t& pos, int len, uint32_t ub, int& err, int& sclen_right);
bool extend_left(char* seq, uint32_t& pos, int len, uint32_t lb, int& err, int& sclen_left);

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
int FilterRead::process_read (	Record* current_record, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
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
int FilterRead::process_read (	Record* current_record1, Record* current_record2, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
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
			current_record1->mr.update_type(NOPROC_MANYHIT);
			return NOPROC_MANYHIT;
		}
		else {
			current_record1->mr.update_type(NOPROC_NOMATCH);
			return NOPROC_NOMATCH;
		}
	}
	if ((forward_best_chain_r1.best_chain_count + backward_best_chain_r1.best_chain_count <= 0) or 
		(forward_best_chain_r2.best_chain_count + backward_best_chain_r2.best_chain_count <= 0)) {
		current_record1->mr.update_type(OEANCH);
		return OEANCH;
	}

	// checking read concordancy
	// any concordancy evidence/explanation ?
	float fc_score_r1 = forward_best_chain_r1.chains[0].score;
	float bc_score_r1 = backward_best_chain_r1.chains[0].score;
	float fc_score_r2 = forward_best_chain_r2.chains[0].score;
	float bc_score_r2 = backward_best_chain_r2.chains[0].score;
	vafprintf(2, stderr, "Scores: fc1=%f, bc1=%f, fc2=%f, bc2=%f\n", fc_score_r1, bc_score_r1, fc_score_r2, bc_score_r2);
	
	int attempt1, attempt2;
	if (fc_score_r1 + bc_score_r2 >= fc_score_r2 + bc_score_r1) {
		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt1 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2, current_record1->mr, true);
		if (attempt1 == CONCRD) {
			return CONCRD;
		}

		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt2 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1, current_record1->mr, false);
		if (attempt2 == CONCRD) {
			return CONCRD;
		}

		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
	else {
		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt1 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1, current_record1->mr, false);
		if (attempt1 == CONCRD) {
			return CONCRD;
		}

		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt2 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2, current_record1->mr, true);
		if (attempt2 == CONCRD) {
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
	sprintf(comment, " %d %s %u %u %d %u %u %c %u %u %d %u %u %c %d %d %d", 
						mr.type, mr.chr.c_str(), 
						mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, r1_dir,
						mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, r2_dir,
						mr.tlen, mr.junc_num, mr.gm_compatible);

	fprintf(temp_fq_r1, "%s%s\n%s%s%s", current_record1->rname, comment, current_record1->seq, current_record1->comment, current_record1->qual);
	fprintf(temp_fq_r2, "%s%s\n%s%s%s", current_record2->rname, comment, current_record2->seq, current_record2->comment, current_record2->qual);

}

void FilterRead::print_mapping (char* rname, const MatchedRead& mr) {
	char r1_dir = (mr.r1_forward) ? '+' : '-';
	char r2_dir = (mr.r2_forward) ? '+' : '-';
	if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ) {
		fprintf(cat_file_pam[mr.type], "%s\t%s\t%u\t%u\t%d\t%u\t%u\t%c\t%u\t%u\t%d\t%u\t%u\t%c\t%d\t%d\t%d\t%d\n", 
										rname, mr.chr.c_str(), 
										mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, r1_dir, 
										mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, r2_dir,
										mr.tlen, mr.junc_num, mr.gm_compatible, mr.type);
	}

	else {
		fprintf(cat_file_pam[mr.type], "%s\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\n", rname);
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

void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer*& frag_l, int& high_hits) {
	int kmer_count = ceil(seq_len / kmer_size);
	int forward_fragment_count, backward_fragment_count;
	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer_size)) - 1;	// considering both overlapping and non-overlapping kmers

	//for (int i = 0; i < max_seg_cnt; i++) {
	//	memset(frag_l[i].junc_dist, 0, FRAGLIM * sizeof(JunctionDist));
	//}

	split_match_hash(read_seq, seq_len, kmer_size, frag_l);
	chain_seeds_sorted_kbest(seq_len, frag_l, best_chain);

	high_hits = 0;
	for (int i = 0; i < max_seg_cnt; i+=2)
		if (((frag_l+i)->frags != NULL) and ((frag_l+i)->frag_count == 0))
			high_hits++;
}

bool extend_chain_left(const chain_t& ch, char* seq, int seq_len, int lb, MatchedMate& mr, int& err) {
	bool left_ok = true;
	int sclen_left = 0;
	int err_left = 0;

	uint32_t lm_pos = ch.frags[0].rpos;
	int remain_beg = ch.frags[0].qpos;

	left_ok = (remain_beg <= 0);
	
	char remain_str_beg[remain_beg+5];
	if (remain_beg > 0) {
		left_ok = extend_left(seq, lm_pos, remain_beg, lb, err_left, sclen_left);
	}
	
	mr.spos = lm_pos;
	mr.matched_len -= (left_ok) ? sclen_left : remain_beg;
	mr.qspos += (left_ok) ? sclen_left : remain_beg;
	mr.sclen_left = sclen_left;

	err = err_left;

	return left_ok;
}

bool extend_chain_right(const chain_t& ch, char* seq, int seq_len, int ub, MatchedMate& mr, int& err) {
	bool right_ok = true;
	int sclen_right = 0;
	int err_right = 0;

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len);

	right_ok = (remain_end <= 0);

	char remain_str_end[remain_end+5];
	if (remain_end > 0) {
		right_ok = extend_right(seq + seq_len - remain_end, rm_pos, remain_end, ub, err_right, sclen_right);
	}

	mr.epos = rm_pos;
	mr.matched_len -= (right_ok)? sclen_right : remain_end;
	mr.qepos -= (right_ok)? sclen_right : remain_end;
	mr.sclen_right = sclen_right;

	err = err_right;

	return right_ok;
}

void update_match_mate_info(bool lok, bool rok, int lerr, int rerr, MatchedMate& mm) {
	mm.left_ok = lok;
	mm.right_ok = rok;
	if (lok and rok and (lerr + rerr <= EDTH)) {
		mm.is_concord = true;
		mm.type = CONCRD;
	}
	else if (lok or rok) {
		mm.type = CANDID;
	}
	else
		mm.type = ORPHAN;
}

void extend_both_mates(const chain_t& lch, const chain_t& rch, char* lseq, char* rseq, int lseq_len, int rseq_len, MatchedMate& lmm, MatchedMate& rmm) {
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

	int llerr;
	int lrerr;
	int rlerr;
	int rrerr;

	if (l_extend) {
		lmm.matched_len = lseq_len;
		lmm.qspos = 1;
		lmm.qepos = lseq_len;
		llok = extend_chain_left(lch, lseq, lseq_len, MINLB, lmm, llerr);
	}

	if (r_extend) {
		rmm.matched_len = rseq_len;
		rmm.qspos = 1;
		rmm.qepos = rseq_len;
		rlok = extend_chain_left(rch, rseq, rseq_len, (l_extend) ? lmm.spos : MINLB, rmm, rlerr);
	}
	
	if (r_extend) {
		rrok = extend_chain_right(rch, rseq, rseq_len, MAXUB, rmm, rrerr);
	}
	
	if (l_extend) {
		lrok = extend_chain_right(lch, lseq, lseq_len, (r_extend) ? rmm.epos : MAXUB, lmm, lrerr);
	}
	
	if (l_extend) {
		update_match_mate_info(llok, lrok, llerr, lrerr, lmm);
	}

	if (r_extend) {
		update_match_mate_info(rlok, rrok, rlerr, rrerr, rmm);
	}
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
	
	char remain_str_beg[remain_beg+5];
	if (remain_beg > 0) {
		left_ok = extend_left(seq, lm_pos, remain_beg, MINLB, err_left, sclen_left);
	}

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len);

	right_ok = (remain_end <= 0);

	char remain_str_end[remain_end+5];
	if (remain_end > 0) {
		right_ok = extend_right(seq + seq_len - remain_end, rm_pos, remain_end, MAXUB, err_right, sclen_right);
	}

	mr.spos = lm_pos;
	mr.epos = rm_pos;
	mr.matched_len = seq_len;
	mr.matched_len -= (left_ok) ? sclen_left : remain_beg;
	mr.matched_len -= (right_ok) ? sclen_right: remain_end;

	mr.qspos = 1 + (left_ok) ? sclen_left : remain_beg;
	mr.qepos = seq_len - (right_ok) ? sclen_right: remain_end;

	mr.dir = dir;
	
	if (left_ok and right_ok and (err_left + err_right <= EDTH)) {
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
					if (!(near_border[contigName[0]-'1'][k] & 2)) {
						same_intron = false;
						break;
					}
				if (same_intron)
					return true;
			}
	
	return false;
}

// sm should start before lm
bool concordant_explanation(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm) {
	if (sm.spos > lm.spos)
		return false;

	int32_t tlen;
	int32_t intron_gap;
	if (sm.exons_spos == NULL or lm.exons_spos == NULL) {
		tlen = lm.spos - sm.epos - 1 + lm.matched_len + sm.matched_len;
		if (tlen <= MAXTLEN)
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
					if (tlen <= MAXTLEN)
						mr.update(sm, lm, chr, shift, tlen, 0, true, CONCRD, r1_sm);
					else
						mr.update(sm, lm, chr, shift, tlen, 0, true, DISCRD, r1_sm);
				}
	}

	if (sm.exons_epos == NULL or lm.exons_spos == NULL) {
		tlen = lm.spos - sm.epos - 1 + sm.matched_len + lm.matched_len;
		if (tlen <= MAXTLEN)
			mr.update(sm, lm, chr, shift, tlen, 0, false, CONCRD, r1_sm);
		else if (tlen <= MAXDISCRDTLEN)
			mr.update(sm, lm, chr, shift, tlen, 0, false, DISCRD, r1_sm);
	}
	else {
		for (int i = 0; i < sm.exons_epos->seg_list.size(); i++)
			for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
				if (sm.exons_epos->seg_list[i].same_exon(lm.exons_spos->seg_list[j])) {
					// end1 and start2 on same exon
					tlen = lm.spos - sm.epos - 1 + sm.matched_len + lm.matched_len;
					if (tlen <= MAXTLEN)
						mr.update(sm, lm, chr, shift, tlen, 0, true, CONCRD, r1_sm);
					else
						mr.update(sm, lm, chr, shift, tlen, 0, true, DISCRD, r1_sm);
				}
				else if (lm.exons_spos->seg_list[j].next_exon(sm.exons_epos->seg_list[i])) {
					// start2 on next exon
					intron_gap = lm.exons_spos->seg_list[j].start - sm.exons_epos->seg_list[i].end;
					tlen = lm.spos - sm.epos - (lm.exons_spos->seg_list[j].start - sm.exons_epos->seg_list[i].end) + sm.matched_len + lm.matched_len;
					if (tlen <= MAXTLEN)
						mr.update(sm, lm, chr, shift, tlen, 1, true, CONCRD, r1_sm);
					else
						mr.update(sm, lm, chr, shift, tlen, 1, true, DISCRD, r1_sm);
				}
				else if (lm.exons_spos->seg_list[j].same_gene(sm.exons_epos->seg_list[i])) {
					tlen = lm.spos - sm.epos - 1 + sm.matched_len + lm.matched_len;
					mr.update(sm, lm, chr, shift, tlen, 2, true, DISCRD, r1_sm);
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
		if ((near_border[contigName[0]-'1'][sm.spos] & 2) and ((near_border[contigName[0]-'1'][lm.spos] & 2)) and 
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

// bool are_chimeric(vector <MatchedMate>& mms, int mms_size, MatchedMate& mm, MatchedRead& mr) {
// 	if (mms_size <= 0 or (mm.dir * mms[0].dir) == 1)	// empty / same orientation
// 		return false;

// 	ConShift con_shift = gtf_parser.get_shift(contigName, mm.spos);
// 	uint32_t shift = con_shift.shift;
	
// 	overlap_to_epos(mm);
// 	overlap_to_spos(mm);

// 	for (int i = 0; i < mms_size; i++) {
// 		MatchedMate omm = mms[i];

// 		overlap_to_spos(omm);
// 		overlap_to_epos(omm);
		
// 		//if (mm.spos <= omm.spos) {
// 		if (mm.dir == 1) {
// 			if (check_bsj(mm, omm, mr, con_shift.contig, shift))
// 				return true;
// 		}
// 		else {
// 			if(check_bsj(omm, mm, mr, con_shift.contig, shift))
// 				return true;
// 		}
// 	}
	
// 	return false;
// }

// bool are_concordant(vector <MatchedMate>& mms, int mms_size, MatchedMate& mm, MatchedRead& mr) {
// 	if (mms_size <= 0 or (mm.dir * mms[0].dir) == 1)	// empty / same orientation
// 		return false;

// 	ConShift con_shift = gtf_parser.get_shift(contigName, mm.spos);
// 	uint32_t shift = con_shift.shift;
	
// 	overlap_to_epos(mm);
// 	overlap_to_spos(mm);

// 	for (int i = 0; i < mms_size; i++) {
// 		MatchedMate omm = mms[i];

// 		overlap_to_spos(omm);
// 		overlap_to_epos(omm);
		
// 		if (mm.dir == 1) {
// 			if (concordant_explanation(mm, omm, mr, con_shift.contig, shift))
// 				return true;
// 			check_chimeric(omm, mm, mr, con_shift.contig, shift);
// 		}


// 		if (mm.dir == -1) {
// 			if (concordant_explanation(omm, mm, mr, con_shift.contig, shift))
// 				return true;
// 			check_chimeric(mm, omm, mr, con_shift.contig, shift);
// 		}
// 	}
	
// 	return false;
// }

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

			if ((forward_exon_list[i] != NULL and reverse_exon_list[j] != NULL and same_gene(forward_exon_list[i], reverse_exon_list[j]))
				or (forward_exon_list[i] != NULL and same_gene(forward_exon_list[i], rs, re))
				or (reverse_exon_list[j] != NULL and same_gene(reverse_exon_list[j], fs, fe))
				or (tlen <= MAXDISCRDTLEN)) {
				//or (forward_exon_list[i] == NULL and reverse_exon_list[j] == NULL and tlen <= MAXTLEN)) {
				MatePair temp;
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
		MatchedMate r1_mm;
		MatchedMate r2_mm;

		r1_mm.dir = 1;
		r2_mm.dir = -1;

		uint32_t forward_start = mate_pairs[i].forward.frags[0].rpos;
		uint32_t reverse_start = mate_pairs[i].reverse.frags[0].rpos;
		uint32_t reverse_end   = mate_pairs[i].reverse.frags[mate_pairs[i].reverse.chain_len-1].rpos + mate_pairs[i].reverse.frags[mate_pairs[i].reverse.chain_len-1].len - 1;
		
		if (forward_start <= reverse_end) {
			extend_both_mates(mate_pairs[i].forward, mate_pairs[i].reverse, forward_rec->seq, backward_rec->rcseq, forward_rec->seq_len, backward_rec->seq_len, r1_mm, r2_mm);
			
			if (r1_mm.type == CONCRD and r2_mm.type == CONCRD) {
				ConShift con_shift = gtf_parser.get_shift(contigName, r1_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				if (concordant_explanation(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward)) {
					return CONCRD;
				}
			}
			
			// potentially back splice junction?
			else if ((r1_mm.type == CANDID and r2_mm.type == CONCRD) or (r1_mm.type == CONCRD and r2_mm.type == CANDID)) {
				ConShift con_shift = gtf_parser.get_shift(contigName, r1_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				check_bsj(r1_mm, r2_mm, mr, con_shift.contig, con_shift.shift, r1_forward);
			}
		}

		if (forward_start > reverse_start) {
			extend_both_mates(mate_pairs[i].reverse, mate_pairs[i].forward, backward_rec->rcseq, forward_rec->seq, backward_rec->seq_len, forward_rec->seq_len, r2_mm, r1_mm);
			
			if (r1_mm.type == CONCRD and r2_mm.type == CONCRD) {
				ConShift con_shift = gtf_parser.get_shift(contigName, r2_mm.spos);
				
				overlap_to_epos(r1_mm);
				overlap_to_spos(r1_mm);

				overlap_to_epos(r2_mm);
				overlap_to_spos(r2_mm);
			
				check_chimeric(r2_mm, r1_mm, mr, con_shift.contig, con_shift.shift, !r1_forward);
			}
			
			// potentially back splice junction?
			else if ((r1_mm.type == CANDID and r2_mm.type == CONCRD) or (r1_mm.type == CONCRD and r2_mm.type == CANDID)) {
				ConShift con_shift = gtf_parser.get_shift(contigName, r2_mm.spos);
				
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

// pos is exclusive
// [ pos+1, pos+len ]
void get_seq_right(char* res_str, char* seq, int seq_len, uint32_t pos, int covered, int remain, uint32_t ub, int& min_ed, int& sclen_best, uint32_t& rmpos_best, bool& consecutive) {
	int len;
	int sclen;
	int indel;
	int edit_dist;
	uint32_t new_rmpos;
	uint32_t exon_remain;

	bool dummy_consec = false;

	// if covered > 0 -> pos = (next_exon_beg - 1) = > search on (pos + 1)
	uint32_t search_pos = (covered == 0) ? pos : pos + 1;
	const IntervalInfo<UniqSeg>* overlapped_exon = gtf_parser.get_location_overlap(search_pos, true);

	if (overlapped_exon == NULL or overlapped_exon->seg_list.size() == 0) { // there is enough distance to the end of exon / intron
		//vafprintf(2, stderr, "Going for %lu - %lu\n", pos + 1, pos + remain);
		
		pac2char(pos + 1, remain, res_str);

		len = covered + remain;
		//edit_dist = alignment.hamming_distance_right(res_str - covered, len, seq, seq_len, sclen);
		edit_dist = alignment.local_alignment_right(res_str - covered, len, seq, seq_len, sclen, indel);
	
		consecutive = true;
		new_rmpos = pos + seq_len - indel;
		
		//vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\n", new_rmpos, seq_len, indel);
		//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist: %d\n", res_str - covered, seq, edit_dist);
		
		if (edit_dist <= EDTH)
			if ((edit_dist < min_ed) or (edit_dist == min_ed and sclen < sclen_best) or (edit_dist == min_ed and sclen == sclen_best and new_rmpos < rmpos_best)) {	// less hamming distance, then less soft clip, then smaller junction distance
				min_ed = edit_dist;
				sclen_best = sclen;
				rmpos_best = new_rmpos;
		}

		return;
	}

	set <GenRegion> trans_extensions;
	//trans_extensions.clear();
	GenRegion new_region;

	for (int i = 0; i < overlapped_exon->seg_list.size(); i++) {
		//vafprintf(2, stderr, "This exon: [%u-%u] -> %u\n", overlapped_exon->seg_list[i].start, overlapped_exon->seg_list[i].end, overlapped_exon->seg_list[i].next_exon_beg);
		if (covered > 0 and overlapped_exon->seg_list[i].start != search_pos)	// when jump to the next exon
			continue;

		exon_remain = overlapped_exon->seg_list[i].end - pos;
		//vafprintf(2, stderr, "Remain on exon: %u\n", exon_remain);
		if (exon_remain >= remain) {	// exonic	
			new_region.set(pos + remain, 0);
			if (trans_extensions.find(new_region) != trans_extensions.end())
				continue;

			trans_extensions.insert(new_region);

			//vafprintf(2, stderr, "Going for %lu - %lu\n", pos + 1, pos + remain);
			
			pac2char(pos + 1, remain, res_str);

			len = covered + remain;

			//edit_dist = alignment.hamming_distance_right(res_str - covered, len, seq, len, sclen);
			edit_dist = alignment.local_alignment_right(res_str - covered, len, seq, seq_len, sclen, indel);
		
			consecutive = true;
			new_rmpos = pos + seq_len - indel;
			
			//vafprintf(2, stderr, "rmpos: %lu\textend len: %d\n", new_rmpos, len);
			//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist: %d\n", res_str - covered, seq, edit_dist);
			
			if (edit_dist <= EDTH)
				if ((edit_dist < min_ed) or (edit_dist == min_ed and sclen < sclen_best) or (edit_dist == min_ed and sclen == sclen_best and new_rmpos < rmpos_best)) {	// less hamming distance, then less soft clip, then smaller junction distance
					min_ed = edit_dist;
					sclen_best = sclen;
					rmpos_best = new_rmpos;
				}
			}
		else {		// junction
			if (overlapped_exon->seg_list[i].next_exon_beg == 0)	// should not be the last exon and partially mappable
				continue;
			
			if (overlapped_exon->seg_list[i].next_exon_beg > ub)	// do not allow junction further than ub
				continue;

			new_region.set(overlapped_exon->seg_list[i].end, overlapped_exon->seg_list[i].next_exon_beg);
			if (trans_extensions.find(new_region) != trans_extensions.end())
				continue;

			trans_extensions.insert(new_region);

			pac2char(pos + 1, exon_remain, res_str);
			get_seq_right(res_str + exon_remain, seq, seq_len, overlapped_exon->seg_list[i].next_exon_beg - 1, covered + exon_remain, remain - exon_remain, ub, min_ed, sclen_best, rmpos_best, dummy_consec);
		}
	}
}

// pos is exclusive
// [ pos+1, pos+len ]
bool extend_right(char* seq, uint32_t& pos, int len, uint32_t ub, int& err, int& sclen) {
	int seq_len = len;
	int ref_len = len + INDELTH;

	char res_str[ref_len+1];
	res_str[ref_len] = '\0';
	
	uint32_t best_rmpos = 0;
	int min_ed = EDTH + 1;
	int sclen_best = SOFTCLIPTH;
	bool consecutive = false;
	int indel;
	err = EDTH + 1;
	
	get_seq_right(res_str, seq, seq_len, pos, 0, ref_len, ub, min_ed, sclen_best, best_rmpos, consecutive);

	if (min_ed <= EDTH) {
		pos = best_rmpos - sclen_best;
		err = min_ed;
		sclen = sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew RM POS: %u\n", min_ed, pos);
		return true;
	}
	
	if (consecutive)
		return false;

	// intron retentaion
	pac2char(pos + 1, ref_len, res_str);	
	//min_ed = alignment.hamming_distance_right(res_str, len, seq, len, sclen_best);
	min_ed = alignment.local_alignment_right(res_str, ref_len, seq, seq_len, sclen_best, indel);
	
	//vafprintf(2, stderr, "Intron Retention:\nrmpos: %lu\textend len: %d\n", pos, len);
	//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist %d\n", res_str, seq, min_ed);

	if (min_ed <= EDTH) {
		pos = pos + seq_len - indel - sclen_best;
		err = min_ed;
		sclen = sclen_best;
		vafprintf(2, stderr, "Intron Retention: Min Edit Dist: %d\tNew RM POS: %u\n", min_ed, pos);
		return true;
	}

	return false;
}

// pos is exclusive
// [ pos-len, pos-1 ]
void get_seq_left(char* res_str, char* seq, int seq_len, uint32_t pos, int covered, int remain, uint32_t lb, int& min_ed, int& sclen_best, uint32_t& lmpos_best, bool& consecutive) {
	int len;
	int sclen;
	int indel;
	int edit_dist;
	uint32_t new_lmpos;
	uint32_t exon_remain;

	bool dummy_consec = false;
	
	// if covered > 0 -> pos = (prev_exon_end + 1) = > search on (pos - 1)
	uint32_t search_pos = (covered == 0) ? pos : pos - 1;
	const IntervalInfo <UniqSeg>* overlapped_exon = gtf_parser.get_location_overlap(search_pos, false);

	if (overlapped_exon == NULL or overlapped_exon->seg_list.size() == 0) {
		//vafprintf(2, stderr, "Going for %lu - %lu\n", pos - remain, pos - 1);

		pac2char(pos - remain, remain, res_str);

		len = covered + remain;
		//edit_dist = alignment.hamming_distance_left(res_str, len, seq, len, sclen);
		edit_dist = alignment.local_alignment_left(res_str, len, seq, seq_len, sclen, indel);
		
		consecutive = true;
		new_lmpos = pos - seq_len + indel;

		//vafprintf(2, stderr, "lmpos: %lu\textend len: %d\t indel: %d\n", new_lmpos, seq_len, indel);
		//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nLeftedit dist: %d\n", res_str, seq, edit_dist);

		if (edit_dist <= EDTH)
			if ((edit_dist < min_ed) or (edit_dist == min_ed and sclen < sclen_best) or (edit_dist == min_ed and sclen == sclen_best and new_lmpos > lmpos_best)) {	// less hamming distance, then less soft clip, then smaller junction distance
				min_ed = edit_dist;
				sclen_best = sclen;
				lmpos_best = new_lmpos;
		}

		return;
	}

	set <GenRegion> trans_extensions;
	//trans_extensions.clear();
	GenRegion new_region;

	for (int i = 0; i < overlapped_exon->seg_list.size(); i++) {
		//vafprintf(2, stderr, "%u <- This exon: [%u-%u]\n", overlapped_exon->seg_list[i].prev_exon_end, overlapped_exon->seg_list[i].start, overlapped_exon->seg_list[i].end);
		if (covered > 0 and overlapped_exon->seg_list[i].end != search_pos)	// when jump to prev exon
			continue;

		exon_remain = pos - overlapped_exon->seg_list[i].start;
		//vafprintf(2, stderr, "Remain on exon: %u\n", exon_remain);
		if (exon_remain >= remain) {
			new_region.set(pos - remain, 0);
			if (trans_extensions.find(new_region) != trans_extensions.end())
				continue;

			trans_extensions.insert(new_region);

			//vafprintf(2, stderr, "Going for %lu - %lu\n", pos - remain, pos - 1);

			pac2char(pos - remain, remain, res_str);

			len = covered + remain;
			//edit_dist = alignment.hamming_distance_left(res_str, len, seq, len, sclen);
			edit_dist = alignment.local_alignment_left(res_str, len, seq, seq_len, sclen, indel);
		
			consecutive = true;
			new_lmpos = pos - seq_len + indel;

			//vafprintf(2, stderr, "lmpos: %lu\textend len: %d\n", new_lmpos, seq_len);
			//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nLeft edit dist: %d\n", res_str, seq, edit_dist);

			if (edit_dist <= EDTH)
				if ((edit_dist < min_ed) or (edit_dist == min_ed and sclen < sclen_best) or (edit_dist == min_ed and sclen == sclen_best and new_lmpos > lmpos_best)) {	// less hamming distance, then less soft clip, then smaller junction distance
					min_ed = edit_dist;
					sclen_best = sclen;
					lmpos_best = new_lmpos;
				}
			}
		else {		// junction
			if (overlapped_exon->seg_list[i].prev_exon_end == 0)	// should not be the first exon and partially mappable
				continue;

			if (overlapped_exon->seg_list[i].prev_exon_end < lb)
				continue;

			new_region.set(overlapped_exon->seg_list[i].start, overlapped_exon->seg_list[i].prev_exon_end);
			if (trans_extensions.find(new_region) != trans_extensions.end())
				continue;

			trans_extensions.insert(new_region);

			pac2char(overlapped_exon->seg_list[i].start, exon_remain, res_str + remain - exon_remain);
			get_seq_left(res_str, seq, seq_len, overlapped_exon->seg_list[i].prev_exon_end + 1, covered + exon_remain, remain - exon_remain, lb, min_ed, sclen_best, lmpos_best, dummy_consec);
		}
	}
}

// pos is exclusive
// [ pos-len, pos-1 ]
bool extend_left(char* seq, uint32_t& pos, int len, uint32_t lb, int& err, int& sclen) {
	int seq_len = len;
	int ref_len = len + INDELTH;
	char res_str[ref_len+1];
	res_str[ref_len] = '\0';
	
	uint32_t lmpos_best = 0;
	int min_ed = EDTH + 1;
	int sclen_best = SOFTCLIPTH;
	bool consecutive = false;
	int indel;
	err = EDTH + 1;

	get_seq_left(res_str, seq, seq_len, pos, 0, ref_len, lb, min_ed, sclen_best, lmpos_best, consecutive);
	
	if (min_ed <= EDTH) {
		pos = lmpos_best + sclen_best;
		err = min_ed;
		sclen = sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\n", min_ed, pos);
		return true;
	}

	if (consecutive)
		return false;

	// intron retentaion
	pac2char(pos - ref_len, ref_len, res_str);
	//min_ed = alignment.hamming_distance_left(res_str, len, seq, len, sclen_best);
	min_ed = alignment.local_alignment_left(res_str, ref_len, seq, seq_len, sclen_best, indel);
	
	//vafprintf(2, stderr, "Intron Retention:\nlmpos: %lu\textend len: %d\n", pos, len);
	//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist %d\n", res_str, seq, min_ed);

	if (min_ed <= EDTH) {
		pos = pos - seq_len + indel + sclen_best;
		err = min_ed;
		sclen = sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\n", min_ed, pos);
		return true;
	}

	return false;
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
