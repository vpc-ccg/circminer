#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "filter.h"
#include "align.h"
#include "common.h"
#include "gene_annotation.h"

extern "C" {
#include "mrsfast/Common.h"
}

void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer*& frag_l);
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedMate& mr, int dir);
int process_mates(const chain_list& forward_chain, const Record* record1, const chain_list& backward_chain, const Record* record2);

bool extend_right(char* seq, uint32_t& pos, int len);
bool extend_left(char* seq, uint32_t& pos, int len);

void overlap_to_epos(MatchedMate& mr);
void overlap_to_spos(MatchedMate& mr);

// updates next fq file to be read if need be (keep file)
FilterRead::FilterRead (char* save_fname, bool pe, char* filter_temp_name, int num_files, char* fq_file1, char* fq_file2) {
	is_pe = pe;
	cat_count = num_files;

	char* output_names[8] = { "concordant", "discordant", "circ_RF", "fusion", "circ_bsj", "keep", "OEA", "orphan" };
	char cat_fname [FILE_NAME_LENGTH];
	
	for (int i = 0; i < num_files; i++) {
		if (is_pe) {
			sprintf(cat_fname, "%s_%s.%s_R1.fastq", save_fname, filter_temp_name, output_names[i]);
			cat_file_r1[i] = open_file(cat_fname, "w");
			sprintf(cat_fname, "%s_%s.%s_R2.fastq", save_fname, filter_temp_name, output_names[i]);
			cat_file_r2[i] = open_file(cat_fname, "w");
		}
		else {
			sprintf(cat_fname, "%s_%s.%s.fastq", save_fname, filter_temp_name, output_names[i]);
			cat_file_r1[i] = open_file(cat_fname, "w");
		}
	}

	// updating fq file to be read in next round
	if (is_pe) {
		sprintf(fq_file1, "%s_%s.%s_R1.fastq", save_fname, filter_temp_name, output_names[DISCRD]);
		sprintf(fq_file2, "%s_%s.%s_R2.fastq", save_fname, filter_temp_name, output_names[DISCRD]);
	}
	else {
		sprintf(fq_file1, "%s_%s.%s.fastq", save_fname, filter_temp_name, output_names[DISCRD]);
	}
}

FilterRead::~FilterRead (void) {
	for (int i = 0; i < CATNUM; i++) {
		close_file(cat_file_r1[i]);
		close_file(cat_file_r2[i]);
	}
}

// SE mode
int FilterRead::process_read (	Record* current_record, int kmer_size, GIMatchedKmer*& fl, GIMatchedKmer*& bl, 
								chain_list& forward_best_chain, chain_list& backward_best_chain) {
	
	vafprintf(1, stderr, "%s\n", current_record->rname);

	MatchedMate mr;
	int ex_ret, min_ret = ORPHAN;

	get_best_chains(current_record->seq, current_record->seq_len, kmer_size, forward_best_chain, fl);
	for (int i = 0; i < forward_best_chain.best_chain_count; i++) {
		ex_ret = extend_chain(forward_best_chain.chains[i], current_record->seq, current_record->seq_len, mr, 1); 
		
		if (ex_ret == CONCRD)
			return CONCRD;

		if (ex_ret < min_ret)
			min_ret = ex_ret;
	}

	get_best_chains(current_record->rcseq, current_record->seq_len, kmer_size, backward_best_chain, bl);
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
	
	// R1
	vafprintf(1, stderr, "R1/%s\n", current_record1->rname);

	get_best_chains(current_record1->seq, current_record1->seq_len, kmer_size, forward_best_chain_r1, fl);
	get_best_chains(current_record1->rcseq, current_record1->seq_len, kmer_size, backward_best_chain_r1, bl);

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

	get_best_chains(current_record2->seq, current_record2->seq_len, kmer_size, forward_best_chain_r2, fl);
	get_best_chains(current_record2->rcseq, current_record2->seq_len, kmer_size, backward_best_chain_r2, bl);

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
	if (forward_best_chain_r1.best_chain_count + backward_best_chain_r1.best_chain_count + forward_best_chain_r2.best_chain_count + backward_best_chain_r2.best_chain_count <= 0)
		return ORPHAN;
	if ((forward_best_chain_r1.best_chain_count + backward_best_chain_r1.best_chain_count <= 0) or (forward_best_chain_r2.best_chain_count + backward_best_chain_r2.best_chain_count <= 0))
		return OEANCH;

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
		attempt1 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2);
		if (attempt1 == CONCRD)
			return CONCRD;

		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt2 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1);
		if (attempt2 == CONCRD)
			return CONCRD;
		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
	else {
		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt1 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1);
		if (attempt1 == CONCRD)
			return CONCRD;

		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt2 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2);
		if (attempt2 == CONCRD)
			return CONCRD;
		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
}

// write reads SE mode
void FilterRead::write_read_category (Record* current_record, int state) {
	state = minM(state, current_record->state);
	int cat = (state >= cat_count) ? DISCRD : state;
	fprintf(cat_file_r1[cat], "%s\n%s%s%d\n%s", current_record->rname, current_record->seq, current_record->comment, state, current_record->qual);
}

// write reads PE mode
void FilterRead::write_read_category (Record* current_record1, Record* current_record2, int state) {
	state = minM(state, current_record1->state);
	int cat = (state >= cat_count) ? DISCRD : state;
	fprintf(cat_file_r1[cat], "%s\n%s%s%d\n%s", current_record1->rname, current_record1->seq, current_record1->comment, state, current_record1->qual);
	fprintf(cat_file_r2[cat], "%s\n%s%s%d\n%s", current_record2->rname, current_record2->seq, current_record2->comment, state, current_record2->qual);
}

void print_mapping(char* rname, const MatchedRead& mr) {
	if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ) {
		fprintf(outputJuncFile, "%s\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%d\t%d\t%d\n", rname, mr.chr.c_str(), 
																			mr.spos_r1, mr.epos_r1, mr.mlen_r1, 
																			mr.spos_r2, mr.epos_r2, mr.mlen_r2, 
																			mr.tlen, mr.junc_num, mr.gm_compatible);
	}
}

bool is_concord(const chain_t& a, int seq_len, MatchedMate& mr) {
	if (a.chain_len < 2) {
		mr.is_concord = false;
	}
	else if ((a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos) >= seq_len) {
		mr.is_concord = true;
		mr.type = CONCRD;
		
		mr.start_pos = a.frags[0].rpos;
		mr.end_pos = a.frags[a.chain_len-1].rpos + a.frags[a.chain_len-1].len - 1;
		mr.matched_len = a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos;
	}
	else {
		mr.is_concord = false;
	}
	return mr.is_concord;
}

void get_best_chains(char* read_seq, int seq_len, int kmer_size, chain_list& best_chain, GIMatchedKmer*& frag_l) {
	int kmer_count = ceil(seq_len / kmer_size);
	int forward_fragment_count, backward_fragment_count;

	split_match_hash(read_seq, seq_len, kmer_size, frag_l);
	chain_seeds_sorted_kbest(seq_len, frag_l, best_chain);
	//chain_seeds_sorted_kbest_old(frag_l, best_chain);
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

	uint32_t lm_pos = ch.frags[0].rpos;
	int remain_beg = ch.frags[0].qpos;
	
	left_ok = (remain_beg <= 0);
	
	char remain_str_beg[remain_beg+5];
	if (remain_beg > 0) {
		left_ok = extend_left(seq, lm_pos, remain_beg);
	}

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len);

	right_ok = (remain_end <= 0);

	char remain_str_end[remain_end+5];
	if (remain_end > 0) {
		right_ok = extend_right(seq + seq_len - remain_end, rm_pos, remain_end);
	}


	mr.start_pos = lm_pos;
	mr.end_pos = rm_pos;
	//mr.matched_len = seq_len - sclen_left - sclen_right;
	mr.dir = dir;
	
	if (left_ok and right_ok) {
		mr.is_concord = true;
		mr.type = CONCRD;
		mr.matched_len = seq_len;
	}
	else if (left_ok or right_ok) {
		mr.type = CANDID;
		mr.matched_len = seq_len - kmer;
	}
	else
		mr.type = ORPHAN;

	return mr.type;
}

// s < l
int32_t calc_tlen(const UniqSeg& s, const UniqSeg& l) {
	int32_t min_tlen;
				
}

// sm should start before lm
bool concordant_explanation(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift) {
	if (sm.start_pos > lm.start_pos)
		return false;

	int32_t tlen;
	int32_t intron_gap;
	if (sm.exons_spos == NULL or lm.exons_spos == NULL) {
		tlen = lm.start_pos - sm.end_pos - 1 + lm.matched_len + sm.matched_len;
		if (tlen <= MAXTLEN)
			mr.update(sm, lm, chr, shift, tlen, 0, false, CONCRD);
	}
	else {
		//fprintf(stderr, "Left Mate [%u-%u] dir=%d, type=%d, Right Mate[%u-%u] dir=%d, type=%d\n", sm.start_pos, sm.end_pos, sm.dir, sm.type, lm.start_pos, lm.end_pos, lm.dir, lm.type);
	// starts on same exon
		for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
			for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
				if (sm.exons_spos->seg_list[i].same_exon(lm.exons_spos->seg_list[j])) {
					// => assume genomic locations
					tlen = lm.start_pos + lm.matched_len - sm.start_pos;
					if (tlen <= MAXTLEN)
						mr.update(sm, lm, chr, shift, tlen, 0, true, CONCRD);
					else
						mr.update(sm, lm, chr, shift, tlen, 0, true, DISCRD);
				}
	}

	if (sm.exons_epos == NULL or lm.exons_spos == NULL) {
		tlen = lm.start_pos - sm.end_pos - 1 + sm.matched_len + lm.matched_len;
		if (tlen <= MAXTLEN)
			mr.update(sm, lm, chr, shift, tlen, 0, false, CONCRD);
	}
	else {
		for (int i = 0; i < sm.exons_epos->seg_list.size(); i++)
			for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
				if (sm.exons_epos->seg_list[i].same_exon(lm.exons_spos->seg_list[j])) {
					// end1 and start2 on same exon
					tlen = lm.start_pos - sm.end_pos - 1 + sm.matched_len + lm.matched_len;
					if (tlen <= MAXTLEN)
						mr.update(sm, lm, chr, shift, tlen, 0, true, CONCRD);
					else
						mr.update(sm, lm, chr, shift, tlen, 0, true, DISCRD);
				}
				else if (lm.exons_spos->seg_list[j].next_exon(sm.exons_epos->seg_list[i])) {
					// start2 on next exon
					intron_gap = lm.exons_spos->seg_list[j].start - sm.exons_epos->seg_list[i].end;
					tlen = lm.start_pos - sm.end_pos - (lm.exons_spos->seg_list[j].start - sm.exons_epos->seg_list[i].end) + sm.matched_len + lm.matched_len;
					if (tlen <= MAXTLEN)
						mr.update(sm, lm, chr, shift, tlen, 1, true, CONCRD);
					else
						mr.update(sm, lm, chr, shift, tlen, 1, true, DISCRD);
				}
	}

	if (mr.type == CONCRD)
		return true;
	else
		return false;
}

bool check_chimeric(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift) {
	if (mr.type == CONCRD)
		return false;

	if (sm.exons_spos == NULL or lm.exons_spos == NULL)
		return false;

	//fprintf(stderr, "Left Mate [%u-%u] dir=%d, type=%d, Right Mate[%u-%u] dir=%d, type=%d\n", sm.start_pos, sm.end_pos, sm.dir, sm.type, lm.start_pos, lm.end_pos, lm.dir, lm.type);
	for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
		for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
			if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j]) and sm.start_pos < lm.start_pos) {
				mr.update(sm, lm, chr, shift, lm.end_pos - sm.start_pos + 1, 0, false, CHIORF);
				return true;
			}
	return false;
}

bool check_bsj(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift) {
	if (mr.type == CONCRD or mr.type == DISCRD)
		return false;

	if (sm.exons_spos == NULL or lm.exons_spos == NULL)
		return false;

	for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
		for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
			if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j])) {
				mr.update(sm, lm, chr, shift, lm.end_pos - sm.start_pos + 1, 0, false, CHIBSJ);
				return true;
			}
	return false;
}

bool are_chimeric(vector <MatchedMate>& mms, int mms_size, MatchedMate& mm, MatchedRead& mr) {
	if (mms_size <= 0 or mm.dir * mms[0].dir == 1)	// empty / same orientation
		return false;

	ConShift con_shift = gtf_parser.get_shift(contigName, mm.start_pos);
	uint32_t shift = con_shift.shift;
	
	//vafprintf(0, stderr, "---%s\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%d\tdir:%d\n", rname, mr.chr, mr.start_pos, mr.end_pos, mr.matched_len, 
	//																		omr.start_pos, omr.end_pos, omr.matched_len, omr.end_pos - mr.start_pos + 1, mr.dir);
		
	overlap_to_epos(mm);
	overlap_to_spos(mm);

	for (int i = 0; i < mms_size; i++) {
		MatchedMate omm = mms[i];

		overlap_to_spos(omm);
		overlap_to_epos(omm);
		
		if (mm.dir == 1) {
			if (check_bsj(mm, omm, mr, con_shift.contig, shift))
				return true;
		}


		if (mm.dir == -1) {
			if(check_bsj(omm, mm, mr, con_shift.contig, shift))
				return true;
		}
	}
	
	return false;
}

bool are_concordant(vector <MatchedMate>& mms, int mms_size, MatchedMate& mm, MatchedRead& mr) {
	if (mms_size <= 0 or mm.dir * mms[0].dir == 1)	// empty / same orientation
		return false;

	ConShift con_shift = gtf_parser.get_shift(contigName, mm.start_pos);
	uint32_t shift = con_shift.shift;
	
	//vafprintf(0, stderr, "---%s\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%d\tdir:%d\n", rname, mr.chr, mr.start_pos, mr.end_pos, mr.matched_len, 
	//																		omr.start_pos, omr.end_pos, omr.matched_len, omr.end_pos - mr.start_pos + 1, mr.dir);
		
	overlap_to_epos(mm);
	overlap_to_spos(mm);

	for (int i = 0; i < mms_size; i++) {
		MatchedMate omm = mms[i];

		overlap_to_spos(omm);
		overlap_to_epos(omm);
		
		if (mm.dir == 1) {
			if (concordant_explanation(mm, omm, mr, con_shift.contig, shift))
				return true;
			check_chimeric(omm, mm, mr, con_shift.contig, shift);
		}


		if (mm.dir == -1) {
			if (concordant_explanation(omm, mm, mr, con_shift.contig, shift))
				return true;
			check_chimeric(mm, omm, mr, con_shift.contig, shift);
		}
	}
	
	return false;
}

int process_mates(const chain_list& forward_chain, const Record* record1, const chain_list& backward_chain, const Record* record2) {
	int fc_size = forward_chain.best_chain_count;
	int bc_size = backward_chain.best_chain_count;

	int max_len = (fc_size >= bc_size) ? fc_size : bc_size;

	vector <MatchedMate> forward_mml(fc_size);
	vector <MatchedMate> backward_mml(bc_size);
	
	vector <MatchedMate> forward_mml_partial(fc_size);
	vector <MatchedMate> backward_mml_partial(bc_size);

	MatchedRead mr;

	int ex_ret;
	int fmml_count = 0;
	int bmml_count = 0;
	int fmmlp_count = 0;
	int bmmlp_count = 0;
	int min_ret1 = ORPHAN;
	int min_ret2 = ORPHAN;

	// concordant?
	for (int i = 0; i < max_len; i++) {
		if (i < fc_size) {
			ex_ret = extend_chain(forward_chain.chains[i], record1->seq, record1->seq_len, forward_mml[fmml_count], 1);
			if (ex_ret == CONCRD) {
				if (are_concordant(backward_mml, bmml_count, forward_mml[fmml_count], mr)) {
					//vafprintf(0, stderr, "%s", record1->rname);
					print_mapping(record1->rname, mr);
					return CONCRD;
				}

				fmml_count++;
			}
			else if (ex_ret == CANDID) {
				forward_mml_partial[fmmlp_count] = forward_mml[fmml_count];
				fmmlp_count++;
			}
			if (ex_ret < min_ret1)	min_ret1 = ex_ret;
		}

		if (i < bc_size) {
			ex_ret = extend_chain(backward_chain.chains[i], record2->rcseq, record2->seq_len, backward_mml[bmml_count], -1); 
			if (ex_ret == CONCRD) {
				if (are_concordant(forward_mml, fmml_count, backward_mml[bmml_count], mr)) {
					//vafprintf(0, stderr, "%s", record1->rname);
					print_mapping(record1->rname, mr);
					return CONCRD;
				}
				bmml_count++;
			}
			else if (ex_ret == CANDID) {
				backward_mml_partial[bmmlp_count] = backward_mml[bmml_count];
				bmmlp_count++;
			}
			if (ex_ret < min_ret2)	min_ret2 = ex_ret;
		}
	}

	// chimeric?
	for (int i = 0; i < fmmlp_count; i++) {
		if (are_chimeric(backward_mml, bmml_count, forward_mml_partial[i], mr)) {
			print_mapping(record1->rname, mr);
			return mr.type;
		}
	}
	for (int i = 0; i < bmmlp_count; i++) {
		if (are_chimeric(forward_mml, fmml_count, backward_mml_partial[i], mr)) {
			print_mapping(record1->rname, mr);
			return mr.type;
		}
	}

	if (mr.type == CONCRD or mr.type == DISCRD or mr.type == CHIORF or mr.type == CHIBSJ) {
		print_mapping(record1->rname, mr);
		return mr.type;
	}

	int ret_val = (((min_ret1 == ORPHAN) and (min_ret2 == CONCRD)) or ((min_ret1 == CONCRD) and (min_ret2 == ORPHAN))) ? OEANCH 
			: ((min_ret1 == ORPHAN) or (min_ret2 == ORPHAN)) ? ORPHAN 
			: CANDID;

	return ret_val;
}

// pos is exclusive
// [ pos+1, pos+len ]
void get_seq_right(char* res_str, char* seq, uint32_t pos, int covered, int remain, int& min_ed, int& sclen_best, uint32_t& rmpos_best) {
	int len;
	int sclen;
	int edit_dist;
	uint32_t new_rmpos;
	uint32_t exon_remain;

	// if covered > 0 -> pos = (next_exon_beg - 1) = > search on (pos + 1)
	uint32_t search_pos = (covered == 0) ? pos : pos + 1;
	vector <UniqSeg> overlapped_exon;
	gtf_parser.get_location_overlap(search_pos, overlapped_exon, true);

	if (overlapped_exon.size() == 0) { // there is enough distance to the end of exon / intron
		//vafprintf(2, stderr, "Going for %lu - %lu\n", pos + 1, pos + remain);
		
		pac2char(pos + 1, remain, res_str);
		new_rmpos = pos + remain;

		len = covered + remain;
		edit_dist = alignment.hamming_distance_right(res_str - covered, len, seq, len, sclen);
		
		//vafprintf(2, stderr, "rmpos: %lu\textend len: %d\n", pos, len);
		//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist: %d\n", res_str - covered, seq, edit_dist);
		
		if ((sclen == sclen_best and edit_dist < min_ed) or (sclen < sclen_best and edit_dist <= EDTH)) {	// less soft clip and then less hamming distance
			min_ed = edit_dist;
			sclen_best = sclen;
			rmpos_best = new_rmpos;
		}
	}

	for (int i = 0; i < overlapped_exon.size(); i++) {
		//vafprintf(2, stderr, "This exon: [%u-%u] -> %u\n", overlapped_exon[i].start, overlapped_exon[i].end, overlapped_exon[i].next_exon_beg);
		if (covered > 0 and overlapped_exon[i].start != search_pos)	// when jump to the next exon
			continue;

		exon_remain = overlapped_exon[i].end - pos;
		//vafprintf(2, stderr, "Remain on exon: %u\n", exon_remain);
		if (exon_remain >= remain) {	// exonic	
			//vafprintf(2, stderr, "Going for %lu - %lu\n", pos + 1, pos + remain);
			
			pac2char(pos + 1, remain, res_str);
			new_rmpos = pos + remain;

			len = covered + remain;
			edit_dist = alignment.hamming_distance_right(res_str - covered, len, seq, len, sclen);
			
			//vafprintf(2, stderr, "rmpos: %lu\textend len: %d\n", pos, len);
			//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist: %d\n", res_str - covered, seq, edit_dist);
			
			if ((sclen == sclen_best and edit_dist < min_ed) or (sclen < sclen_best and edit_dist <= EDTH)) {	// less soft clip and then less hamming distance
				min_ed = edit_dist;
				sclen_best = sclen;
				rmpos_best = new_rmpos;
			}
		}
		else {		// junction
			if (overlapped_exon[i].next_exon_beg == 0)	// should not be the last exon and partially mappable
				continue;
			
			pac2char(pos + 1, exon_remain, res_str);
			get_seq_right(res_str + exon_remain, seq, overlapped_exon[i].next_exon_beg - 1, covered + exon_remain, remain - exon_remain, min_ed, sclen_best, rmpos_best);
		}
	}
}

// pos is exclusive
// [ pos+1, pos+len ]
bool extend_right(char* seq, uint32_t& pos, int len) {
	char res_str[len+5];
	bool right_ok = false;
	uint32_t new_rmpos = pos;
	
	int min_ed = EDTH + 1;
	int sclen_best = SOFTCLIPTH;
	uint32_t best_rmpos = 0;
	
	res_str[len] = '\0';
	get_seq_right(res_str, seq, pos, 0, len, min_ed, sclen_best, best_rmpos);

	if (min_ed <= EDTH) {
		pos = best_rmpos;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew RM POS: %u\n", min_ed, best_rmpos);
		return true;
	}
	
	// intron retentaion
	pac2char(pos + 1, len, res_str);	
	new_rmpos = pos + len;
	right_ok = alignment.hamming_match_right(res_str, len, seq, len);
	//vafprintf(2, stderr, "Intron Retention:\nrmpos: %lu\textend len: %d\n", pos, len);
	//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nright ok? %d\n", res_str, seq, right_ok);

	if (right_ok)
		pos = new_rmpos;

	return right_ok;
}

// pos is exclusive
// [ pos-len, pos-1 ]
bool get_seq_left(char* res_str, char* seq, uint32_t pos, int covered, int remain, int& min_ed, int& sclen_best, uint32_t& lmpos_best) {
	int len;
	int sclen;
	int edit_dist;
	uint32_t new_lmpos;
	uint32_t exon_remain;
	
	// if covered > 0 -> pos = (prev_exon_end + 1) = > search on (pos - 1)
	uint32_t search_pos = (covered == 0) ? pos : pos - 1;
	vector <UniqSeg> overlapped_exon;
	gtf_parser.get_location_overlap(search_pos, overlapped_exon, false);

	if (overlapped_exon.size() == 0) {
		//vafprintf(2, stderr, "Going for %lu - %lu\n", pos - remain, pos - 1);

		new_lmpos = pos - remain;
		pac2char(new_lmpos, remain, res_str);

		len = covered + remain;
		edit_dist = alignment.hamming_distance_left(res_str, len, seq, len, sclen);

		//vafprintf(2, stderr, "lmpos: %lu\textend len: %d\n", new_lmpos, len);
		//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nLeftedit dist: %d\n", res_str, seq, edit_dist);

		if ((sclen == sclen_best and edit_dist < min_ed) or (sclen < sclen_best and edit_dist <= EDTH)) {	// less soft clip and then less hamming distance
			min_ed = edit_dist;
			sclen_best = sclen;
			lmpos_best = new_lmpos;
		}
	}

	for (int i = 0; i < overlapped_exon.size(); i++) {
		//vafprintf(2, stderr, "%u <- This exon: [%u-%u]\n", overlapped_exon[i].prev_exon_end, overlapped_exon[i].start, overlapped_exon[i].end);
		if (covered > 0 and overlapped_exon[i].end != search_pos)	// when jump to prev exon
			continue;

		exon_remain = pos - overlapped_exon[i].start;
		//vafprintf(2, stderr, "Remain on exon: %u\n", exon_remain);
		if (exon_remain >= remain) {
			//vafprintf(2, stderr, "Going for %lu - %lu\n", pos - remain, pos - 1);

			new_lmpos = pos - remain;
			pac2char(new_lmpos, remain, res_str);

			len = covered + remain;
			edit_dist = alignment.hamming_distance_left(res_str, len, seq, len, sclen);

			//vafprintf(2, stderr, "lmpos: %lu\textend len: %d\n", new_lmpos, len);
			//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nLeft edit dist: %d\n", res_str, seq, edit_dist);

			if ((sclen == sclen_best and edit_dist < min_ed) or (sclen < sclen_best and edit_dist <= EDTH)) {	// less soft clip and then less hamming distance
				min_ed = edit_dist;
				sclen_best = sclen;
				lmpos_best = new_lmpos;
			}
		}
		else {		// junction
			if (overlapped_exon[i].prev_exon_end == 0)	// should not be the first exon and partially mappable
				continue;

			pac2char(overlapped_exon[i].start, exon_remain, res_str + remain - exon_remain);
			get_seq_left(res_str, seq, overlapped_exon[i].prev_exon_end + 1, covered + exon_remain, remain - exon_remain, min_ed, sclen_best, lmpos_best);
		}
	}
}

// pos is exclusive
// [ pos-len, pos-1 ]
bool extend_left(char* seq, uint32_t& pos, int len) {
	char res_str[len+5];
	bool left_ok = false;
	uint32_t new_lmpos = pos;
	
	int min_ed = EDTH + 1;
	int sclen_best = SOFTCLIPTH;
	uint32_t lmpos_best = 0;

	res_str[len] = '\0';
	get_seq_left(res_str, seq, pos, 0, len, min_ed, sclen_best, lmpos_best);
	
	if (min_ed <= EDTH) {
		pos = lmpos_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\n", min_ed, pos);
		return true;
	}

	// intron retentaion
	new_lmpos = pos - len;
	pac2char(new_lmpos , len, res_str);	
	left_ok = alignment.hamming_match_left(res_str, len, seq, len);
	
	//vafprintf(2, stderr, "Intron Retention:\nlmpos: %lu\textend len: %d\n", pos, len);
	//vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nleft ok? %d\n", res_str, seq, left_ok);

	if (left_ok)
		pos = new_lmpos;

	return left_ok;
}

void overlap_to_epos(MatchedMate& mr) {
	if (mr.looked_up_epos or mr.exons_epos != NULL)
		return;
	mr.exons_epos = gtf_parser.get_location_overlap(mr.end_pos, false);
	mr.looked_up_epos = true;
	//fprintf(stdout, "End Seg list size: %d\n", mr.exons_epos->seg_list.size());
}

void overlap_to_spos(MatchedMate& mr) {
	if (mr.looked_up_spos or mr.exons_spos != NULL)
		return;
	mr.exons_spos = gtf_parser.get_location_overlap(mr.start_pos, false);
	mr.looked_up_spos = true;
	//fprintf(stdout, "Start Seg list size: %d\n", mr.exons_spos->seg_list.size());
}
