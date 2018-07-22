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
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedRead& mr, int dir);
int process_mates(const chain_list& forward_chain, const Record* record1, const chain_list& backward_chain, const Record* record2);
bool extend_right(char* seq, uint32_t pos, int len);
bool extend_left(char* seq, uint32_t pos, int len);

// updates next fq file to be read if need be (keep file)
FilterRead::FilterRead (char* save_fname, bool pe, char* filter_temp_name, int num_files, char* fq_file1, char* fq_file2) {
	is_pe = pe;
	cat_count = num_files;

	char* output_names[6] = { "ignore", "keep", "OEA", "orphan", "chim_bsj", "chim_fus" };
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
		sprintf(fq_file1, "%s_%s.%s_R1.fastq", save_fname, filter_temp_name, output_names[CANDID]);
		sprintf(fq_file2, "%s_%s.%s_R2.fastq", save_fname, filter_temp_name, output_names[CANDID]);
	}
	else {
		sprintf(fq_file1, "%s_%s.%s.fastq", save_fname, filter_temp_name, output_names[CANDID]);
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
	
	vafprintf(1, stderr, "%s", current_record->rname);

	MatchedRead mr;
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
	vafprintf(1, stderr, "R1/%s", current_record1->rname);

	get_best_chains(current_record1->seq, current_record1->seq_len, kmer_size, forward_best_chain_r1, fl);
	get_best_chains(current_record1->rcseq, current_record1->seq_len, kmer_size, backward_best_chain_r1, bl);

	vafprintf(1, stderr, "R1/%s", current_record1->rname);
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
	vafprintf(1, stderr, "R2/%s", current_record2->rname);

	get_best_chains(current_record2->seq, current_record2->seq_len, kmer_size, forward_best_chain_r2, fl);
	get_best_chains(current_record2->rcseq, current_record2->seq_len, kmer_size, backward_best_chain_r2, bl);

	vafprintf(1, stderr, "R2/%s", current_record2->rname);
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
	int cat = (state >= cat_count) ? CANDID : state;
	fprintf(cat_file_r1[cat], "%s%s%s%d\n%s", current_record->rname, current_record->seq, current_record->comment, state, current_record->qual);
}

// write reads PE mode
void FilterRead::write_read_category (Record* current_record1, Record* current_record2, int state) {
	state = minM(state, current_record1->state);
	int cat = (state >= cat_count) ? CANDID : state;
	fprintf(cat_file_r1[cat], "%s%s%s%d\n%s", current_record1->rname, current_record1->seq, current_record1->comment, state, current_record1->qual);
	fprintf(cat_file_r2[cat], "%s%s%s%d\n%s", current_record2->rname, current_record2->seq, current_record2->comment, state, current_record2->qual);
}



bool is_concord(const chain_t& a, int seq_len, MatchedRead& mr) {
	if (a.chain_len < 2) {
		mr.is_concord = false;
	}
	else if ((a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos) >= seq_len) {
		mr.is_concord = true;
		mr.type = CONCRD;
		
		//mr.chr = chr_name;
		//mr.chr = "1";
		mr.chr = getRefGenomeName();
		mr.start_pos = a.frags[0].rpos;
		mr.end_pos = a.frags[a.chain_len-1].rpos + a.frags[a.chain_len-1].len - 1;
		//int seg_ind = gtf_parser.search_loc(chr_name, true, chr_beg);
		//gtf_parser.get_gene_id(chr, seg_ind, mr.gene_id);
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
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedRead& mr, int dir) {
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

	int match_score = 2;
	
	// mismatch and gap should be equal in order to calculate edit distance
	int mm_pen = -1;
	int gap_pen = -1;
	int ed;

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

	//free(remain_str_beg);
	//free(remain_str_end);

	if (left_ok and right_ok) {
		mr.is_concord = true;
		mr.type = CONCRD;
		mr.chr = getRefGenomeName();
		mr.start_pos = lm_pos;
		mr.end_pos = rm_pos;
		mr.matched_len = rm_pos - lm_pos + 1;
		mr.dir = dir;
	}
	else if (left_ok or right_ok)
		mr.type = CANDID;
	else
		mr.type = ORPHAN;

	return mr.type;
}

bool are_concordant(const vector <MatchedRead>& mrs, int mrs_size, const MatchedRead& mr) {
	if (mrs_size <= 0 or mr.dir * mrs[0].dir == 1)	// empty / same orientation
		return false;

	for (int i = 0; i < mrs_size; i++) {
		vafprintf(2, stderr, "MR[%d]: %s: %d - %d\n", i, mr.chr, mr.start_pos, mr.end_pos);
		MatchedRead omr = mrs[i];
		if (mr.dir == 1 and strcmp(mr.chr, omr.chr) == 0 and mr.start_pos <= omr.start_pos and omr.end_pos <= mr.start_pos + GENETHRESH) {
			vafprintf(2, stderr, "Mapped\n");
			vafprintf(2, stderr, "MR[%d]: %s: %d - %d\n", i, omr.chr, omr.start_pos, omr.end_pos);
			return true;
		}
		if (mr.dir == -1 and strcmp(mr.chr, omr.chr) == 0 and omr.start_pos <= mr.start_pos and mr.end_pos <= omr.start_pos + GENETHRESH) {
			vafprintf(2, stderr, "Mapped\n");
			vafprintf(2, stderr, "MR[%d]: %s: %d - %d\n", i, omr.chr, omr.start_pos, omr.end_pos);
			return true;
		}
	}
	return false;
}

int process_mates(const chain_list& forward_chain, const Record* record1, const chain_list& backward_chain, const Record* record2) {
	int fc_size = forward_chain.best_chain_count;
	int bc_size = backward_chain.best_chain_count;

	int max_len = (fc_size >= bc_size) ? fc_size : bc_size;

	vector <MatchedRead> forward_mrl(fc_size);
	vector <MatchedRead> backward_mrl(bc_size);

	int ex_ret;
	int fmrl_count = 0;
	int bmrl_count = 0;
	int min_ret1 = ORPHAN;
	int min_ret2 = ORPHAN;
	for (int i = 0; i < max_len; i++) {
		if (i < fc_size) {
			ex_ret = extend_chain(forward_chain.chains[i], record1->seq, record1->seq_len, forward_mrl[fmrl_count], 1); 
			if (ex_ret == CONCRD) {
				if (are_concordant(backward_mrl, bmrl_count, forward_mrl[fmrl_count])) {
					return CONCRD;
				}
				fmrl_count++;
			}
			if (ex_ret < min_ret1)	min_ret1 = ex_ret;
		}

		if (i < bc_size) {
			ex_ret = extend_chain(backward_chain.chains[i], record2->rcseq, record2->seq_len, backward_mrl[bmrl_count], -1); 
			if (ex_ret == CONCRD) {
				if (are_concordant(forward_mrl, fmrl_count, backward_mrl[bmrl_count])) {
					return CONCRD;
				}
				bmrl_count++;
			}
			if (ex_ret < min_ret2)	min_ret2 = ex_ret;
		}
	}

	return 	((min_ret1 == ORPHAN) and (min_ret2 == ORPHAN)) ? ORPHAN 
			: (((min_ret1 == ORPHAN) and (min_ret2 == CONCRD)) or ((min_ret1 == CONCRD) and (min_ret2 == ORPHAN))) ? OEANCH 
			: CANDID;
}

// pos is exclusive
// [ pos+1, pos+len ]
bool extend_right(char* seq, uint32_t pos, int len) {
	char res_str[len+5];
	bool right_ok = false;
	
	string chr = getRefGenomeName();
	if (chr == "")
		return false; 

	vafprintf(2, stderr, "Going for %lu - %lu\n", pos+1, pos+len);
	
	vector <UniqSeg> overlapped_exon;
	gtf_parser.get_location_overlap(pos, overlapped_exon);

	UniqSeg seg;
	uint32_t exon_remain;
	for (int i = 0; i < overlapped_exon.size(); i++) {
		res_str[0] = 0;
		exon_remain = overlapped_exon[i].end - pos;
		if (exon_remain >= len) {	
			pac2char(pos + 1, len, res_str);
		}
		else {
			//er = overlapped_exon[i];
			//int covered = 0;
			if (overlapped_exon[i].next_exon_beg == 0)	// should not be the last exon and partially mappable
				continue;
			else {
				pac2char(pos + 1, exon_remain, res_str);
				int remain_len = len - exon_remain;
				//er = ...
				pac2char(overlapped_exon[i].next_exon_beg, remain_len, res_str + exon_remain);
			}
		}
		right_ok = alignment.hamming_match_right(res_str, len, seq, len);
		vafprintf(2, stderr, "rmpos: %lu\textend len: %d\n", pos, len);
		vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nright ok? %d\n", res_str, seq, right_ok);
		if (right_ok)
			return true;
	}

	pac2char(pos + 1, len, res_str);	// intron retentaion
	right_ok = alignment.hamming_match_right(res_str, len, seq, len);
	vafprintf(2, stderr, "Intron Retention:\nrmpos: %lu\textend len: %d\n", pos, len);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nright ok? %d\n", res_str, seq, right_ok);

	return right_ok;
}

// pos is exclusive
// [ pos-len, pos-1 ]
bool extend_left(char* seq, uint32_t pos, int len) {
	char res_str[len+5];
	bool left_ok = false;
	
	string chr = getRefGenomeName();
	if (chr == "")
		return false; 

	vafprintf(2, stderr, "Going for %lu - %lu\n", pos-len, pos-1);
	
	vector <UniqSeg> overlapped_exon;
	gtf_parser.get_location_overlap(pos, overlapped_exon);

	UniqSeg seg;
	uint32_t exon_remain;
	for (int i = 0; i < overlapped_exon.size(); i++) {
		res_str[0] = 0;
		exon_remain = pos - overlapped_exon[i].start;
		if (exon_remain >= len) {	
			pac2char(pos - len, len, res_str);
		}
		else {
			if (overlapped_exon[i].prev_exon_end == 0)	// should not be the first exon and partially mappable
				continue;
			else {
				int remain_len = len - exon_remain;
				pac2char(overlapped_exon[i].prev_exon_end - remain_len + 1, remain_len, res_str);
				pac2char(overlapped_exon[i].start, exon_remain, res_str + remain_len);
			}
		}
		left_ok = alignment.hamming_match_left(res_str, len, seq, len);
		vafprintf(2, stderr, "lmpos: %lu\textend len: %d\n", pos, len);
		vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nleft ok? %d\n", res_str, seq, left_ok);
		if (left_ok)
			return true;
	}

	pac2char(pos - len , len, res_str);	// intron retentaion
	left_ok = alignment.hamming_match_left(res_str, len, seq, len);
	vafprintf(2, stderr, "Intron Retention:\nlmpos: %lu\textend len: %d\n", pos, len);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nleft ok? %d\n", res_str, seq, left_ok);

	return left_ok;
}

