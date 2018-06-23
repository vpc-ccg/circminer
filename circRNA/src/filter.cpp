#include <vector>
#include <cstdlib>
#include <cmath>

#include "filter.h"
#include "align.h"
#include "common.h"
#include "gene_annotation.h"
#include "fragment_list.h"

#define BESTCHAINLIM 30
#define EDTH 3
#define SOFTCLIPTH 7

FilterRead::FilterRead (char* save_fname, bool pe) {
	is_pe = pe;

	char ignore_file1[1000];
	char keep_file1[1000];
	char chimeric_bsj_file1[1000];
	char chimeric_fusion_file1[1000];
	char p_unmap_file1[1000];
	char unmap_file1[1000];

	strcpy (ignore_file1, save_fname);
	strcpy (keep_file1 , save_fname);
	strcpy (chimeric_bsj_file1 , save_fname);
	strcpy (chimeric_fusion_file1 , save_fname);
	strcpy (p_unmap_file1 , save_fname);
	strcpy (unmap_file1 , save_fname);

	strcat (ignore_file1, ".ignore_R1.fastq");
	strcat (keep_file1, ".keep_R1.fastq");
	strcat (chimeric_bsj_file1, ".chimeric.bsj_R1.fastq");
	strcat (chimeric_fusion_file1, ".chimeric.fusion_R1.fastq");
	strcat (p_unmap_file1, ".OEA_R1.fastq");
	strcat (unmap_file1, ".orphan_R1.fastq");

	ignore_r1 = fopen(ignore_file1, "w");
	keep_r1   = fopen(keep_file1, "w");
	chimeric_bsj_r1   = fopen(chimeric_bsj_file1, "w");
	chimeric_fusion_r1   = fopen(chimeric_fusion_file1, "w");
	partly_unmappable_r1   = fopen(p_unmap_file1, "w");
	unmappable_r1   = fopen(unmap_file1, "w");

	if (is_pe) {
		char ignore_file2[1000];
		char keep_file2[1000];
		char chimeric_bsj_file2[1000];
		char chimeric_fusion_file2[1000];
		char p_unmap_file2[1000];
		char unmap_file2[1000];
		strcpy (ignore_file2, save_fname);
		strcpy (keep_file2 , save_fname);
		strcpy (chimeric_bsj_file2 , save_fname);
		strcpy (chimeric_fusion_file2 , save_fname);
		strcpy (p_unmap_file2 , save_fname);
		strcpy (unmap_file2 , save_fname);

		strcat (ignore_file2, ".ignore_R2.fastq");
		strcat (keep_file2, ".keep_R2.fastq");
		strcat (chimeric_bsj_file2, ".chimeric.bsj_R2.fastq");
		strcat (chimeric_fusion_file2, ".chimeric.fusion_R2.fastq");
		strcat (p_unmap_file2, ".OEA_R2.fastq");
		strcat (unmap_file2, ".orphan_R2.fastq");

		ignore_r2 = fopen(ignore_file2, "w");
		keep_r2   = fopen(keep_file2, "w");
		chimeric_bsj_r2   = fopen(chimeric_bsj_file2, "w");
		chimeric_fusion_r2   = fopen(chimeric_fusion_file2, "w");
		partly_unmappable_r2   = fopen(p_unmap_file2, "w");
		unmappable_r2   = fopen(unmap_file2, "w");
	}
}

FilterRead::~FilterRead (void) {
	if (ignore_r1 != NULL)
		fclose(ignore_r1);
	if (keep_r1 != NULL)
		fclose(keep_r1);
	if (chimeric_bsj_r1 != NULL)
		fclose(chimeric_bsj_r1);
	if (chimeric_fusion_r1 != NULL)
		fclose(chimeric_fusion_r1);
	if (partly_unmappable_r1 != NULL)
		fclose(partly_unmappable_r1);
	if (unmappable_r1 != NULL)
		fclose(unmappable_r1);

	if (ignore_r2 != NULL)
		fclose(ignore_r2);
	if (keep_r2 != NULL)
		fclose(keep_r2);
	if (chimeric_bsj_r2 != NULL)
		fclose(chimeric_bsj_r2);
	if (chimeric_fusion_r2 != NULL)
		fclose(chimeric_fusion_r2);
	if (partly_unmappable_r2 != NULL)
		fclose(partly_unmappable_r2);
	if (unmappable_r2 != NULL)
		fclose(unmappable_r2);
}

int FilterRead::process_read (Record* current_record) {
	return find_expanded_positions(current_record->seq, current_record->rcseq, current_record->seq_len);
}

int FilterRead::process_read (Record* current_record1, Record* current_record2, int kmer_size) {
	return check_concordant_mates_expand(current_record1, current_record2, kmer_size);
}


bool is_concord(const chain_t& a, int seq_len, int kmer) {
	if (a.chain_len < 2)
		return false;
	return (a.frags[a.chain_len-1].qpos - a.frags[0].qpos + kmer) >= seq_len;
}

bool is_concord(const chain_t& a, int seq_len, MatchedRead& mr) {
	if (a.chain_len < 2) {
		mr.is_concord = false;
	}
	else if ((a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos) >= seq_len) {
		mr.is_concord = true;
		mr.type = 0;
		
		char* chr_name;
		int32_t chr_len;
		uint32_t chr_beg;
		uint32_t chr_end;

		bwt_get_intv_info((bwtint_t) a.frags[0].rpos, (bwtint_t) (a.frags[a.chain_len-1].rpos + a.frags[a.chain_len-1].len - 1), &chr_name, &chr_len, &chr_beg, &chr_end);
		mr.chr = chr_name;
		mr.start_pos = chr_beg;
		mr.end_pos = chr_end;
		//int seg_ind = gtf_parser.search_loc(chr_name, true, chr_beg);
		//gtf_parser.get_gene_id(chr, seg_ind, mr.gene_id);
		mr.matched_len = a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len - a.frags[0].qpos;
		mr.dir = (a.frags[0].qpos >= 0) ? 1 : -1;
	}
	else {
		mr.is_concord = false;
	}
	return mr.is_concord;
}

void get_best_chain(char* read_seq, int seq_len, int kmer_size, chain_t& forward_best_chain, chain_t& backward_best_chain) {
	int kmer_count = ceil(seq_len / kmer_size);
	int forward_fragment_count, backward_fragment_count;
	vector <fragment_t> forward_fragments(FRAGLIM * kmer_count);
	vector <fragment_t> backward_fragments(FRAGLIM * kmer_count);

	chop_read_match(read_seq, seq_len, kmer_size, 0, true, forward_fragments, forward_fragment_count, backward_fragments, backward_fragment_count);
	//vafprintf(2, stderr, "Forward#: %d\n", forward_fragment_count);
	//for (int i = 0; i < forward_fragment_count; i++)
	//	vafprintf(2, stderr, "rpos: %lu, qpos: %d, len: %lu", forward_fragments[i].rpos, forward_fragments[i].qpos, forward_fragments[i].len);
	//vafprintf(2, stderr, "\n");
	//vafprintf(2, stderr, "Backward#: %d\n", backward_fragment_count);
	//for (int i = 0; i < backward_fragment_count; i++)
	//	vafprintf(2, stderr, "rpos: %lu, qpos: %d, len: %lu", backward_fragments[i].rpos, backward_fragments[i].qpos, backward_fragments[i].len);
	//vafprintf(2, stderr, "\n");

	chain_seeds_n2(forward_fragments, forward_fragment_count, forward_best_chain);
	chain_seeds_n2(backward_fragments, backward_fragment_count, backward_best_chain);
}

void get_best_chains(char* read_seq, int seq_len, int kmer_size, vector <chain_t>& forward_best_chain, vector<chain_t>& backward_best_chain) {
	int kmer_count = ceil(seq_len / kmer_size);
	int forward_fragment_count, backward_fragment_count;

	FragmentList forward_frag_ll;
	FragmentList backward_frag_ll;
	split_match_ll(read_seq, seq_len, kmer_size, forward_frag_ll, backward_frag_ll);
	
	//vector <fragment_t> forward_fragments(FRAGLIM * kmer_count);
	//vector <fragment_t> backward_fragments(FRAGLIM * kmer_count);
	//split_match(read_seq, seq_len, kmer_size, forward_fragments, forward_fragment_count, backward_fragments, backward_fragment_count);
	
	//chop_read_match(read_seq, seq_len, kmer_size, 0, true, forward_fragments, forward_fragment_count, backward_fragments, backward_fragment_count);
	
	//vafprintf(2, stderr, "Forward#: %d\n", forward_fragment_count);
	//for (int i = 0; i < forward_fragment_count; i++)
	//	vafprintf(2, stderr, "rpos: %lu, qpos: %d, len: %lu", forward_fragments[i].rpos, forward_fragments[i].qpos, forward_fragments[i].len);
	//vafprintf(2, stderr, "\n");
	//vafprintf(2, stderr, "Backward#: %d\n", backward_fragment_count);
	//for (int i = 0; i < backward_fragment_count; i++)
	//	vafprintf(2, stderr, "rpos: %lu, qpos: %d, len: %lu", backward_fragments[i].rpos, backward_fragments[i].qpos, backward_fragments[i].len);
	//vafprintf(2, stderr, "\n");

	chain_seeds_n2_kbest(forward_frag_ll, forward_best_chain);
	chain_seeds_n2_kbest(backward_frag_ll, backward_best_chain);

	//chain_seeds_n2_kbest(forward_fragments, forward_fragment_count, forward_best_chain);
	//chain_seeds_n2_kbest(backward_fragments, backward_fragment_count, backward_best_chain);
}

// return:
// 0 is successfully extended
// 3 if not successful from at least one direction
// 5 if orphan (chain length = 0)
int extend_chain(const chain_t& ch, char* seq, int seq_len, MatchedRead& mr) {
	mr.is_concord = false;
	if (ch.chain_len <= 0) {
		mr.type = 5;
		return 5;
	}

	if (is_concord(ch, seq_len, mr))
		return 0;

	char* chr_name;
	int32_t chr_len;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	bool left_ok = true;
	bool right_ok = true;

	int match_score = 2;
	// mismatch and gap should be equal in order to calculate edit distance
	int mm_pen = -1;
	int gap_pen = -1;
	int ed;

	uint32_t lm_pos = ch.frags[0].rpos;
	int remain_beg = (ch.frags[0].qpos >= 0) ? (ch.frags[0].qpos) : (seq_len - 1 + ch.frags[0].qpos);	// based on forward or reverse strand mapping
	
	char* remain_str_beg = (char*) malloc(remain_beg+5);
	if (remain_beg > 0) {
		get_reference_chunk_left(lm_pos, remain_beg, remain_str_beg);
		//get_reference_chunk(lm_pos, -1 * remain_beg, remain_str_beg);
		//vafprintf(1, stderr, "On Ref : %s\n", remain_str_beg);
		//vafprintf(1, stderr, "On Read: %s\n", seq);
		
		//ed = alignment(remain_str_beg, remain_beg, seq, remain_beg, 1, 1);
		
		loc_align align = local_alignment_reverse(remain_str_beg, remain_beg, seq, remain_beg, match_score, mm_pen, gap_pen);
		lm_pos -= align.r_matched;
		ed = alignment(remain_str_beg + remain_beg - align.r_matched, align.r_matched, seq + remain_beg - align.q_matched, align.q_matched, 1, 1);
		//vafprintf(1, stderr, "Score: %d, Edit dist: %d\n", align.score, ed);
		if (remain_beg - align.q_matched > SOFTCLIPTH or ed > EDTH)
		
		//ed = hamming_distance_left(remain_str_beg, remain_beg, seq, remain_beg, SOFTCLIPTH);
		//vafprintf(1, stderr, "Hamming dist: %d\n", ed);
		//if (ed > EDTH)
			left_ok = false;
	}

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = (ch.frags[0].qpos >= 0) ? (seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len)) : (1 - (ch.frags[ch.chain_len-1].len + ch.frags[ch.chain_len-1].qpos));

	char* remain_str_end = (char*) malloc(remain_end+5);
	if (remain_end > 0) {
		get_reference_chunk_right(rm_pos, remain_end, remain_str_end);
		//get_reference_chunk(rm_pos, remain_end, remain_str_end);
		//vafprintf(1, stderr, "On Ref : %s\n", remain_str_end);
		//vafprintf(1, stderr, "On Read: %s\n", seq + seq_len - remain_end);
		//ed = alignment(remain_str_end, remain_end, seq + seq_len - remain_end, remain_end, 1, 1);
		
		loc_align align = local_alignment(remain_str_end, remain_end, seq + seq_len - remain_end, remain_end, match_score, mm_pen, gap_pen);
		rm_pos += align.r_matched;
		ed = alignment(remain_str_end, align.r_matched, seq + seq_len - remain_end, align.q_matched, 1, 1);
		//vafprintf(1, stderr, "Score: %d, Edit dist: %d\n", align.score, ed);
		//vafprintf(1, stderr, "Q matched: %d, R matched: %d\n", align.q_matched, align.r_matched);
		if (remain_end - align.q_matched > SOFTCLIPTH or ed > EDTH)
		
		//ed = hamming_distance_right(remain_str_end, remain_end, seq + seq_len - remain_end, remain_end, SOFTCLIPTH);
		//vafprintf(1, stderr, "Hamming dist: %d\n", ed);
		//if (ed > EDTH)
			right_ok = false;
	}

	free(remain_str_beg);
	free(remain_str_end);

	if (left_ok and right_ok) {
		mr.is_concord = true;
		mr.type = 0;
		bwt_get_intv_info((bwtint_t) lm_pos, (bwtint_t) rm_pos, &chr_name, &chr_len, &chr_beg, &chr_end);
		mr.chr = chr_name;
		mr.start_pos = chr_beg;
		mr.end_pos = chr_end;
		mr.matched_len = rm_pos - lm_pos + 1;
		mr.dir = (ch.frags[0].qpos >= 0) ? 1 : -1;
	}
	else if (left_ok or right_ok)
		mr.type = 3;
	else
		mr.type = 5;

	return mr.type;
}

bool are_concordant(const vector <MatchedRead>& mrs, int mrs_size, const MatchedRead& mr) {
	if (mrs_size <= 0 or mr.dir * mrs[0].dir == 1)	// empty / same orientation
		return false;

	vafprintf(2, stderr, "MR: %s: %d - %d\n", mr.chr, mr.start_pos, mr.end_pos);
	vafprintf(2, stderr, "Comp with %d candidates\n", mrs_size);

	for (int i = 0; i < mrs_size; i++) {
		vafprintf(2, stderr, "MR[%d]: %s: %d - %d\n", i, mr.chr, mr.start_pos, mr.end_pos);
		MatchedRead omr = mrs[i];
		if (mr.dir == 1 and strcmp(mr.chr, omr.chr) == 0 and mr.start_pos <= omr.start_pos and omr.end_pos <= mr.start_pos + GENETHRESH)
			return true;
		if (mr.dir == -1 and strcmp(mr.chr, omr.chr) == 0 and omr.start_pos <= mr.start_pos and mr.end_pos <= omr.start_pos + GENETHRESH)
			return true;
	}
	return false;
}

int process_mates(const vector <chain_t>& forward_chain, const Record* record1, const vector <chain_t>& backward_chain, const Record* record2) {
	int fc_size = forward_chain.size();
	int bc_size = backward_chain.size();
	//if (fc_size == 0 or bc_size == 0)
	//	return 5;

	int max_len = (fc_size >= bc_size) ? fc_size : bc_size;

	vector <MatchedRead> forward_mrl(fc_size);
	vector <MatchedRead> backward_mrl(bc_size);

	int ex_ret;
	int fmrl_count = 0;
	int bmrl_count = 0;
	int min_ret1 = 5;
	int min_ret2 = 5;
	for (int i = 0; i < max_len; i++) {
		if (i < fc_size) {
			ex_ret = extend_chain(forward_chain[i], record1->seq, record1->seq_len, forward_mrl[fmrl_count]); 
			if (ex_ret == 0) {
				if (are_concordant(backward_mrl, bmrl_count, forward_mrl[fmrl_count])) {
					//vafprintf(1, stderr, "----Concordant\n");
					return 0;
				}
				fmrl_count++;
			}
			if (ex_ret < min_ret1)	min_ret1 = ex_ret;
		}

		if (i < bc_size) {
			ex_ret = extend_chain(backward_chain[i], record2->rcseq, record2->seq_len, backward_mrl[bmrl_count]); 
			if (ex_ret == 0) {
				if (are_concordant(forward_mrl, fmrl_count, backward_mrl[bmrl_count])) {
					//vafprintf(1, stderr, "----Concordant\n");
					return 0;
				}
				bmrl_count++;
			}
			if (ex_ret < min_ret2)	min_ret2 = ex_ret;
		}
	}

	//vafprintf(2, stderr, "R1 mappablilty: %d\nR2 mappablity: %d\n", min_ret1, min_ret2);
	
	return 	((min_ret1 == 5) and (min_ret2 == 5)) ? 5 
			: (((min_ret1 == 5) and (min_ret2 == 0)) or ((min_ret1 == 0) and (min_ret2 == 5))) ? 4 
			: 3;
}

int FilterRead::process_read_chain (Record* current_record1, Record* current_record2, int kmer_size) {
	int max_frag_count = current_record1->seq_len / kmer_size + 1;
	bool no_map_r1 = true;
	bool no_map_r2 = true;

	char* chr_name;
	int32_t chr_len;
	uint32_t chr_beg;
	uint32_t chr_end;

	// R1
	vafprintf(1, stderr, "R1/%s", current_record1->rname);
	vector <chain_t> forward_best_chain_r1(BESTCHAINLIM);
	vector <chain_t> backward_best_chain_r1(BESTCHAINLIM);
	for (int i = 0; i < BESTCHAINLIM; i++) {
		forward_best_chain_r1[i].frags.resize(max_frag_count + 1);
		backward_best_chain_r1[i].frags.resize(max_frag_count + 1);
	}

	get_best_chains(current_record1->seq, current_record1->seq_len, kmer_size, forward_best_chain_r1, backward_best_chain_r1);

	//vafprintf(1, stderr, "R1 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r1[0].score, (unsigned long)forward_best_chain_r1.size());
	//for (int j = 0; j < forward_best_chain_r1.size(); j++)
	//	for (int i = 0; i < forward_best_chain_r1[j].chain_len; i++) {
	//		vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, forward_best_chain_r1[j].frags[i].rpos, forward_best_chain_r1[j].frags[i].qpos, forward_best_chain_r1[j].frags[i].len);
	//		bwt_get_intv_info((bwtint_t) forward_best_chain_r1[j].frags[i].rpos, (bwtint_t) (forward_best_chain_r1[j].frags[i].rpos + forward_best_chain_r1[j].frags[i].len - 1), &chr_name, &chr_len, &chr_beg, &chr_end);
	//		vafprintf(1, stderr, "Chr %s: %lu-%lu\n", chr_name, chr_beg, chr_end);
	//	}

	//vafprintf(1, stderr, "R1 Backward score:%.4f,\t len: %lu\n", backward_best_chain_r1[0].score, (unsigned long)backward_best_chain_r1.size());
	//for (int j = 0; j < backward_best_chain_r1.size(); j++)
	//	for (int i = 0; i < backward_best_chain_r1[j].chain_len; i++) {
	//		vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, backward_best_chain_r1[j].frags[i].rpos, backward_best_chain_r1[j].frags[i].qpos, backward_best_chain_r1[j].frags[i].len);
	//		bwt_get_intv_info((bwtint_t) backward_best_chain_r1[j].frags[i].rpos, (bwtint_t) (backward_best_chain_r1[j].frags[i].rpos + backward_best_chain_r1[j].frags[i].len - 1), &chr_name, &chr_len, &chr_beg, &chr_end);
	//		vafprintf(1, stderr, "Chr %s: %lu-%lu\n", chr_name, chr_beg, chr_end);
	//	}

	//vafprintf(1, stderr, "Forward R1\n");
	//vector <MatchedRead> forward_mrl_r1(forward_best_chain_r1.size());
	//int fmrl_r1_count = 0;
	//int ex_ret;
	//for (int j = 0; j < forward_best_chain_r1.size(); j++) {
	//	ex_ret = extend_chain(forward_best_chain_r1[j], current_record1->seq, current_record1->seq_len, forward_mrl_r1[fmrl_r1_count]);
	//	vafprintf(1, stderr, "Extend ret: %d\n", ex_ret);
	//	if (ex_ret == 0)	// only keep concordant ones
	//		fmrl_r1_count++;
	//	if (ex_ret != 5)
	//		no_map_r1 = false;
	//}

	//vafprintf(1, stderr, "Backward R1\n");
	//vector <MatchedRead> backward_mrl_r1(backward_best_chain_r1.size());
	//int bmrl_r1_count = 0;
	//for (int j = 0; j < backward_best_chain_r1.size(); j++) {
	//	ex_ret = extend_chain(backward_best_chain_r1[j], current_record1->rcseq, current_record1->seq_len, backward_mrl_r1[bmrl_r1_count]); 
	//	if (ex_ret == 0)
	//		bmrl_r1_count++;
	//	if (ex_ret != 5)
	//		no_map_r1 = false;
	//}

	// R2
	vafprintf(1, stderr, "R2/%s", current_record2->rname);
	vector <chain_t> forward_best_chain_r2(BESTCHAINLIM);
	vector <chain_t> backward_best_chain_r2(BESTCHAINLIM);
	for (int i = 0; i < BESTCHAINLIM; i++) {
		forward_best_chain_r2[i].frags.resize(max_frag_count + 1);
		backward_best_chain_r2[i].frags.resize(max_frag_count + 1);
	}

	get_best_chains(current_record2->seq, current_record2->seq_len, kmer_size, forward_best_chain_r2, backward_best_chain_r2);

	//vafprintf(1, stderr, "R2 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r2[0].score, (unsigned long)forward_best_chain_r2.size());
	//for (int j = 0; j < forward_best_chain_r2.size(); j++)
	//	for (int i = 0; i < forward_best_chain_r2[j].chain_len; i++) {
	//		vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, forward_best_chain_r2[j].frags[i].rpos, forward_best_chain_r2[j].frags[i].qpos, forward_best_chain_r2[j].frags[i].len);
	//		bwt_get_intv_info((bwtint_t) forward_best_chain_r2[j].frags[i].rpos, (bwtint_t) (forward_best_chain_r2[j].frags[i].rpos + forward_best_chain_r2[j].frags[i].len - 1), &chr_name, &chr_len, &chr_beg, &chr_end);
	//		vafprintf(1, stderr, "Chr %s: %lu-%lu\n", chr_name, chr_beg, chr_end);
	//	}

	//vafprintf(1, stderr, "R2 Backward score:%.4f,\t len: %lu\n", backward_best_chain_r2[0].score, (unsigned long)backward_best_chain_r2.size());
	//for (int j = 0; j < backward_best_chain_r2.size(); j++)
	//	for (int i = 0; i < backward_best_chain_r2[j].chain_len; i++) {
	//		vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, backward_best_chain_r2[j].frags[i].rpos, backward_best_chain_r2[j].frags[i].qpos, backward_best_chain_r2[j].frags[i].len);
	//		bwt_get_intv_info((bwtint_t) backward_best_chain_r2[j].frags[i].rpos, (bwtint_t) (backward_best_chain_r2[j].frags[i].rpos + backward_best_chain_r2[j].frags[i].len - 1), &chr_name, &chr_len, &chr_beg, &chr_end);
	//		vafprintf(1, stderr, "Chr %s: %lu-%lu\n", chr_name, chr_beg, chr_end);
	//	}

	//vafprintf(1, stderr, "Forward R2\n");
	//vector <MatchedRead> forward_mrl_r2(forward_best_chain_r2.size());
	//int fmrl_r2_count = 0;
	//for (int j = 0; j < forward_best_chain_r2.size(); j++) {
	//	ex_ret = extend_chain(forward_best_chain_r2[j], current_record2->seq, current_record2->seq_len, forward_mrl_r2[fmrl_r2_count]);
	//	if (ex_ret == 0)
	//		fmrl_r2_count++;
	//	if (ex_ret != 5)
	//		no_map_r2 = false;
	//}

	//vafprintf(1, stderr, "Backward R2\n");
	//vector <MatchedRead> backward_mrl_r2(backward_best_chain_r2.size());
	//int bmrl_r2_count = 0;
	//for (int j = 0; j < backward_best_chain_r2.size(); j++) {
	//	ex_ret = extend_chain(backward_best_chain_r2[j], current_record2->rcseq, current_record2->seq_len, backward_mrl_r2[bmrl_r2_count]); 
	//	if (ex_ret == 0)
	//		bmrl_r2_count++;
	//	if (ex_ret != 5)
	//		no_map_r2 = false;
	//}

	// Orphan / OEA
	//if (no_map_r1 and no_map_r2)
	//	return 5;
	//if (no_map_r1 or no_map_r2)
	//	return 4;

	if (forward_best_chain_r1.size() + backward_best_chain_r1.size() + forward_best_chain_r2.size() + backward_best_chain_r2.size() <= 0)
		return 5;
	if ((forward_best_chain_r1.size() + backward_best_chain_r1.size() <= 0) or (forward_best_chain_r2.size() + backward_best_chain_r2.size() <= 0))
		return 4;

	// checking read concordancy
	// any concordancy evidence/explanation ?
	float fc_score_r1 = forward_best_chain_r1[0].score;
	float bc_score_r1 = backward_best_chain_r1[0].score;
	float fc_score_r2 = forward_best_chain_r2[0].score;
	float bc_score_r2 = backward_best_chain_r2[0].score;
	vafprintf(2, stderr, "Scores: fc1=%f, bc1=%f, fc2=%f, bc2=%f\n", fc_score_r1, bc_score_r1, fc_score_r2, bc_score_r2);
	int attempt1, attempt2;
	if (fc_score_r1 + bc_score_r2 >= fc_score_r2 + bc_score_r1) {
		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt1 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2);
		if (attempt1 == 0)
			return 0;

		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt2 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1);
		if (attempt2 == 0)
			return 0;
		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
	else {
		vafprintf(1, stderr, "Backward R1 / Forward R2\n");
		attempt1 = process_mates(forward_best_chain_r2, current_record2, backward_best_chain_r1, current_record1);
		if (attempt1 == 0)
			return 0;

		vafprintf(1, stderr, "Forward R1 / Backward R2\n");
		attempt2 = process_mates(forward_best_chain_r1, current_record1, backward_best_chain_r2, current_record2);
		if (attempt2 == 0)
			return 0;
		return (attempt1 < attempt2) ? attempt1 : attempt2;
	}
	
	//
	//MatchedRead mr1, mr2;
	//for (int i = 0; i < fmrl_r1_count; i++) {
	//	for (int j = 0; j < bmrl_r2_count; j++) {
	//		mr1 = forward_mrl_r1[i];
	//		mr2 = backward_mrl_r2[j];
	//		if (strcmp(mr1.chr, mr2.chr) == 0 and mr1.start_pos <= mr2.start_pos and mr2.end_pos <= mr1.start_pos + GENETHRESH)
	//			return 0;
	//	}
	//}
	//for (int i = 0; i < fmrl_r2_count; i++) {
	//	for (int j = 0; j < bmrl_r1_count; j++) {
	//		mr1 = forward_mrl_r2[i];
	//		mr2 = backward_mrl_r1[j];
	//		if (strcmp(mr1.chr, mr2.chr) == 0 and mr1.start_pos <= mr2.start_pos and mr2.end_pos <= mr1.start_pos + GENETHRESH)
	//			return 0;
	//	}
	//}

	//return 3;
}

// write reads SE mode
void FilterRead::write_read (Record* current_record, int is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
	}
	else {
		fprintf(keep_r1, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
	}
}

// write reads PE mode
void FilterRead::write_read (Record* current_record1, Record* current_record2, int is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
}

// write reads PE mode --detailed
void FilterRead::write_read2 (Record* current_record1, Record* current_record2, int is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (is_chimeric == 1) {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else {
		fprintf(chimeric_bsj_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(chimeric_bsj_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
}

// write reads PE mode --detailed / 6 output files
void FilterRead::write_read3 (Record* current_record1, Record* current_record2, int state) {
	if (state == 0 or state == 1) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 2){
		fprintf(chimeric_bsj_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(chimeric_bsj_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 3) {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 4) {
		fprintf(partly_unmappable_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(partly_unmappable_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 5) {
		fprintf(unmappable_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(unmappable_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 6){
		fprintf(chimeric_fusion_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(chimeric_fusion_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
}

