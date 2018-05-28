#include <vector>
#include <cstdlib>

#include "filter.h"
#include "chain.h"
#include "align.h"

#define EDTH 3

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

void get_best_chain(char* read_seq, int seq_len, int kmer_size, chain_t& forward_best_chain, chain_t& backward_best_chain) {
	int forward_fragment_count, backward_fragment_count;
	vector <fragment_t> forward_fragments(FRAGLIM);
	vector <fragment_t> backward_fragments(FRAGLIM);

	chop_read_match(read_seq, seq_len, kmer_size, forward_fragments, forward_fragment_count, backward_fragments, backward_fragment_count);

	chain_seeds_n2(forward_fragments, forward_fragment_count, forward_best_chain);
	chain_seeds_n2(backward_fragments, backward_fragment_count, backward_best_chain);
}

// return:
// 0 is successfully extended
// 3 if not successful from at least one direction
// 5 if orphan (chain length = 0)
int extend_chain(const chain_t& ch, char* seq, int seq_len) {
	if (ch.chain_len <= 0)
		return 5;

	bool left_ok = true;
	bool right_ok = true;

	uint32_t lm_pos = ch.frags[0].rpos;
	int remain_beg = (ch.frags[0].qpos >= 0) ? (ch.frags[0].qpos) : (seq_len - 1 + ch.frags[0].qpos);	// based on forward or reverse strand mapping
	
	char* remain_str_beg = (char*) malloc(remain_beg+5);
	if (remain_beg > 0) {
		get_reference_chunk(lm_pos, -1 * remain_beg, remain_str_beg);
		//fprintf(stderr, "On Ref : %s\n", remain_str_beg);
		//fprintf(stderr, "On Read: %s\n", seq);
		int ed = alignment(remain_str_beg, remain_beg, seq, remain_beg, 1, 1);
		//fprintf(stderr, "Edit dist: %d\n", ed);
		if (ed > EDTH)
			left_ok = false;
	}

	uint32_t rm_pos = ch.frags[ch.chain_len-1].rpos + ch.frags[ch.chain_len-1].len - 1;
	int remain_end = (ch.frags[0].qpos >= 0) ? (seq_len - (ch.frags[ch.chain_len-1].qpos + ch.frags[ch.chain_len-1].len)) : (1 - (ch.frags[ch.chain_len-1].len + ch.frags[ch.chain_len-1].qpos));

	char* remain_str_end = (char*) malloc(remain_end+5);
	if (remain_end > 0) {
		get_reference_chunk(rm_pos, remain_end, remain_str_end);
		//fprintf(stderr, "On Ref : %s\n", remain_str_end);
		//fprintf(stderr, "On Read: %s\n", seq + seq_len - remain_end);
		int ed = alignment(remain_str_end, remain_end, seq + seq_len - remain_end, remain_end, 1, 1);
		//fprintf(stderr, "Edit dist: %d\n", ed);
		if (ed > EDTH)
			right_ok = false;
	}

	free(remain_str_beg);
	free(remain_str_end);
	return (left_ok and right_ok) ? 0 : 3;
}

int process_mates(const chain_t& r1_chain, const chain_t& r2_chain, int seq_len, int kmer_size) {
}

int FilterRead::process_read_chain (Record* current_record1, Record* current_record2, int kmer_size) {
	int max_frag_count = current_record1->seq_len / kmer_size + 1;

	// R1
	//fprintf(stderr, "R1/%s", current_record1->rname);
	chain_t forward_best_chain_r1;
	chain_t backward_best_chain_r1;
	forward_best_chain_r1.frags.resize(max_frag_count + 1);
	backward_best_chain_r1.frags.resize(max_frag_count + 1);

	get_best_chain(current_record1->seq, current_record1->seq_len, kmer_size, forward_best_chain_r1, backward_best_chain_r1);

	//fprintf(stderr, "R1 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r1.score, (unsigned long)forward_best_chain_r1.chain_len);
	//for (int i = 0; i < forward_best_chain_r1.chain_len; i++)
	//	fprintf(stderr, "frag[%d]: %lu\t%d\t%lu\n", i, forward_best_chain_r1.frags[i].rpos, forward_best_chain_r1.frags[i].qpos, forward_best_chain_r1.frags[i].len);
	//fprintf(stderr, "R1 Backward score:%.4f,\t len: %lu\n", backward_best_chain_r1.score, (unsigned long)backward_best_chain_r1.chain_len);
	//for (int i = 0; i < backward_best_chain_r1.chain_len; i++)
	//	fprintf(stderr, "frag[%d]: %lu\t%d\t%lu\n", i, backward_best_chain_r1.frags[i].rpos, backward_best_chain_r1.frags[i].qpos, backward_best_chain_r1.frags[i].len);

	// R2
	//fprintf(stderr, "R2/%s", current_record2->rname);
	chain_t forward_best_chain_r2;
	chain_t backward_best_chain_r2;
	forward_best_chain_r2.frags.resize(max_frag_count + 1);
	backward_best_chain_r2.frags.resize(max_frag_count + 1);

	get_best_chain(current_record2->seq, current_record2->seq_len, kmer_size, forward_best_chain_r2, backward_best_chain_r2);

	//fprintf(stderr, "R2 Forward score:%.4f,\t len: %lu\n", forward_best_chain_r2.score, (unsigned long)forward_best_chain_r2.chain_len);
	//for (int i = 0; i < forward_best_chain_r2.chain_len; i++)
	//	fprintf(stderr, "frag[%d]: %lu\t%d\t%lu\n", i, forward_best_chain_r2.frags[i].rpos, forward_best_chain_r2.frags[i].qpos, forward_best_chain_r2.frags[i].len);
	//fprintf(stderr, "R2 Backward score:%.4f,\t len: %lu\n", backward_best_chain_r2.score, (unsigned long)backward_best_chain_r2.chain_len);
	//for (int i = 0; i < backward_best_chain_r2.chain_len; i++)
	//	fprintf(stderr, "frag[%d]: %lu\t%d\t%lu\n", i, backward_best_chain_r2.frags[i].rpos, backward_best_chain_r2.frags[i].qpos, backward_best_chain_r2.frags[i].len);

	// checking read concordancy
	if (is_concord(forward_best_chain_r1, current_record1->seq_len, kmer_size) and is_concord(backward_best_chain_r2, current_record2->seq_len, kmer_size)) {
		//fprintf(stderr, "%s+++++Concordant+++++\n", current_record1->rname);
		return 0;
	}
	else if (is_concord(forward_best_chain_r2, current_record2->seq_len, kmer_size) and is_concord(backward_best_chain_r1, current_record1->seq_len, kmer_size)) {
		//fprintf(stderr, "+++++Concordant+++++\n");
		return 0;
	}

	int ex_fr1 = extend_chain(forward_best_chain_r1, current_record1->seq, current_record1->seq_len);
	int ex_br1 = extend_chain(backward_best_chain_r1, current_record1->rcseq, current_record1->seq_len);
	int ex_fr2 = extend_chain(forward_best_chain_r2, current_record2->seq, current_record2->seq_len);
	int ex_br2 = extend_chain(backward_best_chain_r2, current_record2->rcseq, current_record2->seq_len);
	
	if ((ex_fr1 == 0 and ex_br2 == 0) or (ex_br1 == 0 and ex_fr2 == 0))
		return 0;

	bool no_map_r1 = false;
	bool no_map_r2 = false;

	if (ex_fr1 == 5 and ex_br1 == 5)
		no_map_r1 = true;
	if (ex_fr2 == 5 and ex_br2 == 5)
		no_map_r2 = true;
	
	if (no_map_r1 and no_map_r2)
		return 5;
	if (no_map_r1 or no_map_r2)
		return 4;

	return 3;
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

