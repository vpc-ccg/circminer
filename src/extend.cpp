#include <set>
#include <cstdlib>
#include <stdint.h>

#include "extend.h"
#include "common.h"
#include "align.h"
#include "match_read.h"
#include "gene_annotation.h"

void get_seq_right(char* res_str, char* seq, int seq_len, uint32_t pos, bool had_junction, int remain, uint32_t ub, AlignRes& best, AlignRes& curr, bool& consecutive);

void get_seq_left(char* res_str, char* seq, int seq_len, uint32_t pos, bool had_junction, int remain, uint32_t lb,  AlignRes& best, AlignRes& curr, bool& consecutive);


// pos is exclusive
// [ pos+1, pos+len ]
bool extend_right(char* seq, uint32_t& pos, int len, uint32_t ub, AlignRes& best_alignment) {
	int seq_len = len;
	int ref_len = len + INDELTH;
	uint32_t orig_pos = pos;

	char res_str[ref_len + 4*INDELTH];
	res_str[ref_len] = '\0';
	
	int indel;
	bool consecutive = false;
	best_alignment.set(pos, EDTH + 1, SOFTCLIPTH + 1, INDELTH + 1, 0);
	AlignRes curr_alignment(ub);
	
	get_seq_right(res_str, seq, seq_len, pos, false, ref_len, ub, best_alignment, curr_alignment, consecutive);

	uint32_t best_rmpos = best_alignment.pos;
	int min_ed = best_alignment.ed;
	int sclen_best = best_alignment.sclen;

	//vafprintf(2, stderr, "Min Edit Dist: %d\tNew RM POS: %u\tCovered len: %d\n", min_ed, best_rmpos, best_alignment.rcovlen);

	if (min_ed <= EDTH) {
		pos = best_rmpos - sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew RM POS: %u\n", min_ed, pos);
		if (best_alignment.qcovlen >= seq_len)
			return true;
	}
	
	if (consecutive)
		return false;

	// intron retention
	if (!pac2char(orig_pos + 1, ref_len, res_str))
		return false;

	min_ed = alignment.local_alignment_right_sc(res_str, ref_len, seq, seq_len, sclen_best, indel);
	
	vafprintf(2, stderr, "Intron Retention:\nrmpos: %lu\textend len: %d\n", orig_pos, len);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist %d\n", res_str, seq, min_ed);

	if (min_ed <= EDTH) {
		curr_alignment.set(orig_pos + seq_len - indel, min_ed, sclen_best, indel, seq_len);
		best_alignment.update_right(curr_alignment);
		pos = orig_pos + seq_len - indel - sclen_best;
		vafprintf(2, stderr, "Intron Retention: Min Edit Dist: %d\tNew RM POS: %u\n", min_ed, pos);
		return true;
	}

	return false;
}

// pos is exclusive
// [ pos-len, pos-1 ]
bool extend_left(char* seq, uint32_t& pos, int len, uint32_t lb, AlignRes& best_alignment) {
	int seq_len = len;
	int ref_len = len + INDELTH;
	uint32_t orig_pos = pos;

	char res_str[ref_len + 4*INDELTH];
	res_str[ref_len] = '\0';
	
	int indel;
	bool consecutive = false;
	best_alignment.set(pos, EDTH + 1, SOFTCLIPTH + 1, INDELTH + 1, 0);
	AlignRes curr_alignment(lb);

	get_seq_left(res_str, seq, seq_len, pos, false, ref_len, lb, best_alignment, curr_alignment, consecutive);
		
	uint32_t lmpos_best = best_alignment.pos;
	int min_ed = best_alignment.ed;
	int sclen_best = best_alignment.sclen;

	//vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\tCovered len: %d\n", min_ed, lmpos_best, best_alignment.rcovlen);

	if (min_ed <= EDTH) {
		pos = lmpos_best + sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\n", min_ed, pos);
		if (best_alignment.qcovlen >= seq_len)
			return true;
	}

	if (consecutive)
		return false;

	// intron retention
	if (!pac2char(orig_pos - ref_len, ref_len, res_str))
		return false;

	min_ed = alignment.local_alignment_left_sc(res_str, ref_len, seq, seq_len, sclen_best, indel);
	
	vafprintf(2, stderr, "Intron Retention:\nlmpos: %lu\textend len: %d\n", orig_pos, len);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist %d\n", res_str, seq, min_ed);

	if (min_ed <= EDTH) {
		curr_alignment.set(orig_pos - seq_len + indel, min_ed, sclen_best, indel, seq_len);
		best_alignment.update_left(curr_alignment);
		pos = orig_pos - seq_len + indel + sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\n", min_ed, pos);
		return true;
	}

	return false;
}


// runs recursivly on consezutive exons
// pos is exclusive
// [ pos+1, pos+len ]
void get_seq_right(char* res_str, char* seq, int seq_len, uint32_t pos, bool had_junction, int remain, uint32_t ub, AlignRes& best, AlignRes& curr, bool& consecutive) {
	if (seq_len <= 0)
		return;
	int sclen;
	int indel;
	int edit_dist;
	uint32_t new_rmpos;
	uint32_t exon_remain;

	bool dummy_consec = false;

	// if we are on an exon after a junction -> pos = (next_exon_beg - 1) = > search on (pos + 1)
	uint32_t search_pos = (! had_junction) ? pos : pos + 1;
	const IntervalInfo<UniqSeg>* overlapped_exon = gtf_parser.get_location_overlap(search_pos, true);

	if (overlapped_exon == NULL or overlapped_exon->seg_list.size() == 0) { // there is enough distance to the end of exon / intron
		vafprintf(2, stderr, "Right 1 Going for %lu - %lu\n", pos + 1, pos + remain);
		
		if (!pac2char(pos + 1, remain, res_str))
			return;

		edit_dist = alignment.local_alignment_right_sc(res_str, remain, seq, seq_len, sclen, indel);
	
		consecutive = true;
		new_rmpos = pos + seq_len - indel;
		
		vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\t sclen: %d\n", new_rmpos, seq_len, indel, edit_dist, sclen);
		vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", res_str, seq);
		
		if (curr.ed + edit_dist <= EDTH) {
			curr.update(edit_dist, sclen, new_rmpos, indel, seq_len);
			best.update_right(curr);
		}

		return;
	}

	set <GenRegion> trans_extensions;
	GenRegion new_region;

	AlignRes orig(curr.pos, curr.ed, curr.sclen, curr.indel, curr.qcovlen);

	for (int i = 0; i < overlapped_exon->seg_list.size(); i++) {
		// reset to orig:
		curr.set(orig.pos, orig.ed, orig.sclen, orig.indel, orig.qcovlen);

		vafprintf(2, stderr, "This exon: [%u-%u] -> %u\n", overlapped_exon->seg_list[i].start, overlapped_exon->seg_list[i].end, overlapped_exon->seg_list[i].next_exon_beg);
		if (had_junction and overlapped_exon->seg_list[i].start != search_pos)	// when jump to the next exon
			continue;

		exon_remain = overlapped_exon->seg_list[i].end - pos;
		vafprintf(2, stderr, "Remain on exon: %u\n", exon_remain);
		if (exon_remain + INDELTH >= remain) {	// exonic	
			int ref_remain = minM(remain, exon_remain);	// will not allow insertions greater than the remaining size on exon
			new_region.set(pos + ref_remain, 0);
			if (trans_extensions.find(new_region) != trans_extensions.end())
				continue;

			trans_extensions.insert(new_region);

			vafprintf(2, stderr, "Right 2 Going for %lu - %lu\n", pos + 1, pos + ref_remain);
			if (!pac2char(pos + 1, ref_remain, res_str))
				return;

			edit_dist = alignment.local_alignment_right_sc(res_str, ref_remain, seq, seq_len, sclen, indel);
		
			consecutive = (ref_remain == remain);
			new_rmpos = pos + seq_len - indel;
			
			vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\tsclen: %d\n", new_rmpos, seq_len, indel, edit_dist, sclen);
			vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", res_str, seq);
			
			if (curr.ed + edit_dist <= EDTH) {
				curr.update(edit_dist, sclen, new_rmpos, indel, seq_len);
				best.update_right(curr);
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

			vafprintf(2, stderr, "Right 3 Going for %lu - %lu\n", pos + 1, pos + exon_remain);
			if (!pac2char(pos + 1, exon_remain, res_str))
				return;

			int seq_remain = minM(exon_remain + INDELTH, seq_len);
			edit_dist = alignment.local_alignment_right(seq, seq_remain, res_str, exon_remain, indel);

			new_rmpos = pos + exon_remain - indel;

			vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\tsclen: %d\n", new_rmpos, exon_remain, indel, edit_dist, sclen);
			vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", res_str, seq);

			if (curr.ed + edit_dist <= EDTH) {
				curr.update(edit_dist, 0, new_rmpos, indel, exon_remain);
				best.update_right(curr);
				get_seq_right(res_str + exon_remain, seq + exon_remain + indel ,seq_len - exon_remain - indel, overlapped_exon->seg_list[i].next_exon_beg - 1, true, remain - exon_remain, ub, best, curr, dummy_consec);
			}
		}
	}
}

// pos is exclusive
// [ pos-len, pos-1 ]
void get_seq_left(char* res_str, char* seq, int seq_len, uint32_t pos, bool had_junction, int remain, uint32_t lb, AlignRes& best, AlignRes& curr, bool& consecutive) {
	if (seq_len <= 0)
		return;

	int sclen;
	int indel;
	int edit_dist;
	uint32_t new_lmpos;
	uint32_t exon_remain;

	bool dummy_consec = false;
	
	// if we are on an exon after a junction -> pos = (prev_exon_end + 1) = > search on (pos - 1)
	uint32_t search_pos = (! had_junction) ? pos : pos - 1;
	const IntervalInfo <UniqSeg>* overlapped_exon = gtf_parser.get_location_overlap(search_pos, false);

	if (overlapped_exon == NULL or overlapped_exon->seg_list.size() == 0) {
		vafprintf(2, stderr, "Left 1 Going for %lu - %lu\n", pos - remain, pos - 1);

		if (!pac2char(pos - remain, remain, res_str))
			return;

		edit_dist = alignment.local_alignment_left_sc(res_str, remain, seq, seq_len, sclen, indel);
		
		consecutive = true; 
		new_lmpos = pos - seq_len + indel;

		vafprintf(2, stderr, "lmpos: %lu\textend len: %d\t indel: %d\tLeft edit dist: %d\n", new_lmpos, seq_len, indel, edit_dist);
		vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", res_str, seq);

		if (curr.ed + edit_dist <= EDTH) {
			curr.update(edit_dist, sclen, new_lmpos, indel, seq_len);
			best.update_left(curr);
		}

		return;
	}

	set <GenRegion> trans_extensions;
	GenRegion new_region;

	AlignRes orig(curr.pos, curr.ed, curr.sclen, curr.indel, curr.qcovlen);

	for (int i = 0; i < overlapped_exon->seg_list.size(); i++) {
		// reset to orig:
		curr.set(orig.pos, orig.ed, orig.sclen, orig.indel, orig.qcovlen);
		vafprintf(2, stderr, "%u <- This exon: [%u-%u]\n", overlapped_exon->seg_list[i].prev_exon_end, overlapped_exon->seg_list[i].start, overlapped_exon->seg_list[i].end);
		if (had_junction and overlapped_exon->seg_list[i].end != search_pos)	// when jump to prev exon
			continue;

		exon_remain = pos - overlapped_exon->seg_list[i].start;
		vafprintf(2, stderr, "Remain on exon: %u\n", exon_remain);
		if (exon_remain + INDELTH >= remain) {
			int ref_remain = minM(remain, exon_remain);	// will not allow insertions greater than the remaining size on exon
			new_region.set(pos - ref_remain, 0);
			if (trans_extensions.find(new_region) != trans_extensions.end())
				continue;

			trans_extensions.insert(new_region);

			vafprintf(2, stderr, "Left 2 Going for %lu - %lu\n", pos - ref_remain, pos - 1);

			if (!pac2char(pos - ref_remain, ref_remain, res_str))
				return;

			edit_dist = alignment.local_alignment_left_sc(res_str, ref_remain, seq, seq_len, sclen, indel);
		
			consecutive = (ref_remain == remain);
			new_lmpos = pos - seq_len + indel;

			vafprintf(2, stderr, "lmpos: %lu\textend len: %d, indel: %d\tLeft edit dist: %d\n", new_lmpos, seq_len, indel, edit_dist);
			vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", res_str, seq);

			if (curr.ed + edit_dist <= EDTH) {
				curr.update(edit_dist, sclen, new_lmpos, indel, seq_len);
				best.update_left(curr);
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

			vafprintf(2, stderr, "Left 3 Going for %lu - %lu\n", overlapped_exon->seg_list[i].start, overlapped_exon->seg_list[i].start + exon_remain - 1);

			if (!pac2char(overlapped_exon->seg_list[i].start, exon_remain, res_str))
				return;

			int seq_remain = minM(exon_remain + INDELTH, seq_len);
			edit_dist = alignment.local_alignment_left(seq + seq_len - seq_remain, seq_remain, res_str, exon_remain, indel);
			new_lmpos = pos - exon_remain + indel;

			vafprintf(2, stderr, "lmpos: %lu\textend len: %d\tLeft edit dist: %d\n", new_lmpos, exon_remain, edit_dist);
			vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", res_str, seq + seq_len - seq_remain);

			if (curr.ed + edit_dist <= EDTH) {
				curr.update(edit_dist, 0, new_lmpos, indel, exon_remain);
				best.update_left(curr);
				get_seq_left(res_str, seq, seq_len - exon_remain - indel, overlapped_exon->seg_list[i].prev_exon_end + 1, true, remain - exon_remain, lb, best, curr, dummy_consec);
			}
		}
	}
}
