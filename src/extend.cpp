#include <set>
#include <cstdlib>
#include <stdint.h>

#include "extend.h"
#include "common.h"
#include "align.h"
#include "match_read.h"
#include "gene_annotation.h"

void get_seq_right(char* res_str, char* seq, int seq_len, uint32_t pos, bool had_junction, int remain, uint32_t ub, AlignRes& best, AlignRes& curr, bool& consecutive);
void get_seq_left (char* res_str, char* seq, int seq_len, uint32_t pos, bool had_junction, int remain, uint32_t lb, AlignRes& best, AlignRes& curr, bool& consecutive);

void extend_right_trans(uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
						uint32_t ub,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res);
void extend_left_trans (uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
						uint32_t lb,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res);


// pos is exclusive
// [ pos+1, pos+len ]
bool extend_right(const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, uint32_t ub, 
					AlignRes& best_alignment) {
	int seq_len = len;
	int ref_len = len + INDELTH;
	uint32_t orig_pos = pos;

	char ref_seq[maxReadLength + 4 * INDELTH];
	ref_seq[ref_len] = '\0';
	
	int indel;
	bool consecutive = false;
	AlignRes curr_alignment(ub);
	best_alignment.set(pos, EDTH + 1, SOFTCLIPTH + 1, INDELTH + 1, 0);

	map <AllCoord, AlignRes> align_res; 
	
	for (int i = 0; i < common_tid.size(); i++) {
		extend_right_trans(common_tid[i], pos, ref_seq, ref_len, seq, seq_len, ub, best_alignment, consecutive, align_res);
		//best_alignment.print();
		// if (best_alignment.qcovlen >= seq_len and best_alignment.ed == 0 and best_alignment.sclen == 0) {
		// 	pos = best_alignment.pos - best_alignment.sclen;
		// 	return true;
		// }
	}

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
	if (!pac2char(orig_pos + 1, ref_len, ref_seq))
		return false;

	min_ed = alignment.local_alignment_right_sc(ref_seq, ref_len, seq, seq_len, sclen_best, indel);
	
	vafprintf(2, stderr, "Intron Retention:\nrmpos: %lu\textend len: %d\n", orig_pos, len);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist %d\n", ref_seq, seq, min_ed);

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
bool extend_left(const vector <uint32_t>& common_tid, char* seq, uint32_t& pos, int len, uint32_t lb, 
					AlignRes& best_alignment) {
	int seq_len = len;
	int ref_len = len + INDELTH;
	uint32_t orig_pos = pos;

	char ref_seq[maxReadLength + 4 * INDELTH];
	ref_seq[ref_len] = '\0';
	
	int indel;
	bool consecutive = false;
	AlignRes curr_alignment(lb);
	best_alignment.set(pos, EDTH + 1, SOFTCLIPTH + 1, INDELTH + 1, 0);

	map <AllCoord, AlignRes> align_res;

	for (int i = 0; i < common_tid.size(); i++) {
		extend_left_trans(common_tid[i], pos, ref_seq, ref_len, seq, seq_len, lb, best_alignment, consecutive, align_res);
		//best_alignment.print();
		// if (best_alignment.qcovlen >= seq_len and best_alignment.ed == 0 and best_alignment.sclen == 0) {
		// 	pos = best_alignment.pos + best_alignment.sclen;
		// 	return true;
		// }
	}
		
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
	if (!pac2char(orig_pos - ref_len, ref_len, ref_seq))
		return false;

	min_ed = alignment.local_alignment_left_sc(ref_seq, ref_len, seq, seq_len, sclen_best, indel);
	
	vafprintf(2, stderr, "Intron Retention:\nlmpos: %lu\textend len: %d\n", orig_pos, len);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\nedit dist %d\n", ref_seq, seq, min_ed);

	if (min_ed <= EDTH) {
		curr_alignment.set(orig_pos - seq_len + indel, min_ed, sclen_best, indel, seq_len);
		best_alignment.update_left(curr_alignment);
		pos = orig_pos - seq_len + indel + sclen_best;
		vafprintf(2, stderr, "Min Edit Dist: %d\tNew LM POS: %u\n", min_ed, pos);
		return true;
	}

	return false;
}

// returns true iff extension was successful
bool extend_right_middle(uint32_t pos, char* ref_seq, uint32_t exon_len, char* qseq, int qseq_len, 
							AlignRes& best, AlignRes& curr, AlignRes& exon_res) {

	vafprintf(2, stderr, "Middle Right Ext Going for %lu - %lu\n", pos + 1, pos + exon_len);
	if (!pac2char(pos + 1, exon_len, ref_seq))
		return false;

	int indel;
	int seq_remain = minM(exon_len + INDELTH, qseq_len);
	int edit_dist = alignment.local_alignment_right(qseq, seq_remain, ref_seq, exon_len, indel);

	uint32_t new_rmpos = pos + exon_len - indel;
	exon_res.set(new_rmpos, edit_dist, 0, indel, exon_len);

	vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\n", 
							new_rmpos, exon_len, indel, edit_dist);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", ref_seq, qseq);

	if (curr.ed + edit_dist <= EDTH) {
		curr.update(edit_dist, 0, new_rmpos, indel, exon_len);
		best.update_right(curr);
		return true;
	}
	return false;
}

void extend_right_end(uint32_t pos, char* ref_seq, uint32_t ref_len, char* qseq, int qseq_len, 
						AlignRes& best, AlignRes& curr, AlignRes& exon_res) {

	vafprintf(2, stderr, "Final Right Ext Going for %lu - %lu\n", pos + 1, pos + ref_len);
	if (!pac2char(pos + 1, ref_len, ref_seq))
		return;

	int sclen, indel;
	int edit_dist = alignment.local_alignment_right_sc(ref_seq, ref_len, qseq, qseq_len, sclen, indel);

	uint32_t new_rmpos = pos + qseq_len - indel;
	exon_res.set(new_rmpos, edit_dist, sclen, indel, qseq_len);
	
	vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\tsclen: %d\n", 
							new_rmpos, qseq_len, indel, edit_dist, sclen);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", ref_seq, qseq);
	
	if (curr.ed + edit_dist <= EDTH) {
		curr.update(edit_dist, sclen, new_rmpos, indel, qseq_len);
		best.update_right(curr);
	}
}

// [pos + 1, pos + len]
void extend_right_trans(uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
						uint32_t ub, AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res) {
	consecutive = false;
	AlignRes curr(ub);
	AlignRes exon_res(ub);

	int it_ind;
	const IntervalInfo<UniqSeg>* it_seg = gtf_parser.get_location_overlap_ind(pos, false, it_ind);
	if (it_seg == NULL) {	// probably wrong chaining to intron 
		return;
		it_seg = gtf_parser.get_interval(it_ind);
		int diff = pos - it_seg->epos;
		qseq -= diff;
		qseq_len += diff;
		ref_len += diff;
		pos = gtf_parser.get_interval(it_ind + 1)->spos - 1;
	}
	
	int it_ind_start = gtf_parser.get_trans_start_ind(contigName, tid);
	int rel_ind = it_ind - it_ind_start;
	int curr_exon_start_ind = it_ind;
	int curr_exon_end_ind = it_ind;
	
	uint32_t rspos = pos;
	int exon_len = it_seg->epos - pos;
	int remain_ref_len = ref_len;
	int covered = 0;
	int indel;

	for (int i = rel_ind + 1; i < gtf_parser.trans2seg[contigName][tid].size(); i++) {
		if (exon_len >= qseq_len - covered)
			break;
		if (gtf_parser.trans2seg[contigName][tid][i] == 1) {
			// go for alignment
			indel = 0;
			if (exon_len > 0) {
				if (rspos + exon_len > ub) {
					return;
				}

				int remain_qseq_len = minM(exon_len + INDELTH, qseq_len - covered);

				// search in map
				AllCoord tmp_coord(rspos, exon_len, covered, remain_qseq_len);
				map <AllCoord, AlignRes>::iterator it;
				it = align_res.find(tmp_coord);
				if (it != align_res.end()) {
					vafprintf(2, stderr, "[Found] Middle Right Ext Going for %lu - %lu\n", rspos + 1, rspos + exon_len);
					vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\n", 
							it->second.pos, it->second.qcovlen, it->second.indel, it->second.ed);

					if (curr.ed + it->second.ed > EDTH)
						return;
					else {
						curr.update(it->second.ed, it->second.sclen, it->second.pos, it->second.indel, it->second.qcovlen);
						best.update_right(curr);
					}

					indel = it->second.indel;
				}
				else {
					bool success = extend_right_middle(rspos, ref_seq, exon_len, qseq + covered, remain_qseq_len, 
													best, curr, exon_res);

					align_res.insert(pair <AllCoord, AlignRes>(tmp_coord, exon_res));
					if (! success)
						return;

					indel = exon_res.indel;
				}
			}
			//

			remain_ref_len -= exon_len;
			covered += exon_len + indel;
			exon_len = 0;
			curr_exon_start_ind = i + it_ind_start;
			curr_exon_end_ind = i + it_ind_start;
			it_seg = gtf_parser.get_interval(i + it_ind_start);
			rspos = it_seg->spos - 1;
		}
		if (gtf_parser.trans2seg[contigName][tid][i] != 0) {
			curr_exon_end_ind = i + it_ind_start;
			it_seg = gtf_parser.get_interval(curr_exon_end_ind);
			exon_len += it_seg->epos - it_seg->spos + 1;
		}
	}

	if (covered >= qseq_len or rspos + qseq_len - covered > ub)
		return;

	consecutive = (rspos == pos);

	remain_ref_len = minM(remain_ref_len, exon_len);
	// search in map
	AllCoord tmp_coord(rspos, remain_ref_len, covered, qseq_len - covered);
	map <AllCoord, AlignRes>::iterator it;
	it = align_res.find(tmp_coord);
	if (it != align_res.end()) {
		vafprintf(2, stderr, "[Found] Final Right Ext Going for %lu - %lu\n", rspos + 1, rspos + remain_ref_len);
		vafprintf(2, stderr, "rmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\tsclen: %d\n", 
							it->second.pos, it->second.qcovlen, it->second.indel, it->second.ed, it->second.sclen);

		if (curr.ed + it->second.ed > EDTH)
			return;
		else {
			curr.update(it->second.ed, it->second.sclen, it->second.pos, it->second.indel, it->second.qcovlen);
			best.update_right(curr);
		}
	}
	else {
		extend_right_end(rspos, ref_seq, remain_ref_len, qseq + covered, qseq_len - covered, best, curr, exon_res);
		align_res.insert(pair <AllCoord, AlignRes>(tmp_coord, exon_res));
	}
}

// returns true iff extension was successful
bool extend_left_middle(uint32_t pos, char* ref_seq, uint32_t exon_len, char* qseq, int qseq_len, 
							AlignRes& best, AlignRes& curr, AlignRes& exon_res) {

	vafprintf(2, stderr, "Middle Left Ext Going for %lu - %lu\n", pos - exon_len, pos - 1);
	if (!pac2char(pos - exon_len, exon_len, ref_seq))
		return false;

	int indel;
	int seq_remain = minM(exon_len + INDELTH, qseq_len);
	int edit_dist = alignment.local_alignment_left(qseq, seq_remain, ref_seq, exon_len, indel);

	uint32_t new_lmpos = pos - exon_len + indel;
	exon_res.set(new_lmpos, edit_dist, 0, indel, exon_len);

	vafprintf(2, stderr, "lmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\n", 
							new_lmpos, exon_len, indel, edit_dist);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", ref_seq, qseq);

	if (curr.ed + edit_dist <= EDTH) {
		curr.update(edit_dist, 0, new_lmpos, indel, exon_len);
		best.update_left(curr);
		return true;
	}
	return false;
}

void extend_left_end(uint32_t pos, char* ref_seq, uint32_t ref_len, char* qseq, int qseq_len, 
						AlignRes& best, AlignRes& curr, AlignRes& exon_res) {

	vafprintf(2, stderr, "Final Left Ext Going for %lu - %lu\n", pos - ref_len, pos - 1);
	if (!pac2char(pos - ref_len, ref_len, ref_seq))
		return;

	int sclen, indel;
	int edit_dist = alignment.local_alignment_left_sc(ref_seq, ref_len, qseq, qseq_len, sclen, indel);

	uint32_t new_lmpos = pos - qseq_len + indel;
	exon_res.set(new_lmpos, edit_dist, sclen, indel, qseq_len);
	
	vafprintf(2, stderr, "lmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\tsclen: %d\n", 
							new_lmpos, qseq_len, indel, edit_dist, sclen);
	vafprintf(2, stderr, "str beg str:  %s\nread beg str: %s\n", ref_seq, qseq);
	
	if (curr.ed + edit_dist <= EDTH) {
		curr.update(edit_dist, sclen, new_lmpos, indel, qseq_len);
		best.update_left(curr);
	}
}

// [pos - len, pos - 1]
void extend_left_trans (uint32_t tid, uint32_t pos, char* ref_seq, int ref_len, char* qseq, int qseq_len, 
						uint32_t lb,  AlignRes& best, bool& consecutive, map <AllCoord, AlignRes>& align_res) {

	consecutive = false;
	AlignRes curr(lb);
	AlignRes exon_res(lb);

	int it_ind;
	const IntervalInfo<UniqSeg>* it_seg = gtf_parser.get_location_overlap_ind(pos, false, it_ind);
	if (it_seg == NULL) {	// probably wrong chaining to intron 
		return;
		it_seg = gtf_parser.get_interval(it_ind + 1);
		int diff = it_seg->spos - pos;
		qseq_len += diff;
		ref_len += diff;
		pos = it_seg->spos + 1;
	}

	int it_ind_start = gtf_parser.get_trans_start_ind(contigName, tid);
	int rel_ind = it_ind - it_ind_start;
	int curr_exon_start_ind = it_ind;
	int curr_exon_end_ind = it_ind;
	
	uint32_t lepos = pos;
	int exon_len = 0;
	int remain_ref_len = ref_len;
	int covered = 0;
	int indel;
	bool first_seg = true;

	for (int i = rel_ind; i >= 0; i--) {
		if (gtf_parser.trans2seg[contigName][tid][i] != 0) {
			curr_exon_start_ind = i + it_ind_start;
			it_seg = gtf_parser.get_interval(curr_exon_start_ind);
			if (first_seg) {
				exon_len = pos - it_seg->spos;
				first_seg = false;
			}
			else {
				if (exon_len == 0) {
					curr_exon_start_ind = i + it_ind_start;
					curr_exon_end_ind = i + it_ind_start;
					lepos = it_seg->epos + 1;
				}
				exon_len += it_seg->epos - it_seg->spos + 1;
			}
		}

		if (exon_len >= qseq_len - covered)
			break;

		if (gtf_parser.trans2seg[contigName][tid][i] == 1) {
			// go for alignment
			indel = 0;
			if (exon_len > 0) {
				if (lepos < lb + exon_len) {
					return;
				}

				int remain_qseq_len = minM(exon_len + INDELTH, qseq_len - covered);

				// search in map
				AllCoord tmp_coord(lepos, exon_len, covered, remain_qseq_len);
				map <AllCoord, AlignRes>::iterator it;
				it = align_res.find(tmp_coord);
				if (it != align_res.end()) {
					vafprintf(2, stderr, "[Found] Middle Left Ext Going for %lu - %lu\n", lepos - exon_len, lepos - 1);
					vafprintf(2, stderr, "lmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\n", 
							it->second.pos, it->second.qcovlen, it->second.indel, it->second.ed);

					if (curr.ed + it->second.ed > EDTH)
						return;
					else {
						curr.update(it->second.ed, it->second.sclen, it->second.pos, it->second.indel, it->second.qcovlen);
						best.update_left(curr);
					}

					indel = it->second.indel;
				}
				else {
					bool success = extend_left_middle(lepos, ref_seq, exon_len, qseq + qseq_len - covered - remain_qseq_len, remain_qseq_len, 
													best, curr, exon_res);

					align_res.insert(pair <AllCoord, AlignRes>(tmp_coord, exon_res));
					if (! success)
						return;

					indel = exon_res.indel;
				}
			}
			//

			remain_ref_len -= exon_len;
			covered += exon_len + indel;
			exon_len = 0;
		}
	}

	if (covered >= qseq_len or lepos < lb + qseq_len - covered)
		return;

	consecutive = (lepos == pos);

	remain_ref_len = minM(remain_ref_len, exon_len);
	// search in map
	AllCoord tmp_coord(lepos, remain_ref_len, covered, qseq_len - covered);
	map <AllCoord, AlignRes>::iterator it;
	it = align_res.find(tmp_coord);
	if (it != align_res.end()) {
		vafprintf(2, stderr, "[Found] Final Left Ext Going for %lu - %lu\n", lepos - remain_ref_len, lepos - 1);
		vafprintf(2, stderr, "lmpos: %lu\textend len: %d\tindel: %d\tedit dist: %d\tsclen: %d\n", 
							it->second.pos, it->second.qcovlen, it->second.indel, it->second.ed, it->second.sclen);

		if (curr.ed + it->second.ed > EDTH)
			return;
		else {
			curr.update(it->second.ed, it->second.sclen, it->second.pos, it->second.indel, it->second.qcovlen);
			best.update_left(curr);
		}
	}
	else {
		extend_left_end(lepos, ref_seq, remain_ref_len, qseq, qseq_len - covered, best, curr, exon_res);
		align_res.insert(pair <AllCoord, AlignRes>(tmp_coord, exon_res));
	}
}

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
