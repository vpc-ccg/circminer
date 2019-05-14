#include "utils.h"
#include "common.h"
#include "align.h"
#include "gene_annotation.h"
#include "match_read.h"


void get_mate_name(char* fq1, char* fq2) {
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
	Alignment alignment;
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

bool is_concord2(const chain_t& a, int seq_len, MatchedMate& mr) {
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
		if (a.frags[0].qpos == 0 or a.frags[a.chain_len-1].qpos + a.frags[a.chain_len-1].len == seq_len) {
			mr.type = CANDID;
		}
	}
	return mr.is_concord;
}


// sm should start before lm
bool concordant_explanation(const MatchedMate& sm, const MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm) {
	if (sm.spos > lm.spos)
		return false;

	int32_t tlen;
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

// check for one valid split mate
// -> one fully mapped mate and one partially mapped
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
		//fprintf(stderr, "R1 start ind: %d\tR2 end ind: %d\n To beg of intron: %d\n", sm.exon_ind_spos, lm.exon_ind_epos, sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos));
		if ((intronic_bs[contigNum][sm.spos]) and ((intronic_bs[contigNum][lm.spos])) and 
			(sm.exon_ind_spos >= 0) and (lm.exon_ind_epos >= 0) and (sm.exon_ind_spos == lm.exon_ind_epos) and 
			(sm.spos - gtf_parser.get_interval_epos(sm.exon_ind_spos) <= LARIAT2BEGTH)) {
			mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHIBSJ, r1_sm);
			return true;
		}
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

// check for two valid split mates
// -> two partially mapped mates
bool check_2bsj(MatchedMate& sm, MatchedMate& lm, MatchedRead& mr, const string& chr, uint32_t shift, bool r1_sm) {
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

	for (int i = 0; i < sm.exons_spos->seg_list.size(); i++)
		for (int j = 0; j < lm.exons_spos->seg_list.size(); j++)
			if (sm.exons_spos->seg_list[i].same_gene(lm.exons_spos->seg_list[j])) {
				mr.update(sm, lm, chr, shift, lm.epos - sm.spos + 1, 0, false, CHI2BSJ, r1_sm);
				return true;
			}
	return false;
}

void intersect_trans(const vector<uint32_t>& tid_l1, const vector<uint32_t>& tid_l2, vector<uint32_t>& common_tid) {
	for (int i = 0; i < tid_l1.size(); i++) {
		for (int j = 0; j < tid_l2.size(); j++) {
			uint32_t tid1 = tid_l1[i];
			uint32_t tid2 = tid_l2[j];
			if (tid1 == tid2) {
				common_tid.push_back(tid1);
				break;
			}
		}
	}
}

bool same_transcript(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r, vector<uint32_t>& common_tid) {
	common_tid.clear();
	if (s == NULL or r == NULL)
		return false;

	vector<uint32_t> seg1_tid;
	vector<uint32_t> seg2_tid;

	for (int i = 0; i < s->seg_list.size(); i++)
		for (int k = 0; k < s->seg_list[i].trans_id.size(); k++)
			seg1_tid.push_back(s->seg_list[i].trans_id[k]);

	for (int j = 0; j < r->seg_list.size(); j++)
		for (int l = 0; l < r->seg_list[j].trans_id.size(); l++)
			seg2_tid.push_back(r->seg_list[j].trans_id[l]);

	intersect_trans(seg1_tid, seg2_tid, common_tid);

	return common_tid.size() != 0;
}

bool same_transcript(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r, const IntervalInfo<UniqSeg>* q, vector<uint32_t>& common_tid) {
	common_tid.clear();
	if (s == NULL or r == NULL or q == NULL)
		return false;
	
	vector <uint32_t> sr_common_tid;
	bool sr_intersect = same_transcript(s, r, sr_common_tid);
	if (! sr_intersect)
		return false;

	vector<uint32_t> seg_tid;

	for (int i = 0; i < s->seg_list.size(); i++)
		for (int k = 0; k < s->seg_list[i].trans_id.size(); k++)
			seg_tid.push_back(s->seg_list[i].trans_id[k]);

	intersect_trans(sr_common_tid, seg_tid, common_tid);

	return common_tid.size() != 0;
}

bool same_transcript(const IntervalInfo<UniqSeg>* s, const IntervalInfo<UniqSeg>* r, 
					 const IntervalInfo<UniqSeg>* q, const IntervalInfo<UniqSeg>* p, vector<uint32_t>& common_tid) {
	
	common_tid.clear();
	if (s == NULL or r == NULL or q == NULL or p == NULL)
		return false;
	
	vector <uint32_t> sr_common_tid;
	bool sr_intersect = same_transcript(s, r, sr_common_tid);
	if (! sr_intersect)
		return false;

	vector <uint32_t> qp_common_tid;
	bool qp_intersect = same_transcript(q, p, qp_common_tid);
	if (! qp_intersect)
		return false;

	intersect_trans(sr_common_tid, qp_common_tid, common_tid);

	return common_tid.size() != 0;
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
