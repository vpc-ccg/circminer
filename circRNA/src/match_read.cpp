#include <cmath>
#include <cassert>
#include "match_read.h"
#include "gene_annotation.h"

extern "C" {
#include "mrsfast/Common.h"
//#include "mrsfast/CommandLineParser.h"
//#include "mrsfast/Reads.h"
//#include "mrsfast/Output.h"
#include "mrsfast/HashTable.h"
//#include "mrsfast/MrsFAST.h"
}

char keep_log1[500];
char keep_log2[500];
static int trig = 0;

void get_mate_name(char* fq1, char* fq2) {
	strcpy(fq2, fq1);
	int i = strlen(fq1) - 1;
	while (fq1[i--] != '.' and i >= 1);
	if (fq1[i] == '1')
		fq2[i] = '2';
	else if (fq1[i] == '2')
		fq2[i] = '1';
	else {
		fprintf(stderr, "PE FASTQ names are not in the correct format\n");
		exit(1);
	}
}

void sort_positions(const bwtint_t& sp, const bwtint_t& ep, const int& len, vector<bwtint_t>& forward_list, bwtint_t& flist_size, vector<bwtint_t>& backward_list, bwtint_t& blist_size) {
	int dir=0;
	bwtint_t tmp;
	bwtint_t f_ind = 0, b_ind = 0;
	forward_list.clear();
	backward_list.clear();
	for (bwtint_t i = sp; i <= ep; i++) {
		tmp = get_pos(&i, len, dir);
		if (dir == 1)
			forward_list[f_ind++] = tmp;
		else
			backward_list[b_ind++] = tmp;

		//char *chr_name;
		//int32_t chr_len;
		//uint32_t chr_beg;
		//uint32_t chr_end;
		//bwt_get_intv_info(tmp, tmp + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		//fprintf(stderr, "Chr: %s,\t%d\n", chr_name, dir);
		//fprintf(stderr, "%"PRIu64"\n", chr_beg);
	}
	
	forward_list[f_ind] = -1;
	backward_list[b_ind] = -1;

	flist_size = f_ind;
	blist_size = b_ind;

	sort(forward_list.begin(), forward_list.end() + f_ind);
		sort(backward_list.begin(), backward_list.end() + b_ind);
}

// assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
bwtint_t binary_search(const vector<bwtint_t>& list, bwtint_t beg, bwtint_t end, const bwtint_t& target) {
	if (end - beg <= 1)
		return end;
	bwtint_t mid = (beg + end) / 2;
	if (target < list[mid])
		return binary_search(list, beg, mid, target);
	else 
		return binary_search(list, mid, end, target);
}

bool frag_comp(const fragment_t& a, const fragment_t& b) {
	return a.rpos < b.rpos;
}

// assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
int frag_binary_search(const vector<fragment_t>& list, int beg, int end, int target) {
	if (end - beg <= 1)
		return end;
	int mid = (beg + end) / 2;
	if (target < list[mid].rpos)
		return frag_binary_search(list, beg, mid, target);
	else 
		return frag_binary_search(list, mid, end, target);
}

// assumes sorted lists
void prune_frag(const ExactMatchRes& large_list,
				const vector<fragment_t>& pre_forward_fragments, int ff_ref_cnt, 
				const vector<fragment_t>& pre_backward_fragments, int bf_ref_cnt, 
				vector<fragment_t>& new_forward_fragments, int& new_forward_fragment_count, 
				vector<fragment_t>& new_backward_fragments, int& new_backward_fragment_count) {

	int valid = 0;
	int ind;
	int dir;
	int end_pos;
	bwtint_t rpos, j;
	for (j = large_list.sp; j <= large_list.ep; j++) {
		rpos = get_pos(&j, large_list.matched_len, dir);
		if (dir == 1) {
			ind = frag_binary_search(pre_forward_fragments, 0, ff_ref_cnt, rpos);
			if (((ind < ff_ref_cnt) and (pre_forward_fragments[ind].rpos - rpos <= FLGENETH)) or ((ind > 0) and (rpos - pre_forward_fragments[ind-1].rpos <= FLGENETH))) {
				valid++;
				if (valid > FRAGLIM) {
					new_forward_fragment_count = 0;
					new_backward_fragment_count = 0;
					return;
				}
				new_forward_fragments[new_forward_fragment_count].qpos = large_list.q_ind;
				new_forward_fragments[new_forward_fragment_count].rpos = rpos;
				new_forward_fragments[new_forward_fragment_count].len = large_list.matched_len;
				new_forward_fragment_count++;
			}
		}
		else {
			ind = frag_binary_search(pre_backward_fragments, 0, bf_ref_cnt, rpos);
			if (((ind < bf_ref_cnt) and (pre_backward_fragments[ind].rpos - rpos <= FLGENETH)) or ((ind > 0) and (rpos - pre_backward_fragments[ind-1].rpos <= FLGENETH))) {
				valid++;
				if (valid > FRAGLIM) {
					new_forward_fragment_count = 0;
					new_backward_fragment_count = 0;
					return;
				}
				end_pos = large_list.q_ind + large_list.matched_len - 1;
				new_backward_fragments[new_backward_fragment_count].qpos = -1 * end_pos;
				new_backward_fragments[new_backward_fragment_count].rpos = rpos;
				new_backward_fragments[new_backward_fragment_count].len = large_list.matched_len;
				new_backward_fragment_count++;
			}
		}
	}
}

int get_exact_locs_hash(char* seq, int32_t qpos, uint32_t len, GIMatchedKmer* mk) {
	vafprintf(2, stderr, "Seq: %s\tHash Val: %d\nWindow Size:%d\n",seq, hashVal(seq), WINDOW_SIZE);
	mk->frag_count = 0;
	mk->frags = NULL;
	mk->qpos = qpos;

	GeneralIndex *it = getCandidates(hashVal(seq));
	uint32_t i, j;

	mk->frag_count = 0;
	if (it == NULL) {
		return 0;
	}

	char* checksum_beg = seq + WINDOW_SIZE;
	uint32_t lb = 1;
	uint32_t ub = it[0].info;
	uint32_t mid;
	uint16_t target = checkSumVal(checksum_beg);
	uint32_t LB = 0;
	uint32_t UB = 0;

	//fprintf(stderr, "E1\n");

	while (lb < ub) {
		mid = (lb + ub) / 2;
		if (target <= it[mid].checksum)
			ub = mid;
		else
			lb = mid + 1;
	}
	
	if (ub < lb || target != it[lb].checksum)
		return 0;
	UB = LB = lb;

	lb = LB;
	ub = it[0].info;
	while (lb < ub) {
		mid = (lb + ub + 1) / 2;
		if (target < it[mid].checksum)
			ub = mid - 1;
		else
			lb = mid;
	}

	if (target == it[lb].checksum)
		UB = lb;
	
	//if (LB > 1)
	//	fprintf(stderr, "<<%d\t", it[LB-1].checksum);
	//fprintf(stderr, "%d\t", it[LB].checksum);
	//fprintf(stderr, "%d\t", it[UB].checksum);
	//if (UB < it[0].info)
	//	fprintf(stderr, "%d>>\n", it[UB+1].checksum);
	//else
	//	fprintf(stderr, "\n");


	//fprintf(stderr, "Size: %d\n", UB-LB+1);

	if (UB-LB+1 > FRAGLIM) {
		mk->frag_count = 0;
		return 0;
	}

	mk->frag_count = UB - LB + 1;
	mk->frags = it + LB;

	//for (i = LB; i <= UB; i++) {
	//	vafprintf(2, stderr, "loc: %lu\tchecksum %d\n", it[i].info, it[i].checksum);
	//	mk->frags[i-LB].rpos = it[i].info;
	//	mk->frags[i-LB].qpos = qpos;
	//	mk->frags[i-LB].len = len;
	//}

	vafprintf(2, stderr, "Occ: %lu\n", UB - LB + 1);

	return UB - LB + 1;
}

int get_exact_locs_hash_linear(char* seq, int32_t qpos, uint32_t len, MatchedKmer* mk) {
	vafprintf(2, stderr, "Seq: %s\tHash Val: %d\nWindow Size:%d\n",seq, hashVal(seq), WINDOW_SIZE);
	GeneralIndex *it = getCandidates(hashVal(seq));
	uint32_t i, j;

	mk->frag_count = 0;
	if (it == NULL) {
		return 0;
	}

	char* checksum_beg = seq + WINDOW_SIZE;
	for (i = 1; i <= it[0].info; i++) {
		if (checkSumVal(checksum_beg) == it[i].checksum) {
			for (j = i; j <= it[0].info; j++) { 
				if (j - i + 1 > FRAGLIM) {
					mk->frag_count = 0;
					return 0;
				}
				if (checkSumVal(checksum_beg) != it[j].checksum)
					break;
				mk->frags[j-i].rpos = it[j].info;
				vafprintf(2, stderr, "loc: %lu\tchecksum %d\n", it[j].info, it[j].checksum);
				mk->frags[j-i].qpos = qpos;
				mk->frags[j-i].len = len;
			}
			mk->frag_count = j - i;
			vafprintf(2, stderr, "Occ: %lu\n", j - i);
			return j - i;
		}
	}
	return 0;
}

void move_on_ll(MatchedKmer*& mk, int steps) {
	while (steps--) {
		assert(mk != NULL);
		mk = mk->next;
	}
}

// return:
// # valid kmers
// is valid if: #fragments > 0 and < FRAGLIM
int kmer_match_skip_hash(char* rseq, int rseq_len, int kmer_size, int shift, int skip, int ll_step, GIMatchedKmer* mk_res, int& em_count) {
	bwtint_t sp, ep;
	int i, occ;
	int sum = 0;
	int dir;
	uint32_t match_len = kmer_size;
	int32_t end_pos;

	//char* chr_name;
	//uint32_t chr_beg;
	//uint32_t chr_end;
	//int32_t chr_len;

	em_count = 0;
	int invalid_kmer = 0;
	for (i = shift; i < rseq_len; i += skip) {
		if (rseq_len - i < kmer_size)
			match_len = rseq_len - i;
		if (match_len < MINKMER) 
			break;

		if (i != shift)
			mk_res += ll_step;
			//move_on_ll(mk_res, ll_step);
		occ = get_exact_locs_hash(rseq + i, i, match_len, mk_res);
		//occ = get_exact_locs_hash1(rseq + i, i, match_len, mk_res[em_count]);
		
		vafprintf(2, stderr, "Occ: %d\tind: %d\tmatch len: %d\n", occ, i, match_len);
		if (occ <= 0) {
			occ = 0;
			invalid_kmer++;
		}

		if (occ > FRAGLIM ) {
			invalid_kmer++;
		}

		em_count++;	
	}

	return em_count - invalid_kmer;
}

// return:
// # valid kmers
// is valid if: #fragments > 0 and < FRAGLIM
int kmer_match_skip(const char* rseq, int rseq_len, int kmer_size, int shift, int skip, vector<ExactMatchRes>& exact_match_res, int& em_count) {
	bwtint_t sp, ep;
	int i, occ;
	int sum = 0;
	int dir;
	uint32_t match_len = kmer_size;
	int32_t end_pos;

	//char* chr_name;
	//uint32_t chr_beg;
	//uint32_t chr_end;
	//int32_t chr_len;

	em_count = 0;
	int invalid_kmer = 0;
	for (i = shift; i < rseq_len; i += skip) {
		if (rseq_len - i < kmer_size)
			match_len = rseq_len - i;
		if (match_len < MINKMER) 
			break;

		occ = get_exact_locs(rseq + i, match_len, &sp, &ep);
		
		vafprintf(2, stderr, "Occ: %d\tind: %d\tmatch len: %d\n", occ, i, match_len);
		if (occ <= 0) {
			occ = 0;
			invalid_kmer++;
		}

		if (occ > FRAGLIM ) {
			invalid_kmer++;
		}

		exact_match_res[em_count].sp = sp;
		exact_match_res[em_count].ep = ep;
		exact_match_res[em_count].q_ind = i;
		exact_match_res[em_count].matched_len = match_len;
		exact_match_res[em_count].occ = occ;

		em_count++;		
	}

	return em_count - invalid_kmer;
}

void swap(fragment_t& a, fragment_t& b, fragment_t& temp) {
	temp = a;
	a = b;
	b = temp;
}

void split_frags_by_strand(const ExactMatchRes& em_res, MatchedKmer* forward, MatchedKmer* backward) {
	bwtint_t j, rpos;
	int dir;
	int32_t end_pos;

	int forward_count = 0;
	int backward_count = 0;

	//vafprintf(2, stderr, "---Occ: %d\tind: %d\n", em_res.occ, em_res.q_ind);
	if (em_res.occ == 0 or em_res.occ > FRAGLIM) {
		return;
	}

	for (j = em_res.sp; j <= em_res.ep; j++) {
		rpos = get_pos(&j, em_res.matched_len, dir);
		
		// do not add a fragment that is comming from two different refrences (concat point is inside it)
		//if (! bwt_uniq_ref(rpos, rpos + exact_match_res[k].matched_len - 1))
		//	continue;

		if (dir == 1) {
			forward->frags[forward->frag_count].qpos = em_res.q_ind;
			forward->frags[forward->frag_count].rpos = rpos;
			forward->frags[forward->frag_count].len = em_res.matched_len;
			forward->frag_count++;
		}
		else {
			end_pos = em_res.q_ind + em_res.matched_len - 1;
			backward->frags[backward->frag_count].qpos = -1 * end_pos;
			backward->frags[backward->frag_count].rpos = rpos;
			backward->frags[backward->frag_count].len = em_res.matched_len;
			backward->frag_count++;
		}
	}
}

// fragment lists will be sorted by qpos
void fill_fragments_ll_hash(	
						vector<ExactMatchHash>& exact_match_res_nonov, int em_count_nonov, 
						vector<ExactMatchHash>& exact_match_res_ov, int em_count_ov,
						FragmentList& forward_fragments) { 

	//forward_fragment_count = 0;
	//backward_fragment_count = 0;

	for (int k = 0; k < em_count_nonov; k++) {
		//vafprintf(2, stderr, "+++Occ: %d\tind: %d\n", exact_match_res_nonov[k].occ, exact_match_res_nonov[k].q_ind);
		if (exact_match_res_nonov[k].occ > 0 and exact_match_res_nonov[k].occ <= FRAGLIM) {
			MatchedKmer* forward = new MatchedKmer();
			//split_frags_by_strand(exact_match_res_nonov[k], forward, backward);
			
			forward->frag_count = exact_match_res_nonov[k].occ;
			for (int i = 0; i < exact_match_res_nonov[k].occ; i++) {
				forward->frags[i].rpos = exact_match_res_nonov[k].locs[i];
				forward->frags[i].qpos = exact_match_res_nonov[k].q_ind;
				forward->frags[i].len = exact_match_res_nonov[k].matched_len;
			}

			forward_fragments.add_back(forward);
		}

		if (k < em_count_ov and exact_match_res_ov[k].occ > 0 and exact_match_res_ov[k].occ <= FRAGLIM) {
			MatchedKmer* forward_ov = new MatchedKmer();
			//split_frags_by_strand(exact_match_res_ov[k], forward_ov, backward_ov);
			
			forward_ov->frag_count = exact_match_res_ov[k].occ;
			for (int i = 0; i < exact_match_res_ov[k].occ; i++) {
				forward_ov->frags[i].rpos = exact_match_res_ov[k].locs[i];
				forward_ov->frags[i].qpos = exact_match_res_ov[k].q_ind;
				forward_ov->frags[i].len = exact_match_res_ov[k].matched_len;
			}

			forward_fragments.add_back(forward_ov);
		}
	}
	//forward_fragments.print();
	//backward_fragments.print();
}

// fragment lists will be sorted by qpos
void fill_fragments_ll(	vector<ExactMatchRes>& exact_match_res_nonov, int em_count_nonov, 
						vector<ExactMatchRes>& exact_match_res_ov, int em_count_ov,
						FragmentList& forward_fragments,  
						FragmentList& backward_fragments) {

	//forward_fragment_count = 0;
	//backward_fragment_count = 0;

	for (int k = 0; k < em_count_nonov; k++) {
		//vafprintf(2, stderr, "+++Occ: %d\tind: %d\n", exact_match_res_nonov[k].occ, exact_match_res_nonov[k].q_ind);
		if (exact_match_res_nonov[k].occ > 0 and exact_match_res_nonov[k].occ <= FRAGLIM) {
			MatchedKmer* forward = new MatchedKmer();
			MatchedKmer* backward = new MatchedKmer();
			split_frags_by_strand(exact_match_res_nonov[k], forward, backward);
			forward_fragments.add_back(forward);
			backward_fragments.add_front(backward);
		}

		if (k < em_count_ov and exact_match_res_ov[k].occ > 0 and exact_match_res_ov[k].occ <= FRAGLIM) {
			MatchedKmer* forward_ov = new MatchedKmer();
			MatchedKmer* backward_ov = new MatchedKmer();
			split_frags_by_strand(exact_match_res_ov[k], forward_ov, backward_ov);
			forward_fragments.add_back(forward_ov);
			backward_fragments.add_front(backward_ov);
		}
	}
	//forward_fragments.print();
	//backward_fragments.print();
}

// fragment lists will be sorted by qpos
void fill_fragments(vector<ExactMatchRes>& exact_match_res, int em_count, vector<fragment_t>& forward_fragments, int& forward_fragment_count, vector<fragment_t>& backward_fragments, int& backward_fragment_count) {
	bwtint_t j, rpos;
	int dir;
	int32_t end_pos;

	forward_fragment_count = 0;
	backward_fragment_count = 0;

	for (int k = 0; k < em_count; k++) {
		//vafprintf(2, stderr, "---Occ: %d\tind: %d\n", exact_match_res[k].occ, exact_match_res[k].q_ind);
		if (exact_match_res[k].occ == 0 or exact_match_res[k].occ > FRAGLIM)
			continue;
		for (j = exact_match_res[k].sp; j <= exact_match_res[k].ep; j++) {
			rpos = get_pos(&j, exact_match_res[k].matched_len, dir);
			
			// do not add a fragment that is comming from two different refrences (concat point is inside it)
			//if (! bwt_uniq_ref(rpos, rpos + exact_match_res[k].matched_len - 1))
			//	continue;

			if (dir == 1) {
				forward_fragments[forward_fragment_count].qpos = exact_match_res[k].q_ind;
				forward_fragments[forward_fragment_count].rpos = rpos;
				forward_fragments[forward_fragment_count].len = exact_match_res[k].matched_len;
				forward_fragment_count++;
			}
			else {
				end_pos = exact_match_res[k].q_ind + exact_match_res[k].matched_len - 1;
				backward_fragments[backward_fragment_count].qpos = -1 * end_pos;
				backward_fragments[backward_fragment_count].rpos = rpos;
				backward_fragments[backward_fragment_count].len = exact_match_res[k].matched_len;
				backward_fragment_count++;
			}
		}
	}
	
	fragment_t temp;
	for (int i = 0; i < backward_fragment_count / 2; i++) {
		swap(backward_fragments[i], backward_fragments[backward_fragment_count - i - 1], temp);
	}
}

// fragment results will be sorted
int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* starting_node) {
	//vector<ExactMatchHash> exact_match_res_nonov(2.0 * rseq_len / kmer_size); 
	//vector<ExactMatchHash> exact_match_res_ov(2.0 * rseq_len / kmer_size); 
	int valid_nonov_kmer = 0;
	int valid_ov_kmer = 0;
	int nonov_em_count;
	int ov_em_count;

	//vector <MatchedKmer*> mk_res_nonov(4);
	//for (int i = 0; i < 4; i++) {
	//	MatchedKmer* tmp = new MatchedKmer();
	//	mk_res_nonov[i] = tmp;
	//}

	valid_nonov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, 0, kmer_size, 2, starting_node, nonov_em_count);
	//for (int i = 0; i < 4; i++) {
	//	forward_fragments.add_back(mk_res_nonov[i]);
	//}

	//if (valid_nonov_kmer < 2)
	//	valid_ov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, kmer_size / 2, kmer_size, exact_match_res_ov, ov_em_count);
	
	vafprintf(1, stderr, "Non-OV valids: %d\nOV valids: %d\n", valid_nonov_kmer, valid_ov_kmer);

	//if (valid_nonov_kmer + valid_ov_kmer > 0)
	//	fill_fragments_ll_hash(exact_match_res_nonov, nonov_em_count, exact_match_res_ov, ov_em_count, forward_fragments);
	return valid_nonov_kmer + valid_ov_kmer;

}

// pos is exclusive
// for Transcriptome
void get_reference_chunk2(uint32_t pos, int len, char* res_str) {
	res_str[0] = 0;
	char* chr_name;
	int32_t chr_len;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	uint32_t avail_len;

	if (len < 0) {	// [ pos-len, pos-1 ]
		len *= -1;
		//char* res_str = (char*) malloc(len+5);

		bwt_get_intv_info((bwtint_t) pos, (bwtint_t) pos, &chr_name, &chr_len, &chr_beg, &chr_end);
		string chr = chr_name;
		if (chr == "")
			return; 

		//fprintf(stderr, "Chrom: %s, pos: %lu\n", chr.c_str(), chr_end+1);	
		
		avail_len = (chr_beg >= len) ? len : chr_beg;
		
		bwt_str_pac2char(pos - avail_len, avail_len, res_str);
		//fprintf(stderr, "Res beg: %s\n", res_str);
	}
	else {		// [ pos+1, pos+len ]
		//char* res_str = (char*) malloc(len+5);

		bwt_get_intv_info((bwtint_t) pos, (bwtint_t) pos, &chr_name, &chr_len, &chr_beg, &chr_end);
		string chr = chr_name;
		if (chr == "")
			return; 

		//fprintf(stderr, "Chrom: %s, pos: %lu\n", chr.c_str(), chr_beg+2);
		
		avail_len = ((chr_len - 1 - pos) >= len) ? len : chr_len - 1 - pos;

		bwt_str_pac2char(pos + 1, avail_len, res_str);
		//fprintf(stderr, "Res end: %s\n", res_str);
	}
}

// pos is exclusive
// [ pos-len, pos-1 ]
void get_reference_chunk_left(uint32_t pos, int len, char* res_str) {
	res_str[0] = 0;
	char* chr_name;
	int32_t chr_len;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	//char* res_str = (char*) malloc(len+5);

	//bwt_get_intv_info((bwtint_t) (pos - len), (bwtint_t) pos, &chr_name, &chr_len, &chr_beg, &chr_end);
	//string chr = chr_name;
	string chr = "1";
	chr_beg = pos - len;
	chr_end = pos;
	chr_len = len;
	if (chr == "")
		return; 
	//fprintf(stderr, "Chrom: %s, pos: %lu\n", chr.c_str(), chr_end+1);	
	int seg_ind = gtf_parser.search_loc(chr, true, chr_end+1);	// chr_end+1: to convert to 1-based coordinate system used by GTF format
	//fprintf(stderr, "Index found: %d\n", seg_ind);
	uint32_t seg_start = gtf_parser.get_start(chr, seg_ind);
	//fprintf(stderr, "Start of exon: %lu\nChr beg: %lu\n", seg_start, chr_beg);
	if (seg_start <= chr_beg+1) {	// 0-based vs. 1-based
		bwt_str_pac2char(pos - len, len, res_str);
		//fprintf(stderr, "Res beg: %s\n", res_str);
		return;
	}
	else {
		int remain_len = seg_start - chr_beg - 1;
		int covered_len = len - remain_len;
		bwt_str_pac2char(pos - covered_len, covered_len, res_str);
		if (seg_ind == 0) {	// reached first exonic region on chromosome
			//fprintf(stderr, "Res beg: %s\n", res_str);
			return;
		}

		uint32_t prev_end = gtf_parser.get_end(chr, seg_ind - 1);
		uint32_t prev_integrated_pos = pos - covered_len - (seg_start - prev_end);
		char* remain_str = (char*) malloc(len+5);
		get_reference_chunk_left(prev_integrated_pos, remain_len, remain_str);
		//fprintf(stderr, "Remain beg: %s\n", remain_str);
		strncat(remain_str, res_str, covered_len);
		strcpy(res_str, remain_str);
		//fprintf(stderr, "Res beg: %s\n", res_str);
		free(remain_str);
		return;
	}
}

// pos is exclusive
// [ pos+1, pos+len ]
void get_reference_chunk_right(uint32_t pos, int len, char* res_str) {
	res_str[0] = 0;
	char* chr_name;
	int32_t chr_len;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	//char* res_str = (char*) malloc(len+5);

	//bwt_get_intv_info((bwtint_t) pos, (bwtint_t) (pos + len), &chr_name, &chr_len, &chr_beg, &chr_end);
	//string chr = chr_name;
	string chr = "1";
	chr_beg = pos;
	chr_end = pos + len;
	chr_len = len;

	if (chr == "")
		return; 
	//fprintf(stderr, "Chrom: %s, pos: %lu\n", chr.c_str(), chr_beg+2);
	int seg_ind = gtf_parser.search_loc(chr, false, chr_beg);	// no need to convert to 1-based coordinate, because of search implementation
	//fprintf(stderr, "Index found: %d\n", seg_ind);
	uint32_t seg_end = gtf_parser.get_end(chr, seg_ind);
	//fprintf(stderr, "End of exon: %lu\nChr end: %lu\n", seg_end, chr_end);
	if (seg_end >= chr_end+1) {	//0-based vs. 1-based
		bwt_str_pac2char(pos + 1, len, res_str);
		//fprintf(stderr, "Res end: %s\n", res_str);
		return;
	}
	else {
		int remain_len = chr_end+1 - seg_end;
		int covered_len = len - remain_len;
		//fprintf(stderr, "Covered len: %d\n", covered_len);
		bwt_str_pac2char(pos + 1, covered_len, res_str);
		//fprintf(stderr, "Here\n");
		if (gtf_parser.is_last_exonic_region(chr, seg_ind)) {	// reached first exonic region on chromosome
			//fprintf(stderr, "Res end: %s\n", res_str);
			return;
		}

		uint32_t next_start = gtf_parser.get_start(chr, seg_ind + 1);
		uint32_t next_integrated_pos = pos + covered_len + (next_start - seg_end);
		char* remain_str = (char*) malloc(len+5);
		get_reference_chunk_right(next_integrated_pos, remain_len, remain_str);
		//fprintf(stderr, "Remain end: %s\n", remain_str);
		strncat(res_str, remain_str, remain_len);
		//fprintf(stderr, "Res end: %s\n", res_str);
		free(remain_str);
		return;
	}
}

// pos is exclusive
void get_reference_chunk(uint32_t pos, int len, char* res_str) {
	res_str[0] = 0;
	char* chr_name;
	int32_t chr_len;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	if (len < 0) {	// [ pos-len, pos-1 ]
		len *= -1;
		//char* res_str = (char*) malloc(len+5);

		bwt_get_intv_info((bwtint_t) (pos - len), (bwtint_t) pos, &chr_name, &chr_len, &chr_beg, &chr_end);
		string chr = chr_name;
		if (chr == "")
			return; 
		//fprintf(stderr, "Chrom: %s, pos: %lu\n", chr.c_str(), chr_end+1);	
		int seg_ind = gtf_parser.search_loc(chr, true, chr_end+1);	// chr_end+1: to convert to 1-based coordinate system used by GTF format
		//fprintf(stderr, "Index found: %d\n", seg_ind);
		uint32_t seg_start = gtf_parser.get_start(chr, seg_ind);
		//fprintf(stderr, "Start of exon: %lu\nChr beg: %lu\n", seg_start, chr_beg);
		if (seg_start <= chr_beg+1) {	// 0-based vs. 1-based
			bwt_str_pac2char(pos - len, len, res_str);
			//fprintf(stderr, "Res beg: %s\n", res_str);
			return;
		}
		else {
			int remain_len = seg_start - chr_beg - 1;
			int covered_len = len - remain_len;
			bwt_str_pac2char(pos - covered_len, covered_len, res_str);
			if (seg_ind == 0) {	// reached first exonic region on chromosome
				//fprintf(stderr, "Res beg: %s\n", res_str);
				return;
			}

			uint32_t prev_end = gtf_parser.get_end(chr, seg_ind - 1);
			uint32_t prev_integrated_pos = pos - covered_len - (seg_start - prev_end);
			char* remain_str = (char*) malloc(len+5);
			get_reference_chunk(prev_integrated_pos, -1 * remain_len, remain_str);
			//fprintf(stderr, "Remain beg: %s\n", remain_str);
			strncat(remain_str, res_str, covered_len);
			strcpy(res_str, remain_str);
			//fprintf(stderr, "Res beg: %s\n", res_str);
			free(remain_str);
			return;
		}
	}
	else {			// [ pos+1, pos+len ]
		//char* res_str = (char*) malloc(len+5);

		bwt_get_intv_info((bwtint_t) pos, (bwtint_t) (pos + len), &chr_name, &chr_len, &chr_beg, &chr_end);
		string chr = chr_name;
		if (chr == "")
			return; 
		//fprintf(stderr, "Chrom: %s, pos: %lu\n", chr.c_str(), chr_beg+2);
		int seg_ind = gtf_parser.search_loc(chr, false, chr_beg);	// no need to convert to 1-based coordinate, because of search implementation
		//fprintf(stderr, "Index found: %d\n", seg_ind);
		uint32_t seg_end = gtf_parser.get_end(chr, seg_ind);
		//fprintf(stderr, "End of exon: %lu\nChr end: %lu\n", seg_end, chr_end);
		if (seg_end >= chr_end+1) {	//0-based vs. 1-based
			bwt_str_pac2char(pos + 1, len, res_str);
			//fprintf(stderr, "Res end: %s\n", res_str);
			return;
		}
		else {
			int remain_len = chr_end+1 - seg_end;
			int covered_len = len - remain_len;
			//fprintf(stderr, "Covered len: %d\n", covered_len);
			bwt_str_pac2char(pos + 1, covered_len, res_str);
			//fprintf(stderr, "Here\n");
			if (gtf_parser.is_last_exonic_region(chr, seg_ind)) {	// reached first exonic region on chromosome
				//fprintf(stderr, "Res end: %s\n", res_str);
				return;
			}

			uint32_t next_start = gtf_parser.get_start(chr, seg_ind + 1);
			uint32_t next_integrated_pos = pos + covered_len + (next_start - seg_end);
			char* remain_str = (char*) malloc(len+5);
			get_reference_chunk(next_integrated_pos, remain_len, remain_str);
			//fprintf(stderr, "Remain end: %s\n", remain_str);
			strncat(res_str, remain_str, remain_len);
			//fprintf(stderr, "Res end: %s\n", res_str);
			free(remain_str);
			return;
		}
	}
}

void print_location_list(int verbosity, const bwtint_t& sp, const bwtint_t& ep, const int& len) {
	if (verbosity < 1)	return;
	bwtint_t tmp;
	int dir;
	fprintf(stderr, "===\nLength of match: %d\n", len);
	for (bwtint_t i = sp; i <= ep; i++) {
		tmp = get_pos(&i, len, dir);
		fprintf(stderr, "strand: %d\t", dir);
		char *chr_name;
		int32_t chr_len;
		uint32_t chr_beg;
		uint32_t chr_end;
		bwt_get_intv_info(tmp, tmp + len - 1, &chr_name, &chr_len, &chr_beg, &chr_end);
		fprintf(stderr, "chr%s: ", chr_name);
		fprintf(stderr, "%"PRIu64"-", chr_beg);
		fprintf(stderr, "%"PRIu64"\t", chr_end);
	}
	fprintf(stderr, "\n===\n");
}
