#include <cmath>
#include "match_read.h"
#include "gene_annotation.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

void print_hits(GIMatchedKmer*, int);

// assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
int frag_binary_search(const vector<fragment_t>& list, int beg, int end, uint32_t target) {
	if (end - beg <= 1)
		return end;
	int mid = (beg + end) / 2;
	if (target < list[mid].rpos)
		return frag_binary_search(list, beg, mid, target);
	else 
		return frag_binary_search(list, mid, end, target);
}

int get_exact_locs_hash(char* seq, int32_t qpos, uint32_t len, GIMatchedKmer* mk) {
	mk->frag_count = 0;
	mk->frags = NULL;
	mk->qpos = qpos;

	int hv = hashVal(seq);
	if (hv < 0)
		return 0;

	char* checksum_beg = seq + WINDOW_SIZE;
	int cv = checkSumVal(checksum_beg);
	if (cv < 0)
		return 0;

	GeneralIndex *it = getCandidates(hv);

	if (it == NULL) {
		return 0;
	}

	uint32_t lb = 1;
	uint32_t ub = it[0].info;
	uint32_t mid;
	int16_t target = cv;
	uint32_t LB = 0;
	uint32_t UB = 0;

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
	
	mk->frag_count = UB - LB + 1;
	mk->frags = it + LB;

	return UB - LB + 1;
}

// reduce number of frags in LargeList (ll) using SmalList (sl)
// ll locations should be smaller than sl
bool reduce_hits_behind(GIMatchedKmer* sl, GIMatchedKmer* ll) {
	if (sl->frag_count == 0 or ll->frag_count == 0)
		return false;

	uint32_t max_dist = ((sl->qpos - ll->qpos) / kmer) * maxIntronLen;

	uint32_t size = 0;
	uint32_t j = 0;
	for (uint32_t i = 0; i < ll->frag_count; i++) {
		while (j < sl->frag_count and sl->frags[j].info <= ll->frags[i].info)
			j++;

		if (j >= sl->frag_count)
			break;

		// now sl->frags[j].info > ll->frags[i].info
		if (ll->frags[i].info + max_dist >= sl->frags[j].info) {
			ll->frags[size].info = ll->frags[i].info;
			size++;
		}
	}

	//fprintf(stderr, "Pre hits size: %d\tCur hits size: %d\n", ll->frag_count, size);

	ll->frag_count = size;
	return (size > 0 and size <= seedLim);
}

// reduce number of frags in LargeList (ll) using SmalList (sl)
// ll locations should be greater than sl
bool reduce_hits_ahead(GIMatchedKmer* sl, GIMatchedKmer* ll) {
	if (sl->frag_count == 0 or ll->frag_count == 0)
		return false;

	uint32_t max_dist = ((ll->qpos - sl->qpos) / kmer) * maxIntronLen;

	uint32_t size = 0;
	uint32_t j = 0;
	for (uint32_t i = 0; i < ll->frag_count; i++) {
		if (j >= sl->frag_count)
			break;

		if (sl->frags[j].info >= ll->frags[i].info)
			continue;

		while (j < sl->frag_count and sl->frags[j].info + max_dist < ll->frags[i].info)
			j++;

		if (j >= sl->frag_count)
			break;

		if (ll->frags[i].info <= max_dist + sl->frags[j].info) {
			ll->frags[size].info = ll->frags[i].info;
			size++;
		}
	}

	//fprintf(stderr, "Pre hits size: %d\tCur hits size: %d\n", ll->frag_count, size);

	ll->frag_count = size;
	return (size > 0 and size <= seedLim);
}

// return:
// # valid kmers
// is valid if: #fragments > 0 and < seedLim
int kmer_match_skip_hash(char* rseq, int rseq_len, int kmer_size, int shift, int skip, int ll_step, GIMatchedKmer* mk_res, int& em_count) {
	int i, j;
	int occ;
	int invalid_kmer = 0;
	int valid_kmer = 0;
	//int valid_kmer_ind = 0;
	uint32_t match_len = kmer_size;

	GIMatchedKmer* cur = mk_res;

	// initialize
	int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer)) - 1;
	for (j = 0; j < max_seg_cnt; j++) {
		(cur+j)->frag_count = 0;
		(cur+j)->frags = NULL;
		(cur+j)->qpos = skip*j;
	}

	em_count = 0;
	j = -1 * ll_step;
	for (i = shift; i < rseq_len; i += skip) {
		j += ll_step;

		if (rseq_len - i < kmer_size)
			match_len = rseq_len - i;
		if (match_len < MINKMER) 
			break;

		if (i != shift)
			cur += ll_step;
		occ = get_exact_locs_hash(rseq + i, i, match_len, cur);
		
		vafprintf(2, stderr, "Occ: %d\tind: %d\tmatch len: %d\n", occ, i, match_len);
		if (occ <= 0) {
			occ = 0;
			invalid_kmer++;
		}

		else if (uint32_t(occ) > seedLim ) {
			invalid_kmer++;
		}

		else {
			//valid_kmer_ind = j;
			valid_kmer++;
		}

		em_count++;	
	}

	valid_kmer = em_count - invalid_kmer;

	// trying to reduce size of hits it there is only one valid kmer
	//if (valid_kmer == 1) {
	//	GIMatchedKmer* small = mk_res + valid_kmer_ind;
	//	for (int i = 0; i < valid_kmer_ind; i += ll_step) {
	//		//fprintf(stderr, "Reducing size of ind = %d\n", i);
	//		if (reduce_hits_behind(small, mk_res + i))
	//			valid_kmer++;
	//	}

	//	for (int i = valid_kmer_ind + ll_step; i < em_count * ll_step; i += ll_step) {
	//		//fprintf(stderr, "Reducing size of ind = %d\n", i);
	//		if (reduce_hits_ahead(small, mk_res + i))
	//			valid_kmer++;
	//	}
	//}
	
	for (int i = 0; i < em_count * ll_step; i += ll_step) {
		if ((mk_res+i)->frag_count > seedLim) {
			(mk_res+i)->frag_count = 0;
			//(mk_res+i)->frags = NULL;
		}
	}

	return valid_kmer;
}

// fragment results will be sorted
int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* starting_node) {
	int valid_nonov_kmer = 0;
	int valid_ov_kmer = 0;
	int nonov_em_count;

	valid_nonov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, 0, kmer_size, 2, starting_node, nonov_em_count);

	//int ov_em_count;
	//if (valid_nonov_kmer < 2)
	//	valid_ov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, kmer_size / 2, kmer_size, 2, starting_node + 1, ov_em_count);
	
	vafprintf(1, stderr, "Non-OV valids: %d\nOV valids: %d\n", valid_nonov_kmer, valid_ov_kmer);

	//print_hits(starting_node, 7);

	return valid_nonov_kmer + valid_ov_kmer;
}

bool pac2char(uint32_t start, int len, char* str) {
	//fprintf(stderr, "extract pos: %d-%d\n", start, start+len-1);
	int ref_len = getRefGenLength();
	if (((int) (start) < 0) or ((int)(start) + len - 1 > ref_len))
		return false;		// do not write ref seq if it goes out of the contig border

	char lookup[8] = { 'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N' };
	
	CompressedSeq* crnext = getCmpRefGenome();
	uint32_t skip = (start - 1) / 21;
	crnext += skip;
	int pass = (start - 1) % 21;

	CompressedSeq crdata = *crnext;
	crdata <<= 3 * pass;

	int i = 0;
	int j = pass;
	int val;

	while (i < len) {
		val = (crdata >> 60) & 7;
		str[i++] = lookup[val];
		if (++j == 21) {
			j = 0;
			crdata = *(++crnext);
		}
		else {
			crdata <<= 3;
		}
	}
	str[i] = '\0';
	return true;
}

//void print_hits(GIMatchedKmer* frag_l, int cnt) {
//	for (int i = 0; i < cnt; i+=2) {
//		fprintf(stderr, "Frag cnt: %"PRIu32"\n", (frag_l+i)->frag_count);
//		for (uint32_t j = 0; j < (frag_l+i)->frag_count; j++) {
//			fprintf(stderr, "%"PRIu32"\t", (frag_l+i)->frags[j].info);
//		}
//		fprintf(stderr, "\n");
//	}
//}
