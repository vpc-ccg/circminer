#include <cmath>
#include <cassert>
#include "match_read.h"
#include "gene_annotation.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

void print_hits(GIMatchedKmer*, int);

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

int get_exact_locs_hash(char* seq, int32_t qpos, uint32_t len, GIMatchedKmer* mk) {
	mk->frag_count = 0;
	mk->frags = NULL;
	mk->qpos = qpos;

	GeneralIndex *it = getCandidates(hashVal(seq));
	uint32_t i, j;

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
	
	if (UB - LB + 1 > FRAGLIM) {
		mk->frag_count = 0;
		return UB - LB + 1;
	}

	mk->frag_count = UB - LB + 1;
	mk->frags = it + LB;

	return UB - LB + 1;
}

// return:
// # valid kmers
// is valid if: #fragments > 0 and < FRAGLIM
int kmer_match_skip_hash(char* rseq, int rseq_len, int kmer_size, int shift, int skip, int ll_step, GIMatchedKmer* mk_res, int& em_count) {
	int i, occ;
	int dir;
	int sum = 0;
	int invalid_kmer = 0;
	uint32_t match_len = kmer_size;
	int32_t end_pos;

	GIMatchedKmer* cur = mk_res;

	em_count = 0;
	for (i = shift; i < rseq_len; i += skip) {
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

		if (occ > FRAGLIM ) {
			invalid_kmer++;
		}

		em_count++;	
	}

	return em_count - invalid_kmer;
}

// fragment results will be sorted
int split_match_hash(char* rseq, int rseq_len, int kmer_size, GIMatchedKmer* starting_node) {
	int valid_nonov_kmer = 0;
	int valid_ov_kmer = 0;
	int nonov_em_count;
	int ov_em_count;

	valid_nonov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, 0, kmer_size, 2, starting_node, nonov_em_count);

	//if (valid_nonov_kmer < 2)
	//	valid_ov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, kmer_size / 2, kmer_size, 2, starting_node + 1, ov_em_count);
	
	vafprintf(1, stderr, "Non-OV valids: %d\nOV valids: %d\n", valid_nonov_kmer, valid_ov_kmer);

	//print_hits(starting_node, 7);

	return valid_nonov_kmer + valid_ov_kmer;
}

void pac2char(uint32_t start, int len, char* str) {
	//fprintf(stderr, "extract pos: %d-%d\n", start, start+len-1);
	int ref_len = getRefGenLength();
	if (start + len - 1 > ref_len)
		return;		// do not write ref seq if it goes out of the contig border

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
	//str[i] = '\0';
}

void print_hits(GIMatchedKmer* frag_l, int cnt) {
	for (int i = 0; i < cnt; i+=2) {
		fprintf(stderr, "Frag cnt: %lu\n", (frag_l+i)->frag_count);
		for (int j = 0; j < (frag_l+i)->frag_count; j++) {
			fprintf(stderr, "%lu\t", (frag_l+i)->frags[j].info);
		}
		fprintf(stderr, "\n");
	}
}
