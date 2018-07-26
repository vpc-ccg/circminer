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
	
	if (UB-LB+1 > FRAGLIM) {
		mk->frag_count = 0;
		return 0;
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
	int sum = 0;
	int dir;
	uint32_t match_len = kmer_size;
	int32_t end_pos;

	em_count = 0;
	int invalid_kmer = 0;
	for (i = shift; i < rseq_len; i += skip) {
		if (rseq_len - i < kmer_size)
			match_len = rseq_len - i;
		if (match_len < MINKMER) 
			break;

		if (i != shift)
			mk_res += ll_step;
		occ = get_exact_locs_hash(rseq + i, i, match_len, mk_res);
		
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
	//	valid_ov_kmer = kmer_match_skip_hash(rseq, rseq_len, kmer_size, kmer_size / 2, kmer_size, exact_match_res_ov, ov_em_count);
	
	vafprintf(1, stderr, "Non-OV valids: %d\nOV valids: %d\n", valid_nonov_kmer, valid_ov_kmer);

//	print_hits(starting_node, 7);

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

// pos is exclusive
// [ pos-len, pos-1 ]
void get_reference_chunk_left(uint32_t pos, int len, char* res_str) {
	res_str[0] = 0;
	char* chr_name;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	string chr = getRefGenomeName();
	chr_beg = pos - len;
	chr_end = pos;
	if (chr == "")
		return; 

	//fprintf(stderr, "Going for %lu - %lu\n", pos-len, pos-1);
	int seg_ind = gtf_parser.search_loc(chr, true, chr_end - 1);
	//fprintf(stderr, "Index found: %d\n", seg_ind);
	uint32_t seg_start = gtf_parser.get_start(chr, seg_ind);
	//fprintf(stderr, "Start of exon: %lu\nChr beg: %lu\n", seg_start, chr_beg);
	if (seg_start <= chr_beg) {	
		pac2char(pos - len, len, res_str);
		return;
	}
	else {
		int remain_len = seg_start - chr_beg;
		int covered_len = len - remain_len;
		pac2char(pos - covered_len, covered_len, res_str);
		if (seg_ind == 0) {	// reached first exonic region on chromosome
			return;
		}

		uint32_t prev_end = gtf_parser.get_end(chr, seg_ind - 1);
		//uint32_t prev_integrated_pos = pos - covered_len - (seg_start - prev_end);
		//char* remain_str = (char*) malloc(len+5);
		char remain_str[len+5];
		get_reference_chunk_left(prev_end+1, remain_len, remain_str);
		strncat(remain_str, res_str, covered_len);
		strcpy(res_str, remain_str);
		//free(remain_str);
		return;
	}
}

// pos is exclusive
// [ pos+1, pos+len ]
void get_reference_chunk_right(uint32_t pos, int len, char* res_str) {
	res_str[0] = 0;
	uint32_t chr_beg;
	uint32_t chr_end;
	
	string chr = getRefGenomeName();
	chr_beg = pos;
	chr_end = pos + len;

	if (chr == "")
		return; 

	vafprintf(2, stderr, "Going for %lu - %lu\n", pos+1, pos+len);
	int seg_ind = gtf_parser.search_loc(chr, false, chr_beg);
	//fprintf(stderr, "Index found: %d\n", seg_ind);
	uint32_t seg_end = gtf_parser.get_end(chr, seg_ind);
	//uint32_t seg_start = gtf_parser.get_start(chr, seg_ind);
	//vafprintf(2, stderr, "Start of exon: %lu\nEnd of exon: %lu\nChr end: %lu\n", seg_start, seg_end, chr_end);
	if (seg_end >= chr_end) {	
		pac2char(pos + 1, len, res_str);
		return;
	}
	else {
		int remain_len = chr_end - seg_end;
		int covered_len = len - remain_len;
		pac2char(pos + 1, covered_len, res_str);
		if (gtf_parser.is_last_exonic_region(chr, seg_ind)) {	// reached first exonic region on chromosome
			return;
		}

		uint32_t next_start = gtf_parser.get_start(chr, seg_ind + 1);
		uint32_t next_integrated_pos = pos + covered_len + (next_start - seg_end);
		//vafprintf(2, stderr, "next start: %lu, %lu\n", next_start, next_integrated_pos);
		char remain_str[len+5];
		get_reference_chunk_right(next_start-1, remain_len, remain_str);
		strncat(res_str, remain_str, remain_len);
		return;
	}
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
