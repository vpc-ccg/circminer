#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <inttypes.h>

#include <string>
#include <vector>
#include <algorithm>

#include "bwt.h"
#include "fastq_parser.h"

using namespace std;

#define GENETHRESH 100000
#define RANGELIM 1000

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

// input interval: [, )
// return: i if target in [i-1, i)
bwtint_t binary_search(vector<bwtint_t>& list, bwtint_t beg, bwtint_t end, bwtint_t& target) {
	if (end - beg <= 1)
		return end;
	bwtint_t mid = (beg + end) / 2;
	if (target < list[mid])
		return binary_search(list, beg, mid, target);
	else 
		return binary_search(list, mid, end, target);
}

bool is_concordant_sorted(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh) {
	if ((ep_b - sp_b + 1) * (ep_f - sp_f + 1) >= RANGELIM)
		return false;

	//fprintf(stderr, "back size: %llu\nfront size: %llu\n", ep_b - sp_b + 1, ep_f - sp_f + 1);

	int dir_front=0, dir_back=0;
	bwtint_t flist_size_f, blist_size_f;
	bwtint_t flist_size_b, blist_size_b;

	bwtint_t spos_front, spos_back;
	vector<bwtint_t> forwardlist_f(ep_f - sp_f + 2);
	vector<bwtint_t> backwardlist_f(ep_f - sp_f + 2);
	vector<bwtint_t> forwardlist_b(ep_b - sp_b + 2);
	vector<bwtint_t> backwardlist_b(ep_b - sp_b + 2);

	sort_positions(sp_f, ep_f, len_f, forwardlist_f, flist_size_f, backwardlist_f, blist_size_f);
	sort_positions(sp_b, ep_b, len_b, forwardlist_b, flist_size_b, backwardlist_b, blist_size_b);
	
	//fprintf(stderr, "forward front size: %llu\nbackward front size: %llu\n", flist_size_f, blist_size_f);
	//fprintf(stderr, "forward back size: %llu\nbackward back size: %llu\n", flist_size_b, blist_size_b);

	//for (int k = 0; k < flist_size_f; k++)
	//	fprintf(stderr, "%llu\t", forwardlist_f[k]);
	//fprintf(stderr, "\n");
	//for (int k = 0; k < blist_size_f; k++)
	//	fprintf(stderr, "%llu\t", backwardlist_f[k]);
	//fprintf(stderr, "\n");
	//for (int k = 0; k < flist_size_b; k++)
	//	fprintf(stderr, "%llu\t", forwardlist_b[k]);
	//fprintf(stderr, "\n");
	//for (int k = 0; k < blist_size_b; k++)
	//	fprintf(stderr, "%llu\t", backwardlist_b[k]);
	//fprintf(stderr, "\n");

	bwtint_t i, j;
	if (flist_size_f > 0 and blist_size_b > 0) {
		if (flist_size_f <= blist_size_b) {
			for (i = 0; i < flist_size_f; i++) {
				if (forwardlist_f[i] < backwardlist_b[0])	// no position less than forwardlist_f[i] available on backwardlist_b
					continue;
				else {
					j = binary_search(backwardlist_b, 0, blist_size_b, forwardlist_f[i]);
				}
				spos_front = forwardlist_f[i];
				spos_back = backwardlist_b[j-1];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				//fprintf(stderr, "HEYYYYY %d\n%llu\n%llu\n", len_f, spos_front, spos_back + len_b + noise_thresh);
				if (spos_front <= spos_back + len_b + noise_thresh)
					return true;
			}
		}
		else {
			for (j = 0; j < blist_size_b; j++) {
				if (backwardlist_b[j] < forwardlist_f[0])
					i = 0;
				else {
					i = binary_search(forwardlist_f, 0, flist_size_f, backwardlist_b[j]);
					if (i >= flist_size_f)	// no position greater than backwardlist_b[j] available on forwardlist_f
						continue;
				}
				spos_front = forwardlist_f[i];
				spos_back = backwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_front <= spos_back + len_b + noise_thresh)
					return true;
			}
		}
	}

	if (blist_size_f > 0 and flist_size_b > 0) {
		if (blist_size_f <= flist_size_b) {
			for (i = 0; i < blist_size_f; i++) {
				if (backwardlist_f[i] < forwardlist_b[0])
					j = 0;
				else {
					j = binary_search(forwardlist_b, 0, flist_size_b, backwardlist_f[i]);
					if (j >= flist_size_b)	// no position greater than backwardlist_f[i] available on backwardlist_f
						continue;
				}
				spos_front = backwardlist_f[i];
				spos_back = forwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_back <= spos_front + len_f + noise_thresh)
					return true;
			}
		}
		else {
			for (j = 0; j < flist_size_b; j++) {
				if (forwardlist_b[j] < backwardlist_f[0])	// no position less than forwardlist_b[j] available on backwardlist_f
					continue;
				else {
					i = binary_search(backwardlist_f, 0, blist_size_f, forwardlist_b[j]);
				}
				spos_front = backwardlist_f[i-1];
				spos_back = forwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_back <= spos_front + len_f + noise_thresh)
					return true;
			}
		}
	}
	return false;
}

bool is_concordant(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh) {
	//fprintf(stderr, "front: %llu\t%llu\n", sp_f, ep_f);
	//fprintf(stderr, "back : %llu\t%llu\n", sp_b, ep_b);
	
	if ((ep_b - sp_b + 1) * (ep_f - sp_f + 1) >= RANGELIM)
		return false;

	//fprintf(stderr, "back size: %llu\nfront size: %llu\n", ep_b - sp_b + 1, ep_f - sp_f + 1);

	int dir_front, dir_back;
	bwtint_t spos_front, spos_back;
	for (bwtint_t i = sp_f; i <= ep_f; i++) {
		for (bwtint_t j = sp_b; j <= ep_b; j++) {
			spos_front = get_pos(&i, len_f, dir_front);
			spos_back = get_pos(&j, len_b, dir_back);
			//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
			
			//char *chr_name;
			//int32_t chr_len;
			//uint32_t chr_beg;
			//uint32_t chr_end;
			//bwt_get_intv_info(spos_front, spos_front + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
			//fprintf(stderr, "Chr: %s\n", chr_name);
			//fprintf(stderr, "%"PRIu64"\n", chr_beg);
			//bwt_get_intv_info(spos_back, spos_back + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
			//fprintf(stderr, "%"PRIu64"\n", chr_beg);

			if (dir_back == 1) {
				if ((dir_front == -1) and (spos_back - spos_front - len_f <= noise_thresh))
					return true;
			}
			else {
				if ((dir_front == 1) and (spos_front - spos_back - len_b <= noise_thresh))
					return true;
			}
		}
	}
	return false;
}

// --return true if potentially chimeric
// 0 : concordant
// 1 : potentially chimeric (keep)
int find_expanded_positions(const char* rseq, const char* rcseq, const int& rseq_len) {
	bwtint_t sp_b, ep_b;
	bwtint_t sp_f, ep_f;
	int exp_len_back = 0;
	int exp_len_front = 0;

	const int noise_thresh = rseq_len / 10;

	exp_len_front = get_expanded_locs(rcseq, rseq_len, &sp_f, &ep_f);
	//fprintf(stderr, "Expanded len from end: %d\n", exp_len_forward);
	if (exp_len_front >= rseq_len - noise_thresh) {	// concordant
		//fprintf(stderr, "[Concordant-e]\tfront matched: %d", exp_len_front);
		return 0;
	}

	exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front, &sp_b, &ep_b);
	//exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front - noise_thresh, &sp_b, &ep_b);
	//fprintf(stderr, "%s\tlen front: %d,\tlen back: %d\n", rseq, exp_len_front, exp_len_back);
	//if (occ > 0)
	//	locate_match(&sp, &ep, window_size);
	
	if (exp_len_front + exp_len_back <= rseq_len - noise_thresh * 2) {
		//fprintf(stderr, "[Chimerici-e]\tfront matched: %d, back matched: %d", exp_len_front, exp_len_back);
		return 1;
	}
	
	//bool consis = check_cosistent_match(&sp_b, &ep_b, exp_len_backward, &sp_f, &ep_f, exp_len_forward);
	bool consis = is_concordant(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh);
	//bool consis = is_concordant_sorted(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh);

	//std::string res = (consis) ? "Concordant" : "Chimeric";
	//fprintf(stderr, "[%s]\tfront matched: %d, back matched: %d", res.c_str(), exp_len_front, exp_len_back);

	if (consis)
		return 0;
	else
		return 1;
}

int find_exact_positions(const char* rseq, int rseq_len, int window_size) {
	bwtint_t sp, ep;
	int occ = 0;

	for (int i = 0; i < rseq_len; i+=window_size) {
		occ = get_exact_locs(rseq + i, window_size, &sp, &ep);
		//fprintf(stderr, "Start pos: %d\n", i);
		//fprintf(stderr, "Number of matches: %d\n", occ);

		if (occ > 0)
			locate_match(&sp, &ep, window_size);
	}

	return occ;
}

int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <FASTA> <FASTQ> <OUT>\n", argv[0]);
		return 1;
	}
	char* ref_file = argv[1];
	char* fq_file1 = argv[2];
	
	char ignore_file[1000], keep_file[1000];
	//char prob_ignore_file[1000];
	strcpy (ignore_file, argv[3]);
	//strcpy (prob_ignore_file, argv[3]);
	strcpy (keep_file, argv[3]);

	strcat (ignore_file, ".ignore.fastq");
	//strcat (prob_ignore_file, ".prob.ignore.fastq");
	strcat (keep_file, ".keep.fastq");

	/*
	if (bwt_index(ref_file)) {
		fprintf(stderr, "Nope! Index failed!\n");
		return 1;
	}
	else
		fprintf(stderr, "Indexed successfully!\n");
	*/

	if (bwt_load(ref_file)) {
		fprintf(stderr, "Index not loaded!\n");
		return 1;
	}
	else
		fprintf(stderr, "Index file successfully loaded!\n");

	//FILE* fqin1 = fopen(fq_file1, "r");
	FILE* ignore_out = fopen(ignore_file, "w");
	//FILE* prob_ignore_out = fopen(prob_ignore_file, "w");
	FILE* keep_out = fopen(keep_file, "w");

	FASTQParser fq_parser(fq_file1);
	Record* current_record;
	fprintf(stderr, "Started reading FASTQ file\n");
	int line = 0;

	while ( fq_parser.has_next() ) { // go line by line on fastq file
		current_record = fq_parser.get_next();
		if (current_record == NULL)	// no new line
			break;

		//fprintf(stderr, "line: %d\n", line);
		line++;

		if (line % 1000000 == 0)
			fprintf(stderr, "[P] %d lines\n", line);

		//fprintf(stderr, "Seq: %s, len: %d\n", current_record->seq, current_record->seq_len);
		//fprintf(stderr, "RC Seq: %s, len: %d\n", current_record->rcseq, current_record->seq_len);
		int find_len = current_record->seq_len;
		//int occ = find_exact_positions(current_record->seq, current_record->seq_len, find_len);
		int is_chimeric = find_expanded_positions(current_record->seq, current_record->rcseq, current_record->seq_len);

		//fprintf(stderr, "\t%s", current_record->rname);

		//if (occ >= 1) {
		if (is_chimeric == 0) {
			fprintf(ignore_out, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
		}
		//else if (is_chimeric == 1) {
		//	fprintf(prob_ignore_out, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
		//}
		else {
			fprintf(keep_out, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
		}
	}

	return 0;
}
