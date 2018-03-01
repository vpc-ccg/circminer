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

// input interval: [, )
// return: i if target in [i-1, i)
bwtint_t binary_search(const vector<bwtint_t>& list, bwtint_t beg, bwtint_t end, const bwtint_t& target) {
	if (end - beg <= 1)
		return end;
	bwtint_t mid = (beg + end) / 2;
	if (target < list[mid])
		return binary_search(list, beg, mid, target);
	else 
		return binary_search(list, mid, end, target);
}

bool is_concordant_sorted(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh) {
	//if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
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
				if (spos_front <= spos_back + len_b + GENETHRESH)
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
				
				if (spos_front <= spos_back + len_b + GENETHRESH)
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
				
				if (spos_back <= spos_front + len_f + GENETHRESH)
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
				
				if (spos_back <= spos_front + len_f + GENETHRESH)
					return true;
			}
		}
	}
	return false;
}

bool is_concordant_sorted2(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh, MatchedRead& mr) {
	//if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
	if ((ep_b - sp_b + 1) * (ep_f - sp_f + 1) >= REGIONSIZELIM) {
		return false;
	}

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
				if (spos_front <= spos_back + len_b + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_back;
					mr.matched_len = spos_front + len_f - spos_back + 1;
					mr.dir = -1;
					return true;
				}
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
				
				if (spos_front <= spos_back + len_b + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_back;
					mr.matched_len = spos_front + len_f - spos_back + 1;
					mr.dir = -1;
					return true;
				}
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
				
				if (spos_back <= spos_front + len_f + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_front;
					mr.matched_len = spos_back + len_b - spos_front + 1;
					mr.dir = 1;
					return true;
				}
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
				
				if (spos_back <= spos_front + len_f + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_front;
					mr.matched_len = spos_back + len_b - spos_front + 1;
					mr.dir = 1;
					return true;
				}
			}
		}
	}
	return false;
}

bool is_concordant_sorted_same_strand(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, const int& noise_thresh, MatchedRead& mr) {
	//if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
	//if ((ep_b - sp_b + 1) * (ep_f - sp_f + 1) >= REGIONSIZELIM) {
	if ((ep_b - sp_b + 1) * (ep_f - sp_f + 1) >= RANGELIM) {
		return false;
	}

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
	if (flist_size_f > 0 and flist_size_b > 0) {
		if (flist_size_f <= flist_size_b) {
			for (i = 0; i < flist_size_f; i++) {
				if (forwardlist_f[i] < forwardlist_b[0])	// no position less than forwardlist_f[i] available on backwardlist_b
					j = 0;
				else {
					j = binary_search(forwardlist_b, 0, flist_size_b, forwardlist_f[i]);
					if (j >= flist_size_b)
						continue;
				}
				spos_front = forwardlist_f[i];
				spos_back = forwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				//fprintf(stderr, "HEYYYYY %d\n%llu\n%llu\n", len_f, spos_front, spos_back + len_b + noise_thresh);
				if (spos_back <= spos_front + len_f + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_front;
					mr.matched_len = spos_back + len_b - spos_front + 1;
					mr.dir = 1;
					return true;
				}
			}
		}
		else {
			for (j = 0; j < flist_size_b; j++) {
				if (forwardlist_b[j] < forwardlist_f[0])
					continue;
				else {
					i = binary_search(forwardlist_f, 0, flist_size_f, forwardlist_b[j]);
				}
				spos_front = forwardlist_f[i-1];
				spos_back = forwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_back <= spos_front + len_f + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_front;
					mr.matched_len = spos_back + len_b - spos_front + 1;
					mr.dir = 1;
					return true;
				}
			}
		}
	}

	if (blist_size_f > 0 and blist_size_b > 0) {
		if (blist_size_f <= blist_size_b) {
			for (i = 0; i < blist_size_f; i++) {
				if (backwardlist_f[i] < backwardlist_b[0])
					continue;
				else {
					j = binary_search(backwardlist_b, 0, blist_size_b, backwardlist_f[i]);
				}
				spos_front = backwardlist_f[i];
				spos_back = backwardlist_b[j-1];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_front <= spos_back + len_b + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_back;
					mr.matched_len = spos_front + len_f - spos_back + 1;
					mr.dir = -1;
					return true;
				}
			}
		}
		else {
			for (j = 0; j < blist_size_b; j++) {
				if (backwardlist_b[j] < backwardlist_f[0])	// no position less than forwardlist_b[j] available on backwardlist_f
					i = 0;
				else {
					i = binary_search(backwardlist_f, 0, blist_size_f, backwardlist_b[j]);
					if (i >= blist_size_f)	// no position greater than backwardlist_f[i] available on backwardlist_f
						continue;
				}
				spos_front = backwardlist_f[i];
				spos_back = backwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_front <= spos_back + len_b + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_back;
					mr.matched_len = spos_front + len_f - spos_back + 1;
					mr.dir = -1;
					return true;
				}
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
				if ((dir_front == -1) and (spos_back - spos_front - len_f <= GENETHRESH))
					return true;
			}
			else {
				if ((dir_front == 1) and (spos_front - spos_back - len_b <= GENETHRESH))
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
	if (exp_len_front >= rseq_len - noise_thresh) {	// concordant
		fprintf(stderr, "[Concordant-e]\tfront matched: %d", exp_len_front);
		return 0;
	}

	exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front, &sp_b, &ep_b);
	//exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front - noise_thresh, &sp_b, &ep_b);
	fprintf(stderr, "%s\tlen front: %d,\tlen back: %d\n", rseq, exp_len_front, exp_len_back);
	//if (occ > 0)
	//	locate_match(&sp, &ep, window_size);
	
	if (exp_len_front + exp_len_back <= rseq_len / 2) {
		fprintf(stderr, "[Chimeric-e]\tfront matched: %d, back matched: %d", exp_len_front, exp_len_back);
		return 1;
	}
	
	//bool consis = check_cosistent_match(&sp_b, &ep_b, exp_len_backward, &sp_f, &ep_f, exp_len_forward);
	bool consis = is_concordant(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh);
	//bool consis = is_concordant_sorted(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh);

	std::string res = (consis) ? "Concordant" : "Chimeric";
	fprintf(stderr, "[%s]\tfront matched: %d, back matched: %d", res.c_str(), exp_len_front, exp_len_back);

	if (consis)
		return 0;
	else
		return 1;
}

// --return true if potentially chimeric
// 0 : concordant
// 1 : potentially chimeric (keep)
int find_expanded_positions2(const char* rseq, const char* rcseq, const int& rseq_len, MatchedRead& mr) {
	bwtint_t sp_b, ep_b;
	bwtint_t sp_f, ep_f;
	int exp_len_back = 0;
	int exp_len_front = 0;

	const int noise_thresh = rseq_len / 10;

	mr.is_concord = false;

	exp_len_front = get_expanded_locs(rcseq, rseq_len, &sp_f, &ep_f);
	if (exp_len_front >= rseq_len - noise_thresh) {	// concordant
		fprintf(stderr, "[Concordant-e]\tfront matched: %d\n", exp_len_front);
		mr.is_concord = true;
		mr.start_pos = get_pos(&sp_f, exp_len_front, mr.dir);
		mr.dir *= -1;
		mr.matched_len = exp_len_front;
		return 0;
	}

	exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front, &sp_b, &ep_b);
	fprintf(stderr, "%s\tlen front: %d,\tlen back: %d\n", rseq, exp_len_front, exp_len_back);
	
	//if (exp_len_front < 23 or  exp_len_back < 23) {
	if (exp_len_front + exp_len_back <= 0.6 * rseq_len) {
		fprintf(stderr, "[Chimeric-e]\tfront matched: %d, back matched: %d\n", exp_len_front, exp_len_back);
		return 1;
	}
	
	//bool consis = check_cosistent_match(&sp_b, &ep_b, exp_len_backward, &sp_f, &ep_f, exp_len_forward);
	//bool consis = is_concordant(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh);
	bool consis = is_concordant_sorted2(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh, mr);

	std::string res = (consis) ? "Concordant" : "Chimeric";
	fprintf(stderr, "[%s]\tfront matched: %d, back matched: %d\n", res.c_str(), exp_len_front, exp_len_back);

	if (consis)
		return 0;
	else
		return 1;
}

int check_concordant_mates(const Record* m1, const Record* m2) {
	MatchedRead mr1, mr2;
	int chimeric1 = find_expanded_positions2(m1->seq, m1->rcseq, m1->seq_len, mr1);
	int chimeric2 = find_expanded_positions2(m2->seq, m2->rcseq, m2->seq_len, mr2);

	if (chimeric1 == 1 and chimeric2 == 1)
		return 1;

	if (chimeric1 == 1 or chimeric2 == 1) {
		//fprintf(stderr, "Half concordant: %s", m1->rname);
		return 2;
	}

	//fprintf(stderr, "%s", m1->rname);
	//fprintf(stderr, "M1: dir:%d,\tstart:%llu,\tend:%llu\n", mr1.dir, mr1.start_pos, mr1.start_pos + mr1.matched_len);
	//fprintf(stderr, "M2: dir:%d,\tstart:%llu,\tend:%llu\n", mr2.dir, mr2.start_pos, mr2.start_pos + mr2.matched_len);

	if (mr1.dir == 1 and mr2.dir == -1 and (mr1.start_pos <= mr2.start_pos))
		return 0;
	if (mr1.dir == -1 and mr2.dir == 1 and (mr2.start_pos <= mr1.start_pos))
		return 0;
	return 3;
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

int find_exact_positions_slide(const char* rseq, int rseq_len, const int& window_size, const int& shift_step, MatchedRead& mr) {
	bwtint_t sp_f, ep_f;
	bwtint_t sp_b, ep_b;
	int occ = 0, i, j;
	MatchedRead mr1, mr2;

	mr1.is_concord = false;
	for (i = 0; i <= rseq_len - 2*window_size; i += shift_step) {
		occ = get_exact_locs(rseq + i, window_size, &sp_f, &ep_f);
		fprintf(stderr, "Start pos: %d\n", i);
		fprintf(stderr, "Number of matches: %d\n", occ);

		if (occ > 0) {
			mr1.is_concord = true;
			mr1.matched_len = window_size;
			break;
		}
	}

	mr2.is_concord = false;
	for (j = rseq_len - window_size; j >= i + window_size; j -= shift_step) {
		occ = get_exact_locs(rseq + j, window_size, &sp_b, &ep_b);
		fprintf(stderr, "Start pos: %d\n", j);
		fprintf(stderr, "Number of matches: %d\n", occ);

		if (occ > 0) {
			mr2.is_concord = true;
			mr2.matched_len = window_size;
			break;
		}
	}

	mr.is_concord = mr1.is_concord && mr2.is_concord;
	if (!mr.is_concord) {
		return 1;
	}
	
	bool consis = is_concordant_sorted_same_strand(sp_f, ep_f, mr1.matched_len, sp_b, ep_b, mr2.matched_len, rseq_len/10, mr);

	std::string res = (consis) ? "Concordant" : "Chimeric";
	fprintf(stderr, "[%s]\tfront matched: %d, back matched: %d\n", res.c_str(), mr1.matched_len, mr2.matched_len);
	fprintf(stderr, "front start rpos: %d, back start rpos: %d\n", i, j);

	if (consis)
		return 0;
	else 
		return 1;
}

// return value:
// 0: concordant
// 1: potentially chimeric
// 2: definitely chimeric
int check_concordant_mates_noexpand(const Record* m1, const Record* m2) {
	MatchedRead mr1, mr2;
	int chimeric1 = find_exact_positions_slide(m1->seq, m1->seq_len, 28, 3, mr1);
	int chimeric2 = find_exact_positions_slide(m2->seq, m2->seq_len, 28, 3, mr2);

	//fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

	if (chimeric1 == 1 and chimeric2 == 1)
		return 1;

	if (chimeric1 == 1 or chimeric2 == 1) {
		//fprintf(stderr, "Half concordant: %s", m1->rname);
		return 1;
	}

	//fprintf(stderr, "%s", m1->rname);
	//fprintf(stderr, "M1: dir:%d,\tstart:%llu,\tend:%llu\n", mr1.dir, mr1.start_pos, mr1.start_pos + mr1.matched_len);
	//fprintf(stderr, "M2: dir:%d,\tstart:%llu,\tend:%llu\n", mr2.dir, mr2.start_pos, mr2.start_pos + mr2.matched_len);

	if (mr1.dir == 1 and mr2.dir == -1 and (mr1.start_pos <= mr2.start_pos))
		return 0;
	if (mr1.dir == -1 and mr2.dir == 1 and (mr2.start_pos <= mr1.start_pos))
		return 0;

	if (mr1.dir == 1 and mr2.dir == -1 and (mr1.start_pos > mr2.start_pos + 8) and (mr1.start_pos <= mr2.start_pos + 1000)) {
		fprintf(stderr, "%s", m1->rname);
		fprintf(stderr, "M1: dir: %d\t%llu\n", mr1.dir, mr1.start_pos);
		fprintf(stderr, "M2: dir: %d\t%llu\n", mr2.dir, mr2.start_pos);

		char *chr_name;
		int32_t chr_len;
		uint32_t chr_beg;
		uint32_t chr_end;
		bwt_get_intv_info(mr1.start_pos, mr1.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		fprintf(stderr, "Chr: %s\t", chr_name);
		fprintf(stderr, "%"PRIu64"\n", chr_beg);
		bwt_get_intv_info(mr2.start_pos, mr2.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		fprintf(stderr, "Chr: %s\t", chr_name);
		fprintf(stderr, "%"PRIu64"\n", chr_beg);

		return 2;
	}
	if (mr1.dir == -1 and mr2.dir == 1 and (mr2.start_pos > mr1.start_pos + 8) and (mr2.start_pos <= mr1.start_pos + 1000)) {
		fprintf(stderr, "%s", m1->rname);
		fprintf(stderr, "M1: dir: %d\t%llu\n", mr1.dir, mr1.start_pos);
		fprintf(stderr, "M2: dir: %d\t%llu\n", mr2.dir, mr2.start_pos);

		char *chr_name;
		int32_t chr_len;
		uint32_t chr_beg;
		uint32_t chr_end;
		bwt_get_intv_info(mr1.start_pos, mr1.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		fprintf(stderr, "Chr: %s\t", chr_name);
		fprintf(stderr, "%"PRIu64"\n", chr_beg);
		bwt_get_intv_info(mr2.start_pos, mr2.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		fprintf(stderr, "Chr: %s\t", chr_name);
		fprintf(stderr, "%"PRIu64"\n", chr_beg);

		return 2;
	}
	return 1;
}

bool is_chimeric_intersect(const vector<bwtint_t>& forwardlist_f, const bwtint_t& flist_size_f, const vector<bwtint_t>& backwardlist_f, const bwtint_t& blist_size_f, const int& len_f, 
							const vector<bwtint_t>& forwardlist_b, const bwtint_t& flist_size_b, const vector<bwtint_t>& backwardlist_b, const bwtint_t& blist_size_b, const int& len_b) {
	//if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
	if (flist_size_f >= REGIONSIZELIM or blist_size_f >= REGIONSIZELIM or flist_size_b >= REGIONSIZELIM or blist_size_b >= REGIONSIZELIM) {
		return false;
	}

	//fprintf(stderr, "back size: %llu\nfront size: %llu\n", ep_b - sp_b + 1, ep_f - sp_f + 1);

	int dir_front=0, dir_back=0;
	bwtint_t spos_front, spos_back;
	
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

	bwtint_t i, j, target;
	if (flist_size_f > 0 and blist_size_b > 0) {
		if (flist_size_f <= blist_size_b) {
			for (i = 0; i < flist_size_f; i++) {
				if (forwardlist_f[i] + len_f < backwardlist_b[0])
					j = 0;
				else {
					target = forwardlist_f[i] + len_f;
					j = binary_search(backwardlist_b, 0, blist_size_b, target);
					if (j >= blist_size_b)
						continue;
				}
				spos_front = forwardlist_f[i];
				spos_back = backwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				//fprintf(stderr, "HEYYYYY %d\n%llu\n", len_f, spos_front);
				if (spos_back <= spos_front + GENETHRESH) {
					return true;
				}
			}
		}
		else {
			for (j = 0; j < blist_size_b; j++) {
				if (backwardlist_b[j] < forwardlist_f[0])
					continue;
				else {
					i = binary_search(forwardlist_f, 0, flist_size_f, backwardlist_b[j]);
				}
				spos_front = forwardlist_f[i-1];
				spos_back = backwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_back <= spos_front + GENETHRESH) {
					return true;
				}
			}
		}
	}

	if (blist_size_f > 0 and flist_size_b > 0) {
		if (blist_size_f <= flist_size_b) {
			for (i = 0; i < blist_size_f; i++) {
				if (backwardlist_f[i] < forwardlist_b[0])
					continue;
				else {
					j = binary_search(forwardlist_b, 0, flist_size_b, backwardlist_f[i]);
				}
				spos_front = backwardlist_f[i];
				spos_back = forwardlist_b[j-1];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_front <= spos_back + GENETHRESH) {
					return true;
				}
			}
		}
		else {
			for (j = 0; j < flist_size_b; j++) {
				if (forwardlist_b[j] + len_b < backwardlist_f[0])	// no position less than forwardlist_b[j] available on backwardlist_f
					i = 0;
				else {
					target = forwardlist_b[j] + len_b;
					i = binary_search(backwardlist_f, 0, blist_size_f, target);
					if (i >= blist_size_f)
						continue;
				}
				spos_front = backwardlist_f[i];
				spos_back = forwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_front <= spos_back + GENETHRESH) {
					return true;
				}
			}
		}
	}
	return false;
}

// return value:
// 1: concordant juinction
// 2: chimeric
// 3: discordant
int intersect(const bwtint_t& sp_f, const bwtint_t& ep_f, const int& len_f, const bwtint_t& sp_b, const bwtint_t& ep_b, const int& len_b, MatchedRead& mr) {
	//if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
	if ((ep_b - sp_b + 1) >= REGIONSIZELIM or (ep_f - sp_f + 1) >= REGIONSIZELIM) {
		return 3;
	}

	//fprintf(stderr, "back size: %llu\nfront size: %llu\n", ep_b - sp_b + 1, ep_f - sp_f + 1);

	mr.is_concord = false;
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

	bwtint_t i, j, target;
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
				
				//fprintf(stderr, "HEYYYYY %d\n%llu\n", len_f, spos_front);
				if (spos_front <= spos_back + len_b + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_back;
					mr.matched_len = spos_front + len_f - spos_back + 1;
					mr.dir = -1;
					return 1;
				}
			}
		}
		else {
			for (j = 0; j < blist_size_b; j++) {
				if (backwardlist_b[j] + len_b < forwardlist_f[0])
					i = 0;
				else {
					target = backwardlist_b[j] + len_b;
					i = binary_search(forwardlist_f, 0, flist_size_f, target);
					if (i >= flist_size_f)	// no position greater than backwardlist_b[j] available on forwardlist_f
						continue;
				}
				spos_front = forwardlist_f[i];
				spos_back = backwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_front <= spos_back + len_b + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_back;
					mr.matched_len = spos_front + len_f - spos_back + 1;
					mr.dir = -1;
					return 1;
				}
			}
		}
	}

	if (blist_size_f > 0 and flist_size_b > 0) {
		if (blist_size_f <= flist_size_b) {
			for (i = 0; i < blist_size_f; i++) {
				if (backwardlist_f[i] + len_f < forwardlist_b[0])
					j = 0;
				else {
					target = backwardlist_f[i] + len_f;
					j = binary_search(forwardlist_b, 0, flist_size_b, target);
					if (j >= flist_size_b)	// no position greater than backwardlist_f[i] available on backwardlist_f
						continue;
				}
				spos_front = backwardlist_f[i];
				spos_back = forwardlist_b[j];
				//fprintf(stderr, "---Start pos front: %llu\n---Start pos back: %llu\n Dir front: %d\n Dir back: %d\n", spos_front, spos_back, dir_front, dir_back);
				
				if (spos_back <= spos_front + len_f + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_front;
					mr.matched_len = spos_back + len_b - spos_front + 1;
					mr.dir = 1;
					return 1;
				}
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
				
				if (spos_back <= spos_front + len_f + GENETHRESH) {
					mr.is_concord = true;
					mr.start_pos = spos_front;
					mr.matched_len = spos_back + len_b - spos_front + 1;
					mr.dir = 1;
					return 1;
				}
			}
		}
	}

	bool is_chimeric = is_chimeric_intersect(forwardlist_f, flist_size_f, backwardlist_f, blist_size_f, len_f, forwardlist_b, flist_size_b, backwardlist_b, blist_size_b, len_b);
	if (is_chimeric)
		return 2;
	return 3;
}

// return value:
// 0 : concordant (exon)
// 1 : concordant (junction)
// 2 : chimeric
// 3 : discordant
// 4 : potentially mappable
// 5 : un-mappable
int find_expanded_sliding_positions(const char* rseq, const char* rcseq, const int& rseq_len, MatchedRead& mr, const int& window_size, const int& step, const int& junction_detect_size_lim) {
	bwtint_t sp_b, ep_b;
	bwtint_t sp_f, ep_f;
	int exp_len_back = 0;
	int exp_len_front = 0;

	mr.is_concord = false;

	int i;
	for (i = 0; i < rseq_len - window_size; i += step) {
		exp_len_front = get_expanded_locs(rcseq, rseq_len - i, &sp_f, &ep_f);
		if (exp_len_front < window_size)
			continue;
		if (i < junction_detect_size_lim and rseq_len - exp_len_front - i < junction_detect_size_lim) {	// concordant (exon)
			//fprintf(stderr, "[Concordant-e]\tfront matched: %d\n", exp_len_front);
			mr.is_concord = true;
			mr.start_pos = get_pos(&sp_f, exp_len_front, mr.dir);
			mr.dir *= -1;
			mr.matched_len = exp_len_front;
			return 0;
		}
		break;
	}

	if (exp_len_front < window_size) {	// un-mappable
		mr.is_concord = false;
		mr.start_pos = -1;
		mr.matched_len = -1;
		return 5;
	}

	const char* remain_rseq = rseq + i + exp_len_front;
	int remain_len = rseq_len - i - exp_len_front;
	if (remain_len < junction_detect_size_lim) {	// won't be able to find an event from the remaining of the read (potentially mappable)
		mr.is_concord = false;
		mr.start_pos = get_pos(&sp_f, exp_len_front, mr.dir);
		mr.dir *= -1;
		mr.matched_len = exp_len_front;
		return 4;
	}

	int j = 0;
	bwtint_t longest_sp_b, longest_ep_b;
	int max_len_back = 0;
	for (j = 0; j < remain_len - exp_len_back; j += step) {
		exp_len_back = get_expanded_locs(remain_rseq, remain_len - j, &sp_b, &ep_b);
		//fprintf(stderr, "%s\tlen front: %d,\tlen back: %d\n", rseq, exp_len_front, exp_len_back);
	
		if (exp_len_back >= window_size)
			break;
		if (exp_len_back > max_len_back) {
			max_len_back = exp_len_back;
			longest_sp_b = sp_b;
			longest_ep_b = ep_b;
		}

		//if (exp_len_front + exp_len_back <= 0.6 * rseq_len) {
		//	fprintf(stderr, "[Chimeric-e]\tfront matched: %d, back matched: %d\n", exp_len_front, exp_len_back);
		//	return 1;
		//}
	}
	
	if (max_len_back <= junction_detect_size_lim) {	// not long enough for detecting junction
		mr.is_concord = false;
		mr.start_pos = get_pos(&sp_f, exp_len_front, mr.dir);
		mr.dir *= -1;
		mr.matched_len = exp_len_front;
		return 4;
	}
	//bool consis = check_cosistent_match(&sp_b, &ep_b, exp_len_backward, &sp_f, &ep_f, exp_len_forward);
	//bool consis = is_concordant(sp_f, ep_f, exp_len_front, sp_b, ep_b, exp_len_back, noise_thresh);
	//bool consis = is_concordant_sorted2(sp_f, ep_f, exp_len_front, longest_sp_b, longest_ep_b, exp_len_back, 10, mr);
	
	int intersect_ret = intersect(sp_f, ep_f, exp_len_front, longest_sp_b, longest_ep_b, exp_len_back, mr);

	//std::string res = (consis) ? "Concordant" : "Chimeric";
	//fprintf(stderr, "[%s]\tfront matched: %d, back matched: %d\n", res.c_str(), exp_len_front, exp_len_back);
	
	return intersect_ret;
}

// return value:
// 0 : concordant (exon)
// 1 : concordant (junction)
// 2 : chimeric
// 3 : discordant
// 4 : potentially mappable
// 5 : un-mappable
int check_concordant_mates_expand(const Record* m1, const Record* m2) {
	int junction_detect_size_lim = 9;
	MatchedRead mr1, mr2;

	//fprintf(stderr, "Read name: %s", m1->rname);
	int mate1_state = find_expanded_sliding_positions(m1->seq, m1->rcseq, m1->seq_len, mr1, 23, 3, junction_detect_size_lim);
	int mate2_state = find_expanded_sliding_positions(m2->seq, m2->rcseq, m2->seq_len, mr2, 23, 3, junction_detect_size_lim);

	//fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

	if (mate1_state == 2 or mate2_state == 2)
		return 2;

	if (mate1_state == 5)
		return mate2_state;
	if (mate2_state == 5)
		return mate1_state;

	if (mate1_state == 3 or mate2_state == 3)
		return 3;

	//fprintf(stderr, "%s", m1->rname);
	//fprintf(stderr, "M1: dir:%d,\tstart:%llu,\tend:%llu\n", mr1.dir, mr1.start_pos, mr1.start_pos + mr1.matched_len);
	//fprintf(stderr, "M2: dir:%d,\tstart:%llu,\tend:%llu\n", mr2.dir, mr2.start_pos, mr2.start_pos + mr2.matched_len);

	if (mr1.dir == 1 and mr2.dir == -1 and (mr1.start_pos <= mr2.start_pos))
		return 0;
	if (mr1.dir == -1 and mr2.dir == 1 and (mr2.start_pos <= mr1.start_pos))
		return 0;

	if (mr1.dir == 1 and mr2.dir == -1 and (mr1.start_pos > mr2.start_pos + 8) and (mr1.start_pos <= mr2.start_pos + GENETHRESH)) {
		//fprintf(stderr, "%s", m1->rname);
		//fprintf(stderr, "M1: dir: %d\t%llu\n", mr1.dir, mr1.start_pos);
		//fprintf(stderr, "M2: dir: %d\t%llu\n", mr2.dir, mr2.start_pos);

		//char *chr_name;
		//int32_t chr_len;
		//uint32_t chr_beg;
		//uint32_t chr_end;
		//bwt_get_intv_info(mr1.start_pos, mr1.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		//fprintf(stderr, "Chr: %s\t", chr_name);
		//fprintf(stderr, "%"PRIu64"\n", chr_beg);
		//bwt_get_intv_info(mr2.start_pos, mr2.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		//fprintf(stderr, "Chr: %s\t", chr_name);
		//fprintf(stderr, "%"PRIu64"\n", chr_beg);

		return 2;
	}
	if (mr1.dir == -1 and mr2.dir == 1 and (mr2.start_pos > mr1.start_pos + 8) and (mr2.start_pos <= mr1.start_pos + GENETHRESH)) {
		//fprintf(stderr, "%s", m1->rname);
		//fprintf(stderr, "M1: dir: %d\t%llu\n", mr1.dir, mr1.start_pos);
		//fprintf(stderr, "M2: dir: %d\t%llu\n", mr2.dir, mr2.start_pos);

		//char *chr_name;
		//int32_t chr_len;
		//uint32_t chr_beg;
		//uint32_t chr_end;
		//bwt_get_intv_info(mr1.start_pos, mr1.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		//fprintf(stderr, "Chr: %s\t", chr_name);
		//fprintf(stderr, "%"PRIu64"\n", chr_beg);
		//bwt_get_intv_info(mr2.start_pos, mr2.start_pos + 20, &chr_name, &chr_len, &chr_beg, &chr_end);
		//fprintf(stderr, "Chr: %s\t", chr_name);
		//fprintf(stderr, "%"PRIu64"\n", chr_beg);

		return 2;
	}
	return 3;
}

