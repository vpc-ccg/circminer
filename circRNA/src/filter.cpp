#include "filter.h"

FilterRead::FilterRead (char* save_fname, bool pe) {
	is_pe = pe;

	char ignore_file1[1000], keep_file1[1000];
	strcpy (ignore_file1, save_fname);
	strcpy (keep_file1  , save_fname);

	strcat (ignore_file1, "_R1.ignore.fastq");
	strcat (keep_file1, "_R1.keep.fastq");

	ignore_r1 = fopen(ignore_file1, "w");
	keep_r1   = fopen(keep_file1, "w");

	if (is_pe) {
		char ignore_file2[1000], keep_file2[1000];
		strcpy (ignore_file2, save_fname);
		strcpy (keep_file2  , save_fname);

		strcat (ignore_file2, "_R2.ignore.fastq");
		strcat (keep_file2, "_R2.keep.fastq");

		ignore_r2 = fopen(ignore_file2, "w");
		keep_r2   = fopen(keep_file2, "w");
	}
}

FilterRead::~FilterRead (void) {
	if (ignore_r1 != NULL)
		fclose(ignore_r1);
	if (keep_r1 != NULL)
		fclose(keep_r1);
	if (ignore_r2 != NULL)
		fclose(ignore_r2);
	if (keep_r2 != NULL)
		fclose(keep_r2);
}

// write reads SE mode
void FilterRead::write_read (Record* current_record, bool is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
	}
	else {
		fprintf(keep_r1, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
	}
}

// write reads PE mode
void FilterRead::write_read (Record* current_record1, Record* current_record2, bool is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
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
	if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
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
	if ((ep_f - sp_f >= REGIONSIZELIM) or (ep_b - sp_b >= REGIONSIZELIM))
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
		//fprintf(stderr, "[Concordant-e]\tfront matched: %d", exp_len_front);
		return 0;
	}

	exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front, &sp_b, &ep_b);
	//exp_len_back = get_expanded_locs(rseq + exp_len_front, rseq_len - exp_len_front - noise_thresh, &sp_b, &ep_b);
	//fprintf(stderr, "%s\tlen front: %d,\tlen back: %d\n", rseq, exp_len_front, exp_len_back);
	//if (occ > 0)
	//	locate_match(&sp, &ep, window_size);
	
	if (exp_len_front + exp_len_back <= rseq_len / 2) {
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
	
	if (exp_len_front < 18 or  exp_len_back < 18) {
		fprintf(stderr, "[Chimerici-e]\tfront matched: %d, back matched: %d\n", exp_len_front, exp_len_back);
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

