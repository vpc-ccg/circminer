#include <vector>
#include <cstring>
#include <algorithm>

#include "process_circ.h"
#include "gene_annotation.h"
#include "hash_table.h"
#include "fastq_parser.h"
#include "match_read.h"
#include "chain.h"
#include "extend.h"
#include "utils.h"

#define BINSIZE 5000
#define MAXHTLISTSIZE 5

ProcessCirc::ProcessCirc (int last_round_num, int ws) {
	sprintf(fq_file1, "%s_%d_remain_R1.fastq", outputFilename, last_round_num);
	sprintf(fq_file2, "%s_%d_remain_R2.fastq", outputFilename, last_round_num);

	sort_fq(fq_file1);
	sort_fq(fq_file2);

	sprintf(fq_file1, "%s_%d_remain_R1.fastq.srt", outputFilename, last_round_num);
	sprintf(fq_file2, "%s_%d_remain_R2.fastq.srt", outputFilename, last_round_num);

	fprintf(stdout, "%s\n",fq_file1 );

	report_file = NULL;

	window_size = ws;
	step = 3;

	pre_contig = -1;

	RegionalHashTable* rht;
	for (int i = 0; i < MAXHTLISTSIZE; ++i) {
		rht = new RegionalHashTable(ws, 0, 0);
		ind2ht.insert(make_pair(i, rht));

		removables.insert(i);

		gid2ind.insert(make_pair(INF-i, i));
		gids.push_back(INF-i);
	}


	int max_kmer_cnt = (maxReadLength - window_size) / step + 1;
	bc.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	for (int i = 0; i < BESTCHAINLIM; i++)
		bc.chains[i].frags = (fragment_t*) malloc(max_kmer_cnt * sizeof(fragment_t));
}

ProcessCirc::~ProcessCirc (void) {
	close_file(report_file);

	for (int i = 0; i < BESTCHAINLIM; i++)
		free(bc.chains[i].frags);

	free(bc.chains);
}

void ProcessCirc::sort_fq(char* fqname) {
	fprintf(stdout, "Sorting remaining read mappings... ");
	if (!system(NULL)) {
		fprintf(stdout, "Failed using command line\n");
		exit(EXIT_FAILURE);
	}

	char command [FILE_NAME_LENGTH];
	sprintf(command, "paste - - - - < %s | sort -k22,22 -k3,3 -k4,4n | tr \"\t\" \"\n\" > %s.srt", fqname, fqname);

	int ret = system(command);
	if (ret == 0)
		fprintf(stdout, "OK\n");
}

void ProcessCirc::do_process (void) {
	/**********************/
	/**Loading Hash Table**/
	/**********************/

	double cputime_start = get_cpu_time();
	double realtime_start = get_real_time();
	double cputime_curr;
	double realtime_curr;

	checkSumLength = (WINDOW_SIZE > kmer) ? 0 : kmer - WINDOW_SIZE;

	char index_file [FILE_NAME_LENGTH];
	strcpy(index_file, referenceFilename);
	strcat(index_file, ".index");

	initCommon();
	
	THREAD_COUNT = threads;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	if (!checkHashTable(index_file))
		return;

	ContigLen* orig_contig_len;
	int contig_cnt;
	if (!initLoadingHashTableMeta(index_file, &orig_contig_len, &contig_cnt))
		return;

	/**********************/
	/**Finised Loading HT**/
	/**********************/

	/*******************/
	/**GTF Parser Init**/
	/*******************/

	if (stage == 1) { 
		gtf_parser.init(gtfFilename, orig_contig_len, contig_cnt);
		if (! gtf_parser.load_gtf()) {
			fprintf(stdout, "Error in reading GTF file.\n");
			exit(1);
		}
		else 
			fprintf(stdout, "GTF file successfully loaded!\n");

		cputime_curr = get_cpu_time();
		realtime_curr = get_real_time();

		fprintf(stdout, "[P] Loaded GTF in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - cputime_start, realtime_curr - realtime_start);

		cputime_start = cputime_curr;
		realtime_start = realtime_curr;
	}

	/*******************/
	/**Finished GTF PI**/
	/*******************/

	double fq_cputime_start = cputime_start;
	double fq_realtime_start = realtime_start;

	contigName = getRefGenomeName();

	bool is_pe = pairedEnd;

	FASTQParser fq_parser1(fq_file1);
	Record* current_record1;

	FASTQParser fq_parser2;
	Record* current_record2;
	if (is_pe) {
		fq_parser2.init(fq_file2);
	}
	
	int line = 0;
	while ( (current_record1 = fq_parser1.get_next()) != NULL ) { // go line by line on fastq file
		if (is_pe)
			current_record2 = fq_parser2.get_next();
		if (current_record1 == NULL)	// no new line
			break;

		line++;
		vafprintf(2, stderr, "Line: %d\n", line);

		if (line % LINELOG == 0) {
			cputime_curr = get_cpu_time();
			realtime_curr = get_real_time();

			fprintf(stdout, "[P] %d reads in %.2lf CPU sec (%.2lf real sec)\t Look ups: %u\n", line, cputime_curr - cputime_start, realtime_curr - realtime_start, lookup_cnt);
			fflush(stdout);

			cputime_start = cputime_curr;
			realtime_start = realtime_curr;
			
			lookup_cnt = 0;
		}

		while (pre_contig != current_record1->mr->contig_num) {
			load_genome();
			pre_contig = atoi(getRefGenomeName()) - 1;
		}

		call_circ(current_record1, current_record2);
	}

	report_events();

	cputime_curr = get_cpu_time();
	realtime_curr = get_real_time();

	fprintf(stdout, "[P] Mapping in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

	finalizeLoadingHashTable();
}

// PE
void ProcessCirc::call_circ(Record* current_record1, Record* current_record2) {
	MatchedRead mr = *(current_record1->mr);
	vafprintf(2, stderr, "%s\n%s\n", current_record1->seq, current_record2->seq);
	vafprintf(2, stderr, "%s\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%d\t%d\t%d\t%d\n", 
									current_record1->rname, 
									mr.chr_r1.c_str(), mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, mr.ed_r1, 
									mr.chr_r2.c_str(), mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, mr.ed_r2,
									mr.tlen, mr.junc_num, mr.gm_compatible, mr.type);

	bool r1_partial = mr.mlen_r1 < mr.mlen_r2;
	char* remain_seq = (r1_partial) ? ((mr.r1_forward) ? current_record1->seq : current_record1->rcseq) : 
									  ((mr.r2_forward) ? current_record2->seq : current_record2->rcseq) ;
	uint32_t qspos = (r1_partial) ? (((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (1) : (mr.qepos_r1 + 1)) :
									(((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (1) : (mr.qepos_r2 + 1)) ;

	uint32_t qepos = (r1_partial) ? (((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (mr.qspos_r1 - 1) : (current_record1->seq_len)) :
									(((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (mr.qspos_r2 - 1) : (current_record2->seq_len)) ;


	int whole_seq_len = (r1_partial) ? current_record1->seq_len : current_record2->seq_len;
	int remain_len = qepos - qspos + 1;
	if (qepos < qspos or remain_len < window_size) {	// it was fully mapped
		return;
	}

	const IntervalInfo<GeneInfo>* gene_info = gtf_parser.get_gene_overlap(mr.spos_r1, false);
	bool found = (gene_info != NULL);
	if (! found) {
		vafprintf(2, stderr, "Gene not found!\n");
		return;
	}
	vafprintf(2, stderr, "# Gene overlaps: %d\n", gene_info->seg_list.size());

	// fill MatchedMate
	MatchedMate mm_r1(mr, 1, current_record1->seq_len);
	MatchedMate mm_r2(mr, 2, current_record2->seq_len);

	check_removables(mr.spos_r1);
	RegionalHashTable* regional_ht;

	for (int i = 0; i < gene_info->seg_list.size(); ++i) {
		uint32_t gene_len = gene_info->seg_list[i].end - gene_info->seg_list[i].start + 1;
		char gene_seq[gene_len + 1];
		gene_seq[gene_len] = '\0';
		pac2char(gene_info->seg_list[i].start, gene_len, gene_seq);
		//vafprintf(2, stderr, "Gene: %s\n", gene_seq);

		regional_ht = get_hash_table(gene_info->seg_list[i], gene_seq);
		
		vafprintf(2, stderr, "R%d partial: [%d-%d]\n", (int) (!r1_partial) + 1, qspos, qepos);
		vafprintf(2, stderr, "%s\n", remain_seq);

		//binning(qspos, qepos, regional_ht, remain_seq, gene_len);
		uint32_t rspos, repos;
		chaining(qspos, qepos, regional_ht, remain_seq, gene_len, gene_info->seg_list[i].start, rspos, repos);

		// find rspos and repos for the best chain
		bool forward = (r1_partial) ? (mr.r1_forward) : (mr.r2_forward);
		int dir = (forward) ? 1 : -1;
		MatchedMate partial_mm;
		find_exact_coord(mm_r1, mm_r2, partial_mm, dir, qspos, remain_seq, remain_len, whole_seq_len);

		if (partial_mm.type == CONCRD) {
			vafprintf(2, stderr, "Coordinates: [%d-%d]\n", partial_mm.spos, partial_mm.epos);
			fprintf(stderr, "%s\t%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t", current_record1->rname, mr.chr_r1.c_str(), 
														partial_mm.spos, partial_mm.epos, partial_mm.qspos, partial_mm.matched_len, partial_mm.dir,
														mm_r1.spos, mm_r1.epos, mm_r1.qspos, mm_r1.matched_len, mm_r1.dir,
														mm_r2.spos, mm_r2.epos, mm_r2.qspos, mm_r2.matched_len, mm_r2.dir);
			
			CircRes cr;
			int type = check_split_map(mm_r1, mm_r2, partial_mm, r1_partial, cr);
			fprintf(stderr, "%d\n", type);
			
			if (type == CR) {
				cr.chr = mr.chr_r1;
				cr.rname = current_record1->rname;
				circ_res.push_back(cr);
				break;	// stop processing next genes
			}
		}
		else {
			vafprintf(2, stderr, "Coordinates: [%d-%d]\n", rspos, repos);
		}
	}
}

void ProcessCirc::binning(uint32_t qspos, uint32_t qepos, RegionalHashTable* regional_ht, char* remain_seq, uint32_t gene_len) {
	int bin_num = gene_len / BINSIZE + 1;
	int bins[bin_num];
	int max_id = 0;
	memset(bins, 0, bin_num * sizeof(int));

	for (int i = qspos - 1; i <= qepos - window_size; i += step) {
		GIMatchedKmer* gl = regional_ht->find_hash(regional_ht->hash_val(remain_seq + i));
		if (gl == NULL) {
			vafprintf(2, stderr, "Hash val not found!!!\n");
		}

		vafprintf(2, stderr, "Occ: %d\n", gl->frag_count);
		
		for (int j = 0; j < gl->frag_count; j++) {
			int bin_id = gl->frags[j].info / BINSIZE;
			bins[bin_id]++;

			if (bins[bin_id] > bins[max_id])
				max_id = bin_id;

			vafprintf(2, stderr, "%d - %d\t", gl->frags[j].info, bins[bin_id]);
		}
		vafprintf(2, stderr, "\n");
	}	
	vafprintf(2, stderr, "Biggest bin: bin[%d][%d - %d] = %d\n", max_id, max_id * BINSIZE, (max_id+1) * BINSIZE - 1, bins[max_id]);

}

void ProcessCirc::chaining(uint32_t qspos, uint32_t qepos, RegionalHashTable* regional_ht, char* remain_seq, 
							uint32_t gene_len, uint32_t shift, uint32_t& rspos, uint32_t& repos) {
	int seq_len = qepos - qspos + 1;
	int kmer_cnt = ((qepos - qspos + 1) - window_size) / step + 1;
	GIMatchedKmer fl[kmer_cnt+1];

	int l = 0;
	for (int i = qspos - 1; i <= qepos - window_size; i += step) {
		GIMatchedKmer* gl = regional_ht->find_hash(regional_ht->hash_val(remain_seq + i));
		if (gl == NULL) {	// has N inside kmer
			vafprintf(2, stderr, "Hash val not found!!!\n");
			continue;
		}
		fl[l] = *gl;
		fl[l].qpos = i;

		// vafprintf(2, stderr, "Occ: %d\n", fl[l].frag_count);
		
		// for (int j = 0; j < fl[l].frag_count; j++) {
		// 	// fl[l]->frags[j].info += shift;
		// 	vafprintf(2, stderr, "%d\t", fl[l].frags[j].info + shift);
		// }
		// vafprintf(2, stderr, "\n");
		l++;
	}

	kmer_cnt = l;
	//printf("kmer cnt: %d\n", kmer_cnt);

	chain_seeds_sorted_kbest2(qepos, fl, bc, window_size, kmer_cnt, shift);

	vafprintf(1, stderr, "Chaining score:%.4f,\t len: %lu\n", bc.chains[0].score, (unsigned long)bc.best_chain_count);

	int allowed_missed_kmers = (qepos - qspos + 1) / 20 * 3 + 1;
	vafprintf(2, stderr, "Allowed missing kmers: %d\n", allowed_missed_kmers);

	rspos = 0;
	repos = 0;
	uint32_t curr_rspos, curr_repos;
	int least_miss = INF;
	int least_miss_ind = -1;
	int missing;
	for (int j = 0; j < bc.best_chain_count; j++) {
		missing = kmer_cnt - bc.chains[j].chain_len;
		vafprintf(2, stderr, "Actual missing: %d\n", missing);
		if (missing > least_miss)
			break;

		curr_rspos = bc.chains[j].frags[0].rpos - bc.chains[j].frags[0].qpos;
		curr_repos = bc.chains[j].frags[bc.chains[j].chain_len - 1].rpos + (seq_len - 1 - bc.chains[j].frags[bc.chains[j].chain_len - 1].qpos);
		uint32_t curr_qlen = bc.chains[j].frags[bc.chains[j].chain_len - 1].qpos - bc.chains[j].frags[0].qpos;
		uint32_t qlen = (least_miss_ind != -1) ? bc.chains[least_miss_ind].frags[bc.chains[least_miss_ind].chain_len - 1].qpos - bc.chains[least_miss_ind].frags[0].qpos : 0;
		if (missing < least_miss or (missing == least_miss and curr_qlen > qlen)) {
			least_miss = missing;
			least_miss_ind = j;
			
			// rspos = curr_rspos + shift;
			// repos = curr_repos + shift;

			rspos = bc.chains[j].frags[0].rpos;
			repos = bc.chains[j].frags[bc.chains[j].chain_len - 1].rpos + bc.chains[j].frags[bc.chains[j].chain_len - 1].len - 1;
		}

		for (int i = 0; i < bc.chains[j].chain_len; i++) {
			vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, bc.chains[j].frags[i].rpos, bc.chains[j].frags[i].qpos, bc.chains[j].frags[i].len);
		}
	}
}

bool ProcessCirc::find_exact_coord(MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, 
									int dir, uint32_t qspos, char* rseq, int rlen, int whole_len) {

	uint32_t partial_spos = bc.chains[0].frags[0].rpos;
	uint32_t partial_epos = bc.chains[0].frags[bc.chains[0].chain_len - 1].rpos + bc.chains[0].frags[bc.chains[0].chain_len - 1].len - 1;
	uint32_t partial_qspos = qspos;
	uint32_t partial_qepos = qspos + rlen - 1;

	partial_mm.set(partial_spos, partial_epos, partial_qspos, partial_qepos, dir);

	--qspos;	// convert to 0-based

	overlap_to_spos(mm_r1);
	overlap_to_spos(mm_r2);
	overlap_to_spos(partial_mm);

	vector<uint32_t> common_tid;
	bool success = same_transcript(mm_r1.exons_spos, mm_r2.exons_spos, partial_mm.exons_spos, common_tid);
	if (!success) {
		// fprintf(stderr, "No common transcript!\n");
		return false;
	}

	partial_mm.middle_ed = calc_middle_ed(bc.chains[0], maxEd, rseq, rlen);
	// fprintf(stderr, "Middle ed: %d\n", partial_mm.middle_ed);
	if (partial_mm.middle_ed > maxEd)
		return false;

	bool extend = true;
	partial_mm.is_concord = false;
	if (bc.chains[0].chain_len <= 0) {
		partial_mm.type = ORPHAN;
		partial_mm.matched_len = 0;
		return false;
	}

	bool lok;
	bool rok;
	int err = partial_mm.middle_ed;
	if (extend) {
		partial_mm.matched_len = rlen;
		lok = extend_chain_left (common_tid, bc.chains[0], rseq + qspos, qspos, MINLB, partial_mm, err);
		if (qspos == 0) {
			rok = extend_chain_right(common_tid, bc.chains[0], rseq, rlen, MAXUB, partial_mm, err);
		}
		else {
			rok = extend_chain_right(common_tid, bc.chains[0], rseq, whole_len, MAXUB, partial_mm, err);
		}
		update_match_mate_info(lok, rok, err, partial_mm);
	}

	return partial_mm.type == CONCRD;
}

void ProcessCirc::check_removables (uint32_t rspos) {
	for (auto it = ind2ht.begin(); it != ind2ht.end(); ++it) {
		if (rspos > it->second->gene_epos) {
			removables.insert(it->first);
			// fprintf(stderr, "Removed gid: %d --removables size: %d\n", gids[it->first], removables.size());
		}
	}
}

RegionalHashTable* ProcessCirc::get_hash_table (const GeneInfo& gene_info, char* gene_seq) {
	uint32_t gid, removable_ind;
	RegionalHashTable* regional_ht;
	RegionalHashTable* new_ht;

	int gene_len = gene_info.end - gene_info.start + 1;
	gid = gene_info.gene_id;
	// fprintf(stderr, "Gene id: %d\n", gid);
	
	if (gid2ind.find(gid) != gid2ind.end()) {
		uint32_t ind = gid2ind[gid];
		regional_ht = ind2ht[ind];
		// fprintf(stderr, "Found gid: %d\n", gid);
	}
	else {
		if (removables.size() == 0) {
			new_ht = new RegionalHashTable(window_size, gene_info.start, gene_info.end);

			uint32_t new_ind = ind2ht.size();
			ind2ht.insert(make_pair(new_ind, new_ht));
			gid2ind.insert(make_pair(gid, new_ind));
			gids.push_back(gid);
			new_ht->create_table(gene_seq, 0, gene_len);
			regional_ht = new_ht;
			// fprintf(stderr, "Allocated new HT gid:\n");
		}
		else {
			removable_ind = *(removables.cbegin());
			removables.erase(removables.cbegin());
			gid2ind.erase(gids[removable_ind]);
			gid2ind.insert(make_pair(gid, removable_ind));
			gids[removable_ind] = gid;

			regional_ht = ind2ht[removable_ind];
			regional_ht->gene_spos = gene_info.start;
			regional_ht->gene_epos = gene_info.end;
			regional_ht->create_table(gene_seq, 0, gene_len);

		}
		// fprintf(stderr, "Added gid: %d\n", gid);
	}

	return regional_ht;
}

int ProcessCirc::check_split_map (MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, bool r1_partial, CircRes& cr) {
	if (r1_partial) {
		if (mm_r1.qspos < partial_mm.qspos)
			return final_check(mm_r2, mm_r1, partial_mm, cr);
		else 
			return final_check(mm_r2, partial_mm, mm_r1, cr);
	}
	else {
		if (mm_r2.qspos < partial_mm.qspos)
			return final_check(mm_r1, mm_r2, partial_mm, cr);
		else 
			return final_check(mm_r1, partial_mm, mm_r2, cr);
	}
	return UD;
}

// full_mm -> not split mate
// split_mm_left -> left hand side of the split read
// split_mm_right -> right hand side of the split read
int ProcessCirc::final_check (MatchedMate& full_mm, MatchedMate& split_mm_left, MatchedMate& split_mm_right, CircRes& cr) {
	if (split_mm_left.epos < split_mm_right.spos) {
		if (full_mm.dir == 1) {
			if (full_mm.spos <= split_mm_left.spos)
				return FR;
			else if (full_mm.epos >= split_mm_right.epos)
				return RF;
		}

		if (full_mm.dir == -1) {
			if (full_mm.epos >= split_mm_right.epos)
				return FR;
			else if (full_mm.spos <= split_mm_left.spos)
				return RF;
		}
	}

	else if (split_mm_right.epos < split_mm_left.spos) {
		if (full_mm.spos >= split_mm_right.spos and full_mm.epos <= split_mm_left.epos) {
			// check splice site
			overlap_to_spos(full_mm);
			overlap_to_epos(full_mm);

			overlap_to_spos(split_mm_right);
			overlap_to_epos(split_mm_right);

			overlap_to_spos(split_mm_left);
			overlap_to_epos(split_mm_left);

			if (split_mm_left.exons_epos == NULL or split_mm_right.exons_spos == NULL)
				return MCR;

			vector <uint32_t> end_tids;
			for (int i = 0; i < split_mm_left.exons_epos->seg_list.size(); ++i) {
				if (split_mm_left.epos == split_mm_left.exons_epos->seg_list[i].end) {
					for (int j = 0; j < split_mm_left.exons_epos->seg_list[i].trans_id.size(); ++j) {
						end_tids.push_back(split_mm_left.exons_epos->seg_list[i].trans_id[j]);
					}
				}
			}

			vector <uint32_t> start_tids;
			for (int i = 0; i < split_mm_right.exons_spos->seg_list.size(); ++i) {
				if (split_mm_right.spos == split_mm_right.exons_spos->seg_list[i].start) {
					for (int j = 0; j < split_mm_right.exons_spos->seg_list[i].trans_id.size(); ++j) {
						start_tids.push_back(split_mm_right.exons_spos->seg_list[i].trans_id[j]);
					}
				}
			}

			for (int i = 0; i < start_tids.size(); ++i)
				for (int j = 0; j < end_tids.size(); ++j)
					if (start_tids[i] == end_tids[j]) {
						cr.set_bp(split_mm_right.spos, split_mm_left.epos);
						return CR;
					}


			if (start_tids.size() > 0 and end_tids.size() > 0)
				return NCR;

			return MCR;
		}
	}
	return UD;
}

void ProcessCirc::report_events (void) {
	open_report_file();
	if (circ_res.size() <= 0)
		return;

	sort(circ_res.begin(), circ_res.end());
	int cnt = 1;
	CircRes last = circ_res[0];
	for (int i = 1; i < circ_res.size(); ++i) {
		if (circ_res[i] == last) {
			cnt++;
		}
		else {
			fprintf(report_file, "%s\t%u\t%u\t%d\n", last.chr.c_str(), last.spos, last.epos, cnt);
			cnt = 1;
			last = circ_res[i];
		}
	}
	fprintf(report_file, "%s\t%u\t%u\t%d\n", last.chr.c_str(), last.spos, last.epos, cnt);
}

void ProcessCirc::open_report_file (void) {
	char temp_fname [FILE_NAME_LENGTH];
	sprintf(temp_fname, "%s.circ_report", outputFilename);
	report_file = open_file(temp_fname, "w");
}

void ProcessCirc::load_genome (void) {
	double cputime_start = get_cpu_time();
	double realtime_start = get_real_time();
	double cputime_curr;
	double realtime_curr;
	double tmpTime;

	fprintf(stdout, "Started loading index...\n");

	int flag = loadCompressedRefGenome ( &tmpTime );  			// Reading a fragment

	cputime_curr = get_cpu_time();
	realtime_curr = get_real_time();

	fprintf(stdout, "[P] Loaded genome index successfully in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - cputime_start, realtime_curr - realtime_start);

	cputime_start = cputime_curr;
	realtime_start = realtime_curr;
}

int ProcessCirc::get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos) {


}