#include <vector>
#include <cstring>
#include <algorithm>
#include <utility>
#include <cstdlib>

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

typedef pair<uint32_t, int> pu32i;

void set_mm(const chain_t& ch, uint32_t qspos, int rlen, int dir, MatchedMate& mm);

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
	bc1.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	for (int i = 0; i < BESTCHAINLIM; i++)
		bc1.chains[i].frags = (fragment_t*) malloc(max_kmer_cnt * sizeof(fragment_t));
	bc2.chains = (chain_t*) malloc(BESTCHAINLIM * sizeof(chain_t));
	for (int i = 0; i < BESTCHAINLIM; i++)
		bc2.chains[i].frags = (fragment_t*) malloc(max_kmer_cnt * sizeof(fragment_t));

	// circ_type.push_back("SingleTranscriptomeCircRNA");
	// circ_type.push_back("MultiTranscriptomeCircRNA");
	// circ_type.push_back("NovelCircRNA");

	circ_type.push_back("STC");
	circ_type.push_back("MTC");
	circ_type.push_back("NC");
}

ProcessCirc::~ProcessCirc (void) {
	close_file(report_file);

	for (int i = 0; i < BESTCHAINLIM; i++)
		free(bc1.chains[i].frags);

	for (int i = 0; i < BESTCHAINLIM; i++)
		free(bc2.chains[i].frags);

	free(bc1.chains);
	free(bc2.chains);
}

void ProcessCirc::sort_fq(char* fqname) {
	fprintf(stdout, "Sorting remaining read mappings... ");
	if (!system(NULL)) {
		fprintf(stdout, "Error: Unable to run system call.\n");
		exit(EXIT_FAILURE);
	}

	char command [FILE_NAME_LENGTH];
	sprintf(command, "sort -k22,22 -k3,3 -k4,4n %s | tr \"\t\" \"\n\" > %s.srt", fqname, fqname);

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
	
	THREAD_COUNT = threadCount;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	if (!checkHashTable(index_file))
		return;

	ContigLen* orig_contig_len;
	int contig_cnt;
	if (!initLoadingCompressedGenomeMeta(index_file, &orig_contig_len, &contig_cnt))
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
			refresh_hash_table_list();
		}

		call_circ(current_record1, current_record2);
	}

	report_events();

	cputime_curr = get_cpu_time();
	realtime_curr = get_real_time();

	fprintf(stdout, "[P] Mapping in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

	finalizeLoadingCompressedGenome();
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

	if (mr.type == CHIBSJ) {
		call_circ_single_split(current_record1, current_record2);
	}
	else if (mr.type == CHI2BSJ) {
		call_circ_double_split(current_record1, current_record2);
	}
}

void ProcessCirc::call_circ_single_split(Record* current_record1, Record* current_record2) {
	MatchedRead mr = *(current_record1->mr);
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

	// convert to position on contig
	gtf_parser.chrloc2conloc(mr.chr_r1, mr.spos_r1, mr.epos_r1);
	gtf_parser.chrloc2conloc(mr.chr_r2, mr.spos_r2, mr.epos_r2);

	const IntervalInfo<GeneInfo>* gene_info = gtf_parser.get_gene_overlap(mr.spos_r1, false);
	bool found = (gene_info != NULL);
	if (! found) {
		vafprintf(2, stderr, "Gene not found!\n");
		return;
	}
	vafprintf(2, stderr, "# Gene overlaps: %d\n", gene_info->seg_list.size());

	// fill MatchedMate
	MatchedMate mm_r1(mr, 1, current_record1->seq_len, r1_partial);
	MatchedMate mm_r2(mr, 2, current_record2->seq_len, !r1_partial);
	ConShift con_shift;

	check_removables(mr.spos_r1);
	RegionalHashTable* regional_ht;

	CircRes best_cr;
	best_cr.type = NF;
	for (int i = 0; i < gene_info->seg_list.size(); ++i) {
		uint32_t gene_len = gene_info->seg_list[i].end - gene_info->seg_list[i].start + 1;
		char gene_seq[gene_len + 1];
		gene_seq[gene_len] = '\0';
		pac2char(gene_info->seg_list[i].start, gene_len, gene_seq);

		regional_ht = get_hash_table(gene_info->seg_list[i], gene_seq);
		
		vafprintf(2, stderr, "R%d partial: [%d-%d]\n", (int) (!r1_partial) + 1, qspos, qepos);
		vafprintf(2, stderr, "%s\n", remain_seq);

		chaining(qspos, qepos, regional_ht, remain_seq, gene_len, gene_info->seg_list[i].start, bc1);
		if (bc1.best_chain_count <= 0)
			continue;

		// find rspos and repos for the best chain
		bool forward = (r1_partial) ? (mr.r1_forward) : (mr.r2_forward);
		int dir = (forward) ? 1 : -1;
		MatchedMate partial_mm;
		find_exact_coord(mm_r1, mm_r2, partial_mm, dir, qspos, remain_seq, remain_len, whole_seq_len, bc1.chains[0]);

		if (partial_mm.type == CONCRD) {
			con_shift = gtf_parser.get_shift(contigNum, mm_r1.spos);
			vafprintf(2, stderr, "Coordinates: [%d-%d]\n", partial_mm.spos - con_shift.shift, partial_mm.epos - con_shift.shift);
			
			print_split_mapping(current_record1->rname, mm_r1, mm_r2, partial_mm, con_shift);
			
			CircRes cr;
			int type = check_split_map(mm_r1, mm_r2, partial_mm, r1_partial, cr);
			fprintf(stderr, "%d\n", type);
			if (type < CR) {
				best_cr.type = type;
				break;
			}
			
			if (type >= CR and type <= MCR and type < best_cr.type) {
				best_cr.chr = con_shift.contig;
				best_cr.rname = current_record1->rname;
				best_cr.spos = cr.spos - con_shift.shift;
				best_cr.epos = cr.epos - con_shift.shift;
				best_cr.type = type;
				if (type == CR)
					break;	// stop processing next genes
			}
		}
	}

	if (best_cr.type >= CR and best_cr.type <= MCR)
		circ_res.push_back(best_cr);
}

void ProcessCirc::call_circ_double_split(Record* current_record1, Record* current_record2) {
	MatchedRead mr = *(current_record1->mr);
	vafprintf(2, stderr, "Double split read...\n");
	char* r1_remain_seq = (mr.r1_forward) ? current_record1->seq : current_record1->rcseq;
	char* r2_remain_seq = (mr.r2_forward) ? current_record2->seq : current_record2->rcseq;

	uint32_t r1_qspos = ((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (1) : (mr.qepos_r1 + 1);
	uint32_t r2_qspos = ((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (1) : (mr.qepos_r2 + 1);

	uint32_t r1_qepos = ((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (mr.qspos_r1 - 1) : (current_record1->seq_len);
	uint32_t r2_qepos = ((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (mr.qspos_r2 - 1) : (current_record2->seq_len);


	int r1_whole_seq_len = current_record1->seq_len;
	int r2_whole_seq_len = current_record2->seq_len;
	
	int r1_remain_len = r1_qepos - r1_qspos + 1;
	int r2_remain_len = r2_qepos - r2_qspos + 1;
	
	// if (qepos < qspos or remain_len < window_size) {	// it was fully mapped
	// 	return;
	// }

	// convert to position on contig
	gtf_parser.chrloc2conloc(mr.chr_r1, mr.spos_r1, mr.epos_r1);
	gtf_parser.chrloc2conloc(mr.chr_r2, mr.spos_r2, mr.epos_r2);

	const IntervalInfo<GeneInfo>* gene_info = gtf_parser.get_gene_overlap(mr.spos_r1, false);
	bool found = (gene_info != NULL);
	if (! found) {
		vafprintf(2, stderr, "Gene not found!\n");
		return;
	}
	vafprintf(2, stderr, "# Gene overlaps: %d\n", gene_info->seg_list.size());

	// fill MatchedMate
	MatchedMate mm_r1(mr, 1, current_record1->seq_len, true);
	MatchedMate mm_r2(mr, 2, current_record2->seq_len, true);
	ConShift con_shift;

	check_removables(mr.spos_r1);
	RegionalHashTable* regional_ht;

	CircRes best_cr;
	best_cr.type = NF;
	for (int i = 0; i < gene_info->seg_list.size(); ++i) {
		uint32_t gene_len = gene_info->seg_list[i].end - gene_info->seg_list[i].start + 1;
		char gene_seq[gene_len + 1];
		gene_seq[gene_len] = '\0';
		pac2char(gene_info->seg_list[i].start, gene_len, gene_seq);

		regional_ht = get_hash_table(gene_info->seg_list[i], gene_seq);
		
		vafprintf(2, stderr, "R1 partial: [%d-%d]\nremain: %s\n", r1_qspos, r1_qepos, r1_remain_seq);
		chaining(r1_qspos, r1_qepos, regional_ht, r1_remain_seq, gene_len, gene_info->seg_list[i].start, bc1);
		vafprintf(2, stderr, "R2 partial: [%d-%d]\nremain: %s\n", r2_qspos, r2_qepos, r2_remain_seq);
		chaining(r2_qspos, r2_qepos, regional_ht, r2_remain_seq, gene_len, gene_info->seg_list[i].start, bc2);
		
		if (bc1.best_chain_count <= 0 or bc2.best_chain_count <= 0)
			continue;

		// find rspos and repos for the best chain
		MatchedMate r1_partial_mm;
		MatchedMate r2_partial_mm;

		set_mm(bc1.chains[0], r1_qspos, r1_remain_len, mm_r1.dir, r1_partial_mm);
		set_mm(bc2.chains[0], r2_qspos, r2_remain_len, mm_r2.dir, r2_partial_mm);

		overlap_to_spos(mm_r1);
		overlap_to_spos(mm_r2);
		overlap_to_spos(r1_partial_mm);
		overlap_to_spos(r2_partial_mm);

		vector <uint32_t> common_tid;
		if (! same_transcript(mm_r1.exons_spos, mm_r2.exons_spos, r1_partial_mm.exons_spos, r2_partial_mm.exons_spos, common_tid))
			continue;

		bool success = false;
		if (bc1.chains[0].frags[0].rpos <= bc2.chains[0].frags[0].rpos) {
			success = extension.extend_both_mates(bc1.chains[0], bc2.chains[0], common_tid, r1_remain_seq, r2_remain_seq, 
						r1_qspos, r2_qspos, r1_qepos, r2_qepos, r1_partial_mm, r2_partial_mm);
		}
		else {
			success = extension.extend_both_mates(bc2.chains[0], bc1.chains[0], common_tid, r2_remain_seq, r1_remain_seq, 
						r2_qspos, r1_qspos, r2_qepos, r1_qepos, r2_partial_mm, r1_partial_mm);
		}

		if (! success)
			continue;

		if (r1_partial_mm.type == CONCRD and r2_partial_mm.type == CONCRD) {
			con_shift = gtf_parser.get_shift(contigNum, mm_r1.spos);
			vafprintf(2, stderr, "R1 Partial Coordinates: [%d-%d]\n", r1_partial_mm.spos - con_shift.shift, r1_partial_mm.epos - con_shift.shift);
			vafprintf(2, stderr, "R2 Partial Coordinates: [%d-%d]\n", r2_partial_mm.spos - con_shift.shift, r2_partial_mm.epos - con_shift.shift);
			
			print_split_mapping(current_record1->rname, mm_r1, mm_r2, r1_partial_mm, r2_partial_mm, con_shift);
			
			CircRes cr;
			int type = check_split_map(mm_r1, mm_r2, r1_partial_mm, r2_partial_mm, cr);
			fprintf(stderr, "%d\n", type);
			if (type < CR) {
				best_cr.type = type;
				break;
			}

			if (type >= CR and type <= MCR and type < best_cr.type) {
				best_cr.chr = con_shift.contig;
				best_cr.rname = current_record1->rname;
				best_cr.spos = cr.spos - con_shift.shift;
				best_cr.epos = cr.epos - con_shift.shift;
				best_cr.type = type;
				if (type == CR)
					break;	// stop processing next genes
			}
		}
	}

	if (best_cr.type >= CR and best_cr.type <= MCR)
		circ_res.push_back(best_cr);
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
							uint32_t gene_len, uint32_t shift, chain_list& bc) {
	int seq_len = qepos - qspos + 1;
	if (seq_len < window_size) {
		bc.best_chain_count = 0;
		return;
	}
	int kmer_cnt = (seq_len - window_size) / step + 1;
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
		if (fl[l].frag_count > seedLim)
			fl[l].frag_count = 0;
		l++;
	}

	kmer_cnt = l;

	chain_seeds_sorted_kbest2(qepos, fl, bc, window_size, kmer_cnt, shift);

	vafprintf(1, stderr, "Chaining score:%.4f,\t len: %lu\n", bc.chains[0].score, (unsigned long)bc.best_chain_count);

	int allowed_missed_kmers = (qepos - qspos + 1) / 20 * 3 + 1;
	vafprintf(2, stderr, "Allowed missing kmers: %d\n", allowed_missed_kmers);

	int missing;
	int least_miss = INF;
	for (int j = 0; j < bc.best_chain_count; j++) {
		missing = kmer_cnt - bc.chains[j].chain_len;
		vafprintf(2, stderr, "Actual missing: %d\n", missing);
		if (missing > least_miss)
			break;

		least_miss = missing;

		for (int i = 0; i < bc.chains[j].chain_len; i++) {
			vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, bc.chains[j].frags[i].rpos, bc.chains[j].frags[i].qpos, bc.chains[j].frags[i].len);
		}
	}
}

bool ProcessCirc::find_exact_coord(MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, 
									int dir, uint32_t qspos, char* rseq, int rlen, int whole_len, const chain_t& bc) {

	set_mm(bc, qspos, rlen, dir, partial_mm);
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

	partial_mm.middle_ed = extension.calc_middle_ed(bc, maxEd, rseq, rlen);
	// fprintf(stderr, "Middle ed: %d\n", partial_mm.middle_ed);
	if (partial_mm.middle_ed > maxEd)
		return false;

	bool extend = true;
	partial_mm.is_concord = false;
	if (bc.chain_len <= 0) {
		partial_mm.type = ORPHAN;
		partial_mm.matched_len = 0;
		return false;
	}

	bool lok;
	bool rok;
	int err = partial_mm.middle_ed;
	if (extend) {
		partial_mm.matched_len = rlen;
		lok = extension.extend_chain_left (common_tid, bc, rseq + qspos, qspos, MINLB, partial_mm, err);
		if (qspos == 0) {
			rok = extension.extend_chain_right(common_tid, bc, rseq, rlen, MAXUB, partial_mm, err);
		}
		else {
			rok = extension.extend_chain_right(common_tid, bc, rseq, whole_len, MAXUB, partial_mm, err);
		}
		update_match_mate_info(lok, rok, err, partial_mm);
	}

	return partial_mm.type == CONCRD;
} 

void ProcessCirc::refresh_hash_table_list (void) {
	for (auto it = ind2ht.begin(); it != ind2ht.end(); ++it) {
		removables.insert(it->first);
	}
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
			// fprintf(stderr, "Allocated new HT gid: %d, %s [%u-%u]\n", gid, mr.chr_r1.c_str(), mr.spos_r1, mr.spos_r2);
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

// for non-overlapping split mates
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

// for overlapping split mates
int ProcessCirc::check_split_map (MatchedMate& mm_r1_1, MatchedMate& mm_r2_1, MatchedMate& mm_r1_2, MatchedMate& mm_r2_2, CircRes& cr) {
	MatchedMate mm_r1_l = (mm_r1_1.spos <= mm_r1_2.spos) ? mm_r1_1 : mm_r1_2;
	MatchedMate mm_r1_r = (mm_r1_1.spos <= mm_r1_2.spos) ? mm_r1_2 : mm_r1_1;

	MatchedMate mm_r2_l = (mm_r2_1.spos <= mm_r2_2.spos) ? mm_r2_1 : mm_r2_2;
	MatchedMate mm_r2_r = (mm_r2_1.spos <= mm_r2_2.spos) ? mm_r2_2 : mm_r2_1;

	bool r1_regular_bsj = (mm_r1_l.qspos < mm_r1_r.qspos);
	bool r2_regular_bsj = (mm_r2_l.qspos < mm_r2_r.qspos);

	// RF or FR check
	if (r1_regular_bsj and r2_regular_bsj) {
		if (mm_r1_l.dir == 1) {
			if (mm_r1_r.spos <= mm_r2_l.spos)
				return FR;
			else if (mm_r1_l.epos >= mm_r2_r.epos)
				return RF;
		}
		if (mm_r1_l.dir == -1) {
			if (mm_r2_r.spos <= mm_r1_l.spos)
				return FR;
			else if (mm_r2_l.epos >= mm_r1_r.epos)
				return RF;
		}
	}

	// single BSJ on R2
	else if (r1_regular_bsj and !r2_regular_bsj) {
		MatchedMate full_mm = mm_r1_l;
		if (!full_mm.merge_to_right(mm_r1_r))
			return UD;
		return final_check(full_mm, mm_r2_l, mm_r2_r, cr);
	}

	// single BSJ on R1
	else if (!r1_regular_bsj and r2_regular_bsj) {
		MatchedMate full_mm = mm_r2_l;
		if (!full_mm.merge_to_right(mm_r2_r))
			return UD;
		return final_check(full_mm, mm_r1_l, mm_r1_r, cr);
	}

	// BSJ on the overlap
	else if (!r1_regular_bsj and !r2_regular_bsj) {
		if (mm_r1_l.spos == mm_r2_l.spos and mm_r1_r.epos == mm_r2_r.epos) {
			overlap_to_spos(mm_r1_l);
			overlap_to_epos(mm_r1_r);

			// if (mm_r1_l.exons_spos == NULL or mm_r1_r.exons_epos == NULL) {
			// 	cr.set_bp(mm_r1_l.spos - mm_r1_l.sclen_left, mm_r1_r.epos + mm_r1_r.sclen_right);
			// 	return MCR;
			// }

			vector <pu32i> end_tids;
			int diff;
			int ind_epos = mm_r1_r.exon_ind_epos;
			const IntervalInfo<UniqSeg>* curr_seg;
			curr_seg = gtf_parser.get_interval(ind_epos);
			while (mm_r1_r.spos < curr_seg->epos) {
				for (int i = 0; i < curr_seg->seg_list.size(); ++i) {
					diff = mm_r1_r.epos + mm_r1_r.sclen_right - curr_seg->seg_list[i].end;
					if (abs(diff) <= BPRES) {
						for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
							end_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
						}
					}
				}
				--ind_epos;
				curr_seg = gtf_parser.get_interval(ind_epos);
			}

			vector <pu32i> start_tids;
			int ind_spos = mm_r1_l.exon_ind_spos;
			curr_seg = gtf_parser.get_interval(ind_spos);
			while (mm_r1_l.epos > curr_seg->spos) {
				for (int i = 0; i < curr_seg->seg_list.size(); ++i) {
					diff = mm_r1_l.spos - mm_r1_l.sclen_left - curr_seg->seg_list[i].start;
					if (abs(diff) <= BPRES) {
						for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
							start_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
						}
					}
				}
				++ind_spos;
				curr_seg = gtf_parser.get_interval(ind_spos);
			}

			int sdiff, ediff;
			for (int i = 0; i < start_tids.size(); ++i) {
				for (int j = 0; j < end_tids.size(); ++j) {
					sdiff = start_tids[i].second;
					ediff = end_tids[j].second;
					if (start_tids[i].first == end_tids[j].first and sdiff == ediff) {
						cr.set_bp(mm_r1_l.spos - mm_r1_l.sclen_left - sdiff, mm_r1_r.epos + mm_r1_r.sclen_right - ediff);
						return CR;
					}
				}
			}

			cr.set_bp(mm_r1_l.spos - mm_r1_l.sclen_left, mm_r1_r.epos + mm_r1_r.sclen_right);
			if (start_tids.size() > 0 and end_tids.size() > 0)
				return NCR;

			return MCR;
		}
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

			// if (split_mm_left.exons_epos == NULL or split_mm_right.exons_spos == NULL) {
			// 	cr.set_bp(split_mm_right.spos - split_mm_right.sclen_left, split_mm_left.epos + split_mm_left.sclen_right);
			// 	return MCR;
			// }

			vector <pu32i> end_tids;
			int diff;
			int ind_epos = split_mm_left.exon_ind_epos;
			const IntervalInfo<UniqSeg>* curr_seg;
			curr_seg = gtf_parser.get_interval(ind_epos);
			while (split_mm_left.spos < curr_seg->epos) {
				for (int i = 0; i < curr_seg->seg_list.size(); ++i) {
					diff = split_mm_left.epos + split_mm_left.sclen_right - curr_seg->seg_list[i].end;
					
					// fprintf(stderr, "\nend transcripts: ");
					// for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j)
					// 	fprintf(stderr, "%s, ", gtf_parser.transcript_ids[contigNum][curr_seg->seg_list[i].trans_id[j]].c_str());
					// fprintf(stderr, "\n\t(epos, right_sc, exon_end, ediff): (%d, %d, %d, %d)\n\n", 
					// 				split_mm_left.epos, split_mm_left.sclen_right, curr_seg->seg_list[i].end, diff);
					
					if (abs(diff) <= BPRES) {
						for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
							end_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
						}
					}
				}
				--ind_epos;
				curr_seg = gtf_parser.get_interval(ind_epos);
			}

			vector <pu32i> start_tids;
			int ind_spos = split_mm_right.exon_ind_spos;
			curr_seg = gtf_parser.get_interval(ind_spos);
			while (split_mm_right.epos > curr_seg->spos) {
				for (int i = 0; i < curr_seg->seg_list.size(); ++i) {
					diff = split_mm_right.spos - split_mm_right.sclen_left - curr_seg->seg_list[i].start;
					
					// fprintf(stderr, "\nstart transcripts: ");
					// for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j)
					// 	fprintf(stderr, "%s, ", gtf_parser.transcript_ids[contigNum][curr_seg->seg_list[i].trans_id[j]].c_str());
					// fprintf(stderr, "\n\t(spos, left_sc, exon_beg, sdiff): (%d, %d, %d, %d)\n\n", 
					// 				split_mm_right.spos, split_mm_right.sclen_left, curr_seg->seg_list[i].start, diff);
					
					if (abs(diff) <= BPRES) {
						for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
							start_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
						}
					}
				}
				++ind_spos;
				curr_seg = gtf_parser.get_interval(ind_spos);
			}

			int sdiff, ediff;
			for (int i = 0; i < start_tids.size(); ++i) {
				for (int j = 0; j < end_tids.size(); ++j) {
					sdiff = start_tids[i].second;
					ediff = end_tids[j].second;
					// if (start_tids[i].first == end_tids[j].first)
						// fprintf(stderr, "tid: %d -> %s - (spos, left_sc, sdiff): (%d, %d, %d) - (epos, right_sc, ediff): (%d, %d, %d)\n", 
										// start_tids[i].first, gtf_parser.transcript_ids[contigNum][start_tids[i].first].c_str(), split_mm_right.spos, split_mm_right.sclen_left, sdiff, split_mm_left.epos, split_mm_left.sclen_right, ediff);
					// if (start_tids[i].first == end_tids[j].first and sdiff + ediff == 0) {
					if (start_tids[i].first == end_tids[j].first and sdiff == ediff) {
						cr.set_bp(split_mm_right.spos - split_mm_right.sclen_left - sdiff, split_mm_left.epos + split_mm_left.sclen_right - ediff);
						return CR;
					}
				}
			}

			cr.set_bp(split_mm_right.spos - split_mm_right.sclen_left, split_mm_left.epos + split_mm_left.sclen_right);
			if (start_tids.size() > 0 and end_tids.size() > 0)
				return NCR;

			return MCR;
		}
	}
	return UD;
}

void ProcessCirc::report_events (void) {
	// fprintf(stdout, "Hash table pool size: %d\n", ind2ht.size());
	// fprintf(stdout, "gid2ind size: %d\n", gid2ind.size());
	open_report_file();
	if (circ_res.size() <= 0)
		return;

	sort(circ_res.begin(), circ_res.end());
	int cnt = 1;
	CircRes last = circ_res[0];
	vector <string> rnames;
	rnames.push_back(circ_res[0].rname);
	for (int i = 1; i < circ_res.size(); ++i) {
		if (circ_res[i] == last) {
			cnt++;
			rnames.push_back(circ_res[i].rname);
		}
		else {
			if (last.type == CR or cnt > 1) {	// won't print novel events with single read support
				fprintf(report_file, "%s\t%u\t%u\t%d\t%s\t", last.chr.c_str(), last.spos, last.epos, cnt, circ_type[last.type-CR].c_str());
				for (int j = 0; j < rnames.size() - 1; ++j)
					fprintf(report_file, "%s,", rnames[j].c_str());
				fprintf(report_file, "%s\n", rnames[rnames.size() - 1].c_str());
			}
			cnt = 1;
			last = circ_res[i];
			rnames.clear();
			rnames.push_back(circ_res[i].rname);
		}
	}
	if (last.type == CR or cnt > 1) {
		fprintf(report_file, "%s\t%u\t%u\t%d\t%s\t", last.chr.c_str(), last.spos, last.epos, cnt, circ_type[last.type-CR].c_str());
		for (int j = 0; j < rnames.size() - 1; ++j)
			fprintf(report_file, "%s,", rnames[j].c_str());
		fprintf(report_file, "%s\n", rnames[rnames.size() - 1].c_str());
	}
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

	pre_contig = atoi(getRefGenomeName()) - 1;
	contigNum = pre_contig;

	cputime_curr = get_cpu_time();
	realtime_curr = get_real_time();

	fprintf(stdout, "[P] Loaded genome index successfully in %.2lf CPU sec (%.2lf real sec)\n\n", cputime_curr - cputime_start, realtime_curr - realtime_start);

	cputime_start = cputime_curr;
	realtime_start = realtime_curr;
}

void ProcessCirc::print_split_mapping (char* rname, MatchedMate& mm_r1, MatchedMate& mm_r2, 
										MatchedMate& partial_mm, ConShift& con_shift) {
	fprintf(stderr, "%s\t%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t", 
					rname, con_shift.contig.c_str(), 
					partial_mm.spos - con_shift.shift, partial_mm.epos - con_shift.shift, 
					partial_mm.qspos, partial_mm.matched_len, partial_mm.dir,
					mm_r1.spos - con_shift.shift, mm_r1.epos - con_shift.shift, 
					mm_r1.qspos, mm_r1.matched_len, mm_r1.dir,
					mm_r2.spos - con_shift.shift, mm_r2.epos - con_shift.shift, 
					mm_r2.qspos, mm_r2.matched_len, mm_r2.dir);
			
}

void ProcessCirc::print_split_mapping (char* rname, MatchedMate& mm_r1, MatchedMate& mm_r2, 
										MatchedMate& r1_partial_mm, MatchedMate& r2_partial_mm, ConShift& con_shift) {
	fprintf(stderr, "%s\t%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t", 
					rname, con_shift.contig.c_str(), 
					r1_partial_mm.spos - con_shift.shift, r1_partial_mm.epos - con_shift.shift, 
					r1_partial_mm.qspos, r1_partial_mm.matched_len, r1_partial_mm.dir,
					r2_partial_mm.spos - con_shift.shift, r2_partial_mm.epos - con_shift.shift, 
					r2_partial_mm.qspos, r2_partial_mm.matched_len, r2_partial_mm.dir,
					mm_r1.spos - con_shift.shift, mm_r1.epos - con_shift.shift, 
					mm_r1.qspos, mm_r1.matched_len, mm_r1.dir,
					mm_r2.spos - con_shift.shift, mm_r2.epos - con_shift.shift, 
					mm_r2.qspos, mm_r2.matched_len, mm_r2.dir);
			
}

int ProcessCirc::get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos) {


}

void set_mm(const chain_t& ch, uint32_t qspos, int rlen, int dir, MatchedMate& mm) {
	uint32_t spos = ch.frags[0].rpos;
	uint32_t epos = ch.frags[ch.chain_len - 1].rpos + ch.frags[ch.chain_len - 1].len - 1;
	uint32_t qepos = qspos + rlen - 1;

	mm.set(spos, epos, qspos, qepos, dir);
}