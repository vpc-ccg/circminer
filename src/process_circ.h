#ifndef __PROCESSCIRC_H__
#define __PROCESSCIRC_H__

#include <cstdio>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include "common.h"
#include "hash_table.h"
#include "fastq_parser.h"

#define FR 0	// Forward Reverse
#define RF 1	// Reverse Forward
#define CR 20	// Circular RNA
#define NCR 21	// Novel Circular RNA
#define MCR 22	// Missed Circular RNA
#define UD 3	// Undefined

#define BPRES 4

using namespace std;

class ProcessCirc {
private:
	char fq_file1[FILE_NAME_LENGTH];
	char fq_file2[FILE_NAME_LENGTH];

	FILE* report_file;

	int window_size;
	int step;
	
	MatchedRead mr;
	chain_list bc;

	int pre_contig;

	unordered_map <uint32_t, RegionalHashTable*> ind2ht;
	unordered_set <uint32_t> removables;
	unordered_map <uint32_t, uint32_t> gid2ind;
	vector <uint32_t> gids;

	vector <CircRes> circ_res;

public:
	ProcessCirc (int last_round_num, int ws);
	~ProcessCirc (void);

	void open_report_file (void);
	void load_genome (void);
	void sort_fq (char* fqname);

	void do_process (void);
	void call_circ (Record* current_record1, Record* current_record2);

	void binning (uint32_t qspos, uint32_t qepos, RegionalHashTable* regional_ht, char* remain_seq, uint32_t gene_len);
	void chaining (uint32_t qspos, uint32_t qepos, RegionalHashTable* regional_ht, char* remain_seq, uint32_t gene_len, uint32_t shift);

	bool find_exact_coord (MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, 
							int dir, uint32_t qspos, char* rseq, int rlen, int whole_len);

	void refresh_hash_table_list (void);
	void check_removables (uint32_t rspos);
	RegionalHashTable* get_hash_table (const GeneInfo& gene_info, char* gene_seq);

	int check_split_map (MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, bool r1_partial, CircRes& cr);
	int final_check (MatchedMate& full_mm, MatchedMate& split_mm_left, MatchedMate& split_mm_right, CircRes& cr);

	void report_events (void);

	int get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos);

	void print_split_mapping (char* rname, MatchedMate& mm_r1, MatchedMate& mm_r2, 
										MatchedMate& partial_mm, ConShift& con_shift);

};

#endif
