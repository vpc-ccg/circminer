#ifndef __PROCESSCIRC_H__
#define __PROCESSCIRC_H__

#include <cstdio>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include "common.h"
#include "hash_table.h"
#include "fastq_parser.h"
#include "extend.h"

#define FR 0	// Forward Reverse
#define RF 1	// Reverse Forward
#define CR 20	// Circular RNA
#define NCR 21	// Novel Circular RNA
#define MCR 22	// Missed Circular RNA
#define UD 30	// Undefined
#define NF 40	// Not found

using namespace std;

class ProcessCirc {
private:
	char fq_file1[FILE_NAME_LENGTH];
	char fq_file2[FILE_NAME_LENGTH];

	FILE* report_file;
	FILE* candid_file;

	int window_size;
	int step;
	
	MatchedRead mr;
	chain_list bc1;
	chain_list bc2;

	int pre_contig;
	string pre_chr;

	unordered_map <uint32_t, RegionalHashTable*> ind2ht;
	unordered_set <uint32_t> removables;
	unordered_map <uint32_t, uint32_t> gid2ind;
	vector <uint32_t> gids;

	vector <CircRes> circ_res;
	vector <string> circ_type;

	TransExtension extension;

public:
	ProcessCirc (int last_round_num, int ws);
	~ProcessCirc (void);

	void open_report_file (void);
	void open_candid_file (void);
	void load_genome (void);
	void sort_fq (char* fqname);

	void do_process (void);
	void call_circ (Record* current_record1, Record* current_record2);

	void call_circ_single_split(Record* current_record1, Record* current_record2);
	void call_circ_double_split(Record* current_record1, Record* current_record2);

	void binning (uint32_t qspos, uint32_t qepos, RegionalHashTable* regional_ht, char* remain_seq, uint32_t gene_len);
	void chaining (uint32_t qspos, uint32_t qepos, RegionalHashTable* regional_ht, char* remain_seq, uint32_t gene_len, uint32_t shift, chain_list& bc);

	bool find_exact_coord (MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, 
							int dir, uint32_t qspos, char* rseq, int rlen, int whole_len, const chain_t& bc);

	void refresh_hash_table_list (void);
	void check_removables (uint32_t rspos);
	RegionalHashTable* get_hash_table (const GeneInfo& gene_info);

	int check_split_map (MatchedMate& mm_r1, MatchedMate& mm_r2, MatchedMate& partial_mm, bool r1_partial, CircRes& cr);
	int check_split_map (MatchedMate& mm_r1_1, MatchedMate& mm_r2_1, MatchedMate& mm_r1_2, MatchedMate& mm_r2_2, CircRes& cr);
	int final_check (MatchedMate& full_mm, MatchedMate& split_mm_left, MatchedMate& split_mm_right, CircRes& cr);

	void report_events (void);

	void print_split_mapping (char* rname, MatchedMate& mm_r1, MatchedMate& mm_r2, 
										MatchedMate& partial_mm, ConShift& con_shift);

	void print_split_mapping (char* rname, MatchedMate& mm_r1, MatchedMate& mm_r2, 
										MatchedMate& r1_partial_mm, MatchedMate& r2_partial_mm, ConShift& con_shift);

};

#endif
