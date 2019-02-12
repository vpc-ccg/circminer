#ifndef __HASHTABLE_H__
#define __HASHTABLE_H__

#include <vector>
#include "common.h"

typedef struct {
	GeneralIndex* locs;
	int cnt;
} GIList;

class RegionalHashTable {
private:
	int size;
	int window_size;
	GIMatchedKmer* table;
	int* kmer_count;

	char nuc_hval[128];

public:
	uint32_t gene_spos;
	uint32_t gene_epos;

	RegionalHashTable ();
	RegionalHashTable (int ws, uint32_t gspos, uint32_t gepos);
	~RegionalHashTable (void);

	void init(int ws, uint32_t gspos, uint32_t gepos);
	void create_table (char* seq, uint32_t start, int len);
	int hash_val(char* seq) const;
	void add_loc (int hv, uint32_t loc);

	GIMatchedKmer* find_hash (int hv) const;

};

#endif //__HASHTABLE_H__
