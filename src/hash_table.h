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
	GIList* table;
	int* kmer_count;

	char nuc_hval[128];

public:
	RegionalHashTable (int ws);
	~RegionalHashTable (void);

	void create_table (char* seq, uint32_t start, int len);
	int hash_val(char* seq);
	void add_loc (int hv, uint32_t loc);

	GeneralIndex* find_hash (int hv);

};

#endif
