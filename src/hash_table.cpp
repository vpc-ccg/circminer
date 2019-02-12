#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "hash_table.h"

#define MAXHIT 1000

RegionalHashTable::RegionalHashTable () {
	
}

RegionalHashTable::RegionalHashTable (int ws, uint32_t gspos, uint32_t gepos) {
	init(ws, gspos, gepos);
}

RegionalHashTable::~RegionalHashTable (void) {
	free(kmer_count);

	for (int i = 0; i < size; i++) {
		if (table[i].frags != NULL)
			free(table[i].frags);
	}
	free(table);
}

void RegionalHashTable::init(int ws, uint32_t gspos, uint32_t gepos) {
	gene_spos = gspos;
	gene_epos = gepos;
	window_size = ws;
	size = 1 << (2 * ws);	// HT size = 4^ws

	table = (GIMatchedKmer*) malloc(size * sizeof(GIMatchedKmer));
	for (int i = 0; i < size; i++)
		table[i].frags = (GeneralIndex*) malloc(MAXHIT * sizeof(GeneralIndex));
	
	kmer_count = (int*) malloc(size * sizeof(int));

	memset(nuc_hval, -1, 128);
	nuc_hval['A'] = 0;
	nuc_hval['C'] = 1;
	nuc_hval['G'] = 2;
	nuc_hval['T'] = 3;
}

void RegionalHashTable::create_table (char* seq, uint32_t start, int len) {
	// memset(table, 0, size * sizeof(GIMatchedKmer));
	memset(kmer_count, 0, size * sizeof(int));

	if (len < window_size)
		return;

	// count kmers
	for (int i = 0; i <= len - window_size; i++) {
		++kmer_count[hash_val(seq + i)];
	}

	// initialize
	for (int i = 0; i < size; i++) {
		table[i].frag_count = 0;
		table[i].qpos = 0;
	}

	// fill hash table
	uint32_t loc = start;
	int hv;
	for (int i = 0; i <= len - window_size; i++) {
		hv = hash_val(seq + i);
		if (kmer_count[hv] <= MAXHIT)
			add_loc(hv, loc);
		++loc;
	}
}

// GIList* RegionalHashTable::find_hash (int hv) const {
// 	if (hv < 0 or hv >= size)
// 		return NULL;

// 	fprintf(stdout, "Hash val: %d \t CNT: %d\n", hv, table[hv].cnt);
// 	return &table[hv];
// }

GIMatchedKmer* RegionalHashTable::find_hash (int hv) const {
	if (hv < 0 or hv >= size)
		return NULL;

	return &table[hv];
}

int RegionalHashTable::hash_val (char* seq) const {
	int i = 0;
	int val = 0;

	while(i < window_size)
	{
		if (nuc_hval[seq[i]] == -1)	// N
			return -1; 
		val = (val << 2) | nuc_hval[seq[i++]]; 
	}
	return val;
}

void RegionalHashTable::add_loc (int hv, uint32_t loc) {
	if (table[hv].frag_count < MAXHIT) {
		table[hv].frags[table[hv].frag_count++].info = loc;
	}
}