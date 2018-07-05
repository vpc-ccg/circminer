#ifndef __LINKED_LIST__
#define __LINKED_LIST__

#include "common.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

class MatchedKmer {
public:
	MatchedKmer() 
				: frag_count(0), next(NULL), prev(NULL) {
					frags = (fragment_t*) malloc(sizeof(fragment_t) * FRAGLIM);
				}

	MatchedKmer(fragment_t* fl, uint32_t fc, MatchedKmer* n = NULL, MatchedKmer* p = NULL) 
				: frags(fl), frag_count(fc), next(n), prev(p) {}

	~MatchedKmer(void) {
		free(frags);
	}

	fragment_t* frags;	// array of fragments
	uint32_t frag_count;
	MatchedKmer* next;
	MatchedKmer* prev;
};

typedef struct {
	GeneralIndex* frags;	// array of locations
	uint32_t frag_count;
	int32_t qpos;
} GIMatchedKmer;

class FragmentList {
private:
	MatchedKmer* head;
	MatchedKmer* tail;
	int size;
	int max_frag_size;

public:
	FragmentList(void);
	FragmentList(int flen);
	~FragmentList(void);

	MatchedKmer* get_head(void);
	MatchedKmer* get_tail(void);

	int get_size(void);
	int get_max_frag_size(void);

	void add_front(fragment_t* frags, uint32_t frag_count);
	void add_front(MatchedKmer* new_mk);

	void add_back(fragment_t* frags, uint32_t frag_count);
	void add_back(MatchedKmer* new_mk);

	void sort_lists(void);

	void print(void);
};

#endif
