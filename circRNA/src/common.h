#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <boost/icl/interval_map.hpp>

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

using namespace std;

#define maxM(A,B) (((A) > (B)) ? (A) : (B))
#define minM(A,B) (((A) < (B)) ? (A) : (B))

#define FILE_NAME_LENGTH 1000
#define LINELOG	100000
#define ASCISIZE 128

//#define GENETHRESH 20000
#define MAXTLEN 500
//#define MAXTLEN 500000
#define MINKMER 15
#define FRAGLIM 5000

#define BESTCHAINLIM 30
#define EDTH 3
#define SOFTCLIPTH 7

// ouput categories
// the order matters:
#define CONCRD 0
#define CANDID 1
#define OEANCH 2
#define ORPHAN 3
#define CHIBSJ 4
#define CHIFUS 5

typedef struct fragment_t{
	uint32_t rpos;
	int32_t qpos;
	uint32_t len;

	bool operator < (const fragment_t& other) const {
		return rpos < other.rpos;
	}
} fragment_t;

typedef struct {
	fragment_t* frags;
	uint32_t chain_len;
	float score;
} chain_t;

typedef struct {
	chain_t* chains;
	int best_chain_count;
} chain_list;

typedef struct {
	GeneralIndex* frags;	// array of locations
	uint32_t frag_count;
	int32_t qpos;
} GIMatchedKmer;

typedef struct UniqSeg {
	string gene_id;
	uint32_t start;
	uint32_t end;
	uint32_t next_exon_beg;
	uint32_t prev_exon_end;

	bool operator < (const UniqSeg& r) const {
		if (start != r.start)
			return start < r.start;
		if (end != r.end)
			return end < r.end;
		if (gene_id != r.gene_id)
			return gene_id < r.gene_id;
		if (next_exon_beg != r.next_exon_beg)
			return next_exon_beg > r.next_exon_beg;		// for backward move after binary search
		return prev_exon_end < r.prev_exon_end;
	}
	
	bool operator == (const UniqSeg& r) const {
		return (start == r.start and end == r.end and gene_id == r.gene_id and next_exon_beg == r.next_exon_beg and prev_exon_end == r.prev_exon_end);
	}

} UniqSeg;

typedef struct UniqSegList {
	vector <UniqSeg> seg_list;
	
	bool operator == (const UniqSegList& r) const {
		if (seg_list.size() != r.seg_list.size())
			return false;

		for (int i = 0; i < seg_list.size(); i++)
			if (!(seg_list[i] == r.seg_list[i]))
				return false;

		return true;
	}

	UniqSegList& operator += (const UniqSeg& r) {
		seg_list.push_back(r);
		return *this;
	}
	
	UniqSegList& operator += (const UniqSegList& r) {
		for (int i = 0; i < r.seg_list.size(); i++)
			seg_list.push_back(r.seg_list[i]);
		return *this;
	}
} UniqSegList;

extern bool pairedEnd;

extern int kmer;
extern int maxReadLength;
extern int verboseMode;

extern char gtfFilename[FILE_NAME_LENGTH];
extern char referenceFilename[FILE_NAME_LENGTH];
extern char fastqFilename[FILE_NAME_LENGTH];
extern char outputFilename[FILE_NAME_LENGTH];
extern char outputDir[FILE_NAME_LENGTH];

extern FILE* outputJuncFile;

extern char* contigName;

extern uint32_t lookup_cnt;
extern bool* is_exon[3];

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

FILE* open_file(char* filename, char* mode);
void close_file(FILE* fp);

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...);

#endif	//__COMMON_H__
