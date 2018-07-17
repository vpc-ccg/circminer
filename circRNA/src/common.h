#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

#define maxM(A,B) (((A) > (B)) ? (A) : (B))
#define minM(A,B) (((A) < (B)) ? (A) : (B))

#define FILE_NAME_LENGTH 1000
#define LINELOG	100000
#define ASCISIZE 128

#define GENETHRESH 20000
//#define GENETHRESH 5000
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

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

FILE* open_file(char* filename, char* mode);
void close_file(FILE* fp);

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...);

#endif	//__COMMON_H__
