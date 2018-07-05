#ifndef __COMMON_H__
#define __COMMON_H__

#define FILE_NAME_LENGTH 1000
#define GENETHRESH 5000
#define FLGENETH 2000000
#define MINKMER 15
#define FRAGLIM 5000

#define BESTCHAINLIM 30
#define EDTH 3
#define SOFTCLIPTH 7

#include <vector>
#include <stdint.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

typedef struct fragment_t{
	uint32_t rpos;
	int32_t qpos;
	uint32_t len;

	bool operator < (const fragment_t& other) const {
		return rpos < other.rpos;
	}
} fragment_t;

typedef struct {
	//std::vector<fragment_t> frags;
	fragment_t* frags;
	uint32_t chain_len;
	float score;
} chain_t;

typedef struct {
	chain_t* chains;
	int best_chain_count;
} chain_list;

extern bool pairedEnd;

extern int kmer;
extern int verboseMode;

extern char gtfFilename[FILE_NAME_LENGTH];
extern char referenceFilename[FILE_NAME_LENGTH];
extern char fastqFilename[FILE_NAME_LENGTH];
extern char outputFilename[FILE_NAME_LENGTH];
extern char outputDir[FILE_NAME_LENGTH];

extern FILE* outputJuncFile;

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...);

#endif
