#ifndef __COMMON_H__
#define __COMMON_H__

#if DEBUG==1
#define DEBUG_MODE
#endif

#include <stdint.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <bitset>

#include "interval_info.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

using namespace std;

#define maxM(A,B) (((A) > (B)) ? (A) : (B))
#define minM(A,B) (((A) < (B)) ? (A) : (B))

#define MAXLINESIZE 400
#define FILE_NAME_LENGTH 1000
#define LINELOG	100000
#define ASCISIZE 128
#define INF 1e9

#define MINKMER 15
#define MAXDISCRDTLEN 20000

#define EDTH 4
#define INDELTH 3
#define SOFTCLIPTH 7

#define MAXTLEN 500
#define FRAGLIM 5000
#define MAXINTRON	2000000
#define BESTCHAINLIM 30

#define LARIAT2BEGTH 1000

// ouput categories
#define CATNUM 11	// number of output categories

// the order matters:
#define CONCRD	0
#define DISCRD	1
#define CHIORF	2
#define CHIBSJ	3
#define CHIFUS	4
#define OEA2	5
#define CANDID	6
#define OEANCH	7
#define ORPHAN	8
#define NOPROC_MANYHIT 9
#define NOPROC_NOMATCH 10

//---------- Global Variables ----------\\

extern bool pairedEnd;

extern int kmer;
extern int maxReadLength;
extern int verboseMode;
extern int scanLevel;
extern int maxEd;
extern int maxSc;
extern int bandWidth;
extern int seedLim;
extern int maxTlen;
extern int maxIntronLen;
extern int threads;
extern int stage;

extern char gtfFilename[FILE_NAME_LENGTH];
extern char referenceFilename[FILE_NAME_LENGTH];
extern char fastqFilename[FILE_NAME_LENGTH];
extern char outputFilename[FILE_NAME_LENGTH];
extern char outputDir[FILE_NAME_LENGTH];

extern char* contigName;
extern int contigNum;

extern uint32_t lookup_cnt;
extern vector <bitset <DEF_CONTIG_MAX_SIZE> > near_border_bs;
extern vector <bitset <DEF_CONTIG_MAX_SIZE> > intronic_bs;

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

//---------- Structures ----------\\

struct fragment_t{
	uint32_t rpos;
	int32_t qpos;
	uint32_t len;

	bool operator < (const fragment_t& other) const {
		return rpos < other.rpos;
	}
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

typedef struct {
	fragment_t* frags;
	uint32_t chain_len;
	float score;
} chain_t;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

typedef struct {
	chain_t* chains;
	int best_chain_count;
} chain_list;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

typedef struct {
	GeneralIndex* frags;	// array of locations

	uint32_t frag_count;
	int32_t qpos;
} GIMatchedKmer;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct GeneInfo {
	uint32_t start;
	uint32_t end;

	friend ostream& operator<<(ostream& os, const GeneInfo& gi);

	bool operator < (const GeneInfo& gi) const;

	uint32_t length() { return end - start + 1; }
};

inline ostream& operator<<(ostream& os, const GeneInfo& gi);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct UniqSeg {
	uint32_t gene_id;
	uint32_t start;
	uint32_t end;
	uint32_t next_exon_beg;
	vector<uint32_t> trans_id;

	friend ostream& operator<<(ostream& os, const UniqSeg& us);

	UniqSeg() : 
			start(0), end(0), next_exon_beg(0), gene_id(0) {}
	UniqSeg(uint32_t gid, uint32_t s, uint32_t e, uint32_t n) : 
			start(s), end(e), next_exon_beg(n), gene_id(gid) {}

	UniqSeg(const UniqSeg& other);

	UniqSeg& operator = (const UniqSeg& other);
	bool operator < (const UniqSeg& r) const;
	bool operator == (const UniqSeg& r) const;

	bool same_gene(const UniqSeg& r) const;
	bool same_exon(const UniqSeg& r) const;
	bool next_exon(const UniqSeg& r) const;
};

// Temporary, just for testing --will be deleted
inline ostream& operator<<(ostream& os, const UniqSeg& us);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct MatchedRead;

struct MatchedMate {
	uint32_t 	spos;
	uint32_t 	epos;
	uint32_t 	qspos;
	uint32_t 	qepos;

	int 		right_ed;
	int 		left_ed;
	int 		middle_ed;

	int			sclen_right;
	int			sclen_left;
	int 		matched_len;
	int 		dir;
	int 		type;

	uint16_t	junc_num;

	bool 		is_concord;

	bool		left_ok;
	bool		right_ok;
	
	bool		looked_up_spos;	// intronic / inter-genic -> looked up but not found (still NULL)
	bool		looked_up_epos;	// intronic / inter-genic -> looked up but not found (still NULL)
	
	bool		looked_up_gene;

	int 		exon_ind_spos;
	int 		exon_ind_epos;

	const IntervalInfo<UniqSeg>* exons_spos;
	const IntervalInfo<UniqSeg>* exons_epos;

	const IntervalInfo<GeneInfo>* gene_info;

	MatchedMate();
	MatchedMate(const MatchedRead& mr, int r1_2, int rlen);

	void operator = (const MatchedMate& mm);

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct MatchedRead {
	// pos on reference
	uint32_t	spos_r1;
	uint32_t	spos_r2;
	uint32_t	epos_r1;
	uint32_t	epos_r2;
	
	// pos on read
	uint32_t 	qspos_r1;
	uint32_t 	qspos_r2;
	uint32_t 	qepos_r1;
	uint32_t 	qepos_r2;

	int 		mlen_r1;
	int 		mlen_r2;

	bool		r1_forward;
	bool		r2_forward;

	int 		ed_r1;
	int 		ed_r2;
	int 		type;
	int32_t 	tlen;
	uint16_t 	junc_num;
	bool		gm_compatible;
	
	string		chr_r1;
	string		chr_r2;

	MatchedRead();
	
	bool update(const MatchedMate& r1, const MatchedMate& r2, const string& chr, uint32_t shift, 
							int32_t tlen, uint16_t jun_between, bool gm_compatible, int type, bool r1_first);

	bool update_type(int type);
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct MatePair {
	float score;
	chain_t forward;
	chain_t reverse;
	vector <uint32_t> common_tid;

	MatePair() : score(-1) {}
	MatePair(const MatePair& other);

	MatePair& operator = (const MatePair& other);
	bool operator < (const MatePair& r) const;
};
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct GenRegion {
	uint32_t last_pos;	// last position on exon
	uint32_t next_pos;	// next position on next exon

	GenRegion () : last_pos(0), next_pos(0) { }
	GenRegion (uint32_t lp, uint32_t np) : last_pos(lp), next_pos(np) { }

	void set (uint32_t lp, uint32_t np);

	bool operator < (const GenRegion& r) const;
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct AllCoord {
	uint32_t rspos;
	uint32_t rlen;
	uint32_t qspos;
	uint32_t qlen;

	AllCoord (uint32_t rs, uint32_t rl, uint32_t qs, uint32_t ql) : rspos(rs), rlen(rl), qspos(qs), qlen(ql) {}

	bool operator < (const AllCoord& r) const;
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

//---------- Functions ----------\\

FILE* open_file(char* filename, char* mode);
void close_file(FILE* fp);

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...);

double get_cpu_time();
double get_real_time();

//------- Implementations --------\\

// verbose-aware fprintf
inline void vafprintf(int verbosity, FILE *stream, const char *format, ...) {
	#ifdef DEBUG_MODE
		if (verbosity > verboseMode)	return;

		va_list args;
		va_start (args, format);
		vfprintf (stream, format, args);
		va_end (args);
	#endif
}

//--------------------------------\\

#endif	//__COMMON_H__
