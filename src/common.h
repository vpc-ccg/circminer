#ifndef __COMMON_H__
#define __COMMON_H__

#if DEBUG==1
#define DEBUG_MODE
#endif

#include <cinttypes>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <bitset>
#include <zlib.h>

#include "interval_info.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

using namespace std;

#define maxM(A,B) (((A) > (B)) ? (A) : (B))
#define minM(A,B) (((A) < (B)) ? (A) : (B))
#define minM3(A,B,C) minM(minM(A, B), C)

#define MAXLINESIZE 600
// #define FILE_NAME_LENGTH 1000
#define LINELOG	100000
#define ASCISIZE 128
#define INF 1e9

#define MINLB 0
#define MAXUB 4294967295	//2^32 - 1

#define MINKMER 15
#define MAXDISCRDTLEN 20000

#define BPRES 5

#define EDTH 4
#define INDELTH 3
#define SOFTCLIPTH 7

#define MAXTLEN 500
#define FRAGLIM 500
#define MAXINTRON	2000000
#define BESTCHAINLIM 30

#define LARIAT2BEGTH 1000

// ouput categories
#define CATNUM 14	// number of output categories

// the order matters:
#define CONCRD	0
#define DISCRD	1
#define CHIORF	2
#define CHIBSJ	3
#define CHI2BSJ	4
#define CONGEN	5
#define CHIFUS	6
#define CONGNM	7
#define OEA2	8
#define CANDID	9
#define OEANCH	10
#define ORPHAN	11
#define NOPROC_MANYHIT 12
#define NOPROC_NOMATCH 13

// mapping file format
#define DISCARDMAPREPORT 0
#define PAMFORMAT 1
#define SAMFORMAT 2

//---------- Global Variables ----------//

extern bool indexMode;
extern bool compactIndex;
extern bool pairedEnd;
extern bool finalCleaning;
extern bool internalSort;

extern uint32_t seedLim;
extern int kmer;
extern int maxReadLength;
extern int verboseMode;
extern int scanLevel;
extern int maxEd;
extern int maxSc;
extern int bandWidth;
extern int maxTlen;
extern int maxIntronLen;
extern int maxChainLen;
extern int threadCount;
extern int stage;
extern int reportMapping;

extern char gtfFilename[FILE_NAME_LENGTH];
extern char referenceFilename[FILE_NAME_LENGTH];
extern char fastqFilename[2][FILE_NAME_LENGTH];
extern char outputFilename[FILE_NAME_LENGTH];
extern char outputDir[FILE_NAME_LENGTH];

extern char* contigName;
extern int contigNum;

extern uint32_t lookup_cnt;
extern vector <bitset <DEF_CONTIG_MAX_SIZE> > near_border_bs;
extern vector <bitset <DEF_CONTIG_MAX_SIZE> > intronic_bs;

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

extern pthread_mutex_t write_lock;
extern pthread_mutex_t pmap_lock;
extern pthread_mutex_t read_lock;


//---------- Structures ----------//

struct fragment_t{
	uint32_t rpos;
	int32_t qpos;
	uint32_t len;

	bool operator < (const fragment_t& other) const {
		return rpos < other.rpos;
	}
};

/**************************************************************************************************/

typedef struct {
	fragment_t* frags;
	uint32_t chain_len;
	float score;
} chain_t;

/**************************************************************************************************/

typedef struct {
	chain_t* chains;
	int best_chain_count;
} chain_list;

/**************************************************************************************************/

typedef struct {
	GeneralIndex* frags;	// array of locations

	uint32_t frag_count;
	int32_t qpos;
} GIMatchedKmer;

/**************************************************************************************************/

struct GeneInfo {
	uint32_t start;
	uint32_t end;
	uint32_t gene_id;

	friend ostream& operator<<(ostream& os, const GeneInfo& gi);

	bool operator < (const GeneInfo& gi) const;

	inline uint32_t length() { return end - start + 1; }
};

inline ostream& operator<<(ostream& os, const GeneInfo& gi);

/**************************************************************************************************/

struct UniqSeg {
	uint32_t start;
	uint32_t end;
	uint32_t next_exon_beg;
	uint32_t gene_id;
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

/**************************************************************************************************/

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
	uint32_t	matched_len;
	int			dir;
	int			type;

	uint16_t	junc_num;

	bool		is_concord;

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

	MatchedMate (void);
	MatchedMate (const MatchedRead& mr, int r1_2, int rlen, bool partial);
	
	void set (uint32_t rs, uint32_t re, uint32_t qs, uint32_t qe, int d);
	bool merge_to_right(const MatchedMate& rmm);

	void operator = (const MatchedMate& mm);

};

/**************************************************************************************************/

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

	uint32_t	mlen_r1;
	uint32_t 	mlen_r2;

	bool		r1_forward;
	bool		r2_forward;

	int 		ed_r1;
	int 		ed_r2;
	int 		type;
	int32_t 	tlen;
	uint16_t 	junc_num;
	bool		gm_compatible;
	int 		contig_num;

	uint64_t	genome_spos;
	
	string		chr_r1;
	string		chr_r2;

	MatchedRead(void);
	MatchedRead(MatchedRead* orig);
	
	bool update(const MatchedMate& r1, const MatchedMate& r2, const string& chr, uint32_t shift, 
							int32_t tlen, uint16_t jun_between, bool gm_compatible, int type, bool r1_first);

	bool update_type(int type);

	inline bool go_for_update(const MatchedMate& r1, const MatchedMate& r2, int32_t tlen, bool gm_compatible, int type);
};

/**************************************************************************************************/

struct MatePair {
	int type;
	float score;
	chain_t forward;
	chain_t reverse;
	vector <uint32_t> common_tid;

	MatePair() : score(-1) {}
	MatePair(const MatePair& other);

	MatePair& operator = (const MatePair& other);
	bool operator < (const MatePair& r) const;
};

/**************************************************************************************************/

typedef struct {
	string contig;
	uint32_t shift;
} ConShift;

/**************************************************************************************************/

struct GenRegion {
	uint32_t last_pos;	// last position on exon
	uint32_t next_pos;	// next position on next exon

	GenRegion () : last_pos(0), next_pos(0) { }
	GenRegion (uint32_t lp, uint32_t np) : last_pos(lp), next_pos(np) { }

	void set (uint32_t lp, uint32_t np);

	bool operator < (const GenRegion& r) const;
};

/**************************************************************************************************/

struct AllCoord {
	uint32_t rspos;
	uint32_t rlen;
	uint32_t qspos;
	uint32_t qlen;

	AllCoord (uint32_t rs, uint32_t rl, uint32_t qs, uint32_t ql) : rspos(rs), rlen(rl), qspos(qs), qlen(ql) {}

	bool operator < (const AllCoord& r) const;
};

/**************************************************************************************************/

struct CircRes {
	string chr;
	string rname;
	uint32_t spos;
	uint32_t epos;
	int type;

	void set_bp (uint32_t sp, uint32_t ep);

	bool operator < (const CircRes& r) const;
	bool operator == (const CircRes& r) const;
};

/**************************************************************************************************/

struct Record {
	char* rname;
	char* seq;
	char* rcseq;
	char* comment;
	char* qual;
	char* rqual;

	uint32_t seq_len;

	MatchedRead* mr;

	Record(void) {
		mr = new MatchedRead;
	}

	~Record(void) {
		delete mr;
	}

	Record(Record* orig);

	void init(void);
	void finalize(void);

	void deep_copy(Record* orig);

	bool operator < (const Record& r) const;
};

struct RecordStr {
	string rname;
	string seq;
	string comment;
	string qual;

	MatchedRead mr;

	RecordStr(void) {
	}

	~RecordStr(void) {
	}

	RecordStr(const RecordStr& other);
	RecordStr(Record* orig);

	void deep_copy(Record* orig);

	bool operator < (const RecordStr& r) const;
};

/**************************************************************************************************/

struct FilterArgs {
	int id;
	int kmer_size;

	GIMatchedKmer* fl;
	GIMatchedKmer* bl;

	chain_list* fbc_r1; 
	chain_list* bbc_r1;
	chain_list* fbc_r2; 
	chain_list* bbc_r2;

	FilterArgs (int k) : id(0), kmer_size(k) {}

	void set(GIMatchedKmer* f, GIMatchedKmer* b, 
			 chain_list* fc1, chain_list* bc1, 
			 chain_list* fc2, chain_list* bc2) {

		fl = f;
		bl = b;

		fbc_r1  = fc1;
		bbc_r1 = bc1;
		fbc_r2  = fc2;
		bbc_r2 = bc2;
	}
};

/**************************************************************************************************/

//---------- Functions ----------//

FILE* open_file(char* filename, char* mode);
gzFile open_gzfile(char* filename, char* mode);
void close_file(FILE* fp);
void close_gzfile(gzFile fp);

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...);

double get_cpu_time();
double get_real_time();

void mutex_lock(pthread_mutex_t* m);
void mutex_unlock(pthread_mutex_t* m);

//------- Implementations --------//

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

//--------------------------------//

#endif	//__COMMON_H__
