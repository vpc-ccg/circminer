#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

#include "interval_info.h"

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

#define INF 1e9

#define GENETHRESH 50000
#define MAXDISCRDTLEN 10000
#define MAXTLEN 500
#define MINKMER 15
#define FRAGLIM 5000
#define MAX_INTRON	2000000

#define BESTCHAINLIM 30
#define EDTH 4
#define INDELTH 3
#define SOFTCLIPTH 7

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
	uint32_t dr;
	uint32_t dl;

	uint32_t range;
	uint32_t max_end;

	bool looked_up;
	bool exonic;
	bool cross_boundry;
} JunctionDist;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

typedef struct {
	GeneralIndex* frags;	// array of locations
	JunctionDist* junc_dist;

	uint32_t frag_count;
	int32_t qpos;
} GIMatchedKmer;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct GeneInfo {
	uint32_t start;
	uint32_t end;

	friend ostream& operator<<(ostream& os, const GeneInfo& gi);

	bool operator < (const GeneInfo& gi) const {
		if (start != gi.start)
			return start < gi.start;
		return end < gi.end;
	}
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

inline ostream& operator<<(ostream& os, const GeneInfo& gi) {
	os << "[" << gi.start << "-" << gi.end << "]";
	return os;
}

struct UniqSeg {
	string gene_id;
	uint32_t start;
	uint32_t end;
	uint32_t next_exon_beg;
	uint32_t prev_exon_end;

	friend ostream& operator<<(ostream& os, const UniqSeg& us);

	UniqSeg() : 
			start(0), end(0), next_exon_beg(0), prev_exon_end(0), gene_id("") {}
	UniqSeg(const string& gid, uint32_t s, uint32_t e, uint32_t n, uint32_t p) : 
			start(s), end(e), next_exon_beg(n), prev_exon_end(p), gene_id(gid) {}

	UniqSeg(const UniqSeg& other) : start(other.start), end(other.end), next_exon_beg(other.next_exon_beg), prev_exon_end(other.prev_exon_end), gene_id(other.gene_id) {}

	UniqSeg& operator = (const UniqSeg& other) {
		if (this == &other)
			return *this;

		start 			= other.start;
		end 			= other.end;
		next_exon_beg 	= other.next_exon_beg;
		prev_exon_end 	= other.prev_exon_end;
		gene_id 		= other.gene_id;

		return *this;
	}

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

	bool same_gene(const UniqSeg& r) const {
		return (gene_id == r.gene_id);
	}

	bool same_exon(const UniqSeg& r) const {
		return (start == r.start and end == r.end);
	}

	bool next_exon(const UniqSeg& r) const {	// is this next exon of r?
		return (r.next_exon_beg == start and prev_exon_end == r.end);
	}
};

// Temporary, just for testing --will be deleted
inline ostream& operator<<(ostream& os, const UniqSeg& us) {
	os << us.prev_exon_end << " [" << us.gene_id << ": " << us.start << "-" << us.end << "] " << us.next_exon_beg;
	return os;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct UniqSegList {
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
	
	bool same_gene(const UniqSegList* r) const {
		if (r == NULL)
			return false;

		for (int i = 0; i < seg_list.size(); i++)
			for (int j = 0; j < r->seg_list.size(); j++)
				if (seg_list[i].gene_id == r->seg_list[j].gene_id)
					return true;
	
		return false;
	}
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct MatchedMate {
	uint32_t 	start_pos;
	uint32_t 	end_pos;
	uint16_t	junc_num;
	int			sclen_right;
	int			sclen_left;
	int 		matched_len;
	int 		dir;
	int 		type;
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

	MatchedMate() : type(ORPHAN), junc_num(0), sclen_right(0), sclen_left(0), left_ok(false), right_ok(false), 
					looked_up_spos(false), looked_up_epos(false), looked_up_gene(false), exons_spos(NULL), exons_epos(NULL), 
					exon_ind_spos(-1), exon_ind_epos(-1), gene_info(NULL) { }

	void operator = (const MatchedMate& mm) {
		start_pos 	= mm.start_pos;
		end_pos		= mm.end_pos;
		junc_num	= mm.junc_num;
		sclen_right	= mm.sclen_right;
		sclen_left	= mm.sclen_left;
		matched_len	= mm.matched_len;
		dir			= mm.dir;
		type		= mm.type;
		is_concord	= mm.is_concord;

		left_ok 	= mm.left_ok;
		right_ok	= mm.right_ok;

		looked_up_spos = mm.looked_up_spos;
		looked_up_epos = mm.looked_up_epos;

		exons_spos	= mm.exons_spos;
		exons_epos	= mm.exons_epos;
	}

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct MatchedRead {
	uint32_t	spos_r1;
	uint32_t	spos_r2;
	uint32_t	epos_r1;
	uint32_t	epos_r2;
	int 		mlen_r1;
	int 		mlen_r2;
	int 		type;
	int32_t 	tlen;
	uint16_t 	junc_num;
	bool		gm_compatible;
	string		chr;

	MatchedRead() : type(NOPROC_NOMATCH), tlen(INF), junc_num(0), gm_compatible(false), chr("-") { }
	
	bool update(const MatchedMate& r1, const MatchedMate& r2, const string& chr, uint32_t shift, int32_t tlen, uint16_t jun_between, bool gm_compatible, int type) {
		if (type > this->type)
			return false;

		if (type == this->type) {
			if (this->gm_compatible and !gm_compatible)
				return false;

			if ((this->gm_compatible == gm_compatible) and (this->tlen < tlen))
				return false;
		}

		this->type = type;
		this->chr = chr;

		spos_r1 = r1.start_pos - shift;
		epos_r1 = r1.end_pos - shift;
		mlen_r1 = r1.matched_len;

		spos_r2 = r2.start_pos - shift;
		epos_r2 = r2.end_pos - shift;
		mlen_r2 = r2.matched_len;

		this->tlen = tlen;
		this->junc_num = jun_between + r1.junc_num + r2.junc_num;
		this->gm_compatible = gm_compatible;

		return true;
	}

	bool update_type(int type) {
		if (type < this->type)
			this->type = type;
	}
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

typedef struct {
	float score;
	chain_t forward;
	chain_t reverse;
} MatePair;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\

struct GenRegion {
	uint32_t last_pos;	// last position on exon
	uint32_t next_pos;	// next position on next exon

	GenRegion () : last_pos(0), next_pos(0) { }
	GenRegion (uint32_t lp, uint32_t np) : last_pos(lp), next_pos(np) { }

	void set(uint32_t lp, uint32_t np) {
		last_pos = lp;
		next_pos = np;
	}
	
	bool operator < (const GenRegion& r) const {
		if (last_pos != r.last_pos)
			return (last_pos < r.last_pos);
		return (next_pos < r.next_pos);
	}
};

//---------- Global Variables ----------\\

extern bool pairedEnd;

extern int kmer;
extern int maxReadLength;
extern int verboseMode;

extern char gtfFilename[FILE_NAME_LENGTH];
extern char referenceFilename[FILE_NAME_LENGTH];
extern char fastqFilename[FILE_NAME_LENGTH];
extern char outputFilename[FILE_NAME_LENGTH];
extern char outputDir[FILE_NAME_LENGTH];

extern char* contigName;

extern uint32_t lookup_cnt;
extern uint8_t* near_border[3];

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

//---------- Functions ----------\\

FILE* open_file(char* filename, char* mode);
void close_file(FILE* fp);

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...);

//--------------------------------\\

#endif	//__COMMON_H__
