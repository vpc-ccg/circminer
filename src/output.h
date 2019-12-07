#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <cinttypes>
#include <cstdio>

#define MAXTOTALTAGLEN 200

struct OptionalTags {
	int 		map_type;	//AT
	int 		ed;			//NM
	uint16_t 	junc_cnt;	//JC
	bool 		gm_compat;	//TC
	char*		all_tags;

	void to_string(void);
};

struct CommonAttr {
	char* 		qname;
	uint8_t 	mapq;
};

struct SeparateAttr {
	uint16_t 		flag;
	char* 			rname;
	uint32_t 		pos;
	char* 			cigar;
	char* 			rnext;
	uint32_t 		pnext;
	int32_t 		tlen;
	char* 			seq;
	char* 			qual;
	OptionalTags 	tags;

};

class SAMOutput {
private:
	FILE* outsam;

	CommonAttr   com_attr;
	SeparateAttr r1_attr;
	SeparateAttr r2_attr;

	uint16_t set_flag_se (MatchedMate* mm);
	uint16_t set_flag_pe (MatchedRead* mr, bool first);

	void set_output_se (Record* rec);
	void set_output_pe (Record* rec1, Record* rec2);

	void print_header (void);

public:
	SAMOutput (void);
	SAMOutput (char* sam_prefix, char* open_mode);
	~SAMOutput (void);

	void init (char* sam_prefix, char* open_mode);

	void finalize (void);

	void write_sam_rec_se (Record* rec);
	void write_pam_rec_se (Record* rec);

	void write_sam_rec_pe (Record* rec1, Record* rec2);
	void write_pam_rec_pe (Record* rec1, Record* rec2);
};

#endif //__OUTPUT_H__
