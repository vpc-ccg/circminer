#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <stdint.h>
#include <cstdio>

struct CommonAttr {
	char* 		qname;
	uint8_t 	mapq;
};

struct SeparateAttr {
	uint16_t 	flag;
	char* 		rname;
	uint32_t 	pos;
	char* 		cigar;
	char* 		rnext;
	uint32_t 	pnext;
	int32_t 	tlen;
	char* 		seq;
	char* 		qual;

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

public:
	SAMOutput (void);
	SAMOutput (char* sam_prefix);
	~SAMOutput (void);

	void init(char* sam_prefix);

	void write_sam_rec_se (Record* rec);
	void write_sam_rec_pe (Record* rec1, Record* rec2);
};

#endif //__OUTPUT_H__
