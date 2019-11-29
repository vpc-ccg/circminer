#include <cstring>

#include "common.h"
#include "output.h"

extern "C" {
#include "mrsfast/HashTable.h"
}

#define PAIRED (1UL << 0)
#define PROPER (1UL << 1)
#define RUNMAP (1UL << 2)
#define MUNMAP (1UL << 3)
#define RREVER (1UL << 4)
#define MREVER (1UL << 5)
#define FIPAIR (1UL << 6)
#define SIPAIR (1UL << 7)

SAMOutput::SAMOutput () {
	outsam = NULL;
}

SAMOutput::SAMOutput (char* sam_prefix) {
	init(sam_prefix);
}

SAMOutput::~SAMOutput (void) {
	close_file(outsam);

	free(r1_attr.rname);
	free(r1_attr.cigar);
	free(r1_attr.rnext);

	free(r2_attr.rname);
	free(r2_attr.cigar);
	free(r2_attr.rnext);
}

void SAMOutput::init (char* sam_prefix) {
	char fname[FILE_NAME_LENGTH];
	sprintf(fname, "%s.mapping.sam", sam_prefix);
	outsam = open_file(fname, "w");

	print_header();

	r1_attr.rname = (char*) malloc(MAXLINESIZE);
	r1_attr.cigar = (char*) malloc(MAXLINESIZE);
	r1_attr.rnext = (char*) malloc(MAXLINESIZE);

	r2_attr.rname = (char*) malloc(MAXLINESIZE);
	r2_attr.cigar = (char*) malloc(MAXLINESIZE);
	r2_attr.rnext = (char*) malloc(MAXLINESIZE);
}

uint16_t SAMOutput::set_flag_se (MatchedMate* mm) {
	uint16_t flag = 0;

	if (mm->type != CONCRD) {
		flag |= RUNMAP;
		return flag;
	}

	if (mm->dir < 0)
		flag |= RREVER;

	return flag;
}

void SAMOutput::set_output_se (Record* rec) {
	com_attr.qname = rec->rname;
	com_attr.mapq  = 255;

	r1_attr.flag = set_flag_pe(rec->mr, true);

	strcpy(r1_attr.cigar, "*");
	strcpy(r1_attr.rname, rec->mr->chr_r1.c_str());
	strcpy(r1_attr.rnext, "*");

	r1_attr.pos   = rec->mr->spos_r1 + 1;
	r1_attr.pnext = 0;
	r1_attr.tlen  = 0;

	if (r1_attr.flag & RREVER) {
		r1_attr.seq  = rec->rcseq;
		r1_attr.qual = rec->rqual;
	}
	else {
		r1_attr.seq  = rec->seq;
		r1_attr.qual = rec->qual;
	}
}

uint16_t SAMOutput::set_flag_pe (MatchedRead* mr, bool first) {
	uint16_t flag = 0;
	flag |= PAIRED;
	
	if (mr->type == CONCRD)
		flag |= PROPER;

	if (!(mr->type <= CHIORF or mr->type == CONGEN or mr->type == CONGNM)) {
		flag |= RUNMAP;
		flag |= MUNMAP;
	}

	if (first) {
		if (!(flag & RUNMAP) and !mr->r1_forward)
			flag |= RREVER;

		if (!(flag & MUNMAP) and !mr->r2_forward)
			flag |= MREVER;

		flag |= FIPAIR;
	}
	else {
		if (!(flag & MUNMAP) and !mr->r1_forward)
			flag |= MREVER;

		if (!(flag & RUNMAP) and !mr->r2_forward)
			flag |= RREVER;

		flag |= SIPAIR;
	}

	return flag;
}

void SAMOutput::set_output_pe (Record* rec1, Record* rec2) {
	com_attr.qname = rec1->rname;
	com_attr.mapq  = 255;

	r1_attr.flag = set_flag_pe(rec1->mr, true);
	r2_attr.flag = set_flag_pe(rec1->mr, false);

	strcpy(r1_attr.cigar, "*");
	strcpy(r2_attr.cigar, "*");

	// tlen
	if (rec1->mr->spos_r1 < rec1->mr->spos_r2) {
		r1_attr.tlen = rec1->mr->tlen;
		r2_attr.tlen = rec1->mr->tlen * -1;
	}
	else {
		r1_attr.tlen = rec1->mr->tlen * -1;
		r2_attr.tlen = rec1->mr->tlen;
	}
	

	if (r1_attr.flag & RUNMAP) {
		strcpy(r1_attr.rname, "*");
		strcpy(r2_attr.rnext, "*");
		r1_attr.pos   = 0;
		r2_attr.pnext = 0;
		r1_attr.tlen  = 0;
		r2_attr.tlen  = 0;
	}
	else {
		strcpy(r1_attr.rname, rec1->mr->chr_r1.c_str());
		strcpy(r2_attr.rnext, (rec1->mr->chr_r1 == rec1->mr->chr_r2)?
								"=" : rec1->mr->chr_r1.c_str());

		r1_attr.pos   = rec1->mr->spos_r1 + 1;
		r2_attr.pnext = rec1->mr->spos_r1 + 1;
	}
	if (r2_attr.flag & RUNMAP) {
		strcpy(r2_attr.rname, "*");
		strcpy(r1_attr.rnext, "*");
		r2_attr.pos   = 0;
		r1_attr.pnext = 0;
		r1_attr.tlen  = 0;
		r2_attr.tlen  = 0;
	}
	else {
		strcpy(r2_attr.rname, rec1->mr->chr_r2.c_str());
		strcpy(r1_attr.rnext, (rec1->mr->chr_r1 == rec1->mr->chr_r2)?
								"=" : rec1->mr->chr_r2.c_str());

		r2_attr.pos   = rec1->mr->spos_r2 + 1;
		r1_attr.pnext = rec1->mr->spos_r2 + 1;
	}

	if (r1_attr.flag & RREVER) {
		r1_attr.seq  = rec1->rcseq;
		r1_attr.qual = rec1->rqual;
	}
	else {
		r1_attr.seq  = rec1->seq;
		r1_attr.qual = rec1->qual;
	}

	if (r2_attr.flag & RREVER) {
		r2_attr.seq  = rec2->rcseq;
		r2_attr.qual = rec2->rqual;
	}
	else {
		r2_attr.seq  = rec2->seq;
		r2_attr.qual = rec2->qual;
	}
}

void SAMOutput::write_sam_rec_se (Record* rec) {
	set_output_se(rec);
	fprintf(outsam, "%s\t%u\t%s\t%u\t%u\t%s\t%s\t%u\t%u\t%s\t%s\n", 
			com_attr.qname, r1_attr.flag, r1_attr.rname, r1_attr.pos, com_attr.mapq,
			r1_attr.cigar, r1_attr.rnext, r1_attr.pnext, r1_attr.tlen, r1_attr.seq, 
			r1_attr.qual);
}

void SAMOutput::write_sam_rec_pe (Record* rec1, Record* rec2) {
	set_output_pe(rec1, rec2);
	fprintf(outsam, "%s\t%u\t%s\t%u\t%u\t%s\t%s\t%u\t%u\t%s\t%s\n", 
			com_attr.qname, r1_attr.flag, r1_attr.rname, r1_attr.pos, com_attr.mapq,
			r1_attr.cigar, r1_attr.rnext, r1_attr.pnext, r1_attr.tlen, r1_attr.seq, 
			r1_attr.qual);

	fprintf(outsam, "%s\t%u\t%s\t%u\t%u\t%s\t%s\t%u\t%u\t%s\t%s\n", 
			com_attr.qname, r2_attr.flag, r2_attr.rname, r2_attr.pos, com_attr.mapq,
			r2_attr.cigar, r2_attr.rnext, r2_attr.pnext, r2_attr.tlen, r2_attr.seq, 
			r2_attr.qual);
}

void SAMOutput::print_header (void) {
	fprintf(outsam, "@HD\tVN:1.4\tSO:unsorted\n");
	int chr_cnt = getChrCnt();
	char** chr_names = getChrNames();
	int* chr_lens = getChrLens();
	for (int i = 0; i < chr_cnt; ++i) {
		fprintf(outsam, "@SQ\tSN:%s\tLN:%d\n", chr_names[i], chr_lens[i]);
	}
}
