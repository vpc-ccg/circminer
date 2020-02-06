#include <sys/time.h>
#include <sys/resource.h>

#include "common.h"

FILE* open_file(char* filename, char* mode) {
	FILE* fp;
	fp = fopen(filename, mode);
	if (fp == NULL) {
		fprintf(stderr, "Error: Could not open file %s\n", filename);
		exit(1);
	}
	return fp;
}

gzFile open_gzfile(char* filename, char* mode) {
	gzFile gzfp;
	gzfp = gzopen(filename, mode);
	if (gzfp == Z_NULL) {
		fprintf(stderr, "Error: Could not open gzip file %s\n", filename);
		exit(1);
	}
	return gzfp;
}

void close_file(FILE* fp) {
	if (fp != NULL)
		fclose(fp);
}

void close_gzfile(gzFile fp) {
	if (fp != Z_NULL)
		gzclose(fp);
}

double get_cpu_time() {
	struct rusage t;
	getrusage(RUSAGE_SELF, &t);
	return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1000000.0 + t.ru_stime.tv_sec + t.ru_stime.tv_usec / 1000000.0;
}

double get_real_time() {
	struct timeval t;
	struct timezone tz;
	gettimeofday(&t, &tz);
	return t.tv_sec + t.tv_usec / 1000000.0;
}

void mutex_lock(pthread_mutex_t* m) {
	if (threadCount > 1)
		pthread_mutex_lock(m);
}

void mutex_unlock(pthread_mutex_t* m) {
	if (threadCount > 1)
		pthread_mutex_unlock(m);
}

//---------- Structures ----------//

bool GeneInfo::operator < (const GeneInfo& gi) const {
	if (start != gi.start)
		return start < gi.start;
	return end < gi.end;
}

inline ostream& operator<<(ostream& os, const GeneInfo& gi) {
	os << "[" << gi.start << "-" << gi.end << "]";
	return os;
}

/**************************************************************************************************/

UniqSeg::UniqSeg(const UniqSeg& other) : start(other.start), end(other.end), next_exon_beg(other.next_exon_beg), gene_id(other.gene_id) {
	trans_id.clear();
	for (unsigned int i = 0; i < other.trans_id.size(); i++)
		trans_id.push_back(other.trans_id[i]);
}

UniqSeg& UniqSeg::operator = (const UniqSeg& other) {
	if (this == &other)
		return *this;

	start 			= other.start;
	end 			= other.end;
	next_exon_beg 	= other.next_exon_beg;
	gene_id 		= other.gene_id;
	
	trans_id.clear();
	for (unsigned int i = 0; i < other.trans_id.size(); i++)
		trans_id.push_back(other.trans_id[i]);

	return *this;
}

bool UniqSeg::operator < (const UniqSeg& r) const {
	if (start != r.start)
		return start < r.start;
	if (end != r.end)
		return end < r.end;
	if (gene_id != r.gene_id)
		return gene_id < r.gene_id;
	return next_exon_beg > r.next_exon_beg;		// for backward move after binary search
}

bool UniqSeg::operator == (const UniqSeg& r) const {
	return (start == r.start and end == r.end and gene_id == r.gene_id and next_exon_beg == r.next_exon_beg);
}

bool UniqSeg::same_gene(const UniqSeg& r) const {
	return (gene_id == r.gene_id);
}

bool UniqSeg::same_exon(const UniqSeg& r) const {
	return (start == r.start and end == r.end);
}

bool UniqSeg::next_exon(const UniqSeg& r) const {	// is this next exon of r?
	return (r.next_exon_beg == start);
}

// Temporary, just for testing --will be deleted
inline ostream& operator<<(ostream& os, const UniqSeg& us) {
	os << "(";
	for (unsigned int i = 1; i < us.trans_id.size(); i++)
		os << us.trans_id[i] << ", ";
	os << ") " << " [" << us.gene_id << ": " << us.start << "-" << us.end << "] " << us.next_exon_beg;
	return os;
}

/**************************************************************************************************/

MatchedMate::MatchedMate(void) : right_ed(maxEd+1), left_ed(maxEd+1), middle_ed(maxEd+1), 
					sclen_right(0), sclen_left(0), type(ORPHAN), junc_num(0), 
					left_ok(false), right_ok(false), 
					looked_up_spos(false), looked_up_epos(false), looked_up_gene(false), 
					exon_ind_spos(-1), exon_ind_epos(-1), exons_spos(NULL), exons_epos(NULL), 
					gene_info(NULL) { }

void MatchedMate::set (uint32_t rs, uint32_t re, uint32_t qs, uint32_t qe, int d) {
	spos = rs; 
	epos = re;
	qspos = qs;
	qepos = qe; 
	matched_len = (qe + 1 >= qs) ? (qe - qs + 1) : 0;
	dir = d;
}

bool MatchedMate::merge_to_right(const MatchedMate& rmm) {
	if (dir != rmm.dir)
		return false;

	epos = rmm.epos;
	qepos = rmm.qepos;

	middle_ed += right_ed + rmm.left_ed;
	right_ed = rmm.right_ed;

	matched_len += rmm.matched_len + sclen_right + rmm.sclen_left;

	middle_ed += sclen_right + rmm.sclen_left;

	sclen_right = rmm.sclen_right;

	// junc_num += rmm.junc_num;
	right_ok = rmm.right_ok;

	looked_up_epos = rmm.looked_up_epos;
	exon_ind_epos = rmm.exon_ind_epos;

	if (left_ed + middle_ed + right_ed > maxEd) {
		return false;
	}
	else 
		return true;
}

// r1_2 = 1 or 2
MatchedMate::MatchedMate(const MatchedRead& mr, int r1_2, int rlen, bool partial) : 
					looked_up_spos(false), looked_up_epos(false), looked_up_gene(false), 
					exon_ind_spos(-1), exon_ind_epos(-1), 
					exons_spos(NULL), exons_epos(NULL), gene_info(NULL) {

	type = mr.type;
	right_ed = 0;
	left_ed = 0;

	if (r1_2 == 1) {
		spos = mr.spos_r1;
		epos = mr.epos_r1;
		qspos = mr.qspos_r1;
		qepos = mr.qepos_r1;
		
		middle_ed = mr.ed_r1;

		matched_len = mr.mlen_r1;
		dir = (mr.r1_forward) ? 1 : -1;
	}
	else {
		spos = mr.spos_r2;
		epos = mr.epos_r2;
		qspos = mr.qspos_r2;
		qepos = mr.qepos_r2;
		
		middle_ed = mr.ed_r2;

		matched_len = mr.mlen_r2;
		dir = (mr.r2_forward) ? 1 : -1;
	}

	if (partial) {
		if ((qspos - 1) > (rlen - qepos)) { // right-side matched
			sclen_left = 0;
			sclen_right = rlen - qepos;
		}
		else { // left-side matched
			sclen_left = qspos - 1;
			sclen_right = 0;
		}
	}
	else {
		sclen_left = qspos - 1;
		sclen_right = rlen - qepos;
	}
}

void MatchedMate::operator = (const MatchedMate& mm) {
	spos 		= mm.spos;
	epos		= mm.epos;
	qspos 		= mm.qspos;
	qepos		= mm.qepos;

	right_ed	= mm.right_ed;
	left_ed		= mm.left_ed;
	middle_ed	= mm.middle_ed;

	sclen_right	= mm.sclen_right;
	sclen_left	= mm.sclen_left;
	matched_len	= mm.matched_len;
	dir			= mm.dir;
	type		= mm.type;
	junc_num	= mm.junc_num;
	is_concord	= mm.is_concord;

	left_ok 	= mm.left_ok;
	right_ok	= mm.right_ok;

	looked_up_spos = mm.looked_up_spos;
	looked_up_epos = mm.looked_up_epos;
	looked_up_gene = mm.looked_up_gene;

	exons_spos	= mm.exons_spos;
	exons_epos	= mm.exons_epos;
}

/**************************************************************************************************/

MatchedRead::MatchedRead(void) : 
		r1_forward(true), r2_forward(true), ed_r1(maxEd+1), ed_r2(maxEd+1), 
		type(NOPROC_NOMATCH), tlen(INF), junc_num(0), gm_compatible(false), 
		contig_num(0), genome_spos(0), chr_r1("-"), chr_r2("-") 
{ }

MatchedRead::MatchedRead(MatchedRead* orig) : 
		spos_r1(orig->spos_r1), spos_r2(orig->spos_r2), 
		epos_r1(orig->epos_r1), epos_r2(orig->epos_r2),
		qspos_r1(orig->qspos_r1), qspos_r2(orig->qspos_r2), 
		qepos_r1(orig->qepos_r1), qepos_r2(orig->qepos_r2),
		mlen_r1(orig->mlen_r1), mlen_r2(orig->mlen_r2),
		r1_forward(orig->r1_forward), r2_forward(orig->r2_forward), 
		ed_r1(orig->ed_r1), ed_r2(orig->ed_r2), 
		type(orig->type), tlen(orig->tlen), junc_num(orig->junc_num), 
		gm_compatible(orig->gm_compatible), contig_num(orig->contig_num), 
		genome_spos(orig->genome_spos), chr_r1(orig->chr_r1), chr_r2(orig->chr_r2) 
{ }
// assuming r1 and r2 are on same chr
bool MatchedRead::update(const MatchedMate& r1, const MatchedMate& r2, const string& chr, uint32_t shift, int32_t tlen, 
							uint16_t jun_between, bool gm_compatible, int type, bool r1_first) {

	if (!go_for_update(r1, r2, tlen, gm_compatible, type))
		return false;

	//int edit_dist = r1.left_ed + r1.middle_ed + r1.right_ed + r2.left_ed + r2.middle_ed + r2.right_ed;

	this->type = type;
	this->chr_r1 = chr;
	this->chr_r2 = chr;

	if (r1_first) {
		spos_r1 = r1.spos - shift;
		epos_r1 = r1.epos - shift;

		qspos_r1 = r1.qspos;
		qepos_r1 = r1.qepos;

		mlen_r1 = r1.matched_len;
		ed_r1 = r1.left_ed + r1.middle_ed + r1.right_ed;


		spos_r2 = r2.spos - shift;
		epos_r2 = r2.epos - shift;

		qspos_r2 = r2.qspos;
		qepos_r2 = r2.qepos;

		mlen_r2 = r2.matched_len;
		ed_r2 = r2.left_ed + r2.middle_ed + r2.right_ed;

		r1_forward = r1.dir > 0;
		r2_forward = r2.dir > 0;
	}
	else {
		spos_r1 = r2.spos - shift;
		epos_r1 = r2.epos - shift;

		qspos_r1 = r2.qspos;
		qepos_r1 = r2.qepos;

		mlen_r1 = r2.matched_len;
		ed_r1 = r2.left_ed + r2.middle_ed + r2.right_ed;


		spos_r2 = r1.spos - shift;
		epos_r2 = r1.epos - shift;

		qspos_r2 = r1.qspos;
		qepos_r2 = r1.qepos;

		mlen_r2 = r1.matched_len;
		ed_r2 = r1.left_ed + r1.middle_ed + r1.right_ed;

		r1_forward = r2.dir > 0;
		r2_forward = r1.dir > 0;
	}

	this->tlen = tlen;
	this->junc_num = jun_between + r1.junc_num + r2.junc_num;
	this->gm_compatible = gm_compatible;

	this->contig_num = contigNum;

	return true;
}

bool MatchedRead::update_type(int type) {
	if (type < this->type) {
		this->type = type;
		return true;
	}
	return false;
}


inline bool MatchedRead::go_for_update(const MatchedMate& r1, const MatchedMate& r2, int32_t tlen, bool gm_compatible, int type) {
// 	if (! ((type == CHIBSJ and this->type == CHI2BSJ) or (type == CHI2BSJ and this->type == CHIBSJ)) ) {
	if (type < this->type)
		return true;
	if (type > this->type)
		return false;

	if (gm_compatible and !this->gm_compatible)
		return true;
	if (!gm_compatible and this->gm_compatible)
		return false;
// 	}

	if (type < CHIBSJ) {
		int edit_dist = r1.left_ed + r1.middle_ed + r1.right_ed + r2.left_ed + r2.middle_ed + r2.right_ed;
		if (this->ed_r1 + this->ed_r2 > edit_dist)
			return true;
		if (this->ed_r1 + this->ed_r2 < edit_dist)
			return false;

		if (this->tlen > tlen)
			return true;
		if (this->tlen < tlen)
			return false;

		if (this->mlen_r1 + this->mlen_r2 < r1.matched_len + r2.matched_len)
			return true;
		if (this->mlen_r1 + this->mlen_r2 > r1.matched_len + r2.matched_len)
			return false;
	}

	else {
		if (this->mlen_r1 + this->mlen_r2 < r1.matched_len + r2.matched_len)
			return true;
		if (this->mlen_r1 + this->mlen_r2 > r1.matched_len + r2.matched_len)
			return false;

		int edit_dist = r1.left_ed + r1.middle_ed + r1.right_ed + r2.left_ed + r2.middle_ed + r2.right_ed;
		if (this->ed_r1 + this->ed_r2 > edit_dist)
			return true;
		if (this->ed_r1 + this->ed_r2 < edit_dist)
			return false;
	}
	return false;
}

/**************************************************************************************************/

MatePair::MatePair(const MatePair& other) : type(other.type), score(other.score), forward(other.forward), reverse(other.reverse) {
	common_tid.clear();
	for (unsigned int i = 0; i < other.common_tid.size(); i++)
		common_tid.push_back(other.common_tid[i]);
}

MatePair& MatePair::operator = (const MatePair& other) {
	type = other.type;
	score = other.score;
	forward = other.forward;
	reverse = other.reverse;

	common_tid.clear();
	for (unsigned int i = 0; i < other.common_tid.size(); i++)
		common_tid.push_back(other.common_tid[i]);

	return *this;
}

bool MatePair::operator < (const MatePair& r) const {
	if (type != r.type)
		return type < r.type;
	return score > r.score;
}

/**************************************************************************************************/

void GenRegion::set (uint32_t lp, uint32_t np) {
	last_pos = lp;
	next_pos = np;
}

bool GenRegion::operator < (const GenRegion& r) const {
	if (last_pos != r.last_pos)
		return (last_pos < r.last_pos);
	return (next_pos < r.next_pos);
}

/**************************************************************************************************/

bool AllCoord::operator < (const AllCoord& r) const {
	if (rspos != r.rspos)
		return rspos < r.rspos;
	if (qspos != r.qspos)
		return qspos < r.qspos;
	if (rlen != r.rlen)
		return rlen < r.rlen;
	return qlen < r.qlen;
}

/**************************************************************************************************/

void CircRes::set_bp (uint32_t sp, uint32_t ep, const string& ssignal, const string& esignal, char* sbref, char* ebref) {
	spos = sp;
	epos = ep;

	start_signal = ssignal;
	end_signal = esignal;

	start_bp_ref = sbref;
	end_bp_ref = ebref;
}

bool CircRes::operator < (const CircRes& r) const {
	if (chr != r.chr) 
		return chr < r.chr;
	if (spos != r.spos)
		return spos < r.spos;
	if (epos != r.epos)
		return epos < r.epos;
	return type < r.type;
}

bool CircRes::operator == (const CircRes& r) const {
	// return (chr == r.chr) and (abs(int(spos - r.spos)) <= BPRES) and (abs(int(epos - r.epos)) <= BPRES);
	return (chr == r.chr) and (spos == r.spos) and (epos == r.epos);
	// return (chr == r.chr) and (spos == r.spos) and (epos == r.epos) and (type == r.type);
}

/**************************************************************************************************/

void Record::init(void) {
	rname	= (char*) malloc(MAXLINESIZE);
	seq		= (char*) malloc(MAXLINESIZE);
	rcseq	= (char*) malloc(MAXLINESIZE);
	comment	= (char*) malloc(MAXLINESIZE);
	qual	= (char*) malloc(MAXLINESIZE);
	rqual	= (char*) malloc(MAXLINESIZE);
}

void Record::finalize(void) {
	free(rname);
	free(seq);
	free(rcseq);
	free(comment);
	free(qual);
	free(rqual);
}

Record::Record(Record* orig) : seq_len(orig->seq_len), mr(orig->mr) {
	strcpy(rname,   orig->rname);
	strcpy(seq,     orig->seq);
	strcpy(rcseq,   orig->rcseq);
	strcpy(comment, orig->comment);
	strcpy(qual,    orig->qual);
	strcpy(rqual,   orig->rqual);
}

void Record::deep_copy(Record* orig) {
	seq_len = orig->seq_len;
	*mr = *(orig->mr);

	strcpy(rname,   orig->rname);
	strcpy(seq,     orig->seq);
	strcpy(rcseq,   orig->rcseq);
	strcpy(comment, orig->comment);
	strcpy(qual,    orig->qual);
	strcpy(rqual,   orig->rqual);
}

bool Record::operator < (const Record& r) const {
	if (mr->contig_num != r.mr->contig_num) 
		return mr->contig_num < r.mr->contig_num;
	if (mr->chr_r1 != r.mr->chr_r1) 
		return mr->chr_r1 < r.mr->chr_r1;
	return mr->spos_r1 < r.mr->spos_r1;
}

/**************************************************************************************************/

RecordStr::RecordStr(const RecordStr& other) : mr(other.mr) {
	rname	= other.rname;
	seq		= other.seq;
	comment	= other.comment;
	qual	= other.qual;

}

RecordStr::RecordStr(Record* orig) : mr(orig->mr) {
	rname	= std::string(orig->rname);
	seq		= std::string(orig->seq);
	comment	= std::string(orig->comment);
	qual	= std::string(orig->qual);
}

void RecordStr::deep_copy(Record* orig) {
	mr = *(orig->mr);

	rname	= std::string(orig->rname);
	seq		= std::string(orig->seq);
	comment	= std::string(orig->comment);
	qual	= std::string(orig->qual);
}

bool RecordStr::operator < (const RecordStr& r) const {
	return mr.genome_spos < r.mr.genome_spos;
}

/**************************************************************************************************/
