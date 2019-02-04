#include <cstring>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <utility>
#include <set>

#include "gene_annotation.h"
#include "common.h"

#define MAXLINESIZE 10000
#define MAXGTFATTR 50
#define INITGTFREC 3e6

GTFParser::GTFParser(void) {
	input = NULL;
}

GTFParser::GTFParser(char* filename, const ContigLen* contig_len, int contig_count) {
	init(filename, contig_len, contig_count);
}

GTFParser::~GTFParser(void) {
	close_file(input);
	free(line);
}

void GTFParser::init(char* filename, const ContigLen* contig_len, int contig_count) {
	input = open_file(filename, "r");

	max_line_size = MAXLINESIZE;
	line = (char*) malloc(max_line_size);

	set_contig_shift(contig_len, contig_count);

	level["gene"] = 1;
	level["transcript"] = 2;
	level["exon"] = 3;
}

bool GTFParser::get_next(void) {
	return has_next() and read_next();
}

bool GTFParser::has_next(void) {
	return !feof(input);
}

bool GTFParser::read_next(void) {
	int len;
	if ((len = getline(&line, &max_line_size, input)) == -1)
		return false;

	// skip header
	if (line[0] == '#')	
		return read_next();

	return true;
}

bool is_delim(char c, const string& delim) {
	for (int i = 0; i < delim.size(); i++)
		if (c == delim[i])
			return true;
	return false;
}

void GTFParser::tokenize(char* line, int len, const string& delim, vector<string>& gtf_fields) {
	while (len && (line[len] == '\r' || line[len] == '\n'))
		line[len--] = 0;
	char *c = line;
	string cur_str = "";
	int cur_field = 0;

	while (*c) {
		if ( is_delim(*c, delim) ) {
			gtf_fields[cur_field] = cur_str;
			if (cur_str != "")	// skipping consecutive delimeters 
				cur_field++;
			cur_str = "";
		}
		else {
			cur_str += *c;
		}
		++c;
	}
	if (cur_str != "")
		gtf_fields[cur_field++] = cur_str;
}

// return true if valid
bool GTFParser::parse_gtf_rec(char* line, int len, GTFRecord* cr) {
	vector<string> gtf_fields (10);
	vector<string> gtf_attr (MAXGTFATTR);
	string major_delim = "\t";
	string minor_delim = " ;\"";
	
	tokenize(line, len, major_delim, gtf_fields);

	if (level.find(gtf_fields[2]) != level.end()) {
		char attr_str[MAXLINESIZE];
		strcpy(attr_str, gtf_fields[8].c_str());
		tokenize(attr_str, gtf_fields[8].length(), minor_delim, gtf_attr);

		cr->chr = gtf_fields[0];
		cr->source = gtf_fields[1];
		cr->type = gtf_fields[2];
		cr->start = atoi(gtf_fields[3].c_str());
		cr->end = atoi(gtf_fields[4].c_str());
		cr->forward_strand = (gtf_fields[6] == "+");

		for (int i = 0; i < gtf_attr.size(); i += 2) {
			if (gtf_attr[i] == "gene_id")
				cr->gene_id = gtf_attr[i+1];
			else if (gtf_attr[i] == "transcript_id")
				cr->trans_id = gtf_attr[i+1];
			else if (gtf_attr[i] == "exon_number") {
				cr->exon_num_int = atoi(gtf_attr[i+1].c_str());
				cr->exon_num = gtf_attr[i+1];
			}
			else if (gtf_attr[i] == "gene_name")
				cr->gene_name = gtf_attr[i+1];
		}
		return true;
	}
	else {
		return false;
	}
}

void copy_seg(GTFRecord* a, GTFRecord* b) {
	a->start = b->start;
	a->end = b->end;
	a->next_start = b->next_start;
	a->prev_end = b->prev_end;

	a->trans_id_int = b->trans_id_int;
	a->exon_num_int = b->exon_num_int;
	
	a->forward_strand = b->forward_strand;

	a->chr = b->chr;
	a->source = b->source;
	a->type = b->type;
	a->gene_id = b->gene_id;
	a->trans_id = b->trans_id;
	a->exon_num = b->exon_num;
	a->gene_name = b->gene_name;
}

void add2merged_exons(map <UniqSeg, string>& mymap, UniqSeg& seg, GTFRecord* rec) {
	map <UniqSeg, string>:: iterator it = mymap.find(seg);
	if (it != mymap.end()) {	// found seg
		string tmp_str = mymap[seg];
		seg = it->first;
		seg.trans_id.push_back(rec->trans_id_int);
		mymap.erase(it);
		mymap[seg] = tmp_str + "\t" + rec->trans_id + "-" + rec->exon_num;
	}
	else {
		seg.trans_id.clear();
		seg.trans_id.push_back(rec->trans_id_int);
		mymap[seg] = rec->trans_id + "-" + rec->exon_num;
	}
}

void GTFParser::chrloc2conloc(string& chr, uint32_t& start, uint32_t& end) {
	if (chr2con.find(chr) != chr2con.end()) {	// chr in genome
		start += chr2con[chr].shift;
		end += chr2con[chr].shift;
		chr = chr2con[chr].contig;
	}
	else 
		chr = "0";
}

bool GTFParser::load_gtf(void) {

	///
	//
	for (int i = 0; i < 3; i++) {
		near_border[i] = (uint8_t*) malloc(1200000000 * sizeof(uint8_t));
		memset(near_border[i], 0, 1200000000 * sizeof(uint8_t));
	}
	//
	///

	bool found;
	UniqSeg seg;
	string tmp_str;

	GTFRecord* current_record = new GTFRecord;
	GTFRecord* prev_record = new GTFRecord;
	
	prev_record->type = "";

	while (has_next()) {
		if (! get_next()) {		// end of file
			break;
		}
 
		found = parse_gtf_rec(line, len, current_record);
		if (!found) continue;

		chrloc2conloc(current_record->chr, current_record->start, current_record->end);

		if (current_record->chr == "0")
			continue;

		if (current_record->type == "gene") {
			int con = current_record->chr[0]-'1';
			if (con >= 0 and con < 3) {
				for (int k = current_record->start; k <= current_record->end; k++)
					near_border[current_record->chr[0]-'1'][k] |= 2;
			}

			GeneInfo tmp = {.start = current_record->start, .end = current_record->end};
			gid2ginfo[current_record->chr][current_record->gene_id] = tmp;

			if (merged_genes[current_record->chr].find(tmp) != merged_genes[current_record->chr].end()) {	// found gene
				merged_genes[current_record->chr][tmp] += "\t" + current_record->gene_id;
			}
			else {
				merged_genes[current_record->chr][tmp] = current_record->gene_id;
			}
		}

		if (current_record->type == "transcript") {
			transcript_ids[current_record->chr].push_back(current_record->trans_id);
		}

		if (current_record->type == "exon") {
			///
			//
			int con = current_record->chr[0]-'1';
			if (con >= 0 and con < 3) {
				for (int k = current_record->start; k <= current_record->end; k++)
					near_border[current_record->chr[0]-'1'][k] &= ~(2);
				
				for (int k = maxM(0, current_record->start - maxReadLength); k < current_record->start; k++)
					near_border[current_record->chr[0]-'1'][k] |= 1;
				for (int k = maxM(0, current_record->end - maxReadLength + 1); k <= current_record->end; k++)
					near_border[current_record->chr[0]-'1'][k] |= 1;

				//for (int k = current_record->start; k < current_record->start + maxReadLength; k++)
				//	near_border[current_record->chr[0]-'1'][k] |= 2;
				//for (int k = current_record->end + 1; k <= current_record->end + maxReadLength; k++)
				//	near_border[current_record->chr[0]-'1'][k] |= 2;
			}
			//
			///
			current_record->trans_id_int = transcript_ids[current_record->chr].size() - 1;
			if (prev_record->type != "exon") {
				//prev_record = current_record;
				copy_seg(prev_record, current_record);
				prev_record->next_start = 0;
				prev_record->prev_end = 0;
				continue;
			}
			else {
				if (prev_record->forward_strand) {
					prev_record->next_start = current_record->start;
					current_record->prev_end = prev_record->end;
				}
				else {
					prev_record->prev_end = current_record->end;
					current_record->next_start = prev_record->start;
				}
				seg.start			= prev_record->start;
				seg.end 			= prev_record->end;
				seg.gene_id			= prev_record->gene_id;
				seg.next_exon_beg	= prev_record->next_start;

				add2merged_exons(merged_exons[prev_record->chr], seg, prev_record);
				
				copy_seg(prev_record, current_record);
			}
		}
		else if (prev_record->type == "exon") {
			if (prev_record->forward_strand) {
				prev_record->next_start = 0;
			}
			else {
				prev_record->prev_end = 0;
			}
			seg.start			= prev_record->start;
			seg.end 			= prev_record->end;
			seg.gene_id			= prev_record->gene_id;
			seg.next_exon_beg	= prev_record->next_start;

			add2merged_exons(merged_exons[prev_record->chr], seg, prev_record);

			prev_record->type = "";
		}
	}

	///////
	
	if (prev_record->type == "exon") {
		if (prev_record->forward_strand) {
			prev_record->next_start = 0;
		}
		else {
			prev_record->prev_end = 0;
		}
		seg.start			= prev_record->start;
		seg.end 			= prev_record->end;
		seg.gene_id			= prev_record->gene_id;
		seg.next_exon_beg	= prev_record->next_start;

		add2merged_exons(merged_exons[prev_record->chr], seg, prev_record);
	}
	//////

	//fprintf(stdout, "162: %s\n", transcript_ids[merged_exons.begin()->first][162].c_str());

	map <string, map <UniqSeg, string> >:: iterator con_it;
	for (con_it = merged_exons.begin(); con_it != merged_exons.end(); con_it++) {
		exons_int_map[con_it->first].build(con_it->second);
		//exons_int_map[con_it->first].print();

		// construct transcript to segment table
		trans2seg[con_it->first].resize(transcript_ids[con_it->first].size());
		exons_int_map[con_it->first].build_trans2seg_table(transcript_ids[con_it->first].size(), trans2seg[con_it->first], trans_start_ind[con_it->first]);
	}

	map <string, map <GeneInfo, string> >:: iterator it;
	for (it = merged_genes.begin(); it != merged_genes.end(); it++) {
		genes_int_map[it->first].build(it->second);
		//genes_int_map[it->first].print();
	}

	int need_lu[3];
	for (int i = 0; i < 3; i++) {
		need_lu[i] = 0;
		for (int j = 0; j < 1200000000; j++)
			if (near_border[i][j] & 1)
				need_lu[i]++;
		fprintf(stdout, "Contig [%d]: Near exon boundaries: %d\n", i+1, need_lu[i]);
	}
	for (int i = 0; i < 3; i++) {
		need_lu[i] = 0;
		for (int j = 0; j < 1200000000; j++)
			if (near_border[i][j] & 2)
				need_lu[i]++;
		fprintf(stdout, "Contig [%d]: Intronic: %d\n", i+1, need_lu[i]);
	}

	delete prev_record;
	delete current_record;
	
	return true;
}

// assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
int GTFParser::binary_search(const vector <ExonSeg>& seg, int beg, int end, bool on_start, uint32_t target) {
	if (end - beg <= 1)
		return end;
	
	int mid = (beg + end) / 2;

	uint32_t to_comp = (on_start)? seg[mid].start : seg[mid].end;
	if (target < to_comp)
		return binary_search(seg, beg, mid, on_start, target);
	else
		return binary_search(seg, mid, end, on_start, target);
}

void GTFParser::print_record(const GTFRecord& r) {
	fprintf(stderr, "%s %s %d %d\n", r.gene_name.c_str(), r.chr.c_str(), r.start, r.end);
}

void GTFParser::set_contig_shift(const ContigLen* chr_info, int contig_count) {
	int curr_contig = 1;
	uint32_t sum_size = 0;
	ConShift con_shift;
	ConShift chr_shift;
	for (int i = 0; i < contig_count; i++) {
		ostringstream con_name;
		con_name << curr_contig;

		con_shift.contig = con_name.str();
		con_shift.shift = sum_size;
		chr2con[chr_info[i].name] = con_shift;

		chr_shift.contig = chr_info[i].name;
		chr_shift.shift = sum_size;
		con2chr[con_shift.contig].push_back(chr_shift);

		sum_size += chr_info[i].len + 1;	// +1 because of N at the end of each chr
		if (sum_size >= MIN_CONTIG_SIZE) {
			sum_size = 0;
			curr_contig++;
		}
	}
}

ConShift GTFParser::get_shift(const string& contig, uint32_t loc) {
	int i;
	for (i = 1; i < con2chr[contig].size(); i++)
		if (loc < con2chr[contig][i].shift)
			return con2chr[contig][i-1];
	return con2chr[contig][i-1];
}

// match an interval:
// [ spos, spos + mlen )
// spos: Start POSition of matched region
// mlen: Matched LENgth
// rlen: the lenght of the read remained to be matched (Rmained LENgth)
uint32_t GTFParser::get_upper_bound_lookup(uint32_t spos, uint32_t mlen, uint32_t rlen, uint32_t& max_end, const IntervalInfo<UniqSeg>*& ol_exons) {
	max_end = 0;

	// fprintf(stderr, "Searching for: [%u-%u], remain len: %u\n", spos, epos, rlen);
	//lookup_cnt++;
	
	int it_ind = -1;
	const IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigName].find_ind(spos, it_ind);

	uint32_t epos = spos + mlen - 1;
	if (ov_res == NULL or ov_res->seg_list.size() == 0) {	// not found => intronic
		ol_exons = NULL;
		max_end = gtf_parser.get_interval(it_ind + 1)->spos - 1;

		if (max_end < epos)	// => crossing the boundry
			return 0;
		else
			return minM(spos + rlen + EDTH, max_end - mlen + 1);
	}

	ol_exons = NULL;
	uint32_t min_end = 1e9;
	uint32_t max_next_exon = 0;
	if (epos > ov_res->epos) {
		for (int i = 0; i < ov_res->seg_list.size(); i++) {
			if (ov_res->seg_list[i].end >= epos) {	// => exonic
				max_end = maxM(max_end, ov_res->seg_list[i].end);
				min_end = minM(min_end, ov_res->seg_list[i].end);
				max_next_exon = maxM(max_next_exon, ov_res->seg_list[i].next_exon_beg);
				// fprintf(stderr, "Min end: %d\n Max end: %d\n  Max next: %d\n", min_end, max_end, max_next_exon);
				// fprintf(stderr, "Exon: [%d-%d]\n", ov_res->seg_list[i].start, ov_res->seg_list[i].end);
			}
		}
	}
	else {
		max_end = ov_res->max_end;
		min_end = ov_res->min_end;
		max_next_exon = ov_res->max_next_exon;
		// fprintf(stderr, "Min end: %d\n Max end: %d\n  Max next: %d\n", min_end, max_end, max_next_exon);
	}

	if (max_end > 0 and max_end >= epos) {	// exonic	
		ol_exons = ov_res;

		// loc excluded
		int32_t min2end = min_end - epos;

		if (min2end < rlen and max_next_exon != 0)	// junction is allowed
			return max_next_exon + mlen - 1;
			//return max_next_exon + rlen - min2end - mlen - 1;
		else
			return max_end - mlen + 1;
	}

	// loc excluded
	// int32_t min2end = min_end - epos;
	// if (max_end > 0 and min2end >= mlen) {	// exonic	
	// 	ol_exons = ov_res;

	// 	if (min2end < rlen and max_next_exon != 0)	// junction is allowed
	// 		return max_next_exon + mlen - 1;
	// 	else
	// 		return max_end - mlen + 1;
	// }

	else {	// on exon boundary
		max_end = 0;
		ol_exons = NULL;
		return 0;
	}
}

// match an interval:
// [ spos, spos + mlen )
// spos: Start POSition of matched region
// mlen: Matched LENgth
// rlen: the lenght of the read remained to be matched (Rmained LENgth)
void GTFParser::get_upper_bound_alu(uint32_t spos, uint32_t mlen, uint32_t rlen, JunctionDist& jd) {
	if (jd.looked_up) {
		//fprintf(stderr, "Looked up before!\n");
		return;
	}

	jd.looked_up = true;

	uint32_t max_beg = 0;
	uint32_t max_end = 0;
	uint32_t min_end = 1e9;
	uint32_t max_next_exon = 0;
	uint32_t epos = spos + mlen - 1;

	//fprintf(stderr, "Searching for: [%u-%u], remain len: %u\n", spos, epos, rlen);
	
	lookup_cnt++;
	IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigName].find(spos);

	if (ov_res == NULL or ov_res->seg_list.size() == 0) {	// not found => intronic
		jd.exonic = false;
		jd.dr = maxReadLength;	// do not consider junction
		jd.dl = maxReadLength;
		
		if (! (near_border[contigName[0]-'1'][spos] & 1)) {	// far from exon start
			jd.cross_boundry = false;
			jd.range = spos + rlen + EDTH;
		}
		else {
			// find end of intron
			// To be modified
			max_end = spos;
			while (near_border[contigName[0]-'1'][max_end] & 1)
				max_end++;
			max_end--;

			if (max_end - spos + 1 < mlen) {	// => crossing the boundry
				jd.cross_boundry = true;
				jd.range = 0;
			}
			else {
				jd.cross_boundry = false;
				jd.range =  max_end - mlen + 1;
				//jd.range = spos + rlen + EDTH;
			}
		}
		
		jd.max_end = max_end;
		return;
	}

	jd.exonic = true;
	for (int i = 0; i < ov_res->seg_list.size(); i++) {
		if (ov_res->seg_list[i].end >= epos) {	// => exonic
			max_beg = maxM(max_beg, ov_res->seg_list[i].start);
			max_end = maxM(max_end, ov_res->seg_list[i].end);
			min_end = minM(min_end, ov_res->seg_list[i].end);
			max_next_exon = maxM(max_next_exon, ov_res->seg_list[i].next_exon_beg);
			//fprintf(stderr, "Min end: %d\n Max end: %d\n  Max next: %d\n", min_end, max_end, max_next_exon);
			//fprintf(stderr, "Exon: [%d-%d]\n", ov_res->seg_list[i].start, ov_res->seg_list[i].end);
		}
	}

	if (max_end > 0) {	// exonic	
		jd.cross_boundry = false;
		
		// loc excluded
		int32_t min2end = min_end - epos;
		int32_t min2beg = spos - max_beg;
		jd.dr = min2end;
		jd.dl = min2beg;
		jd.max_end = max_end;

		if (min2end < rlen and max_next_exon != 0) {	// junction is allowed
			jd.range = max_next_exon + mlen - 1;
		}
		else {
			jd.range = max_end - mlen + 1;
		}
	}

	else {	// on exon boundary
		jd.cross_boundry = true;
		jd.dr = maxReadLength;
		jd.dl = maxReadLength;
		jd.range = 0;
		jd.max_end = 0;
	}
}

// returns intervals overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
const IntervalInfo<UniqSeg>* GTFParser::get_location_overlap(uint32_t loc, bool use_mask) {
	// do not use mask if extending left
	if (use_mask and !(near_border[contigNum][loc] & 1)) {		// intronic
		//fprintf(stderr, "skip lookup\n");
		return NULL;
	}
	
	IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigName].find(loc);

	if (ov_res == NULL or ov_res->seg_list.size() == 0)	// not found => intronic
		return NULL;

	return ov_res;
}

// returns intervals overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
const IntervalInfo<UniqSeg>* GTFParser::get_location_overlap_ind(uint32_t loc, bool use_mask, int& ind) {
	// do not use mask if extending left
	if (use_mask and !(near_border[contigNum][loc] & 1)) {		// intronic
		//fprintf(stderr, "skip lookup\n");
		return NULL;
	}
	
	IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigName].find_ind(loc, ind);

	if (ov_res == NULL or ov_res->seg_list.size() == 0)	// not found => intronic
		return NULL;

	return ov_res;
}

// returns genes overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
const IntervalInfo<GeneInfo>* GTFParser::get_gene_overlap(uint32_t loc, bool use_mask) {
	// do not use mask if extending left
	if (use_mask and !(near_border[contigNum][loc] & 1)) {		// intronic
		//fprintf(stderr, "skip lookup\n");
		return NULL;
	}
	
	IntervalInfo<GeneInfo>* ov_res = genes_int_map[contigName].find(loc);

	if (ov_res == NULL or ov_res->seg_list.size() == 0)	// not found => intronic
		return NULL;

	return ov_res;
}
