#include <cstring>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <utility>
#include <set>

#include "gene_annotation.h"
#include "common.h"

#define MAXGTFLINESIZE 10000
#define MAXGTFATTR 50
#define INITGTFREC 3e6

GTFParser::GTFParser(void) {
	input = NULL;
}

GTFParser::GTFParser(char* filename, const vector <ContigLen>& contig_len) {
	init(filename, contig_len);
}

GTFParser::~GTFParser(void) {
	close_file(input);
	free(line);
}

void GTFParser::init(char* filename, const vector <ContigLen>& contig_len) {
	char* fname = (char*) malloc(FILE_NAME_MAX_LEN);
	char* rmode = (char*) malloc(FILE_NAME_MAX_LEN);

	sprintf(fname, "%s", filename);
	sprintf(rmode, "%c", 'r');

	input = open_file(fname, rmode);

	free(fname);
	free(rmode);

	max_line_size = MAXGTFLINESIZE;
	line = (char*) malloc(max_line_size);

	set_contig_shift(contig_len);

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
	for (unsigned int i = 0; i < delim.size(); i++)
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
		char attr_str[MAXGTFLINESIZE];
		strcpy(attr_str, gtf_fields[8].c_str());
		tokenize(attr_str, gtf_fields[8].length(), minor_delim, gtf_attr);

		cr->chr = gtf_fields[0];
		cr->source = gtf_fields[1];
		cr->type = gtf_fields[2];
		cr->start = atoi(gtf_fields[3].c_str());
		cr->end = atoi(gtf_fields[4].c_str());
		cr->forward_strand = (gtf_fields[6] == "+");

		for (unsigned int i = 0; i < gtf_attr.size(); i += 2) {
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
	a->start 		= b->start;
	a->end 			= b->end;
	a->next_start 	= b->next_start;
	a->prev_end 	= b->prev_end;

	a->gene_id_int 	= b->gene_id_int;
	a->trans_id_int = b->trans_id_int;
	a->exon_num_int = b->exon_num_int;
	a->chr_id 		= b->chr_id;

	a->forward_strand = b->forward_strand;

	a->chr 			= b->chr;
	a->source 		= b->source;
	a->type 		= b->type;
	a->gene_id 		= b->gene_id;
	a->trans_id 	= b->trans_id;
	a->exon_num 	= b->exon_num;
	a->gene_name 	= b->gene_name;
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
	bool found;
	UniqSeg seg;
	string tmp_str;

	GTFRecord* current_record = new GTFRecord;
	GTFRecord* prev_record = new GTFRecord;
	
	prev_record->type = "";

	int tmp_chr;
	while (has_next()) {
		if (! get_next()) {		// end of file
			break;
		}
 
		found = parse_gtf_rec(line, len, current_record);
		if (!found) continue;

		chrloc2conloc(current_record->chr, current_record->start, current_record->end);
		tmp_chr = atoi(current_record->chr.c_str()) - 1;

		if (tmp_chr < 0)	// chr not present in genome index
			continue;

		// if (current_record->chr == "0")	// chr not found in genome index
		// 	continue;

		current_record->chr_id = tmp_chr;

		if (current_record->type == "gene") {
			if (current_record->chr_id >= gene_ids.size()) {
				gene_ids.resize(current_record->chr_id+1);
			}
			gene_ids[current_record->chr_id].push_back(current_record->gene_id);
			
			if (current_record->chr_id >= near_border_bs.size()) {
				near_border_bs.resize(current_record->chr_id+1);
				near_border_bs[current_record->chr_id].reset();
			}
			if (current_record->chr_id >= intronic_bs.size()) {
				intronic_bs.resize(current_record->chr_id+1);
				intronic_bs[current_record->chr_id].reset();
			}

			for (uint32_t k = current_record->start; k <= current_record->end; k++) {
				intronic_bs[current_record->chr_id].set(k, 1);
			}

			if (current_record->chr_id >= gid2ginfo.size()) {
				gid2ginfo.resize(current_record->chr_id+1);
			}
			GeneInfo tmp = {.start = current_record->start, .end = current_record->end, 
							.gene_id = uint32_t(gid2ginfo[current_record->chr_id].size())};
			gid2ginfo[current_record->chr_id].push_back(tmp);

			if (current_record->chr_id >= merged_genes.size()) {
				merged_genes.resize(current_record->chr_id+1);
			}

			if (merged_genes[current_record->chr_id].find(tmp) != merged_genes[current_record->chr_id].end()) {	// found gene
				merged_genes[current_record->chr_id][tmp] += "\t" + current_record->gene_id;
			}
			else {
				merged_genes[current_record->chr_id][tmp] = current_record->gene_id;
			}
		}

		if (current_record->type == "transcript") {
			if (current_record->chr_id >= transcript_ids.size()) {
				transcript_ids.resize(current_record->chr_id+1);
			}
			transcript_ids[current_record->chr_id].push_back(current_record->trans_id);
		}

		if (current_record->type == "exon") {

			// prepare mask
			for (uint32_t k = current_record->start; k <= current_record->end; k++) {
				intronic_bs[current_record->chr_id].set(k, 0);
			}
			
			for (uint32_t k = maxM(0, current_record->start - maxReadLength); k < current_record->start; k++) {
				near_border_bs[current_record->chr_id].set(k, 1);
			}
			for (uint32_t k = maxM(0, current_record->end - maxReadLength + 1); k <= current_record->end; k++) {
				near_border_bs[current_record->chr_id].set(k, 1);
			}
			//

			current_record->trans_id_int = transcript_ids[current_record->chr_id].size() - 1;
			current_record->gene_id_int = gene_ids[current_record->chr_id].size() - 1;
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
				seg.gene_id			= prev_record->gene_id_int;
				seg.next_exon_beg	= prev_record->next_start;

				if (prev_record->chr_id >= merged_exons.size()) {
					merged_exons.resize(prev_record->chr_id+1);
				}
				add2merged_exons(merged_exons[prev_record->chr_id], seg, prev_record);
				
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
			seg.gene_id			= prev_record->gene_id_int;
			seg.next_exon_beg	= prev_record->next_start;

			if (prev_record->chr_id >= merged_exons.size()) {
				merged_exons.resize(prev_record->chr_id+1);
			}
			add2merged_exons(merged_exons[prev_record->chr_id], seg, prev_record);

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
		seg.gene_id			= prev_record->gene_id_int;
		seg.next_exon_beg	= prev_record->next_start;

		if (prev_record->chr_id >= merged_exons.size()) {
			merged_exons.resize(prev_record->chr_id+1);
		}
		add2merged_exons(merged_exons[prev_record->chr_id], seg, prev_record);
	}
	//////

	exons_int_map.resize(merged_exons.size());
	trans2seg.resize(merged_exons.size());
	trans_start_ind.resize(merged_exons.size());

	for (unsigned int con = 0; con < merged_exons.size(); con++) {
		exons_int_map[con].build(merged_exons[con]);
		//exons_int_map[con_it->first].print();

		// construct transcript to segment table
		trans2seg[con].resize(transcript_ids[con].size());
		exons_int_map[con].build_trans2seg_table(transcript_ids[con].size(), trans2seg[con], trans_start_ind[con]);
	}

	genes_int_map.resize(merged_genes.size());
	for (unsigned int con = 0; con < merged_genes.size(); con++) {
		genes_int_map[con].build(merged_genes[con]);
		//genes_int_map[con].print();
	}

	for (unsigned int i = 0; i < near_border_bs.size(); i++) {
		fprintf(stdout, "Contig [%u]: Near exon boundaries: %zu\n", i+1, near_border_bs[i].count());
	}
	for (unsigned int i = 0; i < intronic_bs.size(); i++) {
		fprintf(stdout, "Contig [%u]: Intronic: %zu\n", i+1, intronic_bs[i].count());
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

void GTFParser::set_contig_shift(const vector <ContigLen>& chr_info) {
	unsigned int contig_count = chr_info.size();
	uint32_t curr_contig;

	ConShift con_shift;
	ConShift chr_shift;

	for (unsigned int i = 0; i < contig_count; i++) {
		curr_contig = chr_info[i].contig_id;
		ostringstream con_name;
		con_name << curr_contig;

		con_shift.contig = con_name.str();
		con_shift.shift = chr_info[i].start_pos;
		chr2con[chr_info[i].name] = con_shift;

		chr_shift.contig = chr_info[i].name;
		chr_shift.shift = chr_info[i].start_pos;

		if (curr_contig > con2chr.size()) {
			con2chr.resize(curr_contig);
		}
		con2chr[curr_contig-1].push_back(chr_shift);

	}
}

ConShift GTFParser::get_shift(int contig_id, uint32_t loc) {
	unsigned int i;
	for (i = 1; i < con2chr[contig_id].size(); i++)
		if (loc < con2chr[contig_id][i].shift)
			return con2chr[contig_id][i-1];
	return con2chr[contig_id][i-1];
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
	const IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigNum].find_ind(spos, it_ind);

	uint32_t epos = spos + mlen - 1;
	if (ov_res == NULL or ov_res->seg_list.size() == 0) {	// not found => intronic
		ol_exons = NULL;
		max_end = gtf_parser.get_interval(it_ind + 1)->spos - 1;

		if (max_end < epos)	// => crossing the boundry
			return 0;
		else
			return minM(spos + rlen + maxEd, max_end - mlen + 1);
	}

	ol_exons = NULL;
	uint32_t min_end = 1e9;
	uint32_t max_next_exon = 0;
	if (epos > ov_res->epos) {
		for (unsigned int i = 0; i < ov_res->seg_list.size(); i++) {
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
		//int32_t min2end = min_end - epos;

		if (min_end < rlen + epos and max_next_exon != 0)	// junction is allowed
			return max_next_exon + mlen - 1;
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


// returns intervals overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
const IntervalInfo<UniqSeg>* GTFParser::get_location_overlap(uint32_t loc, bool use_mask) {
	// do not use mask if extending left
	//if (use_mask and !(near_border_bs[contigNum][loc])) {		// intronic
	//	//fprintf(stderr, "skip lookup\n");
	//	return NULL;
	//}
	
	IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigNum].find(loc);

	if (ov_res == NULL or ov_res->seg_list.size() == 0)	// not found => intronic
		return NULL;

	return ov_res;
}

// returns intervals overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
const IntervalInfo<UniqSeg>* GTFParser::get_location_overlap_ind(uint32_t loc, bool use_mask, int& ind) {
	// do not use mask if extending left
	//if (use_mask and !(near_border_bs[contigNum][loc])) {		// intronic
	//	//fprintf(stderr, "skip lookup\n");
	//	return NULL;
	//}
	
	IntervalInfo<UniqSeg>* ov_res = exons_int_map[contigNum].find_ind(loc, ind);

	if (ov_res == NULL or ov_res->seg_list.size() == 0)	// not found => intronic
		return NULL;

	return ov_res;
}

// returns genes overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
const IntervalInfo<GeneInfo>* GTFParser::get_gene_overlap(uint32_t loc, bool use_mask) {
	// do not use mask if extending left
	if (use_mask and !(near_border_bs[contigNum][loc])) {		// intronic
		//fprintf(stderr, "skip lookup\n");
		return NULL;
	}
	
	IntervalInfo<GeneInfo>* ov_res = genes_int_map[contigNum].find(loc);

	if (ov_res == NULL or ov_res->seg_list.size() == 0)	// not found => intronic
		return NULL;

	return ov_res;
}
