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
	delete current_record;
}

void GTFParser::init(char* filename, const ContigLen* contig_len, int contig_count) {
	input = open_file(filename, "r");

	max_line_size = MAXLINESIZE;
	line = (char*) malloc(max_line_size);

	current_record = new GTFRecord;
	//records.reserve(INITGTFREC);
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
bool GTFParser::parse_gtf_rec(char* line, int len) {
	vector<string> gtf_fields (10);
	vector<string> gtf_attr (MAXGTFATTR);
	string major_delim = "\t";
	string minor_delim = " ;\"";
	
	tokenize(line, len, major_delim, gtf_fields);

	if (level.find(gtf_fields[2]) != level.end()) {
		char attr_str[MAXLINESIZE];
		strcpy(attr_str, gtf_fields[8].c_str());
		tokenize(attr_str, gtf_fields[8].length(), minor_delim, gtf_attr);

		current_record->chr = gtf_fields[0];
		current_record->source = gtf_fields[1];
		current_record->type = gtf_fields[2];
		current_record->start = atoi(gtf_fields[3].c_str());
		current_record->end = atoi(gtf_fields[4].c_str());
		current_record->forward_strand = (gtf_fields[6] == "+");

		for (int i = 0; i < gtf_attr.size(); i += 2) {
			if (gtf_attr[i] == "gene_id")
				current_record->gene_id = gtf_attr[i+1];
			else if (gtf_attr[i] == "transcript_id")
				current_record->trans_id = gtf_attr[i+1];
			else if (gtf_attr[i] == "exon_number") {
				current_record->exon_num_int = atoi(gtf_attr[i+1].c_str());
				current_record->exon_num = gtf_attr[i+1];
			}
			else if (gtf_attr[i] == "gene_name")
				current_record->gene_name = gtf_attr[i+1];
		}
		return true;
	}
	else {
//		fprintf(stdout, "Type: %s\n", gtf_fields[2].c_str());
		return false;
	}
}

void copy_seg(GTFRecord* a, GTFRecord* b) {
	a->chr = b->chr;
	a->type = b->type;
	a->start = b->start;
	a->end = b->end;
	a->gene_id = b->gene_id;
	a->next_start = b->next_start;
	a->prev_end = b->prev_end;
	a->trans_id = b->trans_id;
	a->exon_num = b->exon_num;
}

void GTFParser::chrloc2conloc(string& chr, uint32_t& start, uint32_t& end) {
	if (chr2con.find(chr) != chr2con.end()) {	// chr in genome
		start += chr2con[chr].shift;
		end += chr2con[chr].shift;
		chr = chr2con[chr].contig;
	}
}

bool GTFParser::load_gtf(void) {
	bool found;
	UniqSeg seg;
	GTFRecord* prev_record = new GTFRecord;
	prev_record->type = "";

	while (has_next()) {
		if (! get_next()) {		// end of file
			break;
		}
 
		found = parse_gtf_rec(line, len);
		if (!found) continue;

		chrloc2conloc(current_record->chr, current_record->start, current_record->end);
		//if (current_record->chr != "2")		// work on chr 2 for now
		//	continue;

		//fprintf(stdout, "-->%s: %d-%d\n", prev_record->type.c_str(), prev_record->start, prev_record->end);
		//fprintf(stdout, "==>%s: %d-%d\n", current_record->type.c_str(), current_record->start, current_record->end);

		if (current_record->type == "exon") {
			if (prev_record->type != "exon") {
				//prev_record = current_record;
				copy_seg(prev_record, current_record);
				prev_record->next_start = 0;
				prev_record->prev_end = 0;
				continue;
			}
			else {
				if (current_record->forward_strand) {
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
				seg.prev_exon_end	= prev_record->prev_end;

				//fprintf(stdout, "%d [%s: %d-%d] %d\n", seg.prev_exon_end, seg.gene_id.c_str(), seg.start, seg.end, seg.next_exon_beg);

				if (merged_exons[prev_record->chr].find(seg) != merged_exons[prev_record->chr].end()) {	// found seg
					merged_exons[prev_record->chr][seg] += "\t" + prev_record->trans_id + "-" + prev_record->exon_num;
				}
				else {
					merged_exons[prev_record->chr][seg] = prev_record->trans_id + "-" + prev_record->exon_num;
				}
				//prev_record = current_record;
				copy_seg(prev_record, current_record);
				records.push_back(*current_record);
			}
		}
		else if (prev_record->type == "exon") {
			if (current_record->forward_strand) {
				prev_record->next_start = 0;
			}
			else {
				prev_record->prev_end = 0;
			}
			seg.start			= prev_record->start;
			seg.end 			= prev_record->end;
			seg.gene_id			= prev_record->gene_id;
			seg.next_exon_beg	= prev_record->next_start;
			seg.prev_exon_end	= prev_record->prev_end;

			//fprintf(stdout, "%d [%s: %d-%d] %d\n", seg.prev_exon_end, seg.gene_id.c_str(), seg.start, seg.end, seg.next_exon_beg);

			if (merged_exons[prev_record->chr].find(seg) != merged_exons[prev_record->chr].end()) {	// found seg
				merged_exons[prev_record->chr][seg] += "\t" + prev_record->trans_id + "-" + prev_record->exon_num;
			}
			else {
				merged_exons[prev_record->chr][seg] = prev_record->trans_id + "-" + prev_record->exon_num;
			}

			prev_record->type = "";
		}
	}

	///////
	
	if (prev_record->type == "exon") {
		if (current_record->forward_strand) {
			prev_record->next_start = 0;
		}
		else {
			prev_record->prev_end = 0;
		}
		seg.start			= prev_record->start;
		seg.end 			= prev_record->end;
		seg.gene_id			= prev_record->gene_id;
		seg.next_exon_beg	= prev_record->next_start;
		seg.prev_exon_end	= prev_record->prev_end;

		//fprintf(stdout, "%d [%s: %d-%d] %d\n", seg.prev_exon_end, seg.gene_id.c_str(), seg.start, seg.end, seg.next_exon_beg);

		if (merged_exons[prev_record->chr].find(seg) != merged_exons[prev_record->chr].end()) {	// found seg
			merged_exons[prev_record->chr][seg] += "\t" + prev_record->trans_id + "-" + prev_record->exon_num;
		}
		else {
			merged_exons[prev_record->chr][seg] = prev_record->trans_id + "-" + prev_record->exon_num;
		}
		//map <UniqSeg, string>:: iterator it = merged_exons.find(seg);
		//fprintf(stdout, "Next beg: %u\n", (it->first).next_exon_beg);
	}
	//////

	map <string, map <UniqSeg, string> >:: iterator con_it;
	map <UniqSeg, string>:: iterator it;

	for (con_it = merged_exons.begin(); con_it != merged_exons.end(); con_it++) {
		for (it = con_it->second.begin(); it != con_it->second.end(); it++) {
			
			UniqSegList temp_list;
			temp_list += it->first;
			exons_int_map[con_it->first] += make_pair(boost::icl::discrete_interval <uint32_t>::closed (it->first.start, it->first.end), temp_list);
			//if (con_it->first == "1")
			//exons_int_map2 += make_pair(boost::icl::discrete_interval <uint32_t>::closed (it->first.start, it->first.end), temp_list);

			merged_exons_arr[con_it->first].push_back(it->first);
			//fprintf(stderr, "%s -> %10u [%s:%10u%10u] %10u = Val: %s\n", con_it->first.c_str(), it->first.prev_exon_end, it->first.gene_id.c_str(), it->first.start, it->first.end, (it->first).next_exon_beg, it->second.c_str());
			//int n = merged_exons_arr[con_it->first].size() - 1;
			//fprintf(stdout, "%10u [%s:%10u%10u] %10u\n", merged_exons_arr[con_it->first][n].prev_exon_end, merged_exons_arr[con_it->first][n].gene_id.c_str(), merged_exons_arr[con_it->first][n].start, merged_exons_arr[con_it->first][n].end, merged_exons_arr[con_it->first][n].next_exon_beg);
		}
	}

	//for (boost::icl::interval_map<uint32_t, UniqSegList >::iterator it = exons_int_map.begin(); it != exons_int_map.end(); it++) {
	//	fprintf (stdout, "[%10u%10u] %d\n", first(it->first), last(it->first), it->second.seg_list.size());
	//	for (int k = 0; k < it->second.seg_list.size(); k++)
	//		fprintf(stdout, "%---10u [%s:%10u%10u] %10u\n", it->second.seg_list[k].prev_exon_end, it->second.seg_list[k].gene_id.c_str(), it->second.seg_list[k].start, it->second.seg_list[k].end, it->second.seg_list[k].next_exon_beg);
	//}

	set_wild_type();
	return true;
}

void GTFParser::set_wild_type(void) {
	if (records.size() <= 0)
		return;
	sort(records.begin(), records.end());
	//print_records();

	ExonSeg wt_exon;
	wt_exon.gene_id = records[0].gene_id;
	wt_exon.gene_name = records[0].gene_name;
	wt_exon.chr = records[0].chr;
	wt_exon.start = records[0].start;
	wt_exon.end = records[0].end;
	for (int i = 1; i < records.size(); i++) {
		if (records[i].chr != wt_exon.chr or records[i].start > wt_exon.end) {
			if (chr2con.find(wt_exon.chr) != chr2con.end()) {	// chr in genome
				wt_exon.start += chr2con[wt_exon.chr].shift;
				wt_exon.end += chr2con[wt_exon.chr].shift;
				wt_exon.chr = chr2con[wt_exon.chr].contig;
			}
			wt_exons[wt_exon.chr].push_back(wt_exon);

			wt_exon.gene_id = records[i].gene_id;
			wt_exon.gene_name = records[i].gene_name;
			wt_exon.chr = records[i].chr;
			wt_exon.start = records[i].start;
			wt_exon.end = records[i].end;
		}
		else if (records[i].end > wt_exon.end){
			wt_exon.end = records[i].end;
		}
	}
	if (chr2con.find(wt_exon.chr) != chr2con.end()) {	// chr in genome
		wt_exon.start += chr2con[wt_exon.chr].shift;
		wt_exon.end += chr2con[wt_exon.chr].shift;
		wt_exon.chr = chr2con[wt_exon.chr].contig;
	}
	wt_exons[wt_exon.chr].push_back(wt_exon);
}

// assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// => 
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
int GTFParser::binary_search(int beg, int end, uint32_t target) {
	if (end - beg <= 1)
		return end;
	
	int mid = (beg + end) / 2;

	if (target < merged_exons_arr[contigName][mid].start)
		return binary_search(beg, mid, target);
	else
		return binary_search(mid, end, target);
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
	//fprintf(stderr, "beg: %d, mid: %d, seg[mid]: %d\t%d\n", beg, mid, seg[mid].start, seg[mid].end);

	uint32_t to_comp = (on_start)? seg[mid].start : seg[mid].end;
	if (target < to_comp)
		return binary_search(seg, beg, mid, on_start, target);
	else
		return binary_search(seg, mid, end, on_start, target);
}

int GTFParser::search_loc(const string& chr, bool on_start, uint32_t target) {
	//fprintf(stderr, "contig: %s, start? %d, target: %lu\n",chr.c_str(), on_start, target );
	int ret_ind = binary_search(wt_exons[chr], 0, wt_exons[chr].size(), on_start, target);
	return (on_start)? ret_ind-1 : ret_ind;
}

uint32_t GTFParser::get_start(const string& chr, int seg_ind) {
	return wt_exons[chr][seg_ind].start;
}

uint32_t GTFParser::get_end(const string& chr, int seg_ind) {
	return wt_exons[chr][seg_ind].end;
}

void  GTFParser::get_gene_id(const string& chr, int seg_ind, string& gene_id) {
	gene_id = wt_exons[chr][seg_ind].gene_id;
}

bool GTFParser::is_last_exonic_region(const string& chr, int seg_ind) {
	return (seg_ind >= (wt_exons[chr].size() - 1));
}

void GTFParser::print_record(const GTFRecord& r) {
	fprintf(stderr, "%s %s %d %d\n", r.gene_name.c_str(), r.chr.c_str(), r.start, r.end);
}

void GTFParser::print_records(void) {
	for (int i = 0; i < records.size(); i++) {
		fprintf(stderr, "%d: ", i);
		print_record(records[i]);
	}
}

void GTFParser::set_contig_shift(const ContigLen* contig_len, int contig_count) {
	int curr_contig = 1;
	uint32_t sum_size = 0;
	ConShift con_shift;
	for (int i = 0; i < contig_count; i++) {
		ostringstream oss;
		oss << curr_contig;
		con_shift.contig = oss.str();
		con_shift.shift = sum_size;
		chr2con[contig_len[i].name] = con_shift;

		sum_size += contig_len[i].len + 1;	// +1 because of N at the end of each chr
		if (sum_size >= MIN_CONTIG_SIZE) {
			sum_size = 0;
			curr_contig++;
		}
	}
}

uint32_t GTFParser::get_upper_bound(uint32_t loc, int len) {

	uint32_t max_end = 0;
	uint32_t min_end = 1e9;
	uint32_t max_next_exon = 0;
	int i = binary_search(0, merged_exons_arr[contigName].size(), loc) - 1;
	
	while (i >= 0) {
		if (merged_exons_arr[contigName][i].end >= loc) {
			max_end = maxM(max_end, merged_exons_arr[contigName][i].end);
			min_end = minM(min_end, merged_exons_arr[contigName][i].end);
			max_next_exon = maxM(max_next_exon, merged_exons_arr[contigName][i].next_exon_beg);
			//fprintf(stderr, "Min end: %d\n Max end: %d\n  Max next: %d\n", min_end, max_end, max_next_exon);
		}
		i--;
	}
	
	// loc excluded
	int32_t min2end = min_end - loc;

	if (min2end >= len or max_next_exon == 0)
		return max_end;

	else
		return max_next_exon;

}

// match an interval:
// [ spos, spos + mlen )
// spos: Start POSition of matched region
// mlen: Matched LENgth
// rlen: the lenght of the read remained to be matched (Rmained LENgth)
uint32_t GTFParser::get_upper_bound(uint32_t spos, uint32_t mlen, uint32_t rlen, uint32_t& max_end) {
	lookup_cnt++;
	max_end = 0;
	uint32_t min_end = 1e9;
	uint32_t max_next_exon = 0;
	uint32_t epos = spos + mlen - 1;

	//boost::icl::discrete_interval <uint32_t> spos_int = boost::icl::discrete_interval <uint32_t>::closed(spos, spos);
	boost::icl::interval_map<uint32_t, UniqSegList >::const_iterator fit;
	//fit = exons_int_map[contigName].find(spos_int);
	fit = lookup_arr[spos];

	if (fit == exons_int_map[contigName].end())	// not found => intronic
	//if (fit == exons_int_map2.end())	// not found => intronic
		return epos + rlen - mlen + 1;

	for (int i = 0; i < fit->second.seg_list.size(); i++) {
		if (fit->second.seg_list[i].end >= epos) {	// => exonic
			max_end = maxM(max_end, fit->second.seg_list[i].end);
			min_end = minM(min_end, fit->second.seg_list[i].end);
			max_next_exon = maxM(max_next_exon, fit->second.seg_list[i].next_exon_beg);
			//fprintf(stderr, "Min end: %d\n Max end: %d\n  Max next: %d\n", min_end, max_end, max_next_exon);
			//fprintf(stderr, "Exon: [%d-%d]\n", fit->second.seg_list[i].start, fit->second.seg_list[i].end);
		}
	}

	if (max_end > 0) {	// exonic	
		// loc excluded
		int32_t min2end = min_end - epos;

		if (min2end < rlen and max_next_exon != 0)	// junction is allowed
			return max_next_exon + mlen - 1;
		else
			return max_end - mlen + 1;
	}

	else	// on exon boundary
		return 0;
}

// returns intervals overlapping with: loc
// remain lenght does not include loc itself (starting from next location)
void GTFParser::get_location_overlap(uint32_t loc, vector <UniqSeg>& overlap) {
	overlap.clear();

	//boost::icl::discrete_interval <uint32_t> loc_int = boost::icl::discrete_interval <uint32_t>::closed(loc, loc);
	boost::icl::interval_map<uint32_t, UniqSegList >::const_iterator fit;
	//fit = exons_int_map[contigName].find(loc_int);
	fit = lookup_arr[loc];

	if (fit == exons_int_map[contigName].end())	// not found => intronic
	//if (fit == exons_int_map2.end())	// not found => intronic
		return;

	for (int i = 0; i < fit->second.seg_list.size(); i++) {
		overlap.push_back(fit->second.seg_list[i]);
	}

}

// match an interval:
// [ spos, spos + mlen )
// spos: Start POSition of matched region
// mlen: Matched LENgth
// rlen: the lenght of the read remained to be matched (Rmained LENgth)
//uint32_t GTFParser::get_upper_bound(uint32_t spos, uint32_t mlen, uint32_t rlen, uint32_t& max_end) {
//	uint32_t epos = spos + mlen - 1;
//
//	bool partial_overlap = false;
//	max_end = 0;
//	uint32_t min_end = 1e9;
//	uint32_t max_next_exon = 0;
//
//	int i = binary_search(0, merged_exons_arr[contigName].size(), spos) - 1;
//	int j = i + 1;
//
//	while (i >= 0 and j - i < 1000) {
//		if (merged_exons_arr[contigName][i].end >= epos) {
//			max_end = maxM(max_end, merged_exons_arr[contigName][i].end);
//			min_end = minM(min_end, merged_exons_arr[contigName][i].end);
//			max_next_exon = maxM(max_next_exon, merged_exons_arr[contigName][i].next_exon_beg);
//			//fprintf(stderr, "Min end: %d\n Max end: %d\n  Max next: %d\n", min_end, max_end, max_next_exon);
//			//fprintf(stderr, "Exon: [%d-%d]\n", merged_exons_arr[contigName][i].start, merged_exons_arr[contigName][i].end);
//		}
//		else if (merged_exons_arr[contigName][i].end >= spos)
//			partial_overlap = true;
//		i--;
//	}
//
//	if (max_end > 0) {	// exonic	
//		// loc excluded
//		int32_t min2end = min_end - epos;
//
//		if (min2end >= rlen)		// no junction is allowed
//			return max_end - mlen + 1;
//
//		else
//			return max_next_exon + mlen - 1;
//	}
//
//	else {
//		if (j < merged_exons_arr[contigName].size() and merged_exons_arr[contigName][j].start <= epos)
//			partial_overlap = true;
//
//		if (partial_overlap)
//			return 0;
//		else
//			return epos + rlen - mlen + 1;
//	}
//
//}
//
//// returns intervals overlapping with: loc
//// remain lenght does not include loc itself (starting from next location)
//void GTFParser::get_location_overlap(uint32_t loc, vector <UniqSeg>& overlap) {
//	overlap.clear();
//
//	int i = binary_search(0, merged_exons_arr[contigName].size(), loc) - 1;
//	int j = i + 1;
//
//	while (i >= 0 and j - i < 1000) {
//		if (merged_exons_arr[contigName][i].end >= loc) {
//			overlap.push_back(merged_exons_arr[contigName][i]);
//		}
//		i--;
//	}
//}
