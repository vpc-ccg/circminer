#include "gene_annotation.h"
#include <cstring>
#include <algorithm>
#include <iostream>

#define MAXLINESIZE 10000
#define MAXGTFATTR 50
#define INITGTFREC 3e6

GTFParser::GTFParser(void) {
	input = NULL;
}

GTFParser::GTFParser(char* filename) {
	init(filename);
}

GTFParser::~GTFParser(void) {
	if (input != NULL)
		fclose(input);
}

void GTFParser::init(char* filename) {
	input = fopen(filename, "r");
	if (input == NULL) {
		fprintf(stderr, "Could not open %s\n", filename);
		exit(1);
	}

	//fseek(input, 0L, SEEK_END);
	//file_size = ftell(input);
	//fseek(input, 0L, SEEK_SET);

	max_line_size = MAXLINESIZE;
	line = (char*) malloc(max_line_size);

	current_record = new GTFRecord;
	//records.reserve(INITGTFREC);
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

inline bool is_delim(char c, const string& delim) {
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

// return true if it's an exon
bool GTFParser::parse_gtf_rec(char* line, int len) {
	vector<string> gtf_fields (10);
	vector<string> gtf_attr (MAXGTFATTR);
	string major_delim = "\t";
	string minor_delim = " ;\"";
	
	tokenize(line, len, major_delim, gtf_fields);

	if (gtf_fields[2] == "exon") {
		char attr_str[MAXLINESIZE];
		strcpy(attr_str, gtf_fields[8].c_str());
		tokenize(attr_str, gtf_fields[8].length(), minor_delim, gtf_attr);

		//for (int i = 0; i < gtf_fields.size(); i++)
		//	cout << gtf_fields[i] << "--";
		//cout << endl;
		//for (int i = 0; i < gtf_attr.size(); i++)
		//	cout << gtf_attr[i] << "=";
		//cout << endl;

		current_record->chr = gtf_fields[0];
		current_record->source = gtf_fields[1];
		//current_record->type = gtf_fields[2];
		current_record->start = atoi(gtf_fields[3].c_str());
		current_record->end = atoi(gtf_fields[4].c_str());

		for (int i = 0; i < gtf_attr.size(); i += 2) {
			if (gtf_attr[i] == "gene_id")
				current_record->gene_id = gtf_attr[i+1];
			else if (gtf_attr[i] == "transcript_id")
				current_record->trans_id = gtf_attr[i+1];
			else if (gtf_attr[i] == "exon_number")
				current_record->exon_num = atoi(gtf_attr[i+1].c_str());
			else if (gtf_attr[i] == "gene_name")
				current_record->gene_name = gtf_attr[i+1];
		}
		return true;
	}
	else
		return false;
}

bool GTFParser::load_gtf(void) {
	bool is_exon;
	while (has_next()) {
		if (! get_next()) {
			//return false;
			break;
		}
		is_exon = parse_gtf_rec(line, len);
		if (is_exon)
			records.push_back(*current_record);
	}
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
	wt_exons[wt_exon.chr].push_back(wt_exon);
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

int GTFParser::search_loc(const string& chr, bool on_start, uint32_t target) {
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
