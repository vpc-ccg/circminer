#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "fastq_parser.h"
#include "filter.h"

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <FASTA> <FASTQ> <OUT> [--pe]\n", argv[0]);
		return 1;
	}
	char* ref_file = argv[1];
	char* fq_file1 = argv[2];
	char fq_file2[1000];
	bool is_pe = false;
	if (argc >= 5 and strcmp(argv[4], "--pe") == 0) {
		is_pe = true;
		get_mate_name(fq_file1, fq_file2);
	}
	
	/*
	if (bwt_index(ref_file)) {
		fprintf(stderr, "Nope! Index failed!\n");
		return 1;
	}
	else
		fprintf(stderr, "Indexed successfully!\n");
	*/

	if (bwt_load(ref_file)) {
		fprintf(stderr, "Index not loaded!\n");
		return 1;
	}
	else
		fprintf(stderr, "Index file successfully loaded!\n");

	FASTQParser fq_parser1(fq_file1);
	Record* current_record1;

	FASTQParser fq_parser2;
	Record* current_record2;
	if (is_pe) {
		fq_parser2.init(fq_file2);
	}

	FilterRead filter_read(argv[3], is_pe);

	fprintf(stderr, "Started reading FASTQ file...\n");
	int line = 0;

	while ( fq_parser1.has_next() ) { // go line by line on fastq file
		current_record1 = fq_parser1.get_next();
		if (is_pe)
			current_record2 = fq_parser2.get_next();
		if (current_record1 == NULL)	// no new line
			break;

		line++;
		if (line % 1000000 == 0)
			fprintf(stderr, "[P] %d lines\n", line);

		//fprintf(stderr, "Seq: %s, len: %d\n", current_record1->seq, current_record1->seq_len);
		//fprintf(stderr, "RC Seq: %s, len: %d\n", current_record->rcseq, current_record->seq_len);
		//int occ = find_exact_positions(current_record->seq, current_record->seq_len, find_len);

		int is_chimeric;
		if (is_pe) {
			//int is_chimeric2 = find_expanded_positions(current_record2->seq, current_record2->rcseq, current_record2->seq_len);
			is_chimeric = check_concordant_mates(current_record1, current_record2);
			filter_read.write_read(current_record1, current_record2, is_chimeric);
		}
		else {
			is_chimeric = find_expanded_positions(current_record1->seq, current_record1->rcseq, current_record1->seq_len);
			filter_read.write_read(current_record1, is_chimeric);
		}
	}

	return 0;
}
