#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "common.h"
#include "commandline_parser.h"
#include "fastq_parser.h"
#include "filter.h"

char versionNumberMajor[10] = "0";
char versionNumberMinor[10] = "1";

using namespace std;

int main(int argc, char **argv) {
	//if (argc < 4) {
	//	fprintf(stderr, "Usage: %s <FASTA> <FASTQ> <OUT> [--pe]\n", argv[0]);
	//	return 1;
	//}
	
	int exit = parse_command( argc, argv );
	if (exit == 1)
		return 0;

	char* ref_file = referenceFilename;
	char* fq_file1 = fastqFilename;
	char fq_file2[1000];
	bool is_pe = pairedEnd;
	if (pairedEnd) {
		get_mate_name(fq_file1, fq_file2);
	}
	
	if (bwt_load(ref_file)) {
		fprintf(stdout, "Index not found!\n");
		//if (bwt_index(ref_file)) {
		//	fprintf(stderr, "Indexing failed!\n");
		//	return 1;
		//}
		//else
		//	fprintf(stderr, "Indexed successfully!\n");
	}
	else
		fprintf(stdout, "Index file successfully loaded!\n");

	FASTQParser fq_parser1(fq_file1);
	Record* current_record1;

	FASTQParser fq_parser2;
	Record* current_record2;
	if (is_pe) {
		fq_parser2.init(fq_file2);
	}

	FilterRead filter_read(outputFilename, is_pe);

	fprintf(stdout, "Started reading FASTQ file...\n");
	int line = 0;

	while ( fq_parser1.has_next() ) { // go line by line on fastq file
		current_record1 = fq_parser1.get_next();
		if (is_pe)
			current_record2 = fq_parser2.get_next();
		if (current_record1 == NULL)	// no new line
			break;

		line++;
		if (line % 100000 == 0)
			fprintf(stdout, "[P] %d lines\n", line);

		//fprintf(stderr, "Seq: %s, len: %d\n", current_record1->seq, current_record1->seq_len);
		//fprintf(stderr, "RC Seq: %s, len: %d\n", current_record->rcseq, current_record->seq_len);
		//int occ = find_exact_positions(current_record->seq, current_record->seq_len, find_len);

		int is_chimeric;
		if (is_pe) {
			//is_chimeric = check_concordant_mates(current_record1, current_record2);
			is_chimeric = filter_read.process_read(current_record1, current_record2, kmer);
			filter_read.write_read3(current_record1, current_record2, is_chimeric);
		}
		else {
			is_chimeric = filter_read.process_read(current_record1);
			filter_read.write_read(current_record1, is_chimeric);
		}
	}

	return 0;
}
