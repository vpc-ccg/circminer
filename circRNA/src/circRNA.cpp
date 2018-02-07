#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "bwt.h"
#include "fastq_parser.h"

int find_exact_positions(const char* rseq, int rseq_len, int window_size) {
	bwtint_t sp, ep;
	int occ = 0;

	for (int i = 0; i < rseq_len; i+=window_size) {
		occ = get_exact_locs(rseq + i, window_size, &sp, &ep);
		fprintf(stderr, "Start pos: %d\n", i);
		fprintf(stderr, "Number of matches: %d\n", occ);

		if (occ > 0)
			locate_match(&sp, &ep, window_size);
	}

	return occ;
}

int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <FASTA> <FASTQ> <OUT>\n", argv[0]);
		return 1;
	}
	char* ref_file = argv[1];
	char* fq_file1 = argv[2];
	
	char ignore_file[200], keep_file[200];
	strcpy (ignore_file, argv[3]);
	strcpy (keep_file, argv[3]);

	strcat (ignore_file, ".ignore.fastq");
	strcat (keep_file, ".keep.fastq");

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

	//FILE* fqin1 = fopen(fq_file1, "r");
	FILE* ignore_out = fopen(ignore_file, "w");
	FILE* keep_out = fopen(keep_file, "w");

	FASTQParser fq_parser(fq_file1);
	Record* current_record;
	printf("Started reading FASTQ file\n");
	int line = 0;

	while ( fq_parser.has_next() ) { // go line by line on fastq file
		current_record = fq_parser.get_next();
		if (current_record == NULL)	// no new line
			break;

		fprintf(stderr, "line: %d\n", line);
		line += 4;

		fprintf(stderr, "Seq: %s, len: %d\n", current_record->seq, current_record->seq_len);
		int occ = find_exact_positions(current_record->seq, current_record->seq_len, current_record->seq_len / 2);

		if (occ >= 1) {
			fprintf(ignore_out, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
		}
		else {
			fprintf(keep_out, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
		}
	}

	return 0;
}
