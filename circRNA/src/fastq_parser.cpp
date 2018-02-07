#include "fastq_parser.h"

#define MAXLINESIZE 10000

FASTQParser::FASTQParser (char* filename) {
	input = fopen(filename, "r");
	if (input == NULL) {
		fprintf(stderr, "Could not open %s\n", filename);
		exit(1);
	}

	fseek(input, 0L, SEEK_END);
	file_size = ftell(input);
	fseek(input, 0L, SEEK_SET);

	max_line_size = MAXLINESIZE;
	current_record = new Record;

	current_record->rname = (char*) malloc(max_line_size);
	current_record->seq = (char*) malloc(max_line_size);
	current_record->comment = (char*) malloc(max_line_size);
	current_record->qual = (char*) malloc(max_line_size);
}

FASTQParser::~FASTQParser (void) {
	fclose(input);
}

Record* FASTQParser::get_next (void) {
	if (has_next() and read_next())
		return current_record;
	size += 1;
	return NULL;
}

bool FASTQParser::has_next (void) {
	return !feof(input);
}

bool FASTQParser::read_next (void) {
	int len;
	if ((len = getline(&current_record->rname, &max_line_size, input)) == -1)
		return false;
	
	if ((len = getline(&current_record->seq, &max_line_size, input)) == -1)
		return false;

	if ((len = getline(&current_record->comment, &max_line_size, input)) == -1)
		return false;

	if ((len = getline(&current_record->qual, &max_line_size, input)) == -1)
		return false;
	
	current_record->seq_len = strlen(current_record->seq) - 1;	// skipping newline at the end

	assert(current_record->rname[0] == '@');
	assert(current_record->comment[0] == '+');
	
	return true;
}
