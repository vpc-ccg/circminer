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
	set_comp();

	current_record->rname = (char*) malloc(max_line_size);
	current_record->seq = (char*) malloc(max_line_size);
	current_record->rcseq = (char*) malloc(max_line_size);
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
	
	assert(current_record->rname[0] == '@');
	assert(current_record->comment[0] == '+');
	
	current_record->seq_len = strlen(current_record->seq) - 1;	// skipping newline at the end
	//current_record->seq[current_record->seq_len] = '\0';

	set_reverse_comp();

	return true;
}

void FASTQParser::set_comp (void) {
	comp['A'-'A'] = 'T' - 'A';
	comp['C'-'A'] = 'G' - 'C';
	comp['G'-'A'] = 'C' - 'G';
	comp['T'-'A'] = 'A' - 'T';
	comp['N'-'A'] = 'N' - 'N';
}

char FASTQParser::get_comp (char nt) {
	if (nt > 'a')
		return nt + comp[nt - 'a'];
	return nt + comp[nt - 'A'];
}

void FASTQParser::set_reverse_comp (void) {
	if (current_record == NULL) {
		fprintf(stderr, "No read loaded\n");
		return;
	}

	int len = current_record->seq_len;
	char nt;
	for (int i = len-1; i >= 0; --i) {
		nt = get_comp(current_record->seq[i]);
		current_record->rcseq[len-i-1] = nt;
	}
	current_record->rcseq[len] = '\0';
}
