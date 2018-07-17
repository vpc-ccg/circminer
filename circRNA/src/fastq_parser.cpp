#include "fastq_parser.h"

#define MAXLINESIZE 400

FASTQParser::FASTQParser (bool read_state) : read_state(read_state) { 
	input = NULL; 
}

FASTQParser::FASTQParser (char* filename, bool read_state) : read_state(read_state) {
	init(filename);
}

FASTQParser::~FASTQParser (void) {
	close_file(input);
}

void FASTQParser::init (char* filename) {
	input = open_file(filename, "r");

	max_line_size = MAXLINESIZE;
	current_record = new Record;
	set_comp();

	current_record->rname = (char*) malloc(max_line_size);
	current_record->seq = (char*) malloc(max_line_size);
	current_record->rcseq = (char*) malloc(max_line_size);
	current_record->comment = (char*) malloc(max_line_size);
	current_record->qual = (char*) malloc(max_line_size);
}

Record* FASTQParser::get_next (void) {
	if (has_next() and read_next())
		return current_record;
	size += 1;
	return NULL;
}

bool FASTQParser::has_next (void) {
	char c = fgetc(input);
	if (c == EOF)
		return false;
	
	assert(c == '@');	// ensure FASTQ format 
	return true;
}

bool FASTQParser::read_next (void) {
	current_record->rname[0] = '@';	// already read it in has_next()
	char * skip_at = current_record->rname + 1;

	getline(&skip_at, &max_line_size, input);
	getline(&current_record->seq, &max_line_size, input);
	int comment_len = getline(&current_record->comment, &max_line_size, input);
	getline(&current_record->qual, &max_line_size, input);
	
	/** read state from comment **/
	assert(current_record->comment[0] == '+');
	if (read_state and comment_len > 2) {
		current_record->state = current_record->comment[1] - '0';
	}
	else 
		current_record->state = ORPHAN;
	current_record->comment[1] = '\0';
	/****/

	current_record->seq_len = strlen(current_record->seq) - 1;	// skipping newline at the end
	//current_record->seq[current_record->seq_len] = '\0';

	set_reverse_comp();

	return true;
}

void FASTQParser::set_comp (void) {
	comp['A'] = 'T';
	comp['C'] = 'G';
	comp['G'] = 'C';
	comp['T'] = 'A';
	comp['N'] = 'N';

	comp['a'] = 'T';
	comp['c'] = 'G';
	comp['g'] = 'C';
	comp['t'] = 'A';
	comp['n'] = 'N';
}

void FASTQParser::set_reverse_comp (void) {
	if (current_record == NULL) {
		fprintf(stderr, "No read loaded to be reverse completed\n");
		return;
	}

	int len = current_record->seq_len;
	char nt;
	for (int i = len-1; i >= 0; --i) {
		current_record->rcseq[len-i-1] = comp[current_record->seq[i]];
	}
	current_record->rcseq[len] = '\0';
}
