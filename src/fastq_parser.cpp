#include <cstdio> 
#include <cstring>
#include "fastq_parser.h"

FASTQParser::FASTQParser () { 
	input = NULL; 
}

FASTQParser::FASTQParser (char* filename) {
	init(filename);
}

FASTQParser::~FASTQParser (void) {
	if (input != NULL) {
		close_file(input);
		
		free(current_record->rname);
		free(current_record->seq);
		free(current_record->rcseq);
		free(current_record->comment);
		free(current_record->qual);

		delete current_record;
	}
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
	int rname_len = getline(&current_record->rname, &max_line_size, input);
	rname_len = extract_map_info(current_record->rname);

	getline(&current_record->seq, &max_line_size, input);
	current_record->seq_len = strlen(current_record->seq) - 1;	// skipping newline at the end

	getline(&current_record->comment, &max_line_size, input);
	assert(current_record->comment[0] == '+');

	getline(&current_record->qual, &max_line_size, input);
	


	if (current_record->rname[rname_len - 3] == '/')
		current_record->rname[rname_len - 3] = '\0';
	else
		current_record->rname[rname_len - 1] = '\0';

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
	for (int i = len-1; i >= 0; --i) {
		current_record->rcseq[len-i-1] = comp[current_record->seq[i]];
	}
	current_record->rcseq[len] = '\0';
}

int FASTQParser::extract_map_info(char* str) {
	// Returns first token  
	char *token = strtok(str, " "); 
    
	// Keep printing tokens while one of the 
	// delimiters present in str[]. 
	int i = 0;
	while (token != NULL) 
	{ 
		strcpy(tokens[i], token);
		token = strtok(NULL, " "); 
		++i;
	}

	// i == 1 iff there is no comment in the line
	int rname_len = (i == 1) ? strlen(tokens[0]) : strlen(tokens[0]) + 1;
	current_record->rname[rname_len] = '\0';

	fill_map_info(i);
	return rname_len;
}

void FASTQParser::fill_map_info(int cnt) {
	//assert(cnt == 1 or cnt == FQCOMMENTCNT);
	
	if (cnt != FQCOMMENTCNT) {
		current_record->mr->type = NOPROC_NOMATCH;
		current_record->mr->tlen = INF;
		current_record->mr->junc_num = 0;
		current_record->mr->gm_compatible = false;
	}
	else {
		char* stop_string;
		int base = 10;
		current_record->mr->type 	= atoi(tokens[1]);

		if (current_record->mr->type == CONCRD or current_record->mr->type == DISCRD or current_record->mr->type == CHIORF or 
			current_record->mr->type == CHIBSJ or current_record->mr->type == CHI2BSJ) {
			current_record->mr->chr_r1 		= tokens[2];	
			current_record->mr->spos_r1 	= strtoul(tokens[3], &stop_string, base);
			current_record->mr->epos_r1 	= strtoul(tokens[4], &stop_string, base);
			current_record->mr->mlen_r1 	= atoi(tokens[5]);
			current_record->mr->qspos_r1 	= strtoul(tokens[6], &stop_string, base);
			current_record->mr->qepos_r1 	= strtoul(tokens[7], &stop_string, base);
			current_record->mr->r1_forward 	= (tokens[8][0] == '+');
			current_record->mr->ed_r1 		= atoi(tokens[9]);

			current_record->mr->chr_r2 		= tokens[10];
			current_record->mr->spos_r2 	= strtoul(tokens[11], &stop_string, base);
			current_record->mr->epos_r2 	= strtoul(tokens[12], &stop_string, base);
			current_record->mr->mlen_r2 	= atoi(tokens[13]);
			current_record->mr->qspos_r2 	= strtoul(tokens[14], &stop_string, base);
			current_record->mr->qepos_r2 	= strtoul(tokens[15], &stop_string, base);
			current_record->mr->r2_forward 	= (tokens[16][0] == '+');
			current_record->mr->ed_r2 		= atoi(tokens[17]);
			
			current_record->mr->tlen 		= atoi(tokens[18]);
			current_record->mr->junc_num 	= strtoul(tokens[19], &stop_string, base);
			current_record->mr->gm_compatible = (tokens[20][0] == '1');
			current_record->mr->contig_num	= atoi(tokens[21]);
		}
		else {
			current_record->mr->chr_r1 		= "-";	
			current_record->mr->spos_r1 	= 0;
			current_record->mr->epos_r1 	= 0;
			current_record->mr->mlen_r1 	= 0;
			current_record->mr->qspos_r1 	= 0;
			current_record->mr->qepos_r1 	= 0;
			current_record->mr->r1_forward 	= true;
			current_record->mr->ed_r1 		= maxEd+1;

			current_record->mr->chr_r2 		= "-";
			current_record->mr->spos_r2 	= 0;
			current_record->mr->epos_r2 	= 0;
			current_record->mr->mlen_r2 	= 0;
			current_record->mr->qspos_r2 	= 0;
			current_record->mr->qepos_r2 	= 0;
			current_record->mr->r2_forward 	= true;
			current_record->mr->ed_r2 		= maxEd+1;
			
			current_record->mr->tlen 		= INF;
			current_record->mr->junc_num 	= 0;
			current_record->mr->gm_compatible = false;
			current_record->mr->contig_num	= 0;
		}
	}
}
