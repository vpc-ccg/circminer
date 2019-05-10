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
		
		for (int i = 0; i < BLOCKSIZE; ++i) {
			free(current_record[i]->rname);
			free(current_record[i]->seq);
			free(current_record[i]->rcseq);
			free(current_record[i]->comment);
			free(current_record[i]->qual);
			delete current_record[i];
		}
	}
}

void FASTQParser::init (char* filename) {
	input = open_file(filename, "r");

	max_line_size = MAXLINESIZE;
	set_comp();

	for (int i = 0; i < BLOCKSIZE; ++i) {
		current_record[i] = new Record;

		current_record[i]->rname = (char*) malloc(max_line_size);
		current_record[i]->seq = (char*) malloc(max_line_size);
		current_record[i]->rcseq = (char*) malloc(max_line_size);
		current_record[i]->comment = (char*) malloc(max_line_size);
		current_record[i]->qual = (char*) malloc(max_line_size);
	}

	curr_read = 0;
	filled_size = 0;
}

Record* FASTQParser::get_next (void) {
	if (curr_read < filled_size) {
		return current_record[curr_read++];
	}
	else if (read_block()) {
		return current_record[curr_read++];
	}
	return NULL;
}

bool FASTQParser::has_next (void) {
	char c = fgetc(input);
	if (c == EOF)
		return false;
	
	assert(c == '@');	// ensure FASTQ format 
	return true;
}

bool FASTQParser::read_block (void) {
	curr_read = 0;
	filled_size = 0;
	for (int i = 0; i < BLOCKSIZE; ++i) {
		if (!has_next())
			return filled_size > 0;

		int rname_len = getline(&current_record[i]->rname, &max_line_size, input);
		rname_len = extract_map_info(current_record[i]->rname, i);

		getline(&current_record[i]->seq, &max_line_size, input);
		current_record[i]->seq_len = strlen(current_record[i]->seq) - 1;	// skipping newline at the end

		getline(&current_record[i]->comment, &max_line_size, input);
		assert(current_record[i]->comment[0] == '+');

		getline(&current_record[i]->qual, &max_line_size, input);
		


		if (current_record[i]->rname[rname_len - 3] == '/')
			current_record[i]->rname[rname_len - 3] = '\0';
		else
			current_record[i]->rname[rname_len - 1] = '\0';

		//current_record[i]->seq[current_record[i]->seq_len] = '\0';

		set_reverse_comp(i);
		++filled_size;
	}

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

void FASTQParser::set_reverse_comp (int r_ind) {
	if (current_record[r_ind] == NULL) {
		fprintf(stderr, "Error: No read loaded to be reverse completed\n");
		return;
	}

	int len = current_record[r_ind]->seq_len;
	for (int i = len-1; i >= 0; --i) {
		current_record[r_ind]->rcseq[len-i-1] = comp[current_record[r_ind]->seq[i]];
	}
	current_record[r_ind]->rcseq[len] = '\0';
}

int FASTQParser::extract_map_info(char* str, int r_ind) {
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
	current_record[r_ind]->rname[rname_len] = '\0';

	fill_map_info(i, r_ind);
	return rname_len;
}

void FASTQParser::fill_map_info(int cnt, int r_ind) {
	//assert(cnt == 1 or cnt == FQCOMMENTCNT);
	
	if (cnt != FQCOMMENTCNT) {
		current_record[r_ind]->mr->type = NOPROC_NOMATCH;
		current_record[r_ind]->mr->tlen = INF;
		current_record[r_ind]->mr->junc_num = 0;
		current_record[r_ind]->mr->gm_compatible = false;
	}
	else {
		char* stop_string;
		int base = 10;
		current_record[r_ind]->mr->type 	= atoi(tokens[1]);

		if (current_record[r_ind]->mr->type == CONCRD or current_record[r_ind]->mr->type == DISCRD or current_record[r_ind]->mr->type == CHIORF or 
			current_record[r_ind]->mr->type == CHIBSJ or current_record[r_ind]->mr->type == CHI2BSJ) {
			current_record[r_ind]->mr->chr_r1 		= tokens[2];	
			current_record[r_ind]->mr->spos_r1 	= strtoul(tokens[3], &stop_string, base);
			current_record[r_ind]->mr->epos_r1 	= strtoul(tokens[4], &stop_string, base);
			current_record[r_ind]->mr->mlen_r1 	= atoi(tokens[5]);
			current_record[r_ind]->mr->qspos_r1 	= strtoul(tokens[6], &stop_string, base);
			current_record[r_ind]->mr->qepos_r1 	= strtoul(tokens[7], &stop_string, base);
			current_record[r_ind]->mr->r1_forward 	= (tokens[8][0] == '+');
			current_record[r_ind]->mr->ed_r1 		= atoi(tokens[9]);

			current_record[r_ind]->mr->chr_r2 		= tokens[10];
			current_record[r_ind]->mr->spos_r2 	= strtoul(tokens[11], &stop_string, base);
			current_record[r_ind]->mr->epos_r2 	= strtoul(tokens[12], &stop_string, base);
			current_record[r_ind]->mr->mlen_r2 	= atoi(tokens[13]);
			current_record[r_ind]->mr->qspos_r2 	= strtoul(tokens[14], &stop_string, base);
			current_record[r_ind]->mr->qepos_r2 	= strtoul(tokens[15], &stop_string, base);
			current_record[r_ind]->mr->r2_forward 	= (tokens[16][0] == '+');
			current_record[r_ind]->mr->ed_r2 		= atoi(tokens[17]);
			
			current_record[r_ind]->mr->tlen 		= atoi(tokens[18]);
			current_record[r_ind]->mr->junc_num 	= strtoul(tokens[19], &stop_string, base);
			current_record[r_ind]->mr->gm_compatible = (tokens[20][0] == '1');
			current_record[r_ind]->mr->contig_num	= atoi(tokens[21]);
		}
		else {
			current_record[r_ind]->mr->chr_r1 		= "-";	
			current_record[r_ind]->mr->spos_r1 	= 0;
			current_record[r_ind]->mr->epos_r1 	= 0;
			current_record[r_ind]->mr->mlen_r1 	= 0;
			current_record[r_ind]->mr->qspos_r1 	= 0;
			current_record[r_ind]->mr->qepos_r1 	= 0;
			current_record[r_ind]->mr->r1_forward 	= true;
			current_record[r_ind]->mr->ed_r1 		= maxEd+1;

			current_record[r_ind]->mr->chr_r2 		= "-";
			current_record[r_ind]->mr->spos_r2 	= 0;
			current_record[r_ind]->mr->epos_r2 	= 0;
			current_record[r_ind]->mr->mlen_r2 	= 0;
			current_record[r_ind]->mr->qspos_r2 	= 0;
			current_record[r_ind]->mr->qepos_r2 	= 0;
			current_record[r_ind]->mr->r2_forward 	= true;
			current_record[r_ind]->mr->ed_r2 		= maxEd+1;
			
			current_record[r_ind]->mr->tlen 		= INF;
			current_record[r_ind]->mr->junc_num 	= 0;
			current_record[r_ind]->mr->gm_compatible = false;
			current_record[r_ind]->mr->contig_num	= 0;
		}
	}
}
