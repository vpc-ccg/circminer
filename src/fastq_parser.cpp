#include <cstdio> 
#include <cstring>
#include <zlib.h>

#include "fastq_parser.h"

FASTQParser::FASTQParser (void) { 
	init();
}

FASTQParser::FASTQParser (char* filename) {
	init();
	reset(filename);
}

FASTQParser::~FASTQParser (void) {
	free(zbuffer);
	for (int i = 0; i < BLOCKSIZE; ++i) {
		free(current_record[i].rname);
		free(current_record[i].seq);
		free(current_record[i].rcseq);
		free(current_record[i].comment);
		free(current_record[i].qual);
	}
	// free(current_record);
	delete[] current_record;
	
	finalize();
}

void FASTQParser::init (void) {
	input = NULL;
	gzinput = Z_NULL;
	mate_q = NULL;

	max_line_size = MAXLINESIZE;
	set_comp();

	zbuffer = (char*) malloc(BUFFSIZE);
	// current_record = (Record*) malloc(BLOCKSIZE * sizeof(Record));
	current_record = new Record[BLOCKSIZE];
	for (int i = 0; i < BLOCKSIZE; ++i) {
		current_record[i].rname = (char*) malloc(max_line_size);
		current_record[i].seq = (char*) malloc(max_line_size);
		current_record[i].rcseq = (char*) malloc(max_line_size);
		current_record[i].comment = (char*) malloc(max_line_size);
		current_record[i].qual = (char*) malloc(max_line_size);
	}
}

void FASTQParser::reset (char* filename) {
	finalize();

	gzinput = open_gzfile(filename, "r");

	//input = open_file(filename, "r");
	
	buff_pos = 0;
	buff_size = 0;

	curr_read = 0;
	filled_size = 0;
}

void FASTQParser::finalize (void) {
	if (gzinput != NULL) {
		close_gzfile(gzinput);
		gzinput = Z_NULL;
	}
	if (input != NULL) {
		close_file(input);
		input = NULL;
	}
}

void FASTQParser::set_mate(FASTQParser* mq) {
	mate_q = mq;
}

void FASTQParser::read_buffer() {
	buff_size = gzread(gzinput, zbuffer, BUFFSIZE);
	buff_pos = 0;

	if (buff_size == 0 and gzeof(gzinput) == 0) {
		buff_size = -1;
	}
	if (buff_size < 0) {
		int err;
		fprintf(stderr, "gzread error: %s\n", gzerror(gzinput, &err));
		exit(1);
	}
}

uint32_t FASTQParser::read_line(char** seq) {
	char cur;

	uint32_t i = 0;
	while (true) {
		if (buff_pos >= buff_size) {
			read_buffer();
			if (buff_size == 0)
				return 0;
		}

		cur = zbuffer[buff_pos++];
		if (cur == '\n') {
			(*seq)[i] = '\0';
			return i;
		}

		(*seq)[i++] = cur;
	}
}

bool FASTQParser::read_block (void) {
	curr_read = 0;
	filled_size = 0;
	for (int i = 0; i < BLOCKSIZE; ++i) {
		if (has_next()) {
			//int rname_len = getline(&current_record[i].rname, &max_line_size, input);
			int rname_len = read_line(&current_record[i].rname);
			rname_len = extract_map_info(current_record[i].rname, i);

			//getline(&current_record[i].seq, &max_line_size, input);
			current_record[i].seq_len = read_line(&current_record[i].seq);

			//getline(&current_record[i].comment, &max_line_size, input);
			read_line(&current_record[i].comment);
			assert(current_record[i].comment[0] == '+');

			//getline(&current_record[i].qual, &max_line_size, input);
			read_line(&current_record[i].qual);

			set_reverse_comp(i);
			++filled_size;
		}
		else {
			return filled_size > 0;
		}
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
	int len = current_record[r_ind].seq_len;
	for (int i = len-1; i >= 0; --i) {
		current_record[r_ind].rcseq[len-i-1] = comp[current_record[r_ind].seq[i]];
	}
	current_record[r_ind].rcseq[len] = '\0';
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

	int rname_len = strlen(tokens[0]) + 1;
	current_record[r_ind].rname[rname_len] = '\0';

	if (current_record[r_ind].rname[rname_len - 3] == '/')
		current_record[r_ind].rname[rname_len - 3] = '\0';

	fill_map_info(i, r_ind);
	return rname_len;
}

void FASTQParser::fill_map_info(int cnt, int r_ind) {
	//assert(cnt == 1 or cnt == FQCOMMENTCNT);
	
	if (cnt != FQCOMMENTCNT) {
		current_record[r_ind].mr->type = NOPROC_NOMATCH;
		current_record[r_ind].mr->tlen = INF;
		current_record[r_ind].mr->junc_num = 0;
		current_record[r_ind].mr->gm_compatible = false;
	}
	else {
		char* stop_string;
		int base = 10;
		current_record[r_ind].mr->type 	= atoi(tokens[1]);

		if (current_record[r_ind].mr->type == CONCRD or current_record[r_ind].mr->type == DISCRD or current_record[r_ind].mr->type == CHIORF or 
			current_record[r_ind].mr->type == CHIBSJ or current_record[r_ind].mr->type == CHI2BSJ or current_record[r_ind].mr->type == CONGNM  or 
			current_record[r_ind].mr->type == CONGEN) {
			current_record[r_ind].mr->chr_r1 	= tokens[2];	
			current_record[r_ind].mr->spos_r1	= strtoul(tokens[3], &stop_string, base);
			current_record[r_ind].mr->epos_r1	= strtoul(tokens[4], &stop_string, base);
			current_record[r_ind].mr->mlen_r1	= atoi(tokens[5]);
			current_record[r_ind].mr->qspos_r1	= strtoul(tokens[6], &stop_string, base);
			current_record[r_ind].mr->qepos_r1	= strtoul(tokens[7], &stop_string, base);
			current_record[r_ind].mr->r1_forward	= (tokens[8][0] == '+');
			current_record[r_ind].mr->ed_r1		= atoi(tokens[9]);

			current_record[r_ind].mr->chr_r2	= tokens[10];
			current_record[r_ind].mr->spos_r2	= strtoul(tokens[11], &stop_string, base);
			current_record[r_ind].mr->epos_r2	= strtoul(tokens[12], &stop_string, base);
			current_record[r_ind].mr->mlen_r2	= atoi(tokens[13]);
			current_record[r_ind].mr->qspos_r2	= strtoul(tokens[14], &stop_string, base);
			current_record[r_ind].mr->qepos_r2	= strtoul(tokens[15], &stop_string, base);
			current_record[r_ind].mr->r2_forward	= (tokens[16][0] == '+');
			current_record[r_ind].mr->ed_r2	= atoi(tokens[17]);
			
			current_record[r_ind].mr->tlen		= atoi(tokens[18]);
			current_record[r_ind].mr->junc_num	= strtoul(tokens[19], &stop_string, base);
			current_record[r_ind].mr->gm_compatible	= (tokens[20][0] == '1');
			current_record[r_ind].mr->contig_num		= atoi(tokens[21]);
		}
		else {
			current_record[r_ind].mr->chr_r1	= "-";	
			current_record[r_ind].mr->spos_r1	= 0;
			current_record[r_ind].mr->epos_r1	= 0;
			current_record[r_ind].mr->mlen_r1	= 0;
			current_record[r_ind].mr->qspos_r1	= 0;
			current_record[r_ind].mr->qepos_r1	= 0;
			current_record[r_ind].mr->r1_forward	= true;
			current_record[r_ind].mr->ed_r1 		= maxEd+1;

			current_record[r_ind].mr->chr_r2	= "-";
			current_record[r_ind].mr->spos_r2	= 0;
			current_record[r_ind].mr->epos_r2	= 0;
			current_record[r_ind].mr->mlen_r2	= 0;
			current_record[r_ind].mr->qspos_r2	= 0;
			current_record[r_ind].mr->qepos_r2	= 0;
			current_record[r_ind].mr->r2_forward	= true;
			current_record[r_ind].mr->ed_r2	= maxEd+1;
			
			current_record[r_ind].mr->tlen		= INF;
			current_record[r_ind].mr->junc_num	= 0;
			current_record[r_ind].mr->gm_compatible	= false;
			current_record[r_ind].mr->contig_num		= 0;
		}
	}
}
