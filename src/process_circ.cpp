#include <cstring>

#include "process_circ.h"
#include "gene_annotation.h"
#include "hash_table.h"

ProcessCirc::ProcessCirc (char* fname, int ws) {
	circ_file = open_file(fname, "r");
	window_size = ws;

	max_line_size = MAXLINESIZE;
	line = (char*) malloc(max_line_size);
}

ProcessCirc::~ProcessCirc (void) {
	close_file(circ_file);
}

bool ProcessCirc::read_next (void) {
	int ret = getline(&line, &max_line_size, circ_file);
	return ret > 0;
}

void ProcessCirc::tokenize (void) {
	// Returns first token  
	char *token = strtok(line, " "); 
    
	// Keep printing tokens while one of the 
	// delimiters present in str[]. 
	int i = 0;
	while (token != NULL) 
	{ 
		strcpy(tokens[i], token);
		token = strtok(NULL, " "); 
		++i;
	}

	char* stop_string;
	int base = 10;
	mr.type 	= atoi(tokens[0]);
	mr.chr 		= tokens[1];
	mr.spos_r1 	= strtoul(tokens[2], &stop_string, base);
	mr.epos_r1 	= strtoul(tokens[3], &stop_string, base);
	mr.mlen_r1 	= atoi(tokens[4]);
	mr.spos_r2 	= strtoul(tokens[5], &stop_string, base);
	mr.epos_r2 	= strtoul(tokens[6], &stop_string, base);
	mr.mlen_r2 	= atoi(tokens[7]);
	mr.tlen 	= atoi(tokens[8]);
	mr.junc_num = strtoul(tokens[9], &stop_string, base);
	mr.gm_compatible = (tokens[10][0] == '1');
}

void ProcessCirc::do_process (void) {
	while (read_next()) {
		tokenize();

		const IntervalInfo<GeneInfo>* gene_info = gtf_parser.get_gene_overlap(mr.spos_r1, false);

		RegionalHashTable regional_ht(window_size);
		//pac2char();
		//regional_ht.create_table(seq, gene_info.start, gene_info.length());
	}
}