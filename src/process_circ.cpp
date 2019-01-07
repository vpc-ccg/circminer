#include <cstring>

#include "process_circ.h"
#include "gene_annotation.h"
#include "hash_table.h"
#include "fastq_parser.h"
#include "match_read.h"

ProcessCirc::ProcessCirc (int last_round_num, int ws) {
	sprintf(fq_file1, "%s_%d_remain_R1.fastq", outputFilename, last_round_num);
	sprintf(fq_file2, "%s_%d_remain_R2.fastq", outputFilename, last_round_num);

	fprintf(stdout, "%s\n",fq_file1 );
	window_size = ws;
}

ProcessCirc::~ProcessCirc (void) {

}

void ProcessCirc::do_process (void) {
	/**********************/
	/**Loading Hash Table**/
	/**********************/

	int	flag;
	double tmpTime;
	checkSumLength = 5;

	char index_file [FILE_NAME_LENGTH];
	strcpy(index_file, referenceFilename);
	strcat(index_file, ".index");

	initCommon();
	
	THREAD_COUNT = 8;
	fprintf(stdout, "# Threads: %d\n", THREAD_COUNT);
	for (int i = 0; i < 255; i++)
		THREAD_ID[i] = i;

	if (!checkHashTable(index_file))
		return;

	ContigLen* orig_contig_len;
	int contig_cnt;
	if (!initLoadingHashTableMeta(index_file, &orig_contig_len, &contig_cnt))
		return;

	flag = loadHashTable ( &tmpTime );  			// Reading a fragment

	/**********************/
	/**Finised Loading HT**/
	/**********************/

	/*******************/
	/**GTF Parser Init**/
	/*******************/
	time_t pre_time, curr_time;
	double diff_time;
	time(&pre_time);

	gtf_parser.init(gtfFilename, orig_contig_len, contig_cnt);
	if (! gtf_parser.load_gtf()) {
		fprintf(stderr, "Error in reading GTF file.\n");
		exit(1);
	}
	else 
		fprintf(stdout, "GTF file successfully loaded!\n");

	time(&curr_time);
	diff_time = difftime(curr_time, pre_time);
	fprintf(stdout, "[P] Loaded GTF in %.2f sec\n", diff_time);
	pre_time = curr_time;

	/*******************/
	/**Finished GTF PI**/
	/*******************/

	contigName = getRefGenomeName();

	bool is_pe = pairedEnd;

	FASTQParser fq_parser1(fq_file1);
	Record* current_record1;

	FASTQParser fq_parser2;
	Record* current_record2;
	if (is_pe) {
		fq_parser2.init(fq_file2);
	}
	
	int line = 0;
	while ( (current_record1 = fq_parser1.get_next()) != NULL ) { // go line by line on fastq file
		if (is_pe)
			current_record2 = fq_parser2.get_next();
		if (current_record1 == NULL)	// no new line
			break;

		line++;

		fprintf(stdout, "Line: %d\n", line);
		MatchedRead mr = *(current_record1->mr);
		fprintf(stdout, "%s\n%s\n", current_record1->seq, current_record2->seq);
		fprintf(stdout, "%s\t%s\t%u\t%u\t%d\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%d\t%d\t%d\t%d\n", 
										current_record1->rname, mr.chr.c_str(), 
										mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1,  
										mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, 
										mr.tlen, mr.junc_num, mr.gm_compatible, mr.type);

		bool r1_partial = mr.mlen_r1 < mr.mlen_r2;
		char* remain_seq = (r1_partial) ? ((mr.r1_forward) ? current_record1->seq : current_record1->rcseq) : 
										  ((mr.r2_forward) ? current_record2->seq : current_record2->rcseq) ;
		uint32_t qspos = (r1_partial) ? (((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (1) : (mr.qepos_r1 + 1)) :
										(((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (1) : (mr.qepos_r2 + 1)) ;

		uint32_t qepos = (r1_partial) ? (((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (mr.qspos_r1 - 1) : (current_record1->seq_len)) :
										(((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (mr.qspos_r2 - 1) : (current_record2->seq_len)) ;

		const IntervalInfo<GeneInfo>* gene_info = gtf_parser.get_gene_overlap(mr.spos_r1, false);
		fprintf(stderr, "# Gene overlaps: %d\n", gene_info->seg_list.size());

		for (int i = 0; i < gene_info->seg_list.size(); i++) {
			uint32_t gene_len = gene_info->seg_list[i].end - gene_info->seg_list[i].start + 1;
			char gene_seq[gene_len + 1];
			gene_seq[gene_len] = '\0';
			pac2char(gene_info->seg_list[i].start, gene_len, gene_seq);
			fprintf(stderr, "Gene: %s\n", gene_seq);

			RegionalHashTable regional_ht(window_size);
			regional_ht.create_table(gene_seq, 0, gene_len);

			fprintf(stdout, "R%d partial: [%d-%d]\n", (int) (!r1_partial) + 1, qspos, qepos);
			fprintf(stdout, "%s\n", remain_seq);

			for (int i = qspos - 1; i <= qepos - window_size; i += 3) {
				GIList* gl = regional_ht.find_hash(regional_ht.hash_val(remain_seq + i));
				if (gl == NULL) {
					fprintf(stdout, "Hash val not found!!!\n");
				}
				fprintf(stdout, "Occ: %d\n", gl->cnt);
				for (int j = 0; j < gl->cnt; j++) {
					fprintf(stdout, "%d\t", gl->locs[j].info);
				}
				fprintf(stdout, "\n");

			}
		}
	}
}

int ProcessCirc::get_exact_locs_hash (char* seq, uint32_t qspos, uint32_t qepos) {

}