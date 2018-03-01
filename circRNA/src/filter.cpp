#include "filter.h"

FilterRead::FilterRead (char* save_fname, bool pe) {
	is_pe = pe;

	char ignore_file1[1000];
	char keep_file1[1000];
	char chimeric_file1[1000];
	char p_unmap_file1[1000];
	char unmap_file1[1000];
	strcpy (ignore_file1, save_fname);
	strcpy (keep_file1 , save_fname);
	strcpy (chimeric_file1 , save_fname);
	strcpy (p_unmap_file1 , save_fname);
	strcpy (unmap_file1 , save_fname);

	strcat (ignore_file1, ".ignore_R1.fastq");
	strcat (keep_file1, ".keep_R1.fastq");
	strcat (chimeric_file1, ".chimeric_R1.fastq");
	strcat (p_unmap_file1, ".partly_unmappable_R1.fastq");
	strcat (unmap_file1, ".unmappable_R1.fastq");

	ignore_r1 = fopen(ignore_file1, "w");
	keep_r1   = fopen(keep_file1, "w");
	chimeric_r1   = fopen(chimeric_file1, "w");
	partly_unmappable_r1   = fopen(p_unmap_file1, "w");
	unmappable_r1   = fopen(unmap_file1, "w");

	if (is_pe) {
		char ignore_file2[1000];
		char keep_file2[1000];
		char chimeric_file2[1000];
		char p_unmap_file2[1000];
		char unmap_file2[1000];
		strcpy (ignore_file2, save_fname);
		strcpy (keep_file2 , save_fname);
		strcpy (chimeric_file2 , save_fname);
		strcpy (p_unmap_file2 , save_fname);
		strcpy (unmap_file2 , save_fname);

		strcat (ignore_file2, ".ignore_R2.fastq");
		strcat (keep_file2, ".keep_R2.fastq");
		strcat (chimeric_file2, ".chimeric_R2.fastq");
		strcat (p_unmap_file2, ".partly_unmappable_R2.fastq");
		strcat (unmap_file2, ".unmappable_R2.fastq");

		ignore_r2 = fopen(ignore_file2, "w");
		keep_r2   = fopen(keep_file2, "w");
		chimeric_r2   = fopen(chimeric_file2, "w");
		partly_unmappable_r2   = fopen(p_unmap_file2, "w");
		unmappable_r2   = fopen(unmap_file2, "w");
	}
}

FilterRead::~FilterRead (void) {
	if (ignore_r1 != NULL)
		fclose(ignore_r1);
	if (keep_r1 != NULL)
		fclose(keep_r1);
	if (chimeric_r1 != NULL)
		fclose(chimeric_r1);
	if (partly_unmappable_r1 != NULL)
		fclose(partly_unmappable_r1);
	if (unmappable_r1 != NULL)
		fclose(unmappable_r1);

	if (ignore_r2 != NULL)
		fclose(ignore_r2);
	if (keep_r2 != NULL)
		fclose(keep_r2);
	if (chimeric_r2 != NULL)
		fclose(chimeric_r2);
	if (partly_unmappable_r2 != NULL)
		fclose(partly_unmappable_r2);
	if (unmappable_r2 != NULL)
		fclose(unmappable_r2);
}

int FilterRead::process_read (Record* current_record) {
	return find_expanded_positions(current_record->seq, current_record->rcseq, current_record->seq_len);
}

int FilterRead::process_read (Record* current_record1, Record* current_record2) {
	return check_concordant_mates_expand(current_record1, current_record2);
}

// write reads SE mode
void FilterRead::write_read (Record* current_record, int is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
	}
	else {
		fprintf(keep_r1, "%s%s%s%s", current_record->rname, current_record->seq, current_record->comment, current_record->qual);
	}
}

// write reads PE mode
void FilterRead::write_read (Record* current_record1, Record* current_record2, int is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
}

// write reads PE mode --detailed
void FilterRead::write_read2 (Record* current_record1, Record* current_record2, int is_chimeric) {
	if (is_chimeric == 0) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (is_chimeric == 1) {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else {
		fprintf(chimeric_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(chimeric_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
}

// write reads PE mode --detailed / 5 output files
void FilterRead::write_read3 (Record* current_record1, Record* current_record2, int state) {
	if (state == 0 or state == 1) {
		fprintf(ignore_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(ignore_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 2){
		fprintf(chimeric_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(chimeric_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 3) {
		fprintf(keep_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(keep_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 4) {
		fprintf(partly_unmappable_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(partly_unmappable_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
	else if (state == 5) {
		fprintf(unmappable_r1, "%s%s%s%s", current_record1->rname, current_record1->seq, current_record1->comment, current_record1->qual);
		fprintf(unmappable_r2, "%s%s%s%s", current_record2->rname, current_record2->seq, current_record2->comment, current_record2->qual);
	}
}

