//extern "C" {
//#include "mrsfast/RefGenome.h"
//}

#include "genome.h"
#include "common.h"
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

#define CONTIG_SIZE 21

GenomePacker::GenomePacker (void) {

}

GenomePacker::GenomePacker (char* refname) {
	init(refname);
}

GenomePacker::~GenomePacker (void) {
	if (fin.is_open())
		fin.close();

	if (fout.is_open())
		fout.close();

	if (fout_info.is_open())
		fout_info.close();
}

void GenomePacker::init (char* refname) {
	string packed_fname(refname);
	packed_fname += ".packed.fa";

	string info_fname(refname);
	info_fname += ".index.info";

	fin.open(refname, std::ifstream::in);
	fout.open(packed_fname, std::ofstream::out);
	fout_info.open(info_fname, std::ofstream::out);
}

bool GenomePacker::is_next_chr (void) {
	char ch = fin.peek();
	return ch == '>';
}

bool GenomePacker::is_eof (void) {
	char ch = fin.peek();
	return ch == EOF;
}

bool GenomePacker::get_next_chr (string& chr_id, string& content) {
	chr_id = "";
	content = "";
	
	char ch;
	fin.get(ch);
	if (ch != '>')
		return false;
		
	fin >> chr_id;
	fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	string tmp = "";
	while ((! is_eof()) and (! is_next_chr())) {
		fin >> tmp;
		content += tmp;
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	return true;
}

void GenomePacker::pack_genome (void) {
	int contig_num = 0;
	int cur_contig_size = 0;
	string current = "";
	string chr_id, chr_seq;
	string middle_content = "N";

	while (! is_eof()) {
		if (! get_next_chr(chr_id, chr_seq)) {
			exit(EXIT_FAILURE);
		}

		int chr_seq_len = chr_seq.length();

		if (cur_contig_size == 0 or (chr_seq_len + middle_content.length() + cur_contig_size > CONTIG_SIZE)) {
			++contig_num;
			cur_contig_size = 0;

			fout << ">" << contig_num << endl;
			fout << chr_seq << endl;

			fout_info << contig_num << "\t" << cur_contig_size << "\t" << cur_contig_size + chr_seq_len << "\t" << chr_id << endl;

			cur_contig_size += chr_seq_len;
		}
		else {
			fout << middle_content << chr_seq << endl;
			
			fout_info << contig_num << "\t" << cur_contig_size + middle_content.length() << "\t" << 
						cur_contig_size + middle_content.length() + chr_seq_len << "\t" << chr_id << endl;

			cur_contig_size += middle_content.length() + chr_seq_len;
		}
	}
}

