#include <fstream>
#include <iostream>
#include <limits>

#include "genome.h"
#include "common.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

using namespace std;

#define FASTA_LINE_MAX_LEN 50

GenomePacker::GenomePacker(void) {

}

GenomePacker::GenomePacker(char *refname) {
    init(refname);
}

GenomePacker::~GenomePacker(void) {
}

void GenomePacker::init(char *refname) {
    string ref_fname_str(refname);

    ref_fname = ref_fname_str;
    packed_fname = ref_fname_str + ".packed.fa";

    index_fname = packed_fname + ".index";
    index_info_fname = packed_fname + ".index.info";
}

// returns 0 in case of error
int GenomePacker::build_index(bool is_compact_index) {
    fprintf(stdout, "Started packing reference genome...\n");
    pack_genome();
    fprintf(stdout, "Packing done!\n");

    char packed_fname_cstr[FILE_NAME_MAX_LEN];
    strcpy(packed_fname_cstr, packed_fname.c_str());

    char index_fname_cstr[FILE_NAME_MAX_LEN];
    strcpy(index_fname_cstr, index_fname.c_str());

    fprintf(stdout, "Builing index\n");
    if (is_compact_index) {
        if (!generateHashTable(packed_fname_cstr, index_fname_cstr))
            return 0;
    } else {
        if (!generateHashTableOnDisk(packed_fname_cstr, index_fname_cstr))
            return 0;
    }

    return 1;
}

bool GenomePacker::is_next_chr(void) {
    char ch = fin.peek();
    return ch == '>';
}

bool GenomePacker::is_eof(void) {
    char ch = fin.peek();
    return ch == EOF;
}

bool GenomePacker::get_next_chr(string &chr_id, string &content) {
    chr_id = "";
    content = "";

    char ch;
    fin.get(ch);
    if (ch != '>')
        return false;

    fin >> chr_id;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    string tmp = "";
    while ((!is_eof()) and (!is_next_chr())) {
        fin >> tmp;
        content += tmp;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    return true;
}

void GenomePacker::pack_genome(void) {
    // open files
    fin.open(ref_fname, std::ifstream::in);
    fout.open(packed_fname, std::ofstream::out);
    fout_info.open(index_info_fname, std::ofstream::out);

    // pack geonme

    int contig_num = 0;
    int cur_contig_size = 0;
    string current = "";
    string chr_id, chr_seq;
    string middle_content = "N";

    while (!is_eof()) {
        if (!get_next_chr(chr_id, chr_seq)) {
            exit(EXIT_FAILURE);
        }

        int chr_seq_len = chr_seq.length();

        if (cur_contig_size == 0 or (chr_seq_len + middle_content.length() + cur_contig_size > CONTIG_SIZE)) {
            ++contig_num;
            cur_contig_size = 0;

            fout << ">" << contig_num << endl;
            fout << chr_seq << endl;

            fout_info << contig_num << "\t" << cur_contig_size << "\t" << cur_contig_size + chr_seq_len << "\t"
                      << chr_id << endl;

            cur_contig_size += chr_seq_len;
        } else {
            fout << middle_content << chr_seq << endl;

            fout_info << contig_num << "\t" << cur_contig_size + middle_content.length() << "\t" <<
                      cur_contig_size + middle_content.length() + chr_seq_len << "\t" << chr_id << endl;

            cur_contig_size += middle_content.length() + chr_seq_len;
        }
    }

    // close files
    fin.close();
    fout.close();
    fout_info.close();
}

void GenomePacker::load_index_info(vector <ContigLen> &contig_len) {
    fin.open(index_info_fname, std::ifstream::in);

    contig_len.clear();

    string orig_chr;
    uint32_t contig;
    uint32_t start_pos, end_pos;
    while (fin >> contig >> start_pos >> end_pos >> orig_chr) {
        ContigLen tmp = {
                .name = orig_chr,
                .contig_id = contig,
                .start_pos = start_pos,
                .end_pos = end_pos,
                .len = end_pos - start_pos
        };
        contig_len.push_back(tmp);
    }

    fin.close();
}

string GenomePacker::get_index_fname(void) {
    return index_fname;
}

