#ifndef __GENOME_PACKER_H__
#define __GENOME_PACKER_H__

#include <fstream>
#include <vector>

#include "common.h"

using namespace std;

class GenomePacker {
private:
    string ref_fname;
    string packed_fname;
    string index_fname;
    string index_info_fname;

    ifstream fin;
    ofstream fout;
    ofstream fout_info;

    bool is_next_chr(void);
    bool is_eof(void);

    bool get_next_chr(string &chr_id, string &content);

    void pack_genome(void);

public:
    GenomePacker(void);
    GenomePacker(char *refname);
    ~GenomePacker(void);

    void init(char *refname);

    int build_index(bool is_compact_index);

    void load_index_info(vector <ContigLen> &contig_len);

    string get_index_fname(void);
};

extern GenomePacker genome_packer;

#endif //__GENOME_PACKER_H__
