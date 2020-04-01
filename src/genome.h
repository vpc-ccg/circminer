#ifndef __GENOME_PACKER_H__
#define __GENOME_PACKER_H__

#include <fstream>

using namespace std;

class GenomePacker {
private:
	ifstream fin;
	ofstream fout;
	ofstream fout_info;

	bool is_next_chr (void);
	bool is_eof (void);

	bool get_next_chr (string& chr_id, string& content);

public:
	GenomePacker (void);
	GenomePacker (char* refname);
	~GenomePacker (void);

	void init (char* refname);

	void pack_genome (void);

};

extern GenomePacker genome_packer;

#endif //__GENOME_PACKER_H__
