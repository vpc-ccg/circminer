#ifndef __FASTA__
#define __FASTA__
#include <map>
#include <string>


void LoadFasta(const char *fasta_file, std::map<std::string, std::string> &chr_map);
//void LoadFastaShort(char *fasta_file, std::map<std::string, std::string> &chr_map, size_t buffer_size);
//void LoadFasta_Single(char *fasta_file, std::string &chr_name, std::string &chr_seq);

#endif
