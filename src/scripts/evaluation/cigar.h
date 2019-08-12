#ifndef CIGAR_CIGAR
#define CIGAR_CIGAR
#include <string>
#include <iostream>
#include <vector>

#define CIGAR_CHARACTERS "DHIMNPSX="

class cigar{
    std::vector<std::pair<int,char>> cigarray;
    public:
    cigar(const std::string &cgr); 
    decltype(cigarray.begin()) begin();
    decltype(cigarray.begin()) end();
    auto operator [](size_t index) const;

};

#endif
