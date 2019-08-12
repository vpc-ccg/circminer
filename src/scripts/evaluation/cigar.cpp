#include <string>
#include <iostream>
#include <vector>

#include "cigar.h"

using std::string;
using std::vector;

#define CIGAR_CHARACTERS "DHIMNPSX="

cigar::cigar(const string &cgr){
    size_t prev = -1;
    size_t index = cgr.find_first_of(CIGAR_CHARACTERS);
    if(cgr != "*"){
        while(index != string::npos){
            int len = std::stoi(cgr.substr(1+prev,index));
            cigarray.push_back(std::make_pair(len,cgr[index]));
            prev = index;
            index = cgr.find_first_of(CIGAR_CHARACTERS,index + 1);
        }

    }
}
decltype(cigar::cigarray.begin()) cigar::begin(){
    return cigarray.begin();
}

decltype(cigar::cigarray.end()) cigar::end(){
    return cigarray.end();
}

auto cigar::operator [](size_t index) const{
    return cigarray[index];
}

int cigar_test_main(int argc, char **argv){
    cigar cgr("10M1I40M1D1M2I6M1D13M2D7M1D18M1D19M1D64M1D49M1D41M1D5M2I13M1I24M1I12M1D21M1I10M2I2M1I9M1I45M"); 

    for( auto a : cgr){
        std::cout << a.first << a.second << ", ";
    }
    std::cout << "\n";

    std::cout << cgr[0].first << "\n";
    return 0;
}
