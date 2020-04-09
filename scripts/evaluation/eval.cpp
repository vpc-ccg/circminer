

#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <algorithm>
#include "cigar.h"


using std::ifstream;
using std::string;
using std::array;
using std::vector;

enum class cigar_character_type{
    matched, onquery, ontemplate, hardclip, softclip, notcigar,
};

cigar_character_type what_is_this_cigar(char c){
    switch(c){
        case 'M':
        case '=':
        case 'X':
            return cigar_character_type::matched;
        case 'D':
        case 'N':
            return cigar_character_type::ontemplate;
        case 'I':
        case 'P':
            return cigar_character_type::onquery;
        case 'H':
            return cigar_character_type::hardclip;
        case 'S':
            return cigar_character_type::softclip; 
        default:
            return cigar_character_type::notcigar;
    }
}

struct ival{
    int s;
    int e;
    string c;
    ival(int start, int end, string chr) : s(start), e(end), c(chr){}
    ival() : s(0), e(0), c("-1"){}
    friend std::ostream& operator<<(std::ostream& os, const ival& i);
    bool operator==(const ival &other) const{
        return s == other.s && e == other.e && c == other.c;
    }


};

namespace std{
    template <>
        struct hash<ival>
        {
            std::size_t operator()(const ival& k) const
            {
                using std::size_t;
                using std::hash;
                // Compute individual hash values for first,
                // second and third and combine them using XOR
                // and bit shifting:

                return ((hash<int>()(k.s)
                            ^ (hash<int>()(k.e) << 1)) >> 1)
                            ^ (hash<string>()(k.c) << 1);
            }
        };

}

vector<ival> cigar2intervals( string chr, int alignment_start, cigar cig, int max_na){
    int start = alignment_start;
    int end;
    vector<ival> ivals;
    if( cig.begin() == cig.end()){
        ivals.emplace_back(alignment_start,alignment_start+1,chr);
        return ivals;
    }
    for( auto pair : cig){
        int len = pair.first;
        char  c = pair.second;
        cigar_character_type type = what_is_this_cigar(c);
        if( type == cigar_character_type::matched){
            end = start + len;
            ivals.emplace_back(start,end-1,chr);
            start = end;
        }
        else if( type == cigar_character_type::ontemplate){
            start = start + len;
        }
        else if( type == cigar_character_type::onquery){
            // Nothing for now
        }
    }

    if( max_na > 0){
        vector<ival> merged;
        start = ivals[0].s;
        end = ivals[0].e;
        for(auto iter = std::begin(ivals); std::next(iter) != std::end(ivals); ++iter){
            if( std::next(iter)->s - iter->e < max_na){
                end = std::next(iter)->e;
            }
            else{
                merged.emplace_back(start,end,chr);
                start = std::next(iter)->s;
                end = std::next(iter)->e;
            }
        }
        merged.emplace_back(start,end,chr);
        return merged;
    }

    return ivals;
}

struct sam_read{
    string rid;
    string next_rid;
    vector<array<ival,2>> aligs;

    sam_read(ifstream &samin,string prev_rid){
        string new_rid = prev_rid;
        rid = prev_rid;
        array<ival,2> arr = {ival(),ival()};
        while( new_rid == rid){
            for(int i=0; i<2; i++){
                int flag;
                samin >> flag;
                string chr;
                samin >> chr;

                int pos;
                samin >> pos;              
                samin.ignore(32,'\t'); //Ignore qual
                samin.ignore(32,'\t'); //Ignore qual

                string cigar_str;
                samin >> cigar_str;
                cigar c(cigar_str);
                auto ivals = cigar2intervals(chr, pos,c,25000000);
                ival alig = ivals[0];
    
                arr[i].s = alig.s;
                arr[i].e = alig.e;
    
                samin.ignore(100000,'\n'); //ignore rest
                if(samin.eof()){ break;}
                samin >> new_rid;
    
            }

            aligs.push_back(arr);
            if(samin.eof()){ break;}

        }
        next_rid = new_rid;
    

    }
    sam_read(ifstream &samin){
        string new_rid;
        samin >> new_rid;
        sam_read(samin, new_rid);
    }
};

//PAM
//ENST00000000233-100     7       127588557       127589494       100     1       100     -       1       7       127588406       127588505       100     1       100     +       0       251     0       1       0
struct pam_read{
    array<ival,2> alig;
    vector<array<ival,2>> all_aligs; // For sam conversion

    string rid;
    pam_read(ifstream &pamin){

        pamin >> rid;
        
        string chr;
        pamin >> chr;

        int start;
        int end;
        pamin >> start;
        pamin >> end;
        alig[0].s = start;
        alig[0].e = end;
        alig[0].c = chr;


        pamin.ignore(100000,'\t');
        pamin.ignore(100000,'\t'); 
        pamin.ignore(100000,'\t'); 
        pamin.ignore(100000,'\t');
    
        pamin.ignore(100000,'\t');
        pamin.ignore(100000,'\t');
        pamin >> chr;
        pamin >> start;
        pamin >> end;
        alig[1].s = start;
        alig[1].e = end;
        alig[1].c = chr;
        pamin.ignore(10000000,'\n');
    
        std::sort(std::begin(alig),std::end(alig),[](ival &i1, ival &i2) -> bool{
            return i1.s < i2.s;
                });
    }
    pam_read( const sam_read &s){

        alig[0].s = s.aligs[0][0].s;
        alig[0].e = s.aligs[0][0].e;
        alig[1].s = s.aligs[0][1].s;
        alig[1].e = s.aligs[0][1].e;
        std::sort(std::begin(alig),std::end(alig),[](ival &i1, ival &i2) -> bool{
            return i1.s < i2.s;
                });
        for (auto a : s.aligs){
            auto b = a;
            
            std::sort(std::begin(b),std::end(b),[](ival &i1, ival &i2) -> bool{
            return i1.s < i2.s;
                });
            all_aligs.push_back(b);
            
        }

    }
};

int overlap(const pam_read &s1, const pam_read &s2, int permit){
    int sum = 0;
    
    if( s1.all_aligs.size() > 1){
        int max_sum = - 1;
        for( auto a : s1.all_aligs){
            sum = 0;
            for(int i =0; i < 2; i++){
                if(std::abs(a[i].s - s2.alig[i].s) < permit){
                    sum+=1;
                }
                if(std::abs(a[i].e - s2.alig[i].e) < permit){
                    sum+=1;
                }
            }
            if(sum > max_sum){
                max_sum = sum;
            }
        }
        return max_sum;
    }
    else if( s2.all_aligs.size() > 1){
        int max_sum = - 1;
        for( auto a : s2.all_aligs){
            sum = 0;
            for(int i =0; i < 2; i++){
                if(std::abs(s1.alig[i].s - a[i].s) < permit){
                    sum+=1;
                }
                if(std::abs(s1.alig[i].e - a[i].e) < permit){
                    sum+=1;
                }
            }
            if(sum > max_sum){
                max_sum = sum;
            }
        }
        return max_sum;
        
    }
    else{
        for(int i =0; i < 2; i++){
            if(std::abs(s1.alig[i].s - s2.alig[i].s) < permit){
                sum+=1;
            }
            if(std::abs(s1.alig[i].e - s2.alig[i].e) < permit){
                sum+=1;
            }
        }
        return sum;
    }
}

int main(int argc, char **argv){
    if ( argc != 4){
        std::cerr << "./eval art.sam star.sam miner.pam > [output]" << std::endl;
        return -1;
    }
#define PERMIT 10
    ifstream samin(argv[2]);
    ifstream pamin(argv[3]);
    ifstream golin(argv[1]);

    string current_rid;
    string first_tab_star;
    string first_tab_gold;
    golin >>  first_tab_gold;
    while( first_tab_gold[0] == '@'){
        if(golin.eof()){
            std::cerr << "No valid read! for " << argv[1] << std::endl;
            return -1;
        }
        golin.ignore(111111132,'\n'); //Ignore qual
        golin >> first_tab_gold;
    }

    samin >>  first_tab_star;
    while( first_tab_star[0] == '@'){
        if(samin.eof()){
            std::cerr << "No valid read! for " << argv[2] << std::endl;
            return -1;
        }
        samin.ignore(111111132,'\n'); //Ignore qual
        samin >> first_tab_star;
    }


    std::ios::sync_with_stdio(false);
    sam_read r(golin, first_tab_star);
    sam_read g(samin, first_tab_gold);
    pam_read p(pamin);
    int skip_counts[3] ={0};
    while( !samin.eof() && !pamin.eof() && !golin.eof()){  
        
        while(r.rid != p.rid || g.rid!=p.rid){
            auto rids = { p.rid, r.rid, g.rid};
            auto min = std::min_element(std::begin(rids),std::end(rids));
            if( *min == p.rid){
                p = pam_read(pamin);
                skip_counts[0]+=1;
            }
            if( *min == r.rid){
                r = sam_read(samin,r.next_rid);
                skip_counts[1]+=1;
            }
            if( *min == g.rid){
             g = sam_read(golin,g.next_rid);
                skip_counts[2]+=1;
            }
        }

        std::cout << r.rid << "\t" << overlap(p,g,PERMIT) << "\t" << overlap(p,r,PERMIT) << "\t" << overlap(r,g,PERMIT) << "\n";

        p = pam_read(pamin);

        r = sam_read(samin,r.next_rid);

        g = sam_read(golin,g.next_rid);


    }

    if( r.rid == p.rid && p.rid==g.rid){
        std::cout << r.rid << "\t" << overlap(p,g,PERMIT) << "\t" << overlap(p,r,PERMIT) << "\t" << overlap(r,g,PERMIT) << "\n";
    }
    std::cerr << skip_counts[0] << "\t" << skip_counts[1] << "\t" << skip_counts[2] << "\n";
    return 0;
}
