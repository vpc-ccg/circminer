#ifndef __FASTQPARSER_H__
#define __FASTQPARSER_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <zlib.h>

#include "common.h"

#define MAX_THREADS_COUNT 1024
#define FQCOMMENTCNT 23
#define BUFFSIZE 10000000

class FASTQParser {
private:
    FILE *input;
    gzFile gzinput;
    char comp[ASCISIZE];
    char *zbuffer;
    int32_t buff_pos;
    int32_t buff_size;

    Record *current_record;
    size_t max_line_size;

    FASTQParser *mate_q;

    char tokens[FQCOMMENTCNT][100];


    void read_buffer();
    uint32_t read_line(char **seq);

    bool has_next(void);

    void set_comp(void);
    void set_reverse_comp(int r_ind);

    int extract_map_info(char *str, int r_ind);
    void fill_map_info(int cnt, int r_ind);

public:
    FASTQParser(void);
    FASTQParser(char *filename);
    ~FASTQParser(void);

    void init(void);
    void reset(char *filename);
    void finalize(void);

    void set_mate(FASTQParser *mq);

    Record *get_next_read(int thread_id);
};

inline bool FASTQParser::has_next(void) {
    if (buff_pos >= buff_size) {
        read_buffer();
        if (buff_size == 0)
            return false;
    }

    char c = zbuffer[buff_pos++];
    assert(c == '@');    // ensure FASTQ format
    return true;
}


extern FASTQParser fq_parser1;
extern FASTQParser fq_parser2;

#endif    //__FASTQPARSER_H__
