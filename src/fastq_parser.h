#ifndef __FASTQPARSER_H__
#define __FASTQPARSER_H__

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <zlib.h>

#include "common.h"

#define BLOCKSIZE 1000
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
    int curr_read;
    int filled_size;

    FASTQParser *mate_q;

    char tokens[FQCOMMENTCNT][100];


    void read_buffer();
    uint32_t read_line(char **seq);

    bool has_next(void);

    bool read_block(void);

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
    inline Record *get_next(void);
    inline Record *get_next(int rid);
    inline int get_next_rec_id(void);
    Record *get_next_block(void);
    int get_block_size(void);
};


inline int FASTQParser::get_next_rec_id(void) {
    int rid = -1;

    mutex_lock(&read_lock);

    if (curr_read < filled_size) {
        rid = curr_read;
    } else if (read_block()) {
        rid = curr_read;
        if (mate_q != NULL)
            mate_q->read_block();
    }
    ++curr_read;
    mutex_unlock(&read_lock);

    return rid;
}

inline Record *FASTQParser::get_next(int rid) {
    return current_record + rid;
}

inline Record *FASTQParser::get_next(void) {
    Record *r = NULL;

    mutex_lock(&read_lock);

    if (curr_read < filled_size) {
        r = current_record + curr_read;
        ++curr_read;
    } else if (read_block()) {
        r = current_record + curr_read;
        ++curr_read;
    }

    mutex_unlock(&read_lock);

    return r;
}

inline Record *FASTQParser::get_next_block(void) {
    if (read_block())
        return current_record;
    else
        return NULL;
}

inline int FASTQParser::get_block_size(void) {
    return filled_size;
}

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
