#include <cstdio>
#include <cstring>
#include <zlib.h>

#include "fastq_parser.h"

FASTQParser::FASTQParser(void) {
}

FASTQParser::FASTQParser(char *filename) {
    reset(filename);
}

FASTQParser::~FASTQParser(void) {
    free(zbuffer);
    for (int i = 0; i < threadCount; ++i) {
        free(current_record[i].rname);
        free(current_record[i].seq);
        free(current_record[i].rcseq);
        free(current_record[i].comment);
        free(current_record[i].qual);
        free(current_record[i].rqual);
    }
    // free(current_record);
    delete[] current_record;

    finalize();
}

void FASTQParser::init(void) {
    input = NULL;
    gzinput = Z_NULL;
    mate_q = NULL;

    max_line_size = MAXLINESIZE;
    set_comp();

    zbuffer = (char *) malloc(BUFFSIZE);

    current_record = new Record[threadCount];
    for (int i = 0; i < threadCount; ++i) {
        current_record[i].rname = (char *) malloc(max_line_size);
        current_record[i].seq = (char *) malloc(max_line_size);
        current_record[i].rcseq = (char *) malloc(max_line_size);
        current_record[i].comment = (char *) malloc(max_line_size);
        current_record[i].qual = (char *) malloc(max_line_size);
        current_record[i].rqual = (char *) malloc(max_line_size);
    }
}

void FASTQParser::reset(char *filename) {
    finalize();

    char *fname = (char *) malloc(FILE_NAME_MAX_LEN);
    char *rmode = (char *) malloc(FILE_NAME_MAX_LEN);

    sprintf(fname, "%s", filename);
    sprintf(rmode, "%c", 'r');

    gzinput = open_gzfile(fname, rmode);

    //input = open_file(filename, "r");

    free(fname);
    free(rmode);

    buff_pos = 0;
    buff_size = 0;
}

void FASTQParser::finalize(void) {
    if (gzinput != NULL) {
        close_gzfile(gzinput);
        gzinput = Z_NULL;
    }
    if (input != NULL) {
        close_file(input);
        input = NULL;
    }
}

void FASTQParser::set_mate(FASTQParser *mq) {
    mate_q = mq;
}

void FASTQParser::read_buffer() {
    buff_size = gzread(gzinput, zbuffer, BUFFSIZE);
    buff_pos = 0;

    if (buff_size == 0 and gzeof(gzinput) == 0) {
        buff_size = -1;
    }
    if (buff_size < 0) {
        int err;
        fprintf(stderr, "gzread error: %s\n", gzerror(gzinput, &err));
        exit(1);
    }
}

uint32_t FASTQParser::read_line(char **seq) {
    char cur;

    uint32_t i = 0;
    while (true) {
        if (buff_pos >= buff_size) {
            read_buffer();
            if (buff_size == 0)
                return 0;
        }

        cur = zbuffer[buff_pos++];
        if (cur == '\n') {
            (*seq)[i] = '\0';
            return i;
        }

        (*seq)[i++] = cur;
    }
}

Record *FASTQParser::get_next_read(int thread_id) {
    if (has_next()) {
        read_line(&current_record[thread_id].rname);
        extract_map_info(current_record[thread_id].rname, thread_id);

        current_record[thread_id].seq_len = read_line(&current_record[thread_id].seq);

        read_line(&current_record[thread_id].comment);
        assert(current_record[thread_id].comment[0] == '+');
        current_record[thread_id].comment[1] = '\0';

        read_line(&current_record[thread_id].qual);

        set_reverse_comp(thread_id);
        return current_record + thread_id;
    } else {
        return NULL;
    }
}

void FASTQParser::set_comp(void) {
    comp[uint8_t('A')] = 'T';
    comp[uint8_t('C')] = 'G';
    comp[uint8_t('G')] = 'C';
    comp[uint8_t('T')] = 'A';
    comp[uint8_t('N')] = 'N';

    comp[uint8_t('a')] = 'T';
    comp[uint8_t('c')] = 'G';
    comp[uint8_t('g')] = 'C';
    comp[uint8_t('t')] = 'A';
    comp[uint8_t('n')] = 'N';
}

void FASTQParser::set_reverse_comp(int r_ind) {
    uint32_t len = current_record[r_ind].seq_len;
    uint32_t i = len;
    do {
        --i;
        current_record[r_ind].rcseq[len - i - 1] = comp[uint8_t(current_record[r_ind].seq[i])];
    } while (i != 0);
    current_record[r_ind].rcseq[len] = '\0';

    // reverse qual
    uint32_t qual_len = strlen(current_record[r_ind].qual);
    if (qual_len != len) {
        fprintf(stderr, "ERROR: read: %s, length of sequence (%d) does not match with quality (%d)!\nAborting\n",
                current_record[r_ind].rname, len, qual_len);
        exit(1);
    }

    i = len;
    do {
        --i;
        current_record[r_ind].rqual[len - i - 1] = current_record[r_ind].qual[i];
    } while (i != 0);
    current_record[r_ind].rqual[len] = '\0';
}

int FASTQParser::extract_map_info(char *str, int r_ind) {
    // Returns first token
    char *token = strtok(str, " ");

    // Keep printing tokens while one of the
    // delimiters present in str[].
    int i = 0;
    while (token != NULL) {
        strcpy(tokens[i], token);
        token = strtok(NULL, " ");
        ++i;
    }

    int rname_len = strlen(tokens[0]) + 1;
    current_record[r_ind].rname[rname_len] = '\0';

    if (current_record[r_ind].rname[rname_len - 3] == '/')
        current_record[r_ind].rname[rname_len - 3] = '\0';

    fill_map_info(i, r_ind);
    return rname_len;
}

void FASTQParser::fill_map_info(int cnt, int r_ind) {
    //assert(cnt == 1 or cnt == FQCOMMENTCNT);

    if (cnt != FQCOMMENTCNT) {
        current_record[r_ind].mr->type = NOPROC_NOMATCH;
        current_record[r_ind].mr->tlen = INF;
        current_record[r_ind].mr->junc_num = 0;
        current_record[r_ind].mr->gm_compatible = false;
    } else {
        char *stop_string;
        int base = 10;
        current_record[r_ind].mr->type = atoi(tokens[2]);

        if (current_record[r_ind].mr->type == CONCRD or current_record[r_ind].mr->type == DISCRD or
            current_record[r_ind].mr->type == CHIORF or
            current_record[r_ind].mr->type == CHIBSJ or current_record[r_ind].mr->type == CHI2BSJ or
            current_record[r_ind].mr->type == CONGNM or current_record[r_ind].mr->type == CONGEN) {
            current_record[r_ind].mr->genome_spos = strtoull(tokens[1], &stop_string, base);
            current_record[r_ind].mr->chr_r1 = tokens[3];
            current_record[r_ind].mr->spos_r1 = strtoul(tokens[4], &stop_string, base);
            current_record[r_ind].mr->epos_r1 = strtoul(tokens[5], &stop_string, base);
            current_record[r_ind].mr->mlen_r1 = atoi(tokens[6]);
            current_record[r_ind].mr->qspos_r1 = strtoul(tokens[7], &stop_string, base);
            current_record[r_ind].mr->qepos_r1 = strtoul(tokens[8], &stop_string, base);
            current_record[r_ind].mr->r1_forward = (tokens[9][0] == '+');
            current_record[r_ind].mr->ed_r1 = atoi(tokens[10]);

            current_record[r_ind].mr->chr_r2 = tokens[11];
            current_record[r_ind].mr->spos_r2 = strtoul(tokens[12], &stop_string, base);
            current_record[r_ind].mr->epos_r2 = strtoul(tokens[13], &stop_string, base);
            current_record[r_ind].mr->mlen_r2 = atoi(tokens[14]);
            current_record[r_ind].mr->qspos_r2 = strtoul(tokens[15], &stop_string, base);
            current_record[r_ind].mr->qepos_r2 = strtoul(tokens[16], &stop_string, base);
            current_record[r_ind].mr->r2_forward = (tokens[17][0] == '+');
            current_record[r_ind].mr->ed_r2 = atoi(tokens[18]);

            current_record[r_ind].mr->tlen = atoi(tokens[19]);
            current_record[r_ind].mr->junc_num = strtoul(tokens[20], &stop_string, base);
            current_record[r_ind].mr->gm_compatible = (tokens[21][0] == '1');
            current_record[r_ind].mr->contig_num = atoi(tokens[22]);
        } else {
            current_record[r_ind].mr->genome_spos = 0;
            current_record[r_ind].mr->chr_r1 = "-";
            current_record[r_ind].mr->spos_r1 = 0;
            current_record[r_ind].mr->epos_r1 = 0;
            current_record[r_ind].mr->mlen_r1 = 0;
            current_record[r_ind].mr->qspos_r1 = 0;
            current_record[r_ind].mr->qepos_r1 = 0;
            current_record[r_ind].mr->r1_forward = true;
            current_record[r_ind].mr->ed_r1 = maxEd + 1;

            current_record[r_ind].mr->chr_r2 = "-";
            current_record[r_ind].mr->spos_r2 = 0;
            current_record[r_ind].mr->epos_r2 = 0;
            current_record[r_ind].mr->mlen_r2 = 0;
            current_record[r_ind].mr->qspos_r2 = 0;
            current_record[r_ind].mr->qepos_r2 = 0;
            current_record[r_ind].mr->r2_forward = true;
            current_record[r_ind].mr->ed_r2 = maxEd + 1;

            current_record[r_ind].mr->tlen = INF;
            current_record[r_ind].mr->junc_num = 0;
            current_record[r_ind].mr->gm_compatible = false;
            current_record[r_ind].mr->contig_num = 0;
        }
    }
}
