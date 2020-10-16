#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <pthread.h>

#include "common.h"
#include "commandline_parser.h"
#include "fastq_parser.h"
#include "filter.h"
#include "gene_annotation.h"
#include "align.h"
#include "process_circ.h"
#include "hash_table.h"
#include "utils.h"
#include "genome.h"

extern "C" {
#include "mrsfast/Common.h"
#include "mrsfast/HashTable.h"
}

using namespace std;

char versionNumberMajor[10] = "0";
char versionNumberMinor[10] = "4";
char versionNumberPatch[10] = "5";

pthread_mutex_t write_lock;
pthread_mutex_t pmap_lock;
pthread_mutex_t buffer_lock;

GTFParser gtf_parser;
FilterRead filter_read;
ScoreMatrix score_mat;
ScoreMatrix edit_mat;
FASTQParser fq_parser1;
FASTQParser fq_parser2;
GenomeSeeder genome_seeder;
GenomePacker genome_packer;

int mapping(int &last_round_num);
void circ_detect(int last_round_num);
void *map_reads(void *args);

int main(int argc, char **argv) {
    Logger::instance().error.toggle_time().set_prefix("[ERROR] ");
    Logger::instance().error.set_buffer_size(10);
    Logger::instance().info.toggle_time().set_prefix("[INFO] ");
    Logger::instance().info.set_buffer_size(10);
    Logger::instance().debug.toggle_time().set_prefix("[DEBUG] ");

    int exit_c = parse_command(argc, argv);
    if (exit_c == 1)
        return 0;

    genome_packer.init(referenceFilename);

    fq_parser1.init();
    if (pairedEnd)
        fq_parser2.init();

    /****************************************************
     * INDEXING
     ***************************************************/
    if (indexMode) {
        if (!genome_packer.build_index(compactIndex))
            return 1;
    }

        /****************************************************
         * SEARCHING
         ***************************************************/
    else {
        score_mat.init(1, -3, -3, 8);
        edit_mat.init(0, 1, 1, 10000);
        if (pairedEnd)
            fq_parser1.set_mate(&fq_parser2);

        // 0 <= stage <= 2
        int last_round_num = 1;
        if (stage != 1) {
            int map_ret = mapping(last_round_num);
            if (map_ret == 1)
                return 1;
        }

        if (stage != 0) {
            // dirty turnaround --should be removed later
            if (stage == 1)
                last_round_num = 3;
            circ_detect(last_round_num);
        }
    }

    return 0;
}

int mapping(int &last_round_num) {
    double cputime_start = get_cpu_time();
    double realtime_start = get_real_time();
    double cputime_curr;
    double realtime_curr;

    char *fq_file1 = fastqFilename[0];
    char *fq_file2;
    bool is_pe = pairedEnd;
    if (is_pe) {
        fq_file2 = fastqFilename[1];
    }

    if (pthread_mutex_init(&buffer_lock, NULL) != 0) {
        Logger::instance().error("Mutex init has failed\n");
        return 1;
    }
    if (pthread_mutex_init(&pmap_lock, NULL) != 0) {
        Logger::instance().error("Mutex init has failed\n");
        return 1;
    }
    if (pthread_mutex_init(&write_lock, NULL) != 0) {
        Logger::instance().error("Mutex init has failed\n");
        return 1;
    }

    /**********************/
    /**Loading Hash Table**/
    /**********************/

    int flag;
    double tmpTime;

    vector <ContigLen> orig_contig_len;
    genome_packer.load_index_info(orig_contig_len);

    char index_file[FILE_NAME_MAX_LEN];
    strcpy(index_file, genome_packer.get_index_fname().c_str());

    if (!checkHashTable(index_file))
        return 1;

    if (!initLoadingHashTableMeta(index_file))
        return 1;

    int kmer_index = WINDOW_SIZE + checkSumLength;

    Logger::instance().info("Kmer size obtained from index: %d\n", kmer_index);
    if (kmer != kmer_index) {
        Logger::instance().info("Given kmer size (%d) is not equal to the index kmer size (%d).\n",
                                kmer, kmer_index);
        Logger::instance().info("Setting kmer size to %d.\n", kmer_index);
    }

    kmer = kmer_index;
        
    genome_seeder.init(LOADCONTIGSTRINMEM);

    /*********************/
    /**Memory Allocation**/
    /*********************/

    // considering both overlapping and non-overlapping kmers
    int max_seg_cnt = 2 * (ceil(1.0 * maxReadLength / kmer)) - 1;

    vector < GIMatchedKmer * > fl(threadCount);
    vector < GIMatchedKmer * > bl(threadCount);

    vector < chain_list * > fbc_r1(threadCount);
    vector < chain_list * > bbc_r1(threadCount);
    vector < chain_list * > fbc_r2(threadCount);
    vector < chain_list * > bbc_r2(threadCount);

    for (int th = 0; th < threadCount; ++th) {
        fl[th] = (GIMatchedKmer *) malloc(max_seg_cnt * sizeof(GIMatchedKmer));
        bl[th] = (GIMatchedKmer *) malloc(max_seg_cnt * sizeof(GIMatchedKmer));

        fbc_r1[th] = (chain_list *) malloc(1 * sizeof(chain_list));
        bbc_r1[th] = (chain_list *) malloc(1 * sizeof(chain_list));
        fbc_r2[th] = (chain_list *) malloc(1 * sizeof(chain_list));
        bbc_r2[th] = (chain_list *) malloc(1 * sizeof(chain_list));

        fbc_r1[th]->chains = (chain_t *) malloc(maxChainLen * sizeof(chain_t));
        bbc_r1[th]->chains = (chain_t *) malloc(maxChainLen * sizeof(chain_t));
        fbc_r2[th]->chains = (chain_t *) malloc(maxChainLen * sizeof(chain_t));
        bbc_r2[th]->chains = (chain_t *) malloc(maxChainLen * sizeof(chain_t));

        for (int i = 0; i < maxChainLen; i++) {
            fbc_r1[th]->chains[i].frags = (fragment_t *) malloc(max_seg_cnt * sizeof(fragment_t));
            bbc_r1[th]->chains[i].frags = (fragment_t *) malloc(max_seg_cnt * sizeof(fragment_t));
            fbc_r2[th]->chains[i].frags = (fragment_t *) malloc(max_seg_cnt * sizeof(fragment_t));
            bbc_r2[th]->chains[i].frags = (fragment_t *) malloc(max_seg_cnt * sizeof(fragment_t));
        }
    }

    pthread_t *cm_threads = (pthread_t *) malloc(threadCount * sizeof(pthread_t));

    FilterArgs *filter_args[threadCount];
    for (int th = 0; th < threadCount; ++th)
        filter_args[th] = new FilterArgs(kmer);

    /*******************/
    /**GTF Parser Init**/
    /*******************/

    Logger::instance().info("Loading GTF file...\n");

    gtf_parser.init(gtfFilename, orig_contig_len);
    if (!gtf_parser.load_gtf()) {
        fprintf(stdout, "Error in reading GTF file.\n");
        exit(1);
    }

    cputime_curr = get_cpu_time();
    realtime_curr = get_real_time();

    Logger::instance().info("Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n",
                            cputime_curr - cputime_start, realtime_curr - realtime_start);

    cputime_start = cputime_curr;
    realtime_start = realtime_curr;

    /*****************/
    /**Mapping Reads**/
    /*****************/

    bool is_first = true;
    bool is_last = false;

    Logger::instance().info("Genome index type: %s\n", loadFullHashTable ? "Full" : "Compact");
    Logger::instance().info("Starting read extraction\n");
    do {
        Logger::instance().info("+ Loading genome index...\n");

        flag = loadHashTable(&tmpTime);            // Reading a fragment

        cputime_curr = get_cpu_time();
        realtime_curr = get_real_time();

        Logger::instance().info("+ Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n",
                                cputime_curr - cputime_start, realtime_curr - realtime_start);
        Logger::instance().debug("+ kmer size: %d + %d = %d\n", WINDOW_SIZE, checkSumLength,
                                 WINDOW_SIZE + checkSumLength);

        // loading contig sequence from reference index
        Logger::instance().info("+ Loading genome sequence...\n");
        bool success_gload = genome_seeder.pac2char_whole_contig();

        double cputime_loaded_gnome_str = get_cpu_time();
        double realtime_loaded_gnome_str = get_real_time();

        if (!success_gload) {
            fprintf(stderr, "Error: unable to load reference genome to memory\n");
            exit(1);
        }
        Logger::instance().info("+ Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n",
                                cputime_loaded_gnome_str - cputime_curr, realtime_loaded_gnome_str - realtime_curr);

        cputime_curr = cputime_loaded_gnome_str;
        realtime_curr = realtime_loaded_gnome_str;
        // done loading contig sequence

        cputime_start = cputime_curr;
        realtime_start = realtime_curr;

        double fq_cputime_start = cputime_curr;
        double fq_realtime_start = realtime_curr;

        contigName = getRefGenomeName();
        contigNum = atoi(contigName) - 1;
        last_round_num = contigNum + 1;

        Logger::instance().info("+ Starting pseudo-alignment (Round %s)\n", contigName);

        fq_parser1.reset(fq_file1);
        if (is_pe) {
            fq_parser2.reset(fq_file2);
        }

        is_last = !flag;

        if (!is_first) {
            filter_read.finalize();
        }
        filter_read.init(outputFilename, is_pe, contigNum + 1, is_first, is_last, fq_file1, fq_file2, orig_contig_len);
        is_first = false;

        for (int th = 0; th < threadCount; ++th)
            filter_args[th]->set(fl[th], bl[th], fbc_r1[th], bbc_r1[th], fbc_r2[th], bbc_r2[th]);

        lookup_cnt = 0;

        for (int th = 0; th < threadCount; ++th) {
            filter_args[th]->id = th;
            pthread_create(cm_threads + th, NULL, map_reads, filter_args[th]);
        }

        for (int th = 0; th < threadCount; ++th) {
            pthread_join(cm_threads[th], NULL);
        }

        cputime_curr = get_cpu_time();
        realtime_curr = get_real_time();

        Logger::instance().info("+ Completed round %s! (CPU time: %.2lfs; Real time: %.2lfs)\n", contigName,
                                cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

        cputime_start = cputime_curr;
        realtime_start = realtime_curr;

    } while (flag);

    /*************************/
    /**Free Allocated Memory**/
    /*************************/

    filter_read.finalize();
    genome_seeder.finalize();

    finalizeLoadingHashTable();

    for (int th = 0; th < threadCount; ++th) {
        free(fl[th]);
        free(bl[th]);

        for (int i = 0; i < maxChainLen; i++) {
            free(fbc_r1[th]->chains[i].frags);
            free(bbc_r1[th]->chains[i].frags);
            free(fbc_r2[th]->chains[i].frags);
            free(bbc_r2[th]->chains[i].frags);
        }
        free(fbc_r1[th]->chains);
        free(bbc_r1[th]->chains);
        free(fbc_r2[th]->chains);
        free(bbc_r2[th]->chains);

        free(fbc_r1[th]);
        free(bbc_r1[th]);
        free(fbc_r2[th]);
        free(bbc_r2[th]);

        delete filter_args[th];
    }

    free(cm_threads);

    return 0;
}

void circ_detect(int last_round_num) {
    int ws = 8;
    Logger::instance().info("Starting circRNA detection\n");
    ProcessCirc process_circ(last_round_num, ws);
    process_circ.do_process();
}

void *map_reads(void *args) {
    FilterArgs *fa = (struct FilterArgs *) args;
    Logger::instance().debug("--- thread #%d\n", fa->id);

    filter_print_func print_func_ptr;
    if (reportMapping == SAMFORMAT)
        print_func_ptr = &FilterRead::print_sam;
    else if (reportMapping == PAMFORMAT)
        print_func_ptr = &FilterRead::print_pam;
    else
        print_func_ptr = &FilterRead::print_nothing;


    int rid;
    Record *current_record1;
    Record *current_record2;
    int state;
    int is_last = filter_read.get_last_round();
    while (true) { // go line by line on fastq file
        mutex_lock(&buffer_lock);

        current_record1 = fq_parser1.get_next_read(fa->id);
        if (pairedEnd) {
            current_record2 = fq_parser2.get_next_read(fa->id);
        }
        mutex_unlock(&buffer_lock);

        if (current_record1 == NULL)
            break;
        if (pairedEnd) {
            state = filter_read.process_read(fa->id, current_record1, current_record2, fa->kmer_size, fa->fl, fa->bl,
                                             *(fa->fbc_r1), *(fa->bbc_r1), *(fa->fbc_r2), *(fa->bbc_r2));
            bool skip = (scanLevel == 0 and state == CONCRD) or
                        (scanLevel == 1 and state == CONCRD and current_record1->mr->gm_compatible and
                        (current_record1->mr->ed_r1 + current_record1->mr->ed_r2 == 0) and
                        (current_record1->mr->mlen_r1 + current_record1->mr->mlen_r2 ==
                         current_record1->seq_len + current_record2->seq_len));

            if (skip or is_last) {
                (filter_read.*(print_func_ptr))(current_record1, current_record2);
            }
            if ((!is_last and !skip) or
                (is_last and (current_record1->mr->type == CHIBSJ or current_record1->mr->type == CHI2BSJ)))
                filter_read.write_read_category(current_record1, current_record2, *(current_record1->mr));
        } else {
            state = filter_read.process_read(fa->id, current_record1, fa->kmer_size, fa->fl, fa->bl,
                                             *(fa->fbc_r1), *(fa->bbc_r1));
            filter_read.write_read_category(current_record1, state);
        }
    }

    return NULL;
}
