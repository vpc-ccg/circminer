#include <vector>
#include <cstring>
#include <algorithm>
#include <utility>
#include <cstdlib>

#include "process_circ.h"
#include "gene_annotation.h"
#include "hash_table.h"
#include "fastq_parser.h"
#include "match_read.h"
#include "chain.h"
#include "extend.h"
#include "utils.h"
#include "genome.h"

#define BINSIZE 5000
#define MAXHTLISTSIZE 0
#define TOPCHAIN 10

typedef pair<uint32_t, int> pu32i;

void set_mm(chain_t &ch, uint32_t qspos, int rlen, int dir, MatchedMate &mm);

ProcessCirc::ProcessCirc(int last_round_num, int ws) : extension(0, EDIT_ALIGNMENT) {
    Logger::instance().info.set_prefix("[INFO] + ");

    sprintf(fq_file1, "%s_%d_remain_R1.fastq", outputFilename, last_round_num);
    sprintf(fq_file2, "%s_%d_remain_R2.fastq", outputFilename, last_round_num);

    double cpu_time = get_cpu_time();
    double real_time = get_real_time();

    if (internalSort) {
        Logger::instance().info("Sorting remaining read mappings internally... \n");
        sort_fq_internal(fq_file1);
        sort_fq_internal(fq_file2);
    } else {
        Logger::instance().info("Sorting remaining read mappings using GNU sort... \n");
        int ret1 = sort_fq(fq_file1);
        int ret2 = sort_fq(fq_file2);
//		if (ret1 and ret2)
//		    fprintf(stdout, "OK\n");
    }

    double final_cpu_time = get_cpu_time();
    double final_real_time = get_real_time();

    Logger::instance().info("Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n",
                            final_cpu_time - cpu_time, final_real_time - real_time);
    fflush(stdout);

    sprintf(fq_file1, "%s_%d_remain_R1.fastq.srt", outputFilename, last_round_num);
    sprintf(fq_file2, "%s_%d_remain_R2.fastq.srt", outputFilename, last_round_num);

    report_file = NULL;
    candid_file = NULL;

    window_size = ws;
    step = 3;

    pre_contig = -1;
    pre_chr = "-";

    RegionalHashTable *rht;
    for (int i = 0; i < MAXHTLISTSIZE; ++i) {
        rht = new RegionalHashTable(ws, 0, 0);
        ind2ht.insert(make_pair(i, rht));

        removables.insert(i);

        gid2ind.insert(make_pair(INF - i, i));
        gids.push_back(INF - i);
    }


    int max_kmer_cnt = (maxReadLength - window_size) / step + 1;
    bc1.chains = (chain_t *) malloc(maxChainLen * sizeof(chain_t));
    for (int i = 0; i < maxChainLen; i++)
        bc1.chains[i].frags = (fragment_t *) malloc(max_kmer_cnt * sizeof(fragment_t));
    bc2.chains = (chain_t *) malloc(maxChainLen * sizeof(chain_t));
    for (int i = 0; i < maxChainLen; i++)
        bc2.chains[i].frags = (fragment_t *) malloc(max_kmer_cnt * sizeof(fragment_t));

    // circ_type.push_back("SingleTranscriptomeCircRNA");
    // circ_type.push_back("MultiTranscriptomeCircRNA");
    // circ_type.push_back("NovelCircRNA");

    circ_type.push_back("STC");
    circ_type.push_back("MTC");
    circ_type.push_back("NC");
}

ProcessCirc::~ProcessCirc(void) {
    check_removables(MAXUB);

    close_file(report_file);
    close_file(candid_file);

    for (int i = 0; i < maxChainLen; i++)
        free(bc1.chains[i].frags);

    for (int i = 0; i < maxChainLen; i++)
        free(bc2.chains[i].frags);

    free(bc1.chains);
    free(bc2.chains);

    // Remove intermediate files
    if (finalCleaning) {
        char command[FILE_NAME_MAX_LEN + 100];
        sprintf(command, "rm %s*_remain_R*.fastq{,.srt}", outputFilename);

        system(command);
    }
}

bool record_ptr_compare(Record *r, Record *o) {
    if (r->mr->contig_num != o->mr->contig_num)
        return r->mr->contig_num < o->mr->contig_num;
    if (r->mr->chr_r1 != o->mr->chr_r1)
        return r->mr->chr_r1 < o->mr->chr_r1;
    return r->mr->spos_r1 < o->mr->spos_r1;

}

void ProcessCirc::sort_fq_internal(char *fqname) {
    Logger::instance().debug("Filename: %s\n", fqname);
    fflush(stdout);

    char *sorted_file = (char *) malloc(FILE_NAME_MAX_LEN);
    char *wmode = (char *) malloc(FILE_NAME_MAX_LEN);

    sprintf(sorted_file, "%s.srt", fqname);
    sprintf(wmode, "%c", 'w');

    FASTQParser fq_parser;
    fq_parser.init();
    fq_parser.reset(fqname);
    vector <RecordStr> all_records;
    Record *current_record;

    while ((current_record = fq_parser.get_next_read(0)) != NULL) {
        all_records.push_back(RecordStr(current_record));
    }

    sort(all_records.begin(), all_records.end());

    FILE *sorted_fout = open_file(sorted_file, wmode);

    char sep = '\n';
    MatchedRead *mr;
    char comment[400];
    char r1_dir, r2_dir;
    for (unsigned int i = 0; i < all_records.size(); ++i) {
        mr = &(all_records[i].mr);

        r1_dir = (mr->r1_forward) ? '+' : '-';
        r2_dir = (mr->r2_forward) ? '+' : '-';

		sprintf(comment, " %" PRId64 " %d %s %u %u %d %u %u %c %d %s %u %u %d %u %u %c %d %d %d %d %d",
                mr->genome_spos, mr->type,
                mr->chr_r1.c_str(), mr->spos_r1, mr->epos_r1, mr->mlen_r1,
                mr->qspos_r1, mr->qepos_r1, r1_dir, mr->ed_r1,
                mr->chr_r2.c_str(), mr->spos_r2, mr->epos_r2, mr->mlen_r2,
                mr->qspos_r2, mr->qepos_r2, r2_dir, mr->ed_r2,
                mr->tlen, mr->junc_num, mr->gm_compatible, mr->contig_num);

        fprintf(sorted_fout, "@%s%s%c%s%c%s%c%s\n", all_records[i].rname.c_str(), comment, sep,
                all_records[i].seq.c_str(), sep, all_records[i].comment.c_str(), sep, all_records[i].qual.c_str());
    }

    close_file(sorted_fout);

    free(sorted_file);
    free(wmode);
}

int ProcessCirc::sort_fq(char *fqname) {
    Logger::instance().debug("Filename: %s\n", fqname);
    fflush(stdout);
    if (!system(NULL)) {
        fprintf(stdout, "Error: Unable to run system call.\n");
        exit(EXIT_FAILURE);
    }

    char command[FILE_NAME_MAX_LEN + 500];
    sprintf(command, "cat %s | paste - - - - | sort --parallel=%d -S 8G -k2,2n | tr \"\t\" \"\n\" > %s.srt", fqname,
            threadCount, fqname);

    int ret = system(command);
    return ret;
}

void ProcessCirc::do_process(void) {

    /**********************/
    /**Loading Hash Table**/
    /**********************/

    double cputime_start = get_cpu_time();
    double realtime_start = get_real_time();
    double cputime_curr;
    double realtime_curr;

    vector <ContigLen> orig_contig_len;
    genome_packer.load_index_info(orig_contig_len);

    char index_file[FILE_NAME_MAX_LEN];
    strcpy(index_file, genome_packer.get_index_fname().c_str());

    initCommon();

    THREAD_COUNT = threadCount;
    for (int i = 0; i < 255; i++)
        THREAD_ID[i] = i;

    if (!checkHashTable(index_file))
        return;

    if (!initLoadingCompressedGenomeMeta(index_file))
        return;

    kmer = WINDOW_SIZE + checkSumLength;

    genome_seeder.init(LOADCONTIGSTRINMEM);

    /**********************/
    /**Finised Loading HT**/
    /**********************/

    /*******************/
    /**GTF Parser Init**/
    /*******************/

    if (stage == 1) {
        Logger::instance().info("Loading GTF file...\n");

        gtf_parser.init(gtfFilename, orig_contig_len);
        if (!gtf_parser.load_gtf()) {
            fprintf(stdout, "Error in reading GTF file.\n");
            exit(1);
        }

        cputime_curr = get_cpu_time();
        realtime_curr = get_real_time();

        Logger::instance().info("Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n", cputime_curr - cputime_start,
                                realtime_curr - realtime_start);

        cputime_start = cputime_curr;
        realtime_start = realtime_curr;
    }

    /*******************/
    /**Finished GTF PI**/
    /*******************/

    double fq_cputime_start = cputime_start;
    double fq_realtime_start = realtime_start;

    bool is_pe = pairedEnd;

    FASTQParser fq_parser1;
    fq_parser1.init();
    fq_parser1.reset(fq_file1);
    Record *current_record1 = NULL;

    FASTQParser fq_parser2;
    fq_parser2.init();
    Record *current_record2 = NULL;
    if (is_pe) {
        fq_parser2.reset(fq_file2);
    }

    open_candid_file();

    int line = 0;
    //while ( (current_record1 = fq_parser1.get_next()) != NULL ) { // go line by line on fastq file
    while (true) {

        current_record1 = fq_parser1.get_next_read(0);

        if (current_record1 == NULL)
            break;
        if (is_pe)
            current_record2 = fq_parser2.get_next_read(0);

        line++;
        vafprintf(2, stderr, "Line: %d\n", line);

        if (line % LINELOG == 0) {
            cputime_curr = get_cpu_time();
            realtime_curr = get_real_time();

//			fprintf(stdout, "[P] %d reads in %.2lf CPU sec (%.2lf real sec)\t Look ups: %u\n", line, cputime_curr - cputime_start, realtime_curr - realtime_start, lookup_cnt);
//			fflush(stdout);

            cputime_start = cputime_curr;
            realtime_start = realtime_curr;

            lookup_cnt = 0;
        }

        while (pre_contig != current_record1->mr->contig_num) {
            load_genome();
            refresh_hash_table_list();
        }

        //if (pre_chr != current_record1->mr->chr_r1)
        //	refresh_hash_table_list();
        //pre_chr = current_record1->mr->chr_r1;

        if (is_pe)
            call_circ(current_record1, current_record2);
    }

    refresh_hash_table_list();

    report_events();

    cputime_curr = get_cpu_time();
    realtime_curr = get_real_time();

    Logger::instance().info.set_prefix("[INFO] ");
    Logger::instance().info("CircRNA detection Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n",
                            cputime_curr - fq_cputime_start, realtime_curr - fq_realtime_start);

    finalizeLoadingCompressedGenome();
    genome_seeder.finalize();
}

// PE
void ProcessCirc::call_circ(Record *current_record1, Record *current_record2) {
    fullmap_seq = NULL;
    remain_seq = NULL;
    r1_seq = NULL;
    r2_seq = NULL;

    fullmap_seq_len = 0;
    remain_seq_len = 0;
    r1_seq_len = 0;
    r2_seq_len = 0;

    MatchedRead mr = *(current_record1->mr);
    vafprintf(2, stderr, "%s\n%s\n", current_record1->seq, current_record2->seq);
    vafprintf(2, stderr, "%s\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%s\t%u\t%u\t%d\t%u\t%u\t%d\t%d\t%d\t%d\t%d\n",
              current_record1->rname,
              mr.chr_r1.c_str(), mr.spos_r1, mr.epos_r1, mr.mlen_r1, mr.qspos_r1, mr.qepos_r1, mr.ed_r1,
              mr.chr_r2.c_str(), mr.spos_r2, mr.epos_r2, mr.mlen_r2, mr.qspos_r2, mr.qepos_r2, mr.ed_r2,
              mr.tlen, mr.junc_num, mr.gm_compatible, mr.type);

    if (mr.type == CHIBSJ) {
        call_circ_single_split(current_record1, current_record2);
    } else if (mr.type == CHI2BSJ) {
        call_circ_double_split(current_record1, current_record2);
    }
}

void ProcessCirc::call_circ_single_split(Record *current_record1, Record *current_record2) {
    MatchedRead mr = *(current_record1->mr);
    bool r1_partial = mr.mlen_r1 < mr.mlen_r2;
    remain_seq = (r1_partial) ? ((mr.r1_forward) ? current_record1->seq : current_record1->rcseq) :
                 ((mr.r2_forward) ? current_record2->seq : current_record2->rcseq);
    fullmap_seq = (!r1_partial) ? ((mr.r1_forward) ? current_record1->seq : current_record1->rcseq) :
                  ((mr.r2_forward) ? current_record2->seq : current_record2->rcseq);

    remain_seq_len = (r1_partial) ? current_record1->seq_len : current_record2->seq_len;
    fullmap_seq_len = (!r1_partial) ? current_record1->seq_len : current_record2->seq_len;

    // convert to position on contig
    gtf_parser.chrloc2conloc(mr.chr_r1, mr.spos_r1, mr.epos_r1);
    gtf_parser.chrloc2conloc(mr.chr_r2, mr.spos_r2, mr.epos_r2);

    // fill MatchedMate
    MatchedMate mm_r1(mr, 1, current_record1->seq_len, r1_partial);
    MatchedMate mm_r2(mr, 2, current_record2->seq_len, !r1_partial);

    //if (r1_partial)
    //	remove_side_introns(mm_r1, current_record1->seq_len);
    //else
    //	remove_side_introns(mm_r2, current_record2->seq_len);

    uint32_t qspos = (r1_partial) ? (((mm_r1.qspos - 1) > (current_record1->seq_len - mm_r1.qepos)) ? (1) :
                     (mm_r1.qepos + 1)) : (((mm_r2.qspos - 1) > (current_record2->seq_len - mm_r2.qepos)) ? (1) :
                     (mm_r2.qepos + 1));

    uint32_t qepos = (r1_partial) ? (((mm_r1.qspos - 1) > (current_record1->seq_len - mm_r1.qepos)) ? (mm_r1.qspos - 1)
                                  : (current_record1->seq_len))
                                  : (((mm_r2.qspos - 1) > (current_record2->seq_len - mm_r2.qepos)) ? (mm_r2.qspos - 1)
                                  : (current_record2->seq_len));


    int whole_seq_len = (r1_partial) ? current_record1->seq_len : current_record2->seq_len;
    int remain_len = qepos - qspos + 1;
    if (qepos < qspos or remain_len < window_size) {    // it was fully mapped
        return;
    }

    const IntervalInfo<GeneInfo> *gene_info = gtf_parser.get_gene_overlap(mm_r1.spos, false);
    bool found = (gene_info != NULL);
    if (!found) {
        vafprintf(2, stderr, "Gene not found!\n");
        return;
    }
    vafprintf(2, stderr, "# Gene overlaps: %d\n", gene_info->seg_list.size());


    ConShift con_shift;

    check_removables(mm_r1.spos);
    RegionalHashTable *regional_ht;

    CircRes best_cr;
    best_cr.type = NF;
    for (unsigned int i = 0; i < gene_info->seg_list.size(); ++i) {
        uint32_t gene_len = gene_info->seg_list[i].end - gene_info->seg_list[i].start + 1;

        regional_ht = get_hash_table_smart(gene_info->seg_list[i]);

        vafprintf(2, stderr, "R%d partial: [%d-%d]\n", (int) (!r1_partial) + 1, qspos, qepos);
        vafprintf(2, stderr, "%s\n", remain_seq);

        chaining(qspos, qepos, regional_ht, remain_seq, gene_len, gene_info->seg_list[i].start, bc1);

        if (bc1.best_chain_count <= 0)
            continue;

        // find rspos and repos for the best chain
        bool forward = (r1_partial) ? (mr.r1_forward) : (mr.r2_forward);
        int dir = (forward) ? 1 : -1;

        for (int j = 0; j < minM(bc1.best_chain_count, TOPCHAIN); ++j) {
// 			fprintf(stderr, "Trying chain #%d\n\n\n", j);
            MatchedMate partial_mm;
            find_exact_coord(mm_r1, mm_r2, partial_mm, dir, qspos, remain_seq, remain_len, whole_seq_len,
                             bc1.chains[j]);

            if (partial_mm.type == CONCRD) {
                con_shift = gtf_parser.get_shift(contigNum, mm_r1.spos);
                vafprintf(2, stderr, "Coordinates: [%d-%d]\n", partial_mm.spos - con_shift.shift,
                          partial_mm.epos - con_shift.shift);


                CircRes cr;

// 				fprintf(stderr, "%s\t", current_record1->rname);

                int type = check_split_map(mm_r1, mm_r2, partial_mm, r1_partial, cr);
                print_split_mapping(current_record1->rname, mm_r1, mm_r2, partial_mm, con_shift);

// 				fprintf(stderr, "\n");

                fprintf(candid_file, "%d\n", type);

                if (type < CR) {
                    best_cr.type = type;
                    return;
                }

                if (type >= CR and type <= MCR and type < best_cr.type) {
                    best_cr.chr = con_shift.contig;
                    best_cr.rname = current_record1->rname;
                    best_cr.spos = cr.spos - con_shift.shift;
                    best_cr.epos = cr.epos - con_shift.shift;
                    best_cr.type = type;
                    best_cr.start_signal = cr.start_signal;
                    best_cr.end_signal = cr.end_signal;
                    best_cr.start_bp_ref = cr.start_bp_ref;
                    best_cr.end_bp_ref = cr.end_bp_ref;

                    if (type == CR) {
                        circ_res.push_back(best_cr);
                        return;
                    }
                }
            }
        }
    }

    if (best_cr.type >= CR and best_cr.type <= MCR)
        circ_res.push_back(best_cr);
}

void ProcessCirc::call_circ_double_split(Record *current_record1, Record *current_record2) {
    MatchedRead mr = *(current_record1->mr);
    vafprintf(2, stderr, "Double split read...\n");
    char *r1_remain_seq = (mr.r1_forward) ? current_record1->seq : current_record1->rcseq;
    char *r2_remain_seq = (mr.r2_forward) ? current_record2->seq : current_record2->rcseq;

    r1_seq = r1_remain_seq;
    r2_seq = r2_remain_seq;
    r1_seq_len = current_record1->seq_len;
    r2_seq_len = current_record2->seq_len;

    uint32_t r1_qspos = ((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (1) : (mr.qepos_r1 + 1);
    uint32_t r2_qspos = ((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (1) : (mr.qepos_r2 + 1);

    uint32_t r1_qepos = ((mr.qspos_r1 - 1) > (current_record1->seq_len - mr.qepos_r1)) ? (mr.qspos_r1 - 1)
                        : (current_record1->seq_len);
    uint32_t r2_qepos = ((mr.qspos_r2 - 1) > (current_record2->seq_len - mr.qepos_r2)) ? (mr.qspos_r2 - 1)
                        : (current_record2->seq_len);


    int r1_remain_len = r1_qepos - r1_qspos + 1;
    int r2_remain_len = r2_qepos - r2_qspos + 1;

    // it was fully mapped
    if (r1_remain_len < window_size and r2_remain_len < window_size) {
        return;
    }

    // single split
    if (r1_remain_len < window_size or r2_remain_len < window_size) {
        call_circ_single_split(current_record1, current_record2);
    }

    // convert to position on contig
    gtf_parser.chrloc2conloc(mr.chr_r1, mr.spos_r1, mr.epos_r1);
    gtf_parser.chrloc2conloc(mr.chr_r2, mr.spos_r2, mr.epos_r2);

    const IntervalInfo<GeneInfo> *gene_info = gtf_parser.get_gene_overlap(mr.spos_r1, false);
    bool found = (gene_info != NULL);
    if (!found) {
        vafprintf(2, stderr, "Gene not found!\n");
        return;
    }
    vafprintf(2, stderr, "# Gene overlaps: %d\n", gene_info->seg_list.size());

    // fill MatchedMate
    MatchedMate mm_r1(mr, 1, current_record1->seq_len, true);
    MatchedMate mm_r2(mr, 2, current_record2->seq_len, true);
    ConShift con_shift;

    check_removables(mr.spos_r1);
    RegionalHashTable *regional_ht;

    CircRes best_cr;
    best_cr.type = NF;
    for (unsigned int i = 0; i < gene_info->seg_list.size(); ++i) {
        uint32_t gene_len = gene_info->seg_list[i].end - gene_info->seg_list[i].start + 1;

        regional_ht = get_hash_table_smart(gene_info->seg_list[i]);

        vafprintf(2, stderr, "R1 partial: [%d-%d]\nremain: %s\n", r1_qspos, r1_qepos, r1_remain_seq);
        chaining(r1_qspos, r1_qepos, regional_ht, r1_remain_seq, gene_len, gene_info->seg_list[i].start, bc1);
        vafprintf(2, stderr, "R2 partial: [%d-%d]\nremain: %s\n", r2_qspos, r2_qepos, r2_remain_seq);
        chaining(r2_qspos, r2_qepos, regional_ht, r2_remain_seq, gene_len, gene_info->seg_list[i].start, bc2);

        if (bc1.best_chain_count <= 0 and bc2.best_chain_count <= 0) {
// 			fprintf(stderr, "LostRead: Not chained!\n");
            continue;
        }

        if (bc1.best_chain_count <= 0 or bc2.best_chain_count <= 0) {
// 			fprintf(stderr, "Doing one single split check instead of double!\n");
            call_circ_single_split(current_record1, current_record2);
            continue;
        }

        for (int j = 0; j < minM(bc1.best_chain_count, TOPCHAIN); ++j) {
            for (int k = 0; k < minM(bc2.best_chain_count, TOPCHAIN); ++k) {
                // find rspos and repos for the best chain
                MatchedMate r1_partial_mm;
                MatchedMate r2_partial_mm;

                set_mm(bc1.chains[j], r1_qspos, r1_remain_len, mm_r1.dir, r1_partial_mm);
                set_mm(bc2.chains[k], r2_qspos, r2_remain_len, mm_r2.dir, r2_partial_mm);

                overlap_to_spos(mm_r1);
                overlap_to_spos(mm_r2);
                overlap_to_spos(r1_partial_mm);
                overlap_to_spos(r2_partial_mm);

                vector <uint32_t> common_tid;
// 				if (! same_transcript(mm_r1.exons_spos, mm_r2.exons_spos, r1_partial_mm.exons_spos, r2_partial_mm.exons_spos, common_tid))
// 					continue;

                vector <MatchedMate> segments;
                segments.push_back(mm_r1);
                segments.push_back(mm_r2);
                segments.push_back(r1_partial_mm);
                segments.push_back(r2_partial_mm);
                if (!same_transcript(segments, 4, common_tid))
                    continue;

                bool success = false;
                if (bc1.chains[j].frags[0].rpos <= bc2.chains[k].frags[0].rpos) {
                    success = extension.extend_both_mates(bc1.chains[j], bc2.chains[k], common_tid, r1_remain_seq,
                                                          r2_remain_seq, r1_qspos, r2_qspos, r1_qepos, r2_qepos,
                                                          r1_partial_mm, r2_partial_mm);
                } else {
                    success = extension.extend_both_mates(bc2.chains[k], bc1.chains[j], common_tid, r2_remain_seq,
                                                          r1_remain_seq, r2_qspos, r1_qspos, r2_qepos, r1_qepos,
                                                          r2_partial_mm, r1_partial_mm);
                }

                if (!success) {
// 					fprintf(stderr, "LostRead: Not extended!\n");
                    continue;
                }

                if (r1_partial_mm.type == CONCRD and r2_partial_mm.type == CONCRD) {
                    con_shift = gtf_parser.get_shift(contigNum, mm_r1.spos);
                    vafprintf(2, stderr, "R1 Partial Coordinates: [%d-%d]\n", r1_partial_mm.spos - con_shift.shift,
                              r1_partial_mm.epos - con_shift.shift);
                    vafprintf(2, stderr, "R2 Partial Coordinates: [%d-%d]\n", r2_partial_mm.spos - con_shift.shift,
                              r2_partial_mm.epos - con_shift.shift);

                    CircRes cr;
                    int type = check_split_map(mm_r1, mm_r2, r1_partial_mm, r2_partial_mm, cr);
                    print_split_mapping(current_record1->rname, mm_r1, mm_r2, r1_partial_mm, r2_partial_mm, con_shift);
                    fprintf(candid_file, "%d\n", type);
                    if (type < CR) {
                        best_cr.type = type;
                        return;
                    }

                    if (type >= CR and type <= MCR and type < best_cr.type) {
                        best_cr.chr = con_shift.contig;
                        best_cr.rname = current_record1->rname;
                        best_cr.spos = cr.spos - con_shift.shift;
                        best_cr.epos = cr.epos - con_shift.shift;
                        best_cr.type = type;
                        best_cr.start_signal = cr.start_signal;
                        best_cr.end_signal = cr.end_signal;
                        best_cr.start_bp_ref = cr.start_bp_ref;
                        best_cr.end_bp_ref = cr.end_bp_ref;

                        if (type == CR) {
                            circ_res.push_back(best_cr);
                            return;
                        }
                    }
                }
            }
        }
    }

    if (best_cr.type >= CR and best_cr.type <= MCR) {
        circ_res.push_back(best_cr);
    } else {
        call_circ_single_split(current_record1, current_record2);
    }
}

void ProcessCirc::binning(uint32_t qspos, uint32_t qepos, RegionalHashTable *regional_ht, char *remain_seq,
                          uint32_t gene_len) {
    int bin_num = gene_len / BINSIZE + 1;
    int bins[bin_num];
    int max_id = 0;
    memset(bins, 0, bin_num * sizeof(int));

    for (uint32_t i = qspos - 1; i <= qepos - window_size; i += step) {
        GIMatchedKmer *gl = regional_ht->find_hash(regional_ht->hash_val(remain_seq + i));
        if (gl == NULL) {
            vafprintf(2, stderr, "Hash val not found!!!\n");
        }

        vafprintf(2, stderr, "Occ: %d\n", gl->frag_count);

        for (uint32_t j = 0; j < gl->frag_count; j++) {
            int bin_id = gl->frags[j].info / BINSIZE;
            bins[bin_id]++;

            if (bins[bin_id] > bins[max_id])
                max_id = bin_id;

            vafprintf(2, stderr, "%d - %d\t", gl->frags[j].info, bins[bin_id]);
        }
        vafprintf(2, stderr, "\n");
    }
    vafprintf(2, stderr, "Biggest bin: bin[%d][%d - %d] = %d\n", max_id, max_id * BINSIZE, (max_id + 1) * BINSIZE - 1,
              bins[max_id]);

}

void ProcessCirc::chaining(uint32_t qspos, uint32_t qepos, RegionalHashTable *regional_ht, char *remain_seq,
                           uint32_t gene_len, uint32_t shift, chain_list &bc) {
    int seq_len = qepos - qspos + 1;
    if (seq_len < window_size) {
        bc.best_chain_count = 0;
        return;
    }
    int kmer_cnt = (seq_len - window_size) / step + 1;
    GIMatchedKmer fl[kmer_cnt + 1];

    int l = 0;
    for (uint32_t i = qspos - 1; i <= qepos - window_size; i += step) {
        GIMatchedKmer *gl = regional_ht->find_hash(regional_ht->hash_val(remain_seq + i));
        if (gl == NULL) {    // has N inside kmer
            vafprintf(2, stderr, "Hash val not found!!!\n");
            continue;
        }
        fl[l] = *gl;
        fl[l].qpos = i;

        //vafprintf(2, stderr, "Occ: %d\n", fl[l].frag_count);

        //for (int j = 0; j < gl->frag_count; j++) {
        //	vafprintf(2, stderr, "%d - %d\t", gl->frags[j].info, fl[l].frags[j].info);
        //}
        //vafprintf(2, stderr, "\n");

        if (fl[l].frag_count > seedLim)
            fl[l].frag_count = 0;
        l++;
    }

    kmer_cnt = l;

    chain_seeds_sorted_kbest2(qepos, fl, bc, window_size, kmer_cnt, shift);

    vafprintf(1, stderr, "Chaining score:%.4f,\t len: %lu\n", bc.chains[0].score, (unsigned long) bc.best_chain_count);

    int allowed_missed_kmers = (qepos - qspos + 1) / 20 * 3 + 1;
    vafprintf(2, stderr, "Allowed missing kmers: %d\n", allowed_missed_kmers);

    int missing;
    int least_miss = INF;
    for (int j = 0; j < bc.best_chain_count; j++) {
        missing = kmer_cnt - bc.chains[j].chain_len;
        vafprintf(2, stderr, "Actual missing: %d\n", missing);
        if (missing > least_miss) {
            // this one does not count
            bc.best_chain_count = j;
            break;
        }

        least_miss = missing;

        for (uint32_t i = 0; i < bc.chains[j].chain_len; i++) {
            vafprintf(1, stderr, "#%d\tfrag[%d]: %lu\t%d\t%d\n", j, i, bc.chains[j].frags[i].rpos - shift,
                      bc.chains[j].frags[i].qpos, bc.chains[j].frags[i].len);
        }
    }
}

bool ProcessCirc::find_exact_coord(MatchedMate &mm_r1, MatchedMate &mm_r2, MatchedMate &partial_mm,
                                   int dir, uint32_t qspos, char *rseq, int rlen, int whole_len, chain_t &bc) {

    set_mm(bc, qspos, rlen, dir, partial_mm);
    --qspos;    // convert to 0-based

    overlap_to_spos(mm_r1);
    overlap_to_spos(mm_r2);
    overlap_to_spos(partial_mm);

    vector <uint32_t> common_tid;
// 	bool success = same_transcript(mm_r1.exons_spos, mm_r2.exons_spos, partial_mm.exons_spos, common_tid);

    vector <MatchedMate> segments;
    segments.push_back(mm_r1);
    segments.push_back(mm_r2);
    segments.push_back(partial_mm);
    bool success = same_transcript(segments, 3, common_tid);
    if (!success) {
// 		fprintf(stderr, "No common transcript!\n");
        return false;
    }

    partial_mm.middle_ed = extension.calc_middle_ed(bc, maxEd, rseq, rlen);
    if (partial_mm.middle_ed > maxEd)
        return false;

    bool extend = true;
    partial_mm.is_concord = false;
    if (bc.chain_len <= 0) {
        partial_mm.type = ORPHAN;
        partial_mm.matched_len = 0;
        return false;
    }

    bool lok;
    bool rok;
    int err = partial_mm.middle_ed;
    if (extend) {
        partial_mm.matched_len = rlen;
        lok = extension.extend_chain_left(common_tid, bc, rseq + qspos, qspos, MINLB, partial_mm, err);
        if (qspos == 0) {
            rok = extension.extend_chain_right(common_tid, bc, rseq, rlen, MAXUB, partial_mm, err);
        } else {
            rok = extension.extend_chain_right(common_tid, bc, rseq, whole_len, MAXUB, partial_mm, err);
        }
        update_match_mate_info(lok, rok, err, partial_mm);
    }

    return partial_mm.type == CONCRD;
}

void ProcessCirc::refresh_hash_table_list(void) {
    for (auto it = ind2ht.begin(); it != ind2ht.end(); ++it) {
        delete it->second;
        //removables.insert(it->first);

    }

    removables.clear();
    gids.clear();
    gid2ind.clear();
    ind2ht.clear();

}

void ProcessCirc::check_removables(uint32_t rspos) {
    for (auto it = ind2ht.begin(); it != ind2ht.end(); ++it) {
        if (rspos > it->second->gene_epos) {
            removables.insert(it->first);
            // fprintf(stderr, "Removed gid: %d --removables size: %d\n", gids[it->first], removables.size());
        }
    }
}

RegionalHashTable *ProcessCirc::get_hash_table(const GeneInfo &gene_info) {
    RegionalHashTable *regional_ht;
    RegionalHashTable *new_ht;

    int gene_len = gene_info.end - gene_info.start + 1;

    char gene_seq[gene_len + 1];
    gene_seq[gene_len] = '\0';
    genome_seeder.pac2char_otf(gene_info.start, gene_len, gene_seq);

    new_ht = new RegionalHashTable(window_size, gene_info.start, gene_info.end);

    new_ht->create_table(gene_seq, 0, gene_len);
    regional_ht = new_ht;

    return regional_ht;
}

RegionalHashTable *ProcessCirc::get_hash_table_smart(const GeneInfo &gene_info) {
    uint32_t gid, removable_ind;
    RegionalHashTable *regional_ht;
    RegionalHashTable *new_ht;

    int gene_len = gene_info.end - gene_info.start + 1;
    gid = gene_info.gene_id;
    //fprintf(stderr, "Gene id: %d\n", gid);

    if (gid2ind.find(gid) != gid2ind.end()) {
        uint32_t ind = gid2ind[gid];
        regional_ht = ind2ht[ind];
        //fprintf(stderr, "Found gid: %d\t---Skipping HT create\n", gid);
    } else {
        char gene_seq[gene_len + 1];
        gene_seq[gene_len] = '\0';
        //genome_seeder.pac2char(gene_info.start, gene_len, gene_seq);
        genome_seeder.pac2char_otf(gene_info.start, gene_len, gene_seq);

        if (removables.size() == 0) {
            new_ht = new RegionalHashTable(window_size, gene_info.start, gene_info.end);

            uint32_t new_ind = ind2ht.size();
            ind2ht.insert(make_pair(new_ind, new_ht));
            gid2ind.insert(make_pair(gid, new_ind));
            gids.push_back(gid);
            new_ht->create_table(gene_seq, 0, gene_len);
            regional_ht = new_ht;
            //fprintf(stderr, "Allocated new HT gid: %d, %s [%u-%u]\n", gid, mr.chr_r1.c_str(), mr.spos_r1, mr.spos_r2);
        } else {
            removable_ind = *(removables.cbegin());
            removables.erase(removables.cbegin());
            gid2ind.erase(gids[removable_ind]);
            gid2ind.insert(make_pair(gid, removable_ind));

            //fprintf(stderr, "Removed gid: %d\tTotal: %d\n", gids[removable_ind], gid2ind.size());

            gids[removable_ind] = gid;

            regional_ht = ind2ht[removable_ind];

            regional_ht->gene_spos = gene_info.start;
            regional_ht->gene_epos = gene_info.end;
            regional_ht->create_table(gene_seq, 0, gene_len);

        }
        //fprintf(stderr, "Added gid: %d\tTotal: %d\n", gid, gid2ind.size());
        //fprintf(stderr, "----------------------------------------------------------\n");
        //fprintf(stderr, "Map Info \n");

        //for (auto it = ind2ht.begin(); it != ind2ht.end(); ++it) {
        //	fprintf(stderr, "ind: %d,\tgid: %d\n", it->first, gids[it->first]);
        //}
        //fprintf(stderr, "----------------------------------------------------------\n");
    }

    return regional_ht;
}

// for non-overlapping split mates
int ProcessCirc::check_split_map(MatchedMate &mm_r1, MatchedMate &mm_r2, MatchedMate &partial_mm, bool r1_partial,
                                 CircRes &cr) {
    int valid = NF;
    int split_read_ed = 0;

    if (r1_partial) {
        split_read_ed = mm_r1.right_ed + mm_r1.left_ed + mm_r1.middle_ed +
                        partial_mm.right_ed + partial_mm.left_ed + partial_mm.middle_ed;

        if (mm_r1.qspos < partial_mm.qspos)
            valid = final_check(mm_r2, mm_r1, partial_mm, cr);
        else
            valid = final_check(mm_r2, partial_mm, mm_r1, cr);
    } else {
        split_read_ed = mm_r2.right_ed + mm_r2.left_ed + mm_r2.middle_ed +
                        partial_mm.right_ed + partial_mm.left_ed + partial_mm.middle_ed;
        if (mm_r2.qspos < partial_mm.qspos)
            valid = final_check(mm_r1, mm_r2, partial_mm, cr);
        else
            valid = final_check(mm_r1, partial_mm, mm_r2, cr);
    }

    if (split_read_ed > maxEd) {
        //fprintf(stderr, "non-overlapping edit distance exceeded! ed: %d\n", split_read_ed);
        valid = UD;
    }
    return valid;
}

// for overlapping split mates
int ProcessCirc::check_split_map(MatchedMate &mm_r1_1, MatchedMate &mm_r2_1, MatchedMate &mm_r1_2, MatchedMate &mm_r2_2,
                                 CircRes &cr) {

    int r1_ed = mm_r1_1.right_ed + mm_r1_1.left_ed + mm_r1_1.middle_ed +
                mm_r1_2.right_ed + mm_r1_2.left_ed + mm_r1_2.middle_ed;
    int r2_ed = mm_r2_1.right_ed + mm_r2_1.left_ed + mm_r2_1.middle_ed +
                mm_r2_2.right_ed + mm_r2_2.left_ed + mm_r2_2.middle_ed;

    if (r1_ed > maxEd or r2_ed > maxEd) {
        //fprintf(stderr, "overlapping edit distance exceeded! ed1: %d, ed2: %d\n", r1_ed, r2_ed);
        return UD;
    }

    MatchedMate mm_r1_l = (mm_r1_1.spos <= mm_r1_2.spos) ? mm_r1_1 : mm_r1_2;
    MatchedMate mm_r1_r = (mm_r1_1.spos <= mm_r1_2.spos) ? mm_r1_2 : mm_r1_1;

    MatchedMate mm_r2_l = (mm_r2_1.spos <= mm_r2_2.spos) ? mm_r2_1 : mm_r2_2;
    MatchedMate mm_r2_r = (mm_r2_1.spos <= mm_r2_2.spos) ? mm_r2_2 : mm_r2_1;

    bool r1_regular_bsj = (mm_r1_l.qspos < mm_r1_r.qspos);
    bool r2_regular_bsj = (mm_r2_l.qspos < mm_r2_r.qspos);

    string ssignal_final;
    string esignal_final;

    string ssignal1;
    string esignal1;
    string ssignal2;
    string esignal2;

    // RF or FR check
    if (r1_regular_bsj and r2_regular_bsj) {
        if (mm_r1_l.dir == 1) {
            if (mm_r1_r.spos <= mm_r2_l.spos)
                return FR;
            else if (mm_r1_l.epos >= mm_r2_r.epos)
                return RF;
        }
        if (mm_r1_l.dir == -1) {
            if (mm_r2_r.spos <= mm_r1_l.spos)
                return FR;
            else if (mm_r2_l.epos >= mm_r1_r.epos)
                return RF;
        }
    }

        // single BSJ on R2
    else if (r1_regular_bsj and !r2_regular_bsj) {
        MatchedMate full_mm = mm_r1_l;
        if (!full_mm.merge_to_right(mm_r1_r))
            return UD;
        remain_seq = r2_seq;
        remain_seq_len = r2_seq_len;
        return final_check(full_mm, mm_r2_l, mm_r2_r, cr);
    }

        // single BSJ on R1
    else if (!r1_regular_bsj and r2_regular_bsj) {
        MatchedMate full_mm = mm_r2_l;
        if (!full_mm.merge_to_right(mm_r2_r))
            return UD;
        remain_seq = r1_seq;
        remain_seq_len = r1_seq_len;
        return final_check(full_mm, mm_r1_l, mm_r1_r, cr);
    }

        // BSJ on the overlap
    else if (!r1_regular_bsj and !r2_regular_bsj) {
        if (mm_r1_l.spos == mm_r2_l.spos and mm_r1_r.epos == mm_r2_r.epos) {
            overlap_to_spos(mm_r1_l);
            overlap_to_epos(mm_r1_r);

            // if (mm_r1_l.exons_spos == NULL or mm_r1_r.exons_epos == NULL) {
            // 	cr.set_bp(mm_r1_l.spos - mm_r1_l.sclen_left, mm_r1_r.epos + mm_r1_r.sclen_right);
            // 	return MCR;
            // }

            vector <pu32i> end_tids;
            int diff;
            int ind_epos = mm_r1_r.exon_ind_epos;
            const IntervalInfo<UniqSeg> *curr_seg;
            curr_seg = gtf_parser.get_interval(ind_epos);
            while (mm_r1_r.spos < curr_seg->epos) {
                for (unsigned int i = 0; i < curr_seg->seg_list.size(); ++i) {
                    diff = mm_r1_r.epos + mm_r1_r.sclen_right - curr_seg->seg_list[i].end;
                    if (abs(diff) <= BPRES) {
                        for (unsigned int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
                            end_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
                        }
                    }
                }
                --ind_epos;
                curr_seg = gtf_parser.get_interval(ind_epos);
            }

            vector <pu32i> start_tids;
            int ind_spos = mm_r1_l.exon_ind_spos;
            curr_seg = gtf_parser.get_interval(ind_spos);
            while (mm_r1_l.epos > curr_seg->spos) {
                for (unsigned int i = 0; i < curr_seg->seg_list.size(); ++i) {
                    diff = mm_r1_l.spos - mm_r1_l.sclen_left - curr_seg->seg_list[i].start;
                    if (abs(diff) <= BPRES) {
                        for (unsigned int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
                            start_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
                        }
                    }
                }
                ++ind_spos;
                curr_seg = gtf_parser.get_interval(ind_spos);
            }

            int sdiff, ediff;
            int best_ed1 = maxEd + 1;
            int best_ed2 = maxEd + 1;
            uint32_t qcutpos;
            uint32_t beg_bp;
            uint32_t end_bp;
            vector <uint32_t> common_tid;
            for (unsigned int i = 0; i < start_tids.size(); ++i) {
                for (unsigned int j = 0; j < end_tids.size(); ++j) {
                    sdiff = start_tids[i].second;
                    ediff = end_tids[j].second;
                    if (start_tids[i].first == end_tids[j].first and sdiff == ediff) {
                        common_tid.clear();
                        common_tid.push_back(start_tids[i].first);
                        beg_bp = mm_r1_l.spos - mm_r1_l.sclen_left - sdiff;
                        end_bp = mm_r1_r.epos + mm_r1_r.sclen_right - ediff;

                        qcutpos = mm_r1_r.qepos + mm_r1_r.sclen_right - ediff;
                        int ed1 = split_realignment(qcutpos, beg_bp, end_bp, r1_seq, r1_seq_len, common_tid, mm_r1_r,
                                                    mm_r1_l);

                        if (qcutpos < 2 or qcutpos + 2 > r1_seq_len) {
                            esignal1 = "";
                            ssignal1 = "";
                        } else {
                            esignal1 = string() + r1_seq[qcutpos - 2] + r1_seq[qcutpos - 1];
                            ssignal1 = string() + r1_seq[qcutpos] + r1_seq[qcutpos + 1];
                        }

                        qcutpos = mm_r2_r.qepos + mm_r2_r.sclen_right - ediff;
                        int ed2 = split_realignment(qcutpos, beg_bp, end_bp, r2_seq, r2_seq_len, common_tid, mm_r2_r,
                                                    mm_r2_l);

                        if (qcutpos < 2 or qcutpos + 2 > r2_seq_len) {
                            ssignal2 = "";
                            esignal2 = "";
                        } else {
                            esignal2 = string() + r2_seq[qcutpos - 2] + r2_seq[qcutpos - 1];
                            ssignal2 = string() + r2_seq[qcutpos] + r2_seq[qcutpos + 1];
                        }

                        if (ed1 < best_ed1 and ed2 < best_ed2) {
                            char near_start_bp[10];
                            genome_seeder.pac2char_otf(beg_bp, 2, near_start_bp);

                            char near_end_bp[10];
                            genome_seeder.pac2char_otf(end_bp - 1, 2, near_end_bp);

                            if (ssignal1 == "") {
                                cr.set_bp(beg_bp, end_bp, ssignal2, esignal2, near_start_bp, near_end_bp);
                            } else if (ssignal2 == "") {
                                cr.set_bp(beg_bp, end_bp, ssignal1, esignal1, near_start_bp, near_end_bp);
                            } else {
                                ssignal_final = get_consensus(ssignal1, ssignal2);
                                esignal_final = get_consensus(esignal1, esignal2);

                                cr.set_bp(beg_bp, end_bp, ssignal_final, esignal_final, near_start_bp, near_end_bp);
                            }
                            best_ed1 = ed1;
                            best_ed2 = ed2;
                        }
                    }
                }
            }

            if (best_ed1 <= maxEd and best_ed2 <= maxEd)
                return CR;

            qcutpos = mm_r1_r.qepos + mm_r1_r.sclen_right;
            beg_bp = mm_r1_l.spos - mm_r1_l.sclen_left;
            end_bp = mm_r1_r.epos + mm_r1_r.sclen_right;

            if (qcutpos < 2 or qcutpos > (r1_seq_len - 2) or qcutpos > (r2_seq_len - 2))
                return MCR;

            esignal1 = string() + r1_seq[qcutpos - 2] + r1_seq[qcutpos - 1];
            ssignal1 = string() + r1_seq[qcutpos] + r1_seq[qcutpos + 1];
            esignal2 = string() + r2_seq[qcutpos - 2] + r2_seq[qcutpos - 1];
            ssignal2 = string() + r2_seq[qcutpos] + r2_seq[qcutpos + 1];

            ssignal_final = get_consensus(ssignal1, ssignal2);
            esignal_final = get_consensus(esignal1, esignal2);

            char near_start_bp[10];
            genome_seeder.pac2char_otf(beg_bp, 2, near_start_bp);

            char near_end_bp[10];
            genome_seeder.pac2char_otf(end_bp - 1, 2, near_end_bp);

            cr.set_bp(beg_bp, end_bp, ssignal_final, esignal_final, near_start_bp, near_end_bp);
            if (start_tids.size() > 0 and end_tids.size() > 0)
                return NCR;

            return MCR;
        }
    }
    return UD;
}

// full_mm -> not split mate
// split_mm_left -> left hand side of the split read
// split_mm_right -> right hand side of the split read
int
ProcessCirc::final_check(MatchedMate &full_mm, MatchedMate &split_mm_left, MatchedMate &split_mm_right, CircRes &cr) {
    string ssignal;
    string esignal;

    if (split_mm_left.epos < split_mm_right.spos) {
        if (full_mm.dir == 1) {
            if (full_mm.spos <= split_mm_left.spos)
                return FR;
            else if (full_mm.epos >= split_mm_right.epos)
                return RF;
        }

        if (full_mm.dir == -1) {
            if (full_mm.epos >= split_mm_right.epos)
                return FR;
            else if (full_mm.spos <= split_mm_left.spos)
                return RF;
        }
    }

        //else if (split_mm_right.epos < split_mm_left.spos) {
        // allow detection of short circRNAs
    else if (split_mm_right.spos <= split_mm_left.spos and split_mm_left.epos >= split_mm_right.epos) {
        if (full_mm.spos < split_mm_right.spos) {
            int off = split_mm_right.spos - full_mm.spos;
            int sc_remained = maxSc - full_mm.sclen_left;
            if (off <= sc_remained) {
                full_mm.spos = split_mm_right.spos;
                full_mm.sclen_left += off;
                full_mm.qspos += off;
                full_mm.matched_len -= off;
            }

        }
        if (full_mm.epos > split_mm_left.epos) {
            int off = full_mm.epos - split_mm_left.epos;
            int sc_remained = maxSc - full_mm.sclen_right;
            if (off <= sc_remained) {
                full_mm.epos = split_mm_left.epos;
                full_mm.sclen_right += off;
                full_mm.qepos -= off;
                full_mm.matched_len -= off;
            }
        }
        if (full_mm.spos >= split_mm_right.spos and full_mm.epos <= split_mm_left.epos) {
            // check splice site
            overlap_to_spos(full_mm);
            overlap_to_epos(full_mm);

            overlap_to_spos(split_mm_right);
            overlap_to_epos(split_mm_right);

            overlap_to_spos(split_mm_left);
            overlap_to_epos(split_mm_left);

            // if (split_mm_left.exons_epos == NULL or split_mm_right.exons_spos == NULL) {
            // 	cr.set_bp(split_mm_right.spos - split_mm_right.sclen_left, split_mm_left.epos + split_mm_left.sclen_right);
            // 	return MCR;
            // }

            vector <pu32i> end_tids;
            int diff;
            int ind_epos = split_mm_left.exon_ind_epos;
            const IntervalInfo<UniqSeg> *curr_seg;
            curr_seg = gtf_parser.get_interval(ind_epos);
            while (split_mm_left.spos < curr_seg->epos) {
                for (unsigned int i = 0; i < curr_seg->seg_list.size(); ++i) {
                    diff = split_mm_left.epos + split_mm_left.sclen_right - curr_seg->seg_list[i].end;

                    // fprintf(stderr, "\nend transcripts: ");
                    // for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j)
                    // 	fprintf(stderr, "%s, ", gtf_parser.transcript_ids[contigNum][curr_seg->seg_list[i].trans_id[j]].c_str());
                    // fprintf(stderr, "\n\t(epos, right_sc, exon_end, ediff): (%d, %d, %d, %d)\n\n",
                    // 				split_mm_left.epos, split_mm_left.sclen_right, curr_seg->seg_list[i].end, diff);

                    if (abs(diff) <= BPRES) {
                        for (unsigned int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
                            end_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
                        }
                    }
                }
                --ind_epos;
                curr_seg = gtf_parser.get_interval(ind_epos);
            }

            vector <pu32i> start_tids;
            int ind_spos = split_mm_right.exon_ind_spos;
            curr_seg = gtf_parser.get_interval(ind_spos);
            while (split_mm_right.epos > curr_seg->spos) {
                for (unsigned int i = 0; i < curr_seg->seg_list.size(); ++i) {
                    diff = split_mm_right.spos - split_mm_right.sclen_left - curr_seg->seg_list[i].start;

                    // fprintf(stderr, "\nstart transcripts: ");
                    // for (int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j)
                    // 	fprintf(stderr, "%s, ", gtf_parser.transcript_ids[contigNum][curr_seg->seg_list[i].trans_id[j]].c_str());
                    // fprintf(stderr, "\n\t(spos, left_sc, exon_beg, sdiff): (%d, %d, %d, %d)\n\n",
                    // 				split_mm_right.spos, split_mm_right.sclen_left, curr_seg->seg_list[i].start, diff);

                    if (abs(diff) <= BPRES) {
                        for (unsigned int j = 0; j < curr_seg->seg_list[i].trans_id.size(); ++j) {
                            start_tids.push_back(make_pair(curr_seg->seg_list[i].trans_id[j], diff));
                        }
                    }
                }
                ++ind_spos;
                curr_seg = gtf_parser.get_interval(ind_spos);
            }

            int sdiff, ediff;
            int best_ed = maxEd + 1;
            vector <uint32_t> common_tid;
            for (unsigned int i = 0; i < start_tids.size(); ++i) {
                for (unsigned int j = 0; j < end_tids.size(); ++j) {
                    sdiff = start_tids[i].second;
                    ediff = end_tids[j].second;
                    // if (start_tids[i].first == end_tids[j].first)
                    // fprintf(stderr, "tid: %d -> %s - (spos, left_sc, sdiff): (%d, %d, %d) - (epos, right_sc, ediff): (%d, %d, %d)\n",
                    // start_tids[i].first, gtf_parser.transcript_ids[contigNum][start_tids[i].first].c_str(), split_mm_right.spos, split_mm_right.sclen_left, sdiff, split_mm_left.epos, split_mm_left.sclen_right, ediff);
                    // if (start_tids[i].first == end_tids[j].first and sdiff + ediff == 0) {
                    if (start_tids[i].first == end_tids[j].first and sdiff == ediff) {
                        common_tid.clear();
                        common_tid.push_back(start_tids[i].first);
                        uint32_t qcutpos = split_mm_left.qepos + split_mm_left.sclen_right - ediff;
                        uint32_t beg_bp = split_mm_right.spos - split_mm_right.sclen_left - sdiff;
                        uint32_t end_bp = split_mm_left.epos + split_mm_left.sclen_right - ediff;

                        if (full_mm.sclen_right > 0) {
                            if (full_mm.epos + full_mm.sclen_right > end_bp) {
                                // other mate is crossing end_bp hence split realignment
                                uint32_t fm_qcutpos = full_mm.qepos + (end_bp - full_mm.epos);
                                MatchedMate full_mm_right;
                                int fm_ed = split_realignment(fm_qcutpos, beg_bp, end_bp, fullmap_seq, fullmap_seq_len,
                                                              common_tid, full_mm, full_mm_right);
                                if (fm_ed > maxEd)
                                    continue;
                            } else if (full_mm.sclen_right > maxSc)
                                continue;
                        }

                        if (full_mm.sclen_left > 0) {
                            if (full_mm.spos - full_mm.sclen_left < beg_bp) {
                                // other mate is crossing end_bp hence split realignment
                                uint32_t fm_qcutpos = full_mm.sclen_left + (full_mm.spos - beg_bp);
                                MatchedMate full_mm_left;
                                int fm_ed = split_realignment(fm_qcutpos, beg_bp, end_bp, fullmap_seq, fullmap_seq_len,
                                                              common_tid, full_mm_left, full_mm);
                                if (fm_ed > maxEd)
                                    continue;
                            } else if (full_mm.sclen_left > maxSc)
                                continue;
                        }


                        int ed = split_realignment(qcutpos, beg_bp, end_bp, remain_seq, remain_seq_len, common_tid,
                                                   split_mm_left, split_mm_right);

                        if (ed < best_ed) {
                            esignal = string() + remain_seq[qcutpos - 2] + remain_seq[qcutpos - 1];
                            ssignal = string() + remain_seq[qcutpos] + remain_seq[qcutpos + 1];

                            char near_start_bp[10];
                            genome_seeder.pac2char_otf(beg_bp, 2, near_start_bp);

                            char near_end_bp[10];
                            genome_seeder.pac2char_otf(end_bp - 1, 2, near_end_bp);

                            cr.set_bp(beg_bp, end_bp, ssignal, esignal, near_start_bp, near_end_bp);

                            if (ed == 0)
                                return CR;

                            best_ed = ed;
                        }
                    }
                }
            }

            if (best_ed <= maxEd)
                return CR;

            uint32_t qcutpos = split_mm_left.qepos + split_mm_left.sclen_right;
            uint32_t beg_bp = split_mm_right.spos - split_mm_right.sclen_left;
            uint32_t end_bp = split_mm_left.epos + split_mm_left.sclen_right;

            if (qcutpos < 2 or qcutpos > (remain_seq_len - 2))
                return MCR;

            ssignal = string() + remain_seq[qcutpos - 2] + remain_seq[qcutpos - 1];
            esignal = string() + remain_seq[qcutpos] + remain_seq[qcutpos + 1];

            char near_start_bp[10];
            genome_seeder.pac2char_otf(beg_bp, 2, near_start_bp);

            char near_end_bp[10];
            genome_seeder.pac2char_otf(end_bp - 1, 2, near_end_bp);

            cr.set_bp(beg_bp, end_bp, ssignal, esignal, near_start_bp, near_end_bp);
            if (start_tids.size() > 0 and end_tids.size() > 0)
                return NCR;

            return MCR;
        }
    }
    return rescue_overlapping_bsj(full_mm, split_mm_left, split_mm_right, cr);
    //return UD;
}

int ProcessCirc::split_realignment(uint32_t qcutpos, uint32_t beg_bp, uint32_t end_bp, char *seq, uint32_t seq_len,
                                   const vector <uint32_t> &common_tid, MatchedMate &split_mm_left,
                                   MatchedMate &split_mm_right) {

// 	fprintf(stderr, "Qcutpos: %u, seq_len: %u\n", qcutpos, seq_len);
    if (qcutpos <= 0 or qcutpos >= seq_len)
        return maxEd + 1;

// 	fprintf(stderr, "beg_bp: %u, end_bp: %u\n", beg_bp, end_bp);

    char single_bp[10];
    genome_seeder.pac2char_otf(end_bp, 1, single_bp);
    int last_bp_err = (seq[qcutpos - 1] == single_bp[0]) ? 0 : 1;
// 	fprintf(stderr, "last_on_read: %c, on_ref: %c\n", seq[qcutpos-1], single_bp[0]);

    genome_seeder.pac2char_otf(beg_bp, 1, single_bp);
    int first_bp_err = (seq[qcutpos] == single_bp[0]) ? 0 : 1;
// 	fprintf(stderr, "last_on_read: %c, on_ref: %c\n", seq[qcutpos], single_bp[0]);

    uint32_t lm_pos = end_bp;
    uint32_t rm_pos = beg_bp;

    uint32_t lb = beg_bp;
    uint32_t ub = end_bp;

    AlignRes best_alignment_left(lb);
    AlignRes best_alignment_right(ub);

    extension.set_query_seq(seq);
    extension.set_query_seq_len(seq_len);
    extension.set_query_spos(0);

    bool lok = extension.extend_left(common_tid, seq, lm_pos, qcutpos - 1, maxEd - last_bp_err, lb,
                                     best_alignment_left);

    extension.set_query_spos(qcutpos + 1);
    bool rok = extension.extend_right(common_tid, seq + qcutpos + 1, rm_pos, seq_len - qcutpos - 1,
                                      maxEd - first_bp_err, ub, best_alignment_right);

    best_alignment_left.ed += last_bp_err;
    best_alignment_right.ed += first_bp_err;

// 	fprintf(stderr, "lok: %d, rok: %d\nleft ed: %d\nright ed: %d\n", lok, rok, best_alignment_left.ed, best_alignment_right.ed);


    if (lok and rok and (best_alignment_left.ed + best_alignment_right.ed) <= maxEd)
        return best_alignment_left.ed + best_alignment_right.ed;
    else
        return maxEd + 1;
}

int ProcessCirc::split_realignment(uint32_t qcutpos, MatchedMate &full_mm, MatchedMate &split_mm_left,
                                   MatchedMate &split_mm_right, CircRes &cr) {

    if (qcutpos <= 0 or qcutpos >= fullmap_seq_len)
        return UD;

    qcutpos += full_mm.qspos - 1;

    if (qcutpos <= 0 or qcutpos >= fullmap_seq_len)
        return UD;

    overlap_to_spos(split_mm_left);
    overlap_to_epos(split_mm_left);
    overlap_to_spos(split_mm_right);
    overlap_to_epos(split_mm_right);

    vector <uint32_t> common_tid;
// 	bool success = same_transcript(split_mm_left.exons_spos, split_mm_left.exons_epos, 
// 					split_mm_right.exons_spos, split_mm_right.exons_epos, common_tid);

    vector <MatchedMate> segments;
    segments.push_back(split_mm_left);
    segments.push_back(split_mm_right);
    bool success = same_transcript(segments, 2, common_tid);
    if (!success)
        return UD;

    char single_bp[10];
    genome_seeder.pac2char_otf(split_mm_left.epos, 1, single_bp);
    int last_bp_err = (fullmap_seq[qcutpos - 1] == single_bp[0]) ? 0 : 1;
    genome_seeder.pac2char_otf(split_mm_right.spos, 1, single_bp);
    int first_bp_err = (fullmap_seq[qcutpos] == single_bp[0]) ? 0 : 1;

    uint32_t lm_pos = split_mm_left.epos;
    uint32_t rm_pos = split_mm_right.spos;

    uint32_t lb = split_mm_right.spos;
    uint32_t ub = split_mm_left.epos;

    AlignRes best_alignment_left(lb);
    AlignRes best_alignment_right(ub);

    extension.set_query_seq(fullmap_seq);
    extension.set_query_seq_len(fullmap_seq_len);
    extension.set_query_spos(0);

    bool lok = extension.extend_left(common_tid, fullmap_seq, lm_pos, qcutpos - 1, maxEd - last_bp_err, lb,
                                     best_alignment_left);

    extension.set_query_spos(qcutpos + 1);
    bool rok = extension.extend_right(common_tid, fullmap_seq + qcutpos + 1, rm_pos, fullmap_seq_len - qcutpos - 1,
                                      maxEd - first_bp_err, ub, best_alignment_right);

    best_alignment_left.ed += last_bp_err;
    best_alignment_right.ed += first_bp_err;

    if (!lok or !rok or (best_alignment_left.ed + best_alignment_right.ed) > maxEd)
        return UD;

    MatchedMate new_split_left;
    new_split_left.spos = lm_pos;
    new_split_left.epos = split_mm_left.epos;
    new_split_left.qspos = best_alignment_left.sclen;
    new_split_left.qepos = qcutpos;
    new_split_left.dir = full_mm.dir;
    new_split_left.matched_len = qcutpos - best_alignment_left.sclen;
    new_split_left.sclen_left = best_alignment_left.sclen;
    new_split_left.sclen_right = 0;
    new_split_left.left_ed = best_alignment_left.ed;
    new_split_left.right_ed = 0;
    new_split_left.middle_ed = 0;
    new_split_left.left_ok = true;
    new_split_left.right_ok = true;

    MatchedMate new_split_right;
    new_split_right.spos = split_mm_right.spos;
    new_split_right.epos = rm_pos;
    new_split_right.qspos = qcutpos + 1;
    new_split_right.qepos = fullmap_seq_len - best_alignment_right.sclen;
    new_split_right.dir = full_mm.dir;
    new_split_right.matched_len = fullmap_seq_len - qcutpos - best_alignment_right.sclen;
    new_split_right.sclen_left = 0;
    new_split_right.sclen_right = best_alignment_right.sclen;
    new_split_right.left_ed = 0;
    new_split_right.right_ed = best_alignment_right.ed;
    new_split_right.middle_ed = 0;
    new_split_right.left_ok = true;
    new_split_right.right_ok = true;

    r1_seq = remain_seq;
    r2_seq = fullmap_seq;
    r1_seq_len = remain_seq_len;
    r2_seq_len = fullmap_seq_len;

    return check_split_map(split_mm_right, new_split_right, split_mm_left, new_split_left, cr);
}

int ProcessCirc::rescue_overlapping_bsj(MatchedMate &full_mm, MatchedMate &split_mm_left, MatchedMate &split_mm_right,
                                        CircRes &cr) {
    ConShift con_shift = gtf_parser.get_shift(contigNum, full_mm.spos);

    // start of BP
    if (split_mm_right.spos <= full_mm.epos and split_mm_right.spos > full_mm.spos) {
        // fprintf(stderr, "On the start BP:\n%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t",
        // 			con_shift.contig.c_str(),
        // 			full_mm.spos - con_shift.shift, full_mm.epos - con_shift.shift,
        // 			full_mm.qspos, full_mm.matched_len, full_mm.dir,
        // 			split_mm_right.spos - con_shift.shift, split_mm_right.epos - con_shift.shift,
        // 			split_mm_right.qspos, split_mm_right.matched_len, split_mm_right.dir,
        // 			split_mm_left.spos - con_shift.shift, split_mm_left.epos - con_shift.shift,
        // 			split_mm_left.qspos, split_mm_left.matched_len, split_mm_left.dir);

        get_junctions(full_mm);
// 		full_mm.junc_info.print();

        uint32_t qcutpos = 0;
        for (unsigned int i = 0; i < full_mm.junc_info.count; ++i) {
            if (full_mm.junc_info.junc_info[i].end == split_mm_right.spos)
                qcutpos = full_mm.junc_info.junc_info[i].bp_matched;
        }
        // check for intron retention
        if (qcutpos == 0) {
            qcutpos = split_mm_right.spos - full_mm.spos;
// 			fprintf(stderr, "Rescue intron retention: qcutpos: %d\n", qcutpos);
        }
        if (split_realignment(qcutpos, full_mm, split_mm_left, split_mm_right, cr) == CR)
            return CR;
    }

    // end of BP
    if (split_mm_left.epos >= full_mm.spos and split_mm_left.epos < full_mm.epos) {
        // fprintf(stderr, "On the end BP:\n%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t",
        // 			con_shift.contig.c_str(),
        // 			full_mm.spos - con_shift.shift, full_mm.epos - con_shift.shift,
        // 			full_mm.qspos, full_mm.matched_len, full_mm.dir,
        // 			split_mm_right.spos - con_shift.shift, split_mm_right.epos - con_shift.shift,
        // 			split_mm_right.qspos, split_mm_right.matched_len, split_mm_right.dir,
        // 			split_mm_left.spos - con_shift.shift, split_mm_left.epos - con_shift.shift,
        // 			split_mm_left.qspos, split_mm_left.matched_len, split_mm_left.dir);

        get_junctions(full_mm);
// 		full_mm.junc_info.print();

        uint32_t qcutpos = 0;
        for (unsigned int i = 0; i < full_mm.junc_info.count; ++i) {
            if (full_mm.junc_info.junc_info[i].beg == split_mm_left.epos)
                qcutpos = full_mm.junc_info.junc_info[i].bp_matched;
        }
        // check for intron retention
        if (qcutpos == 0) {
            qcutpos = full_mm.matched_len - (full_mm.epos - split_mm_left.epos);
// 			fprintf(stderr, "Rescue intron retention: qcutpos: %d\n", qcutpos);
        }
        if (split_realignment(qcutpos, full_mm, split_mm_left, split_mm_right, cr) == CR)
            return CR;
    }

    return UD;
}

void ProcessCirc::both_side_consensus(const vector <CircRes> &bsj_reads, string &ssignal, string &esignal) {
    ssignal = "";
    esignal = "";

    vector <string> ss;
    vector <string> es;

    for (unsigned int i = 0; i < bsj_reads.size(); ++i) {
        ss.push_back(bsj_reads[i].start_signal);
        es.push_back(bsj_reads[i].end_signal);
    }

    ssignal = get_consensus(ss);
    esignal = get_consensus(es);
}

void ProcessCirc::report_events(void) {
    // fprintf(stdout, "Hash table pool size: %d\n", ind2ht.size());
    // fprintf(stdout, "gid2ind size: %d\n", gid2ind.size());
    open_report_file();
    if (circ_res.size() <= 0)
        return;

    //for (int i = 0; i < circ_res.size(); ++i) {
    //	fprintf(stdout, "%s\t%u\t%u\t%s", circ_res[i].chr.c_str(), circ_res[i].spos, circ_res[i].epos, circ_res[i].rname.c_str());
    //}

    sort(circ_res.begin(), circ_res.end());
    int cnt = 1;
    CircRes last = circ_res[0];
    vector <CircRes> bsj_reads;
    bsj_reads.push_back(circ_res[0]);
    string ss_con;
    string es_con;
    for (unsigned int i = 1; i < circ_res.size(); ++i) {
        if (circ_res[i] == last) {
            cnt++;
            bsj_reads.push_back(circ_res[i]);
        } else {
            //if (last.type == CR or cnt > 1) {	// won't print novel events with single read support
            if (last.type == CR) {    // won't print novel events
                both_side_consensus(bsj_reads, ss_con, es_con);

                string correct_bp = ((ss_con.compare(last.start_bp_ref) == 0) and
                                     (es_con.compare(last.end_bp_ref) == 0)) ?
                                    "Pass" : "Fail";

                fprintf(report_file, "%s\t%u\t%u\t%d\t%s\t%s-%s\t%s-%s\t%s\t",
                        last.chr.c_str(), last.spos, last.epos, cnt, circ_type[last.type - CR].c_str(),
                        ss_con.c_str(), es_con.c_str(), last.start_bp_ref.c_str(), last.end_bp_ref.c_str(),
                        correct_bp.c_str());

                for (unsigned int j = 0; j < bsj_reads.size() - 1; ++j)
                    fprintf(report_file, "%s,", bsj_reads[j].rname.c_str());
                fprintf(report_file, "%s\n", bsj_reads[bsj_reads.size() - 1].rname.c_str());
            }
            cnt = 1;
            last = circ_res[i];
            bsj_reads.clear();
            bsj_reads.push_back(circ_res[i]);
        }
    }
    //if (last.type == CR or cnt > 1) {
    if (last.type == CR) {
        both_side_consensus(bsj_reads, ss_con, es_con);

        string correct_bp = ((ss_con.compare(last.start_bp_ref) == 0) and (es_con.compare(last.end_bp_ref) == 0)) ?
                            "Pass" : "Fail";

        fprintf(report_file, "%s\t%u\t%u\t%d\t%s\t%s-%s\t%s-%s\t%s\t",
                last.chr.c_str(), last.spos, last.epos, cnt, circ_type[last.type - CR].c_str(),
                ss_con.c_str(), es_con.c_str(), last.start_bp_ref.c_str(), last.end_bp_ref.c_str(), correct_bp.c_str());

        for (unsigned int j = 0; j < bsj_reads.size() - 1; ++j)
            fprintf(report_file, "%s,", bsj_reads[j].rname.c_str());
        fprintf(report_file, "%s\n", bsj_reads[bsj_reads.size() - 1].rname.c_str());
    }
}

void ProcessCirc::open_report_file(void) {
    char *temp_fname = (char *) malloc(FILE_NAME_MAX_LEN);
    char *mode = (char *) malloc(FILE_NAME_MAX_LEN);

    sprintf(temp_fname, "%s.circ_report", outputFilename);
    sprintf(mode, "%c", 'w');

    report_file = open_file(temp_fname, mode);

    free(temp_fname);
    free(mode);
}

void ProcessCirc::open_candid_file(void) {
    char *temp_fname = (char *) malloc(FILE_NAME_MAX_LEN);
    char *mode = (char *) malloc(FILE_NAME_MAX_LEN);

    sprintf(temp_fname, "%s.candidates.pam", outputFilename);
    sprintf(mode, "%c", 'w');

    candid_file = open_file(temp_fname, mode);

    free(temp_fname);
    free(mode);
}

void ProcessCirc::load_genome(void) {
    double cputime_start = get_cpu_time();
    double realtime_start = get_real_time();
    double cputime_curr;
    double realtime_curr;
    double tmpTime;

    Logger::instance().info("Loading genome sequence...\n");

    loadCompressedRefGenome(&tmpTime);            // Reading a fragment

    pre_contig = atoi(getRefGenomeName()) - 1;
    contigNum = pre_contig;

    genome_seeder.pac2char_whole_contig();

    cputime_curr = get_cpu_time();
    realtime_curr = get_real_time();

    Logger::instance().info("Completed! (CPU time: %.2lfs; Real time: %.2lfs)\n",
                            cputime_curr - cputime_start, realtime_curr - realtime_start);

    cputime_start = cputime_curr;
    realtime_start = realtime_curr;
}

void ProcessCirc::print_split_mapping(char *rname, MatchedMate &mm_r1, MatchedMate &mm_r2,
                                      MatchedMate &partial_mm, ConShift &con_shift) {
    fprintf(candid_file, "%s\t%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t",
            rname, con_shift.contig.c_str(),
            partial_mm.spos - con_shift.shift, partial_mm.epos - con_shift.shift,
            partial_mm.qspos, partial_mm.matched_len, partial_mm.dir,
            mm_r1.spos - con_shift.shift, mm_r1.epos - con_shift.shift,
            mm_r1.qspos, mm_r1.matched_len, mm_r1.dir,
            mm_r2.spos - con_shift.shift, mm_r2.epos - con_shift.shift,
            mm_r2.qspos, mm_r2.matched_len, mm_r2.dir);

}

void ProcessCirc::print_split_mapping(char *rname, MatchedMate &mm_r1, MatchedMate &mm_r2,
                                      MatchedMate &r1_partial_mm, MatchedMate &r2_partial_mm, ConShift &con_shift) {
    fprintf(candid_file, "%s\t%s\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\t",
            rname, con_shift.contig.c_str(),
            r1_partial_mm.spos - con_shift.shift, r1_partial_mm.epos - con_shift.shift,
            r1_partial_mm.qspos, r1_partial_mm.matched_len, r1_partial_mm.dir,
            r2_partial_mm.spos - con_shift.shift, r2_partial_mm.epos - con_shift.shift,
            r2_partial_mm.qspos, r2_partial_mm.matched_len, r2_partial_mm.dir,
            mm_r1.spos - con_shift.shift, mm_r1.epos - con_shift.shift,
            mm_r1.qspos, mm_r1.matched_len, mm_r1.dir,
            mm_r2.spos - con_shift.shift, mm_r2.epos - con_shift.shift,
            mm_r2.qspos, mm_r2.matched_len, mm_r2.dir);

}

void set_mm(chain_t &ch, uint32_t qspos, int rlen, int dir, MatchedMate &mm) {
    uint32_t spos = ch.frags[0].rpos;
    uint32_t epos = ch.frags[ch.chain_len - 1].rpos + ch.frags[ch.chain_len - 1].len - 1;
    uint32_t qepos = qspos + rlen - 1;

    //// remove chain into intron on the left
    //int spos_seg_ind;
    //const IntervalInfo<UniqSeg>* it_seg_left = gtf_parser.get_location_overlap_ind(spos, false, spos_seg_ind);
    //if (it_seg_left == NULL) {	// into intron
    //	it_seg_left = gtf_parser.get_interval(spos_seg_ind + 1);
    //	int diff = 0;
    //	if (it_seg_left != NULL)
    //		diff = it_seg_left->spos - spos;
    //	if (diff > 0 and diff < kmer) {
    //		spos = it_seg_left->spos;
    //		qspos += diff;
    //		ch.frags[0].rpos += diff;
    //		ch.frags[0].qpos += diff;
    //		ch.frags[0].len  -= diff;
    //	}
    //}

    //
    //// remove chain into intron on the right
    //int epos_seg_ind;
    //const IntervalInfo<UniqSeg>* it_seg_right = gtf_parser.get_location_overlap_ind(epos, false, epos_seg_ind);
    //if (it_seg_right == NULL) {	// into intron
    //	it_seg_right = gtf_parser.get_interval(epos_seg_ind);
    //	int diff = 0;
    //	if (it_seg_right != NULL)
    //		diff = epos - it_seg_right->epos;
    //	if (diff > 0 and diff < kmer) {
    //		epos = it_seg_right->epos;
    //		qepos -= diff;
    //		ch.frags[ch.chain_len - 1].len -= diff;
    //	}
    //}

    mm.set(spos, epos, qspos, qepos, dir);
}
