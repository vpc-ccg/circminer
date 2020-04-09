#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>
#include <set>

#include "chain.h"
#include "gene_annotation.h"

#define REWARD_COEF        2e4
#define PENALTY_COEF    0.1

inline double score_beta(int distr, int distt, int frag_len) {
    int maxd = distr < distt ? distt : distr;
    int mind = distr < distt ? distr : distt;

    return PENALTY_COEF * (maxd - mind);
}

inline double score_alpha(int distr, int distl, int frag_len) {
    return REWARD_COEF * frag_len;
}

bool compare_frag(fragment_t a, fragment_t b) {
    return a.qpos < b.qpos;
}

bool check_junction(uint32_t s1, uint32_t s2, const IntervalInfo<UniqSeg> *ol_exons, int kmer, int read_dist,
                    int &trans_dist) {
    trans_dist = INF;
    if (ol_exons == NULL)
        return false;

    uint32_t e1 = s1 + kmer - 1;
    if (s2 <= e1)    // overlapping kmers
        return false;

    int e12end, beg2s2;
    int trans_dist2intron = -1;
    for (unsigned int i = 0; i < ol_exons->seg_list.size(); i++) {
        e12end = ol_exons->seg_list[i].end - e1;
        beg2s2 = s2 - ol_exons->seg_list[i].next_exon_beg;

        // enforcing 2nd kmer to be completely in the middle of immediate intron
        if (e12end >= 0 and e12end < read_dist and beg2s2 + kmer < 0)
            trans_dist2intron = s2 - e1 - 1;

        if (e12end < 0 or beg2s2 < 0)
            continue;

        trans_dist = e12end + beg2s2;
        if (abs(trans_dist - read_dist) <= maxEd)
            return true;
    }

    //fprintf(stderr, "S1: %d\tS2: %d\tTrans Dist2intron: %d\n", s1, s2, trans_dist2intron);
    if (trans_dist2intron != -1) {
        trans_dist = trans_dist2intron;
        return true;
    } else {
        trans_dist = INF;
        return false;
    }
}

// Assumption: fragment list sorted by reference position
//
// f(i) = max{max_{j>i}{f(j) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
// w_i = kmer size
void chain_seeds_sorted_kbest(int seq_len, GIMatchedKmer *fragment_list, chain_list &best_chain) {
    best_chain.best_chain_count = 0;

    int kmer_cnt = 2 * ceil(1.0 * seq_len / kmer) - 1;
    int max_frag_cnt = seedLim;
    chain_cell dp[kmer_cnt][max_frag_cnt + 1];

    uint32_t max_best = maxChainLen;
    uint32_t best_count = 0;
    double best_score = -1;

    int distr, distt;
    int genome_dist;
    int trans_dist;
    int read_dist;

    double temp_score;

    int ii, jj;
    uint32_t i, j;
    uint32_t lb_ind[kmer_cnt];
    uint32_t max_lpos_lim;
    uint32_t read_remain;
    uint32_t seg_start;
    uint32_t seg_end;
    uint32_t max_exon_end;
    GIMatchedKmer *cur_mk;
    GIMatchedKmer *pc_mk;    // previously calculated matched kmer
    const IntervalInfo<UniqSeg> *ol_exons;

    map<double, chain_cell_list> score2chain;

    chain_cell tmp_cell;
    chain_cell_list empty_chain_cell_list;
    empty_chain_cell_list.count = 0;

    // Ignore empty fragment list at the back
    //fprintf(stderr, "Got here before frag list in chaining\n");
    //fprintf(stderr, "frag list kmer count: %d\n", (fragment_list + kmer_cnt - 1)->frag_count);
    while ((kmer_cnt >= 1) and (fragment_list + kmer_cnt - 1)->frag_count <= 0)
        kmer_cnt--;

    if (kmer_cnt <= 0)
        return;

    // Initialize dp array
    for (ii = kmer_cnt - 1; ii >= 0; ii--) {
        cur_mk = fragment_list + ii;
        for (i = 0; i < cur_mk->frag_count; i++) {
            dp[ii][i].score = kmer;
            dp[ii][i].prev_list = -1;
            dp[ii][i].prev_ind = -1;
        }
    }

    // Updating stage
    for (ii = kmer_cnt - 2; ii >= 0; ii--) {
        cur_mk = fragment_list + ii;
        read_remain = seq_len - cur_mk->qpos - kmer;
        //memset(lb_ind, 0, kmer_cnt * sizeof(int));
        for (int k = 0; k < kmer_cnt; k++)
            lb_ind[k] = 0;

        for (i = 0; i < cur_mk->frag_count; i++) {
            seg_start = cur_mk->frags[i].info;
            seg_end = cur_mk->frags[i].info + kmer - 1;

            //max_lpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end);	// = 0 means not found
            max_lpos_lim = MAXUB;

            for (jj = ii + 1; jj < kmer_cnt; jj++) {
                pc_mk = fragment_list + jj;
                if (pc_mk->frag_count <= 0 or lb_ind[jj] >= pc_mk->frag_count)    // no fragment left in pc_mk to chain
                    continue;

                // will not chain to any fragment in this list
                if (cur_mk->frags[i].info + maxIntronLen < pc_mk->frags[lb_ind[jj]].info)
                    continue;

                // skip fragments starting before target
                while ((lb_ind[jj] < pc_mk->frag_count) and (pc_mk->frags[lb_ind[jj]].info <= cur_mk->frags[i].info))
                {
                    lb_ind[jj]++;
                }

                // no more fragment from pre are in range for the remaining of the t
                if (lb_ind[jj] >= pc_mk->frag_count)
                    continue;

                if (max_lpos_lim == MAXUB) {
                    // = 0 means not found
                    max_lpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end, ol_exons);
                    //vafprintf(2, stderr, "[%d-%d] -> upper bound: %u\n", seg_start, seg_end, max_lpos_lim);
                }

                distr = pc_mk->qpos - cur_mk->qpos - kmer;
                read_dist = distr;

                j = lb_ind[jj];
                while ((j < pc_mk->frag_count) and (pc_mk->frags[j].info <= max_lpos_lim)) {
                    if (max_exon_end == 0 or (pc_mk->frags[j].info + kmer - 1) <= max_exon_end) {
                        // allowed to put on genome
                        genome_dist = pc_mk->frags[j].info - seg_end - 1;
                    } else {
                        genome_dist = INF;
                    }

                    if (abs(genome_dist - read_dist) <= maxEd) {
                        distt = genome_dist;
                    } else if (check_junction(seg_start, pc_mk->frags[j].info, ol_exons, kmer, read_dist, trans_dist)) {
                        distt = trans_dist;
                    } else {
                        j++;
                        continue;
                    }

                    temp_score = dp[jj][j].score + score_alpha(distr, distt, kmer) - score_beta(distr, distt, kmer);

                    if (temp_score > dp[ii][i].score) {
                        dp[ii][i].score = temp_score;
                        dp[ii][i].prev_list = jj;
                        dp[ii][i].prev_ind = j;

                        if (score2chain.find(temp_score) == score2chain.end()) {
                            score2chain[temp_score] = empty_chain_cell_list;
                        }

                        if (score2chain[temp_score].count < max_best) {
                            tmp_cell.score = dp[ii][i].score;
                            tmp_cell.prev_list = ii;
                            tmp_cell.prev_ind = i;
                            score2chain[temp_score].chain_list[score2chain[temp_score].count++] = tmp_cell;
                        }

                    }
                    j++;
                }
            }
        }
    }

    // Finding best score

    // chain_cell best_indices[max_best];
    //for (ii = kmer_cnt - 1; ii >= 0; ii--) {
    //	cur_mk = fragment_list + ii;
    //	for (i = 0; i < cur_mk->frag_count; i++) {
    //		if (dp[ii][i].score > best_score) {
    //			best_score = dp[ii][i].score;
    //			best_count = 1;
    //			best_indices[0].prev_list = ii;
    //			best_indices[0].prev_ind = i;
    //		}
    //		else if (dp[ii][i].score == best_score and best_count < max_best) {
    //			best_indices[best_count].prev_list = ii;
    //			best_indices[best_count].prev_ind = i;
    //			best_count++;
    //		}
    //	}
    //}

    // Back-tracking
    chain_cell best_index;
    int tmp_list_ind;

    //best_chain.best_chain_count = best_count;
    //for (j = 0; j < best_count; j++) {
    //	best_index = best_indices[j];

    uint32_t spos;
    best_count = 0;
    best_score = (score2chain.rbegin() != score2chain.rend()) ? score2chain.rbegin()->first : kmer;

    set <uint32_t> repeats;
    map<double, chain_cell_list>::reverse_iterator it;

    for (it = score2chain.rbegin(); it != score2chain.rend(); ++it) {
        for (uint32_t l = 0; l < it->second.count; l++) {
            if (best_count >= max_best)
                break;

            best_index = it->second.chain_list[l];
            spos = (fragment_list + best_index.prev_list)->frags[best_index.prev_ind].info;
            if (best_index.score < best_score and repeats.find(spos) != repeats.end())
                continue;

            i = 0;
            j = best_count++;
            while (best_index.prev_list != -1) {
                best_chain.chains[j].frags[i].rpos = (fragment_list +
                                                      best_index.prev_list)->frags[best_index.prev_ind].info;
                best_chain.chains[j].frags[i].qpos = (fragment_list + best_index.prev_list)->qpos;
                best_chain.chains[j].frags[i].len = kmer;

                if (i != 0) {
                    repeats.insert(best_chain.chains[j].frags[i].rpos);
                }

                tmp_list_ind = best_index.prev_list;
                best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
                best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;

                i++;
            }

            best_chain.chains[j].score = best_index.score;
            best_chain.chains[j].chain_len = i;
        }
    }

    // add chains of size 1
    if (best_count == 0) {
        for (ii = kmer_cnt - 1; ii >= 0; ii--) {
            cur_mk = fragment_list + ii;
            for (i = 0; i < cur_mk->frag_count; i++) {
                if (best_count >= max_best)
                    break;
                j = best_count++;
                best_chain.chains[j].frags[0].rpos = cur_mk->frags[i].info;
                best_chain.chains[j].frags[0].qpos = cur_mk->qpos;
                best_chain.chains[j].frags[0].len = kmer;
                best_chain.chains[j].score = dp[ii][i].score;
                best_chain.chains[j].chain_len = 1;
            }
        }
    }

    best_chain.best_chain_count = best_count;
}

// Assumption: fragment list sorted by reference position
//
// f(i) = max{max_{j>i}{f(j) + a(i, j) - b(i, j)}, w_i}
// a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
// b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
// b(i, j) = gap_cost
// w_i = kmer size
void chain_seeds_sorted_kbest2(int seq_len, GIMatchedKmer *fragment_list, chain_list &best_chain,
                               int kmer, int kmer_cnt, int shift) {
    best_chain.best_chain_count = 0;

    if (kmer_cnt <= 0)
        return;

    int max_frag_cnt = seedLim;
    chain_cell dp[kmer_cnt][max_frag_cnt + 1];

    uint32_t max_best = maxChainLen;
    uint32_t best_count = 0;
    double best_score = -1;

    int distr, distt;
    int genome_dist;
    int trans_dist;
    int read_dist;

    double temp_score;

    int ii, jj;
    uint32_t i, j;
    uint32_t lb_ind[kmer_cnt];
    uint32_t max_lpos_lim;
    uint32_t read_remain;
    uint32_t seg_start;
    uint32_t seg_end;
    uint32_t max_exon_end;
    GIMatchedKmer *cur_mk;
    GIMatchedKmer *pc_mk;    // previously calculated matched kmer
    const IntervalInfo<UniqSeg> *ol_exons;

    map<double, chain_cell_list> score2chain;

    chain_cell tmp_cell;
    chain_cell_list empty_chain_cell_list;
    empty_chain_cell_list.count = 0;

    // Ignore empty fragment list at the back
    while ((kmer_cnt >= 1) and (fragment_list + kmer_cnt - 1)->frag_count <= 0)
        kmer_cnt--;

    if (kmer_cnt <= 0)
        return;

    // Initialize dp array
    for (ii = kmer_cnt - 1; ii >= 0; ii--) {
        cur_mk = fragment_list + ii;
        for (i = 0; i < cur_mk->frag_count; i++) {
            dp[ii][i].score = kmer;
            dp[ii][i].prev_list = -1;
            dp[ii][i].prev_ind = -1;
        }
    }

    // Updating stage
    for (ii = kmer_cnt - 2; ii >= 0; ii--) {
        cur_mk = fragment_list + ii;
        read_remain = seq_len - cur_mk->qpos - kmer;
        //memset(lb_ind, 0, kmer_cnt * sizeof(int));
        for (int k = 0; k < kmer_cnt; k++)
            lb_ind[k] = 0;

        for (i = 0; i < cur_mk->frag_count; i++) {
            seg_start = cur_mk->frags[i].info;
            seg_end = cur_mk->frags[i].info + kmer - 1;

            //max_lpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end);	// = 0 means not found
            max_lpos_lim = MAXUB;

            for (jj = ii + 1; jj < kmer_cnt; jj++) {
                pc_mk = fragment_list + jj;
                if (pc_mk->frag_count <= 0 or lb_ind[jj] >= pc_mk->frag_count)    // no fragment left in pc_mk to chain
                    continue;

                // will not chain to any fragment in this list
                if (cur_mk->frags[i].info + maxIntronLen < pc_mk->frags[lb_ind[jj]].info)
                    continue;

                // skip fragments starting before target
                while ((lb_ind[jj] < pc_mk->frag_count) and (pc_mk->frags[lb_ind[jj]].info <= cur_mk->frags[i].info))
                {
                    lb_ind[jj]++;
                }

                // no more fragment from pre are in range for the remaining of the t
                if (lb_ind[jj] >= pc_mk->frag_count)
                    continue;

                if (max_lpos_lim == MAXUB) {
                    // = 0 means not found
                    max_lpos_lim = gtf_parser.get_upper_bound(seg_start, kmer, read_remain, max_exon_end, ol_exons);
                    //vafprintf(2, stderr, "[%d-%d] -> upper bound: %u\n", seg_start, seg_end, max_lpos_lim);
                }

                distr = pc_mk->qpos - cur_mk->qpos - kmer;
                read_dist = distr;

                j = lb_ind[jj];
                while ((j < pc_mk->frag_count) and (pc_mk->frags[j].info <= max_lpos_lim)) {
                    if (max_exon_end == 0 or (pc_mk->frags[j].info + kmer - 1) <= max_exon_end)
                        // allowed to put on genome
                        genome_dist = pc_mk->frags[j].info - seg_end - 1;
                    else
                        genome_dist = INF;

                    //fprintf(stderr, "Genome dist: %d\nRead dist: %d\n\n", genome_dist, read_dist);
                    if (abs(genome_dist - read_dist) <= maxEd) {
                        distt = genome_dist;
                    } else if (check_junction(seg_start, pc_mk->frags[j].info, ol_exons, kmer, read_dist, trans_dist)) {
                        distt = trans_dist;
                    } else {
                        j++;
                        continue;
                    }

                    temp_score = dp[jj][j].score + score_alpha(distr, distt, kmer) - score_beta(distr, distt, kmer);

                    if (temp_score > dp[ii][i].score) {
                        dp[ii][i].score = temp_score;
                        dp[ii][i].prev_list = jj;
                        dp[ii][i].prev_ind = j;

                        if (score2chain.find(temp_score) == score2chain.end()) {
                            score2chain[temp_score] = empty_chain_cell_list;
                        }

                        if (score2chain[temp_score].count < max_best) {
                            tmp_cell.score = dp[ii][i].score;
                            tmp_cell.prev_list = ii;
                            tmp_cell.prev_ind = i;
                            score2chain[temp_score].chain_list[score2chain[temp_score].count++] = tmp_cell;
                        }

                    }
                    j++;
                }
            }
        }
    }

    // Finding best score

    // chain_cell best_indices[max_best];
    //for (ii = kmer_cnt - 1; ii >= 0; ii--) {
    //	cur_mk = fragment_list + ii;
    //	for (i = 0; i < cur_mk->frag_count; i++) {
    //		if (dp[ii][i].score > best_score) {
    //			best_score = dp[ii][i].score;
    //			best_count = 1;
    //			best_indices[0].prev_list = ii;
    //			best_indices[0].prev_ind = i;
    //		}
    //		else if (dp[ii][i].score == best_score and best_count < max_best) {
    //			best_indices[best_count].prev_list = ii;
    //			best_indices[best_count].prev_ind = i;
    //			best_count++;
    //		}
    //	}
    //}

    // Back-tracking
    chain_cell best_index;
    int tmp_list_ind;

    //best_chain.best_chain_count = best_count;
    //for (j = 0; j < best_count; j++) {
    //	best_index = best_indices[j];

    uint32_t spos;
    best_count = 0;
    best_score = (score2chain.rbegin() != score2chain.rend()) ? score2chain.rbegin()->first : kmer;

    set <uint32_t> repeats;
    map<double, chain_cell_list>::reverse_iterator it;

    for (it = score2chain.rbegin(); it != score2chain.rend(); ++it) {
        for (uint32_t l = 0; l < it->second.count; l++) {
            if (best_count >= max_best)
                break;

            best_index = it->second.chain_list[l];
            spos = (fragment_list + best_index.prev_list)->frags[best_index.prev_ind].info;
            if (best_index.score < best_score and repeats.find(spos) != repeats.end())
                continue;

            i = 0;
            j = best_count++;
            while (best_index.prev_list != -1) {
                best_chain.chains[j].frags[i].rpos =
                        shift + (fragment_list + best_index.prev_list)->frags[best_index.prev_ind].info;
                best_chain.chains[j].frags[i].qpos = (fragment_list + best_index.prev_list)->qpos;
                best_chain.chains[j].frags[i].len = kmer;

                if (i != 0) {
                    repeats.insert(best_chain.chains[j].frags[i].rpos);
                }

                tmp_list_ind = best_index.prev_list;
                best_index.prev_list = dp[best_index.prev_list][best_index.prev_ind].prev_list;
                best_index.prev_ind = dp[tmp_list_ind][best_index.prev_ind].prev_ind;

                i++;
            }

            best_chain.chains[j].score = best_index.score;
            best_chain.chains[j].chain_len = i;
        }
    }

    // add chains of size 1
    if (best_count == 0) {
        for (ii = kmer_cnt - 1; ii >= 0; ii--) {
            cur_mk = fragment_list + ii;
            for (i = 0; i < cur_mk->frag_count; i++) {
                if (best_count >= max_best)
                    break;
                j = best_count++;
                best_chain.chains[j].frags[0].rpos = shift + cur_mk->frags[i].info;
                best_chain.chains[j].frags[0].qpos = cur_mk->qpos;
                best_chain.chains[j].frags[0].len = kmer;
                best_chain.chains[j].score = dp[ii][i].score;
                best_chain.chains[j].chain_len = 1;
            }
        }
    }

    best_chain.best_chain_count = best_count;
}
