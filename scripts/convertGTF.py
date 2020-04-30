import sys

def get_min_max(lines):
    min_st = 1000000000
    max_en = 0
    for l in lines:
        ll = l.strip().split()
        if ll[2] == 'exon':
            min_st = min(min_st, int(ll[3]))
            max_en = max(max_en, int(ll[4]))

    return min_st, max_en

def get_gid(l):
    return l[l.index('gene_id') + 1]

def get_tid(l):
    return l[l.index('transcript_id') + 1]

def process_trans(lines, out_gtf):
    if len(lines) <= 0:
        return

    ll = lines[0].strip().split()
    rec_type = ll[2]
    # transcript record is already there
    if rec_type == 'transcript':
        out_gtf.write('{}\n'.format(lines[0]))
        del lines[0]
        return

    min_st, max_en = get_min_max(lines)

    gid = get_gid(ll)
    tid = get_tid(ll)
    out_gtf.write('{}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id {} transcript_id {}\n'.format(ll[0], ll[1], min_st, max_en, ll[6], gid, tid))

def process_gene(lines, out_gtf):
    if len(lines) <= 0:
        return

    min_st, max_en = get_min_max(lines)

    ll = lines[0].strip().split()
    gid = get_gid(ll)

    rec_type = ll[2]
    # gene record is not there
    if rec_type != 'gene':
        out_gtf.write('{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tgene_id {}\n'.format(ll[0], ll[1], min_st, max_en, ll[6], gid))

    myl = []
    pre_tid = ''
    exon_cnt = 0
    for l in lines:
        ll = l.strip().split()
        tid = get_tid(ll)

        if tid == pre_tid:
            myl.append(l)
            if ll[2] == 'exon':
                exon_cnt += 1
        else:
            if exon_cnt > 0:
                process_trans(myl, out_gtf)
                strand = myl[0].split()[6]
                if strand == '+':
                    for e in myl:
                        out_gtf.write('{}\n'.format(e))
                else:
                    for e in reversed(myl):
                        out_gtf.write('{}\n'.format(e))

            myl = []
            myl.append(l)
            exon_cnt = 1 if ll[2] == 'exon' else 0

        pre_tid = tid

    if exon_cnt > 0 and len(myl) > 0:
        process_trans(myl, out_gtf)
        strand = myl[0].split()[6]
        if strand == '+':
            for e in myl:
                out_gtf.write('{}\n'.format(e))
        else:
            for e in reversed(myl):
                out_gtf.write('{}\n'.format(e))

def usage():
    print('Usage: python {} INPUT_GTF OUTPUT_GTF'.format(sys.argv[0]))

def main():
    args = sys.argv[1:]
    if len(args) != 2:
        usage()
        exit(1)
    
    igtf = args[0]
    ogtf = args[1]

    out_gtf = open(ogtf, 'w')

    with open(igtf) as gf:
        lines = []
        gid = ''
        pre_gid = ''
        exon_cnt = 0
        for l in gf:
            ll = l.strip().split()
            gid = get_gid(ll)

            if gid == pre_gid:
                lines.append(l.strip())
                if ll[2] == 'exon':
                    exon_cnt += 1
            else:
                if exon_cnt > 0:
                    process_gene(lines, out_gtf)
                    lines = []
                lines.append(l.strip())
                exon_cnt = 1 if ll[2] == 'exon' else 0

            pre_gid = gid

        if exon_cnt > 0 and len(lines) > 0:
            process_gene(lines, out_gtf)


    out_gtf.close()

if __name__ == '__main__':
    main()
