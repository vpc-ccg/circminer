import sys

igtf = sys.argv[1]
ogtf = sys.argv[2]

out_gtf = open(ogtf, 'w')

def get_min_max(lines):
	min_st = 1000000000
	max_en = 0
	for l in lines:
		ll = l.strip().split()
		if ll[2] == 'exon':
			min_st = min(min_st, int(ll[3]))
			max_en = max(max_en, int(ll[4]))
	
	return min_st, max_en

def process_trans(lines):
	if len(lines) <= 0:
		return

	min_st, max_en = get_min_max(lines)

	ll = lines[0].strip().split()
	gid = ll[9]
	tid = ll[11]
	print >>out_gtf, '{}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id {} transcript_id {}'.format(ll[0], ll[1], min_st, max_en, ll[6], gid, tid)

def process_gene(lines):
	if len(lines) <= 0:
		return

	min_st, max_en = get_min_max(lines)

	ll = lines[0].strip().split()
	gid = ll[9]
	print >>out_gtf, '{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tgene_id {}'.format(ll[0], ll[1], min_st, max_en, ll[6], gid)

	myl = []
	pre_tid = ''
	exon_cnt = 0
	for l in lines:
		ll = l.strip().split()
		tid = ll[11]

		if tid == pre_tid:
			myl.append(l)
			if ll[2] == 'exon':
				exon_cnt += 1
		else:
			if exon_cnt > 0:
				process_trans(myl)
				strand = myl[0].split()[6]
				if strand == '+':
					for l in myl:
						print >>out_gtf, l
				else:
					for l in reversed(myl):
						print >>out_gtf, l

				myl = []
				
			myl.append(l)
			exon_cnt = 1 if ll[2] == 'exon' else 0

		pre_tid = tid

	if exon_cnt > 0 and len(myl) > 0:
		process_trans(myl)
		strand = myl[0].split()[6]
		if strand == '+':
			for l in myl:
				print >>out_gtf, l
		else:
			for l in reversed(myl):
				print >>out_gtf, l

with open(igtf) as gf:
	lines = []
	gid = ''
	pre_gid = ''
	exon_cnt = 0
	for l in gf:
		ll = l.strip().split()
		gid = ll[9]
		
		if gid == pre_gid:
			lines.append(l.strip())
			if ll[2] == 'exon':
				exon_cnt += 1
		else:
			if exon_cnt > 0:
				process_gene(lines)
				lines = []
			lines.append(l.strip())
			exon_cnt = 1 if ll[2] == 'exon' else 0

		pre_gid = gid

	if exon_cnt > 0 and len(lines) > 0:
		process_gene(lines)


out_gtf.close()
