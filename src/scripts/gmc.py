import sys

gtf_fname = sys.argv[1]
junc_fname = sys.argv[2]

th = 5
exons = [ ]
starts = { }
ends = { }

def valid_start_bsj(exons, chrom, start, end):
	for exon in exons:
		if chrom != exon[0]:
			continue
		if abs(exon[1] - start) <= th:
			return True
	return False

def valid_end_bsj(exons, chrom, start, end):
	for exon in exons:
		if chrom != exon[0]:
			continue
		if abs(exon[2] - end) <= th:
			return True
	return False

def find_target(mlist, beg, end, target):
	if end - beg <= 0:
		return False
	mid = (beg + end) / 2
	if target > mlist[mid] + th:
		return find_target(mlist, mid+1, end, target)
	if target < mlist[mid] - th:
		return find_target(mlist, beg, mid, target)

	# abs(mlist[mid] - target) < th
	return True

with open(gtf_fname) as gf:
	for l in gf:
		if l[0] == '#':
			continue
		
		l = l.strip().split()
		type = l[2]
		if type != 'exon':
			continue

		chrom = l[0]
		if chrom not in starts.keys():
			starts[chrom] = []
		if chrom not in ends.keys():
			ends[chrom] = []

		start = int(l[3])
		end = int(l[4])
		exons.append((chrom, start, end))
		starts[chrom].append(start)
		ends[chrom].append(end)

#print exons
for ch in starts.keys():
	starts[ch].sort()
	ends[ch].sort()

with open(junc_fname) as jf:
	line = 0
	for ll in jf:
		if line % 1000 == 0:
			print >>sys.stderr, line
		line += 1
		ll = ll.strip()
		l = ll.split()
		if len(l) < 9:
			continue
		chrom = l[1]
		start = int(l[2])
		end = int(l[3])
		if not find_target(ends[chrom], 0, len(ends[chrom]), end):
			continue

		chrom = l[5]
		start = int(l[6])
		end = int(l[7])
		#if valid_start_bsj(exons, chrom, start, end):
		if find_target(starts[chrom], 0, len(starts[chrom]), start):
			print ll

