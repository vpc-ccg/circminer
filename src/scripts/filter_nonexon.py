import sys
import re

th = 9

def intersect(chrom, beg, end, echrom, ebeg, eend):
	if chrom != echrom:
		if beg < ebeg:
			return 1
		else:
			return -1
	if end < ebeg:
		return -1
	if beg > eend:
		return 1
	
	if beg < ebeg:
		if ebeg - beg < th:
			return 0
		else:
			return -1
	if end > eend:
		if end - eend < th:
			return 0
		else:
			return 1
	return 0

def check_record(chrom, beg, end, last_chrom, last_beg, last_end, gtff):
	int_res = intersect(chrom, beg, end, last_chrom, last_beg, last_end)
	#print 'LAST check:', int_res, chrom, beg, end, last_chrom, last_beg, last_end
	if int_res == -1:
		return False, (last_chrom, last_beg, last_end)

	if int_res == 0:
		return True, (last_chrom, last_beg, last_end)

	for l in gtff:
		if l[0] == '#':
			continue
		ll = l.strip().split()
		if ll[2] != 'exon':
			continue
	
		echrom = ll[0]
		ebeg = int(ll[3])
		eend = int(ll[4])
	
		int_res = intersect(chrom, beg, end, echrom, ebeg, eend)
		#print 'XXXX:', int_res, chrom, beg, end, echrom, ebeg, eend

		if int_res == -1:
			return False, (echrom, ebeg, eend)

		if int_res == 0:
			return True, (echrom, ebeg, eend)

	return False, ('-1', -1, -1)

def get_len(cigar):
	mysum = 0
	ind = 0
	lens = re.split("M|D|N|=|X|I|S|H|P", cigar)
	if lens[-1] == '':
		del lens[-1]

	#print cigar
	#print lens
	for x in lens:
		ind += len(x)
		if cigar[ind] in 'MDN=X':
			mysum += int(x)
		ind += 1
	return mysum

chimeric_sam = sys.argv[1]
gtf_file = sys.argv[2]
out_file = sys.argv[3]

nonexon_rname = open(out_file + '.rname', 'w')	# sorted
nonexon_sam = open(out_file + '.sam', 'w')	# sorted

with open (gtf_file) as gtff:
	with open(chimeric_sam) as sf:
		last_chrom = '0'
		last_beg = 1000000000
		last_end = 1000000001
		for l in sf:
			if l[0] == '@':
				continue

			ll = l.strip().split()
			rname = ll[0]
			chrom = ll[2]
			spos = int(ll[3])
			cigar = ll[5]
			rlen = get_len(cigar)
			#print rlen

			res, (ch, b, e) = check_record(chrom, spos, spos + rlen - 1, last_chrom, last_beg, last_end, gtff)
			if not res:
				print >> nonexon_rname, rname
				print 'NO', rname
				print >> nonexon_sam, l.strip()

			else:
				print 'YES', rname

			last_chrom = ch
			last_beg = b
			last_end = e
