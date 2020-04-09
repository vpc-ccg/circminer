import sys
import re

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
bed_out = sys.argv[2]

bedf = open(bed_out, 'w')

with open(chimeric_sam) as sf:
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

		print >>bedf, '{}\t{}\t{}\t{}'.format(chrom, spos-2, spos+rlen, rname)
