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

chimeric_sam = sys.argv[1]	# sorted
out_file = sys.argv[2]

inv_rname = open(out_file + '.rname', 'w')
inv_sam = open(out_file + '.inv.sam', 'w')
bsj_sam = open(out_file + '.bsj.sam', 'w')

with open(chimeric_sam) as sf:
	mates = { }
	for l in sf:
		if l[0] == '@':
			continue

		ll = l.strip().split()
		rname = ll[0]
		flag = int(ll[1])
		chrom = ll[2]
		spos = int(ll[3])
		cigar = ll[5]
		seq = ll[9]

		rlen = get_len(cigar)
		#print rlen

		if rname not in mates.keys():
			mates[rname] = [ (flag, chrom, seq, l.strip()) ]
		else:
			mates[rname].append( (flag, chrom, seq, l.strip()) )

for rname in mates.keys():
	primaries = []
	for x in mates[rname]:
		if x[0] < 256:
			primaries.append(x)

	inver = True
	for x in mates[rname]:
		if x[0] < 256:
			continue

		for pri in primaries:
			if x[2] == pri[2]:
				inver = False
		
		if inver:
			print 'INV', rname
			print >> inv_rname, rname
			for y in mates[rname]:
				print >> inv_sam, y[3]

		else:
			print 'ORD', rname
			print >> inv_rname, rname
			for y in mates[rname]:
				print >> inv_sam, y[3]

