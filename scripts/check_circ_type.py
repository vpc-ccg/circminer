import sys
import re

igtf = sys.argv[1]
circf = sys.argv[2]

tid2exons = []
tids = []

with open(igtf) as gf:
	pre_tid = ''
	i = -1
	for l in gf:
		ll = l.strip().split()
		if l[0] == '#' or ll[2] != 'exon':
			continue
		ch = ll[0]
		st = ll[3]
		end = ll[4]
		tid = ll[13][1:-2]
		#print ch,st,end,tid

		if tid != pre_tid:
			i += 1
			tid2exons.append([])
			tid2exons[i].append((ch, st, end))
			tids.append(tid)
		else:
			tid2exons[i].append((ch, st, end))

		pre_tid = tid

print >>sys.stderr, 'Finished loading GTF, loaded {} transcripts. {}'.format(len(tid2exons), len(tids))

with open(circf) as cf:
	for l in cf:
		ll = l.strip().split()
		ch = ll[0]
		st = '{}'.format(int(ll[1]) + 1)
		end = ll[2]
		tid = ll[6].split('.')[0]
		
		ind = -1
		for i in range(len(tid2exons)):
			if tids[i] == tid:
				ind = i

		if ind < 0:
			print >>sys.stderr, 'tid {} not found!'.format(tid)

		st_found = False
		en_found = False
		for e in tid2exons[ind]:
			(ech, est, een) = e
			if ech == ch and est == st:
				st_found = True
			if ech == ch and een == end:
				en_found = True

		if st_found and en_found:
			typ = "BothIn"
		elif st_found or en_found:
			typ = "OneEnd"
		else:
			typ = "None"
		
		print '{}\t{}\t{}\t{}'.format(ch, st, end, typ)


