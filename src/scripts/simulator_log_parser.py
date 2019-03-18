import sys
import re

events_file = sys.argv[1]
mapping_file = sys.argv[2]
ciri_log_file = sys.argv[3]

ev = {}
with open(events_file) as ef:
	for l in ef:
		ll = l.strip().split()
		ev[(ll[0], ll[1], ll[2])] = ll[3]

#print ev

def find_event(ch, ac):
	for e in ev.keys():
		(c, x, y) = e
		if c == ch and x in ac and y in ac:
			return e
	return None


reads = {}
with open(mapping_file) as mf:
	for l in mf:
		ll = l.strip().split()
		rname = re.split(':|/', ll[0])[1]
		ch = ll[1]
		all_coord = [ ll[2], ll[3], ll[7], ll[8], ll[12], ll[13] ]
		typ = ll[-1]
		evt = find_event(ch, all_coord)
		if typ == "20" and evt != None:
			if evt in reads.keys():
				reads[evt].append(rname)
			else:
				reads[evt] = [rname]

#print '---'
#print reads


read2tid = {}

for e in reads.keys():
	(c, x, y) = e
	#print '{}\t{}\t{}\t{}'.format(c, x, y, len(reads[e]))

iso = ''
gid = ''
tid = ''
cid = 0
rnames = []
split = []
circRNAs = {}
with open(ciri_log_file) as cf:
	for l in cf:
		ll = l.strip().split()
		if l[0] == '!':
			continue
		elif l[0] == 'i':
			iso = ll[0]
		elif ll[0] == '>':
			rnames.append(ll[2])
		elif ll[0] == '**':
			split.append(rnames[-1])
		else:
			circRNAs[cid] = (gid, tid, rnames, split)
			#print gid, tid, iso
			#print rnames
			#print split
			#print '---'
			rnames = []
			split = []
			gid = ll[2]
			tid = ll[3]
			cid += 1
	circRNAs[iso] = (gid, tid, rnames, split)

#print circRNAs

for iso in sorted(circRNAs.keys()):
	(gid, tid, rnames, split) = circRNAs[iso]
	#print gid, tid, iso
	#print rnames
	#print split
	#print '---'

def search_read(r):
	for iso in circRNAs.keys():
		(gid, tid, rnames, split) = circRNAs[iso]
		if r in split:
			return tid
	return None

tids = {}
def classify():
	for e in reads.keys():
		(c, x, y) = e
		real_split = 0
		tids = {}
		for r in reads[e]:
			tid = search_read(r)
			if tid != None:
				real_split += 1
				if tid in tids.keys():
					tids[tid] += 1
				else:
					tids[tid] = 1
		print '{}\t{}\t{}\t{}\t{}\t{}'.format(c, x, y, len(reads[e]), real_split, len(tids.keys())) , 
		sorted_keys = sorted(tids, key=tids.get, reverse=True)
		for t in sorted_keys:
			print '{}\t{}'.format(t, tids[t]) ,
		print 


classify()
