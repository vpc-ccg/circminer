import sys

rname_fname = sys.argv[1]

rnames = [ ]
with open(rname_fname) as rf:
	for l in rf:
		rnames.append(l.strip())

def bin_search(mlist, beg, end, target):
	if end - beg <= 0:
		return False
	mid = (beg + end) / 2
	if target < mlist[mid]:
		return bin_search(mlist, beg, mid, target)
	elif target > mlist[mid]:
		return bin_search(mlist, mid+1, end, target)
	return True

for l in sys.stdin:
	if l[0] == '@':
		continue
	l = l.strip()
	ll = l.split()
	rname = ll[0].split('.')[1]
	if bin_search(rnames, 0, len(rnames), rname):
		print l
