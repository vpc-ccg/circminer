import sys

sam_fname = sys.argv[1]
myout_fname = sys.argv[2]

SAMSIZE = 40000000
LOCDIFTH = 10

sam_locs = [ 0, 0 ]
sam_locs[0] = [ ]
sam_locs[1] = [ ]

exact_locs = [ ]
scat = [ ]

for i in range(SAMSIZE):
	sam_locs[0].append([])
	sam_locs[1].append([])
	exact_locs.append(-1)
	scat.append(0)

is_r1 = 0
read = 0

print "Loading sam file..."
with open(sam_fname) as sf:
	for l in sf:
		ll = l.strip().split()
		rname = int(ll[0].split('.')[1])
		chrom = ll[2]
		if chrom == "*":
			continue
		
		read += 1
		if (read % 200000) == 0:
			print "# {} lines".format(read)
		flag = int(ll[1])
		spos = int(ll[3])
		
		cat = ll[-1]
		scat[rname] = int(cat[-1])

		is_r1 = (flag & 64) / 64
		sam_locs[is_r1][rname].append( (chrom, spos) )

#print "R1: ", sam_locs[1]
#print "R2: ", sam_locs[0]

print "Size: ", len(sam_locs[0]) - sam_locs[0].count([]), len(sam_locs[1]) - sam_locs[1].count([])

def closest_vacinity(chrom, p, l):
	return min( [ abs (x - p) for (ch, x) in l if ch == chrom ] )

# score of each mate:
# 4: exact match
# 1: in vacinity
# 0: not matched
def set_score(chrom, p, l):
	cl = 2 * LOCDIFTH
	
	try:
		cl = closest_vacinity(chrom, p, l)
	except ValueError:
		cl = LOCDIFTH + 1

	if cl == 0:		#Exact match
		return 4
	elif cl <= LOCDIFTH:
		return 1
	else:
		return 0

def check_loc(rname, chrom, p1, p2):
	s1 = 0
	s2 = 0
	
	s1 += set_score(chrom, p1, sam_locs[1][rname])
	s1 += set_score(chrom, p2, sam_locs[0][rname])

	s2 += set_score(chrom, p2, sam_locs[1][rname])
	s2 += set_score(chrom, p1, sam_locs[0][rname])

	return max(s1, s2)

read = 0
print "Loading my output file..."
with open(myout_fname) as mf:
	for l in mf:
		read += 1
		if (read % 200000) == 0:
			print "# {} reads".format(read)

		ll = l.strip().split()
		rname = int(ll[0].split('.')[1])
		chrom = ll[1]
		pos1 = int(ll[2])
		pos2 = int(ll[10])
		mtype = int(ll[-1])
		
		exact_locs[rname] = max(exact_locs[rname], check_loc(rname, chrom, pos1, pos2))

print "Printing results..."
for rname in range(len(exact_locs)):
	if rname % 100000 == 0:
		print "# {} reads".format(rname)
	if exact_locs[rname] != -1:
		print>> sys.stderr, rname, exact_locs[rname], "XC:i:{}".format(scat[rname])
