import sys
import re
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#import numpy as np

def parse_cigar(cigar):
	mysum = 0
	ind = 0
	lens = re.split("M|D|N|=|X|I|S|H|P", cigar)
	if lens[-1] == '':
		del lens[-1]

	types = [ c for c in cigar if not c.isdigit() ]

	#print cigar
	#print types
	#print lens

	count = { }
	for t in "MDN=XISHP":
		count[t] = 0
	for i in range(len(types)):
		count[types[i]] += int(lens[i])

	#print count
	return count['S'], count['I'] + count['D']

#def plot(dist):
#	fig = plt.figure()
#	ax = fig.add_subplot(111, projection='3d')
#	xs = range(76)
#	for d in dist:
#		ax.bar(xs, d)
#	
#	ax.set_xlabel('CS')
#	ax.set_ylabel('ERR')
#	ax.set_zlabel('CNT')
#
#	plt.show()

sam_fname = sys.argv[1]

dist = [[[ 0 for x in range(76) ] for y in range(76) ] for z in range(76) ]

line = 0
mate_sc = 0
mate_err = 0
mate_mm = 0
mate_indel = 0

with open(sam_fname, 'r') as sf:
	for l in sf:
		if l[0] == '@':
			continue

		ll = l.strip().split()
		
		flag = int(ll[1])
		if flag > 255:
			continue

		tlen = int(ll[8])
		if not(tlen > 20000 or tlen < -20000):
		#if tlen > 20000 or tlen < -20000:
			continue
		
		line += 1
		cigar = ll[5]
		if cigar == "*":
			print "bad cigar!!!"
			continue

		nm = ll[15]
		#print cigar, nm
		sc, indel = parse_cigar(cigar)
		mm = int(nm.split(':')[2])
		#print sc, indel, mm
		err = mm + indel

		if line % 2 == 1:
			mate_err = err
			mate_sc = sc
			mate_mm = mm
			mate_indel = indel
		else:
			dist[indel+mate_indel][mm+mate_mm][sc+mate_sc] += 1
			#if (indel + mate_indel == 0) and (mm + mate_mm < 3) and (sc + mate_sc < 3):
			if (indel + mate_indel == 0) and (mm + mate_mm > 10):
				print>> sys.stderr, l.strip()

sums = [ ]
for i in range(76):
	sums.append(sum(sum(dist[i], [])))
print "indel > 0: {}".format( sum(sums[1:]) )

for d in dist[0][:21]:
	for e in d[:21]:
		print e,
	print

print sum(sum(dist[0],[]))
