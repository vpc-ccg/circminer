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

dist = [[ 0 for x in range(76) ] for y in range(76) ]

line = 0
mate_sc = 0
mate_err = 0
with open(sam_fname, 'r') as sf:
	for l in sf:
		if l[0] == '@':
			continue

		line += 1
		ll = l.strip().split()
		cigar = ll[5]
		if cigar == "*":
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
		else:
			dist[sc+mate_sc][err+mate_err] += 1
			if sc + mate_sc + err + mate_err == 0:
				print l.strip()

#for d in dist[:17]:
#	for e in d[:17]:
#		print e,
#	print

#print sum(sum(dist,[]))
