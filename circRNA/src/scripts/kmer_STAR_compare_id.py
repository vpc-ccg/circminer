#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import os, sys, operator, math, time


def usage():
	print '\nUsage: python ' + sys.argv[0] + ' kmer_fastq STAR_chimera_sam output'
	sys.exit(-1)

####################
def main():
	args = sys.argv[1:]
	if len(args) != 3:
		usage()

	kmer_fname = sys.argv[1]
	STAR_fname = sys.argv[2]

	STAR_id = [ ]
	k_id = [ ]

	sr_kmer = open( sys.argv[1], 'r')
	sr_STAR = open( sys.argv[2], 'r')
	skip = False

	print "Reading " + sr_STAR
	with open(STAR_fname, 'r') as sr_STAR:
	for line in sr_STAR:
		if ( "@" != line[0]):
			STAR_id.add(line.strip().split()[0])

	print "Reading " + 
 	for line in sr_kmer:
		if skip:
			skip = False
			continue
		elif ( "+" == line[0]):
			skip = True
			continue
		elif ( "@" == line[0]):
			continue
		else:
			k_sequences.add(line[:-2])				

	print "Comparing " + sys.argv[2] + " and " + sys.argv[1]

	common_seq = list(STAR_sequences & k_sequences)
	konly_seq = list(k_sequences - STAR_sequences)
	STARonly_seq = list(STAR_sequences - k_sequences)
	common_tot = len(common_seq)
	konly_tot = len(konly_seq)
	STARonly_tot = len(STARonly_seq)

	print "Found " + str(common_tot) + " common sequences"
	print "Found " + str(konly_tot) + " sequences in only " + sys.argv[1]
	print "Found " + str(STARonly_tot) + " sequences in only " + sys.argv[2]

	sw = open( sys.argv[3] + ".intersection", 'w')
	sw.write( "Common Sequences")
	sw.write('\n'.join([str(s) for s in common_seq]))
	sw.close()

	sw = open( sys.argv[3] + ".onlySTAR", 'w')
	sw.write( "Only in " + sys.argv[2] + ":\n")
	sw.write('\n'.join([str(s) for s in konly_seq]))
	sw.close()

	sw = open (sys.argv[3] + ".onlyKmer", 'w')
	sw.write( "Only in " + sys.argv[1] + ":\n")
	sw.write('\n'.join([str(s) for s in STARonly_seq]))
	sw.close()	
	
	#print("--- Writing %s seconds ---" % (time.time() - start_time))		



#############################################################################################
if __name__ == "__main__":
    sys.exit(main())
