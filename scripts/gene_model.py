#!/usr/bin/env python
#-*- encoding:utf-8 -*-
#


import sys, os, math
from operator import itemgetter

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def usage():
	print '\nUsage: python ' + sys.argv[0] + ' GTF'
	sys.exit(-1)

####################
def overlap_interval( s1, e1, s2, e2):
	ol = 0
	if ( s1 <= s2 and e2 <= s2 ):
		ol = e2 - s2 + 1
	elif ( s2 <= s1 and e1 <= e2 ):
		ol = e1 - s1 + 1
	elif ( s2 <= s1 and s1 <= e2 ):
		ol = e2 - s1 + 1
	elif ( s1 <= s2 and s2 <= e1 ):
		ol = e1 - s2 + 1
	return ol

########## Load GTF with versioon >=75 ##########
##### 
# chr -> list of sorted genes: start, end, list_of_isoforms
# an isoform is a list of exons
def load_GTF( gtf_file ):
	gene_dict = {} # from chromosome to a list of sorted gene record
	gene_list = [] # sorted exon record for a gene

	# tmp struc for exon for a gene
	coor_dict = {}
	coor_list = []
	iso_list  = [] # for promotor
	
	gene_name = ""
	ref = ""
	start = 0
	end = 0
	
	sr = open( gtf_file, 'r' )
	for line in sr:
		if ("#" == line[0]):	
			continue
		z_list = line.split("\t")
		# ignoring all pacthed records 
		if (  2 < len(z_list[0])):
			continue

		if ( ref != z_list[0] ):	# storing the list of all genes into a ref entry
			if 0 < len(coor_list):
 				iso_list.append( [s, e, coor_list ])
			if 0 < len(iso_list):
				iso_list.sort( key = itemgetter(0,1) )
				gene_list.append( [ ref, start, end, gene_name, iso_list ])
			if ( 0 < len(gene_list) ):
				gene_list.sort( key=itemgetter(1,2) )
				gene_dict[ ref ] = gene_list
				print "Loading Chromosome " + ref + " with " + str(len(gene_list))  + " genes"
			
			gene_list = []
			iso_list  = []
			coor_list = []
			start = 0
			end = 0
			ref = z_list[0]
			print "Starting Chromosome " +ref

		if ("gene" == z_list[2]):
			if 0 < len(coor_list):
				iso_list.append( [ s, e, coor_list ] )
			if 0 < len(iso_list):
				iso_list.sort( key = itemgetter(0,1) )
				gene_list.append( [ ref, start, end, gene_name, iso_list ])
			coor_list = []
			iso_list  = []
			coor_dict = {}

			gene_name = line.split("\"")[5]
			start = int(z_list[3])
			end   = int(z_list[4])
		
		if ("transcript" == z_list[2]):
			if 0 < len(coor_list):
				iso_list.append( [ s, e, coor_list ] )
			s = int(z_list[3])
			e = int(z_list[4])
			coor_list = []
		elif ("exon" == z_list[2]):
			coor_list.append( [ int(z_list[3]), int(z_list[4]) ] )	

	sr.close()

	# Last exon
	if 0 < len(coor_list):
		iso_list.append( [ s, e, coor_list ] )
	if 0 < len(iso_list):
		iso_list.sort( key = itemgetter(0,1) )
		#iso_list.sort(key = itemgetter(1,2) )
		gene_list.append( [ ref, start, end, gene_name, iso_list ])
	if ( 0 < len(gene_list) ):
		gene_list.sort( key=itemgetter(1,2) )
		gene_dict[ ref ] = gene_list
		print "Loading Chromosome " + ref + " with " + str(len(gene_list))  + " genes"

	return gene_dict
#########################	
def main():
	args = sys.argv[1:]
	if len(args) != 1:
		usage()

	gene_model = load_GTF( sys.argv[1] )
		
#########################	
if __name__ == "__main__":
	sys.exit(main())
