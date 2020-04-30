#!/usr/bin/env python
#-*- encoding:utf-8 -*-

import sys, os, math
from operator import itemgetter
import gene_model as GM

diff_th = 10

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def usage():
	print '\nUsage: python ' + sys.argv[0] + ' GTF STAR_Mappings all.junc Prefix'
	sys.exit(-1)

####################
def overlap_interval( s1, e1, s2, e2):
	ol = 0
	if ( s1 <= s2 and s2 <= e1 ):
		ol = e1 - s2 + 1
		if e2 < e1:
			ol = e2 - s2 + 1
	elif ( s2 <= s1 and s1 <= e2 ):
		ol = e2 - s1 + 1
		if e1 < e2:
			ol = e1 - s1 + 1
	return ol
			
# Assumption: READ is encoded as LIBRARY.ID
#########################
def get_readid( raw_str ):
	return int( raw_str.split(".")[1] )
# rid -> [ class, line ]
#########################
def load_Prediction( inp ):
	read_dict = {}
	sr = open(inp, 'r')
	for line in sr:
		t_list = line.split("\t")

		rid = get_readid( t_list[0] ) 
		z_class = int(t_list[-1].strip() )
		
		if rid not in read_dict:
			read_dict[ rid ] = [ z_class, line ]
		elif z_class < read_dict[rid ][0]:
			read_dict[ rid ] = [ z_class, line ]
	sr.close()
	return read_dict

#########################	
def parse_cigar( pos, cigar ):
	tmp  = 0
	shift= 0
	ml   = 0
	s    = 0
	e    = 0
	z_list = []

	for ch in cigar:
		if ('N' == ch):
			e = s + shift - 1
			ml += shift
			z_list.append( [ pos+s, pos+e] )
			s += shift
			s += tmp
			shift = 0
			tmp   = 0
		elif ('S' == ch or 'H' == ch):
			if 0 < shift:#tmp:
				#shift += tmp
				e = s + shift - 1
				ml += shift
				z_list.append( [ pos+s, pos+e] )
				s  += shift
			shift = 0
			tmp   = 0
		elif ('I' == ch):
			tmp = 0
		elif( ch.isalpha() ):
			shift += tmp
			tmp = 0
		else:
			tmp = tmp*10 + ord(ch)-ord('0')
	if 0 < shift:
		e = s + shift - 1
		z_list.append( [ pos+s, pos+e] )
		ml += shift
	
	return pos+e, ml, z_list
#########################	
def parse_mapping( loci_list ):
	ref1, ref2 = "", ""
	s1, e1, l1, s2, e2, l2 = 0, 0, 0, 0, 0, 0
	z1 = []
	z2 = []
	for item in loci_list:
		if 256 < item[3]:
			continue
		#if item[3] & 64:
		if 0 < item[4]:# & 64:
			ref1 = item[0]
			s1 = item[1]
			e1, l1, z1 = parse_cigar( item[1], item[2])
		else:
			ref2 = item[0]
			s2 = item[1]
			e2, l2, z2 = parse_cigar( item[1], item[2])
	return ref1, s1, e1, l1, z1, ref2, s2, e2, l2, z2

#########################	
def compare_isoform( seg, exon_list ):
	best_iso = -1
	idx = 0 
	match = []
	for z in seg:
		match.append(0)
		
	exon_num = len( exon_list )
	for i in range( exon_num):
		#print "-----> comp with exon"
		#print seg
		#print exon_list[i]
		if( i < exon_num and ( exon_list[i][1] < seg[idx][0] ) ):
			continue
		while( idx < len(match)-1 and ( exon_list[i][0] > seg[idx][1] ) ):
			idx +=1
		if 0 < overlap_interval( seg[idx][0], seg[idx][1], exon_list[i][0], exon_list[i][1]):
			#print 
			#print overlap_interval( seg[idx][0], seg[idx][1], exon_list[i][0], exon_list[i][1])
			best_iso = 1
			match[idx] = 1
			if 1 == len(match) and exon_list[i][0] <= seg[idx][0] and seg[idx][1] <= exon_list[i][1]:
				match[idx] = 2
			if 0 == idx and seg[idx][1] == exon_list[i][1]:
				match[idx] = 2
			if idx == len(match) and seg[idx][0] == exon_list[i][0]:
				match[idx] = 2
			
		if seg[-1][1] < exon_list[i][0]:
			break
	
	flag = 2
	for item in match:
		if 1 == item :
			flag = 1
		elif 0 == item:
			flag = 0
	#print "===" + str(flag)
	#print seg
	#print match
	if -1 == best_iso:
		flag = -1
	return flag

#########################	
# 0 : intergenic
# 1 : intron
# 2 : cross-bd
# 3 : exon 
def compare_gene( s, e, ml, seg, gene_list ):
	flag = 0
	isoflag = 0
	best_gene = -1
	#best_iso  = -1
	for i in range( len(gene_list) ):
		gene = gene_list[i]
		if e < gene[1]:
			break
		elif gene[2] < s :
			continue
		elif 0 < overlap_interval( s, e, gene[1], gene[2]):
			#print "comparing with "
			#print gene[1]
			#print gene[2]
			#flag = 1
			for isoform in gene[4]:
				if e < isoform[0]:
					break
				elif isoform[1] < s:
					continue
				elif 0 < overlap_interval( s, e, isoform[0], isoform[1]):
					#print "-->comparing with isoforms"
					#print gene[1]
					#print gene[2]
					# at least intronic
					isoflag = compare_isoform( seg, isoform[2])
					# 2: cDNA
					# 1: cross bd
					# 0: intron
					if isoflag+1 > flag:
						flag = isoflag+1
						best_gene = i
	return flag, best_gene					




#########################	
def dump_stat( out, c_list, total):
	sw = open(out, 'w')
	#sw.write( "Total\tTotal\tTotal\tTotal\t" + c_list[0] + "\n")
	sw.write( "Analyzed\t\t\t\t" + str(c_list[0]) + "\n")
	v1 = ["cont", "junc"]
	v2 = ["none", "one-side", "two-side"]
	v3 = ["inter", "intron", "cross", "exon", "fusion"]
	c = 1
	for i in v1:
		for j in v2:
			for k in v3:
				sw.write(i + "\t" + j + "\t" + k + "\t" + str(c_list[c]) + "\n")
				c+=1
	sw.write("Orphan\t\t\t\t" + str(c_list[-1]) + "\n")
	miss = total - c_list[0]
	if 0 == total:
		miss = 0
		total = c_list[0]
	sw.write("Missing\t\t\t\t" + str(miss)  + "\n")
	sw.write("Total\t\t\t\t"   + str(total) + "\n")
	sw.close()

def compare_loc(a, b):
	return abs(a - b) <= diff_th

#########################	
def main():
	args = sys.argv[1:]
	if len(args) != 4:
		usage()

	gene_model= GM.load_GTF( sys.argv[1]  )
	print "Loading Gene Models with " + str( len(gene_model) ) + " chromosomes"
	#exit(0)

	total = [0, 0, 0, 0, 0]
	predict_dict = load_Prediction( sys.argv[3] )
	for x,y in predict_dict.iteritems():
		z_class = y[0]#int(y.strip().split("\t")[-1])
		total[z_class] +=1
	print total
	count = 0
	sam_buffer = []
	loci_list  = []
	prev_read  = ""
	flag       = 0 
	prefix = sys.argv[4]

	c0 = [0]*32
	c1 = [0]*32
	c2 = [0]*32
	c4 = [0]*32
	c5 = [0]*32
	sw0 = open( prefix + ".con",  'w')
	sw1 = open( prefix + ".dis",  'w')
	sw2 = open( prefix + ".rf",  'w')
	sw4 = open( prefix + ".split",  'w')
	sw5 = open( prefix + ".remain",  'w')
	
	sr = open(sys.argv[2], 'r')
	for line in sr:
		if ("@" == line[0]):
			continue
		t_list = line.split("\t")
		if t_list[0] != prev_read:
			if "" != prev_read:
				#for rc in sam_buffer:
				#	print rc
				rid = int( prev_read.split(".")[1] )
				if 12 == flag & 12:
					if rid in predict_dict:
						z_tmp = predict_dict[rid][1]
						z_class = predict_dict[rid][0]
						z_sam = ""
						for rc in sam_buffer:
							z_sam += rc.strip() + "\tXC:i:4\tXD:i:0\tXV:i:0\tXY:Z:00\n"
						if 0 == z_class:
							sw0.write(z_sam)
							c0[0] += 1
							c0[-1] += 1
						elif 1 == z_class:
							c1[0] += 1
							c1[-1] += 1
							sw1.write(z_sam)
						elif 2 == z_class:
							c2[0] += 1
							c2[-1] += 1
							sw2.write(z_sam)
						elif 4 == z_class:
							c4[0] += 1
							c4[-1] += 1
							sw4.write(z_sam)
					
				elif 12 != flag &12:
					#print sam_buffer
					#print flag&12
					for item in loci_list:
						#print item
						if 0 <= item[3]:
							ref1 = item[0]
							s1 = item[1]
							e1, m1 ,seg1 = parse_cigar( item[1], item[2]) 
							#print seg1
							gtf_flag1, best_gene1 = compare_gene( s1, e1, m1, seg1, gene_model[item[0]] )
						elif 0 > item[3]:
							ref2 = item[0]
							s2 = item[1]
							e2, m2 ,seg2 = parse_cigar( item[1], item[2]) 
							#print seg2
							gtf_flag2, best_gene2 = compare_gene( s2, e2, m2, seg2, gene_model[item[0]] )


					#print str(gtf_flag1) + " " + str(best_gene1)
					#print str(gtf_flag2) + " " + str(best_gene2)
					# exon or junction
					c_flag = 1
					if 1 < len(seg1) or 1 < len(seg2):
						c_flag = 2
					#print "--->c_flag"
					#print c_flag
					
					# 0 for interginc, 1 for inron, 2 for cross-b, 3 for cDNA, 4 for fusion
					d_flag = gtf_flag1
					if d_flag > gtf_flag2:
						d_flag = gtf_flag2
					if -1 < best_gene1 and -1 < best_gene2  and best_gene1 != best_gene2:
						d_flag = 4

					if rid not in predict_dict:
						z_sam = ""
						for rc in sam_buffer:
							z_sam += rc.strip() + "\tXC:i:" + str(c_flag) + "\tXD:i:" + str(d_flag) + "\tXV:i:0\tXY:Z:" + str(gtf_flag1) + str(gtf_flag2) + "\n"
						sw5.write(z_sam)
						c5[0]+=1
						c5[( (c_flag-1)*12) + (d_flag+1) + (cmp_flag*4)]+=1
					else:
						p_list = predict_dict[ rid ][1].split("\t")
						z_class = predict_dict[ rid ][0]
						cmp_flag = 0
						# 2: both-side identical
						# 1: one-side
						# 0: different
						if p_list[1] == ref1:
							if compare_loc(int(p_list[2]), s1) and compare_loc(int(p_list[3]), e1):
								cmp_flag += 1
								if compare_loc(int(p_list[5]), s2) and compare_loc(int(p_list[6]), e2):
									cmp_flag += 1
							elif compare_loc(int(p_list[5]), s1) and compare_loc(int(p_list[6]), e1):
								cmp_flag += 1
								if compare_loc(int(p_list[2]), s2) and compare_loc(int(p_list[3]), e2):
									cmp_flag += 1
							elif compare_loc(int(p_list[2]), s2) and compare_loc(int(p_list[3]), e2):
								cmp_flag += 1
							elif compare_loc(int(p_list[5]), s2) and compare_loc(int(p_list[6]), e2):
								cmp_flag += 1
					#for rc in sam_buffer:
					#	print rc
					#print str(gtf_flag1) + " " + str(gtf_flag2)
						#if 3 != gtf_flag1 or 3!= gtf_flag2:
						#	for rc in sam_buffer:
						#		print rc
						#	print "Not Compatible"
						#	print "left " + str(gtf_flag1) 
						#	print "right " + str(gtf_flag2) 
						z_sam = ""
						for rc in sam_buffer:
							z_sam += rc.strip() + "\tXC:i:" + str(c_flag) + "\tXD:i:" + str(d_flag) + "\tXV:i:" + str(cmp_flag) + "\tXY:Z:" + str(gtf_flag1) + str(gtf_flag2) + "\n"
						if 0 == z_class:
							sw0.write(z_sam)
							c0[0]+=1
							c0[( (c_flag-1)*12) + (d_flag+1) + (cmp_flag*4)]+=1
						elif 1 == z_class:
							sw1.write(z_sam)
							c1[0]+=1
							c1[( (c_flag-1)*12) + (d_flag+1) + (cmp_flag*4)]+=1
						elif 2 == z_class:
							sw2.write(z_sam)
							c2[0]+=1
							c2[( (c_flag-1)*12) + (d_flag+1) + (cmp_flag*4)]+=1
						elif 4 == z_class:
							sw4.write(z_sam)
							c4[0]+=1
							c4[( (c_flag-1)*12) + (d_flag+1) + (cmp_flag*4)]+=1

			prev_read = t_list[0]
			sam_buffer = []
			loci_list  = []
		sam_buffer.append( line )
		count += 1
		flag = int(t_list[1])
		if ( 256 > flag ):
			loci_list.append( [ t_list[2], int(t_list[3]), t_list[5], int(t_list[8] ), flag  ] ) 

		if 0 == count%1000000:
			print  >> sys.stderr, str(count) + " line parsed"
	sr.close()
	sw0.close()
	sw1.close()
	sw2.close()
	sw4.close()
	sw5.close()
	dump_stat( "concord.stat", c0, total[0])
	dump_stat( "disconcord.stat", c1, total[1] )
	dump_stat( "circ_RF.stat", c2, total[2])
	dump_stat( "circ_BSJ.stat", c4, total[4])
	dump_stat( "other.stat", c5, 0)
		
if __name__ == "__main__":
	sys.exit(main())
