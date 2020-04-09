#!/usr/bin/env python
#-*- encoding:utf-8 -*-

import sys, os, math
class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def usage():
	print '\nUsage: python ' + sys.argv[0] + ' fastq/pam_file output STAR_Mapings Chimeric_Mappings Chimeric_RID'
	sys.exit(-1)

def get_read_id( fname ):
	r_dict = {}
	count = 0
	sr = open( fname, 'r')

	if fname[-3:] == 'pam':
		for line in sr:
			rname = line.strip().split()[0]
			r_dict[ rname.split(".")[1] ] = 0
	else:
		for line in sr:
			if 0 == count % 4:
				r_dict[ line.strip().split(".")[1] ] = 0
			count += 1

	sr.close()
	return r_dict
	
def main():
	args = sys.argv[1:]
	if len(args) != 5:
		usage()


	r_dict = get_read_id( sys.argv[1] )
	print "Predicted " + str( len( r_dict) ) + " reads "

	bs_dict = {}
	sr      = open( sys.argv[5], 'r')
	for line in sr:
		bs_dict[ line.split()[0]] = 1
	sr.close()
	
	total_count = [0, 0, 0, 0, 0, 0, 0, 0]
	# exon: 1
	# junction: 2
	# oea: 3
	# unmpped: 4
	count = 0
	sw = open( sys.argv[2], 'w' )
	sr = open( sys.argv[3], 'r' )
	old_id = ""
	r_flag = 1
	tmp_list = []
	for line in sr:
		if "@" == line[0]:
			sw.write(line)
			continue
		t_list = line.strip().split()
		
		if ( t_list[0] != old_id ):
			if "" != old_id:
				gn = old_id.split(".")[1]
				if gn in r_dict:
					r_dict[gn] = r_flag
					for item in tmp_list:
						sw.write( item + " XC:i:" + str(r_flag) + "\n")
				total_count[r_flag]+=1
			old_id = t_list[0]
			del tmp_list[:]
			r_flag = 1
		
		tmp_list.append(line.strip())
		flag = int(t_list[1])
		if ( 4 == (4&flag) and 8 == (8&flag)):
			r_flag = 4
		elif ( 4 == (4&flag) or 8 == (8&flag)):
			r_flag = 3
		else:
			if ( -1 < t_list[5].find("N")):
				r_flag = 2

		count +=1
		if 0 == count%1000000:
			print "Parsing " +str(count) + " reads"

	sr.close()

	#print old_id
	gn = old_id.split(".")[1]
	if gn in r_dict:
		r_dict[gn] = r_flag
		for item in tmp_list:
			sw.write( line.strip() + " XC:i:" + str(r_flag) + "\n")
		total_count[r_flag]+=1
	del tmp_list[:]
	r_flag = 6
	#exit(0)
	
	# fusion = 5
	# back   = 6
	del tmp_list[:]	
	sw = open( sys.argv[2], 'a' )
	sr = open( sys.argv[4], 'r' )
	for line in sr:
		if "@" == line[0]:
			continue
		
		t_list = line.strip().split()
		
		#if ( old_id != "" and t_list[0] != old_id ):
		if (  t_list[0] != old_id ):
			if "" != old_id:
				gn = old_id.split(".")[1]
				if gn in r_dict:
					if gn in bs_dict:
						r_flag = 7
					r_dict[gn] = r_flag
					for item in tmp_list:
						sw.write( item + " XC:i:" + str(r_flag) + "\n")
				total_count[r_flag]+=1
			old_id = t_list[0]
			del tmp_list[:]
			r_flag = 6
		
		tmp_list.append(line.strip())
		if ( ("="!=t_list[6]) and ( t_list[2] != t_list[6] )) :
			r_flag = 5

	sr.close()
	gn = old_id.split(".")[1]
	if gn in r_dict:
		if gn in bs_dict:
			r_flag = 7
		r_dict[gn] = r_flag
		for item in tmp_list:
			sw.write( line.strip() + " XC:i:" + str(r_flag) + "\n")
		total_count[r_flag]+=1
	#old_id = t_list[0]
	#del tmp_list[:]
	
	count_list= [0,0,0,0,0,0,0,0]

	sw1 = open( sys.argv[2] + ".summary", "w")
	for x,y in r_dict.iteritems():
		sw1.write(x + "\t" + str(y) + "\n")
		count_list[ y ] +=1
	sw1.close()
	print "Total Prediction:\t" + str(len(r_dict))
	print "Not reported:\t" + str(count_list[0]) 
	print "Exon Mappings:\t" + str(count_list[1]) + " / " + str(total_count[1])
	print "Junc Mappings:\t" + str(count_list[2])+ " / " + str(total_count[2])
	print "OEA :\t" + str(count_list[3])+ " / " + str(total_count[3])
	print "Orphan:\t" + str(count_list[4])+ " / " + str(total_count[4])
	print "Fusion :\t" + str(count_list[5])+ " / " + str(total_count[5])
	print "Chimeric:\t" + str(count_list[6])+ " / " + str(total_count[6])
	print "CIRC:\t" + str(count_list[7])+ " / " + str( len( bs_dict) )
		
if __name__ == "__main__":
	sys.exit(main())
