import sys
from collections import defaultdict

def load_gtf(gtffilename):
	print('Started loading GTF file...')

	tr_exons = defaultdict(list)
	exons_tr_beg = {}
	exons_tr_end = {}
	with open(gtffilename) as gin:
		cur_tr_id = ''
		for l in gin:
			if l[0] == '#':
				continue
			ltok = l.strip().split()
			if ltok[2] == 'transcript':
				cur_tr_id = ltok[13][1:-2]
			elif ltok[2] == 'exon':
				ch = ltok[0]
				beg = int(ltok[3])
				end = int(ltok[4])
				gname = ltok[19][1:-2]
				exon_num = ltok[17][1:-2]
				
				if ch not in exons_tr_beg.keys():
					exons_tr_beg[ch] = {}

				if ch not in exons_tr_end.keys():
					exons_tr_end[ch] = {}

				if beg not in exons_tr_beg[ch].keys():
					exons_tr_beg[ch][beg] = []

				if end not in exons_tr_end[ch].keys():
					exons_tr_end[ch][end] = []

				exons_tr_beg[ch][beg].append(('{}({})'.format(cur_tr_id,gname), exon_num))
				exons_tr_end[ch][end].append(('{}({})'.format(cur_tr_id, gname), exon_num))

				tr_exons[cur_tr_id].append((ch, beg, end))

	print('Successfully loaded GTF file...')

	return tr_exons, exons_tr_beg, exons_tr_end

def intersection(l1, l2):
	l1_vals = [ v for (v, e) in l1 ]
	l2_vals = [ v for (v, e) in l2 ]
	comm_vals = [ val for val in l1_vals if val in l2_vals ]
	beg_enums = [ enum for (v, enum) in l1 if v in comm_vals ]
	end_enums = [ enum for (v, enum) in l2 if v in comm_vals ]

	return comm_vals, beg_enums, end_enums

def annotate(crfilename, exons_tr_beg, exons_tr_end, outfilename):
	print('Started annotating circRNA report file...')

	fout = open(outfilename, 'w')

	with open(crfilename) as fin:
		for l in fin:
			ltok = l.strip().split()
			if ltok[0] == 'chr':
				fout.write('{}\ttranscripts\n'.format(l.strip()))
				continue
			
			ch = ltok[0]
			beg = int(ltok[1])
			end = int(ltok[2])

			trans = ''
			if ch not in exons_tr_beg.keys() or beg not in exons_tr_beg[ch].keys() or ch not in exons_tr_end.keys() or end not in exons_tr_end[ch].keys():
				trans = 'NA'
			else:
				trans_list, beg_enums, end_enums = intersection(exons_tr_beg[ch][beg], exons_tr_end[ch][end])
				for i in range(len(trans_list)):
					trans_list[i] = '{}[{}-{}]'.format(trans_list[i], beg_enums[i], end_enums[i])

				if len(trans_list) == 0:
					trans = 'NA'
				else:
					trans = ', '.join(trans_list)
			
			fout.write('{}\t{}\n'.format(l.strip(), trans))

	fout.close()

	print('Annotating finished.')

def usage():
	print('\nUsage:\npython {} circ_report_file gtf_file output_file'.format(sys.argv[0]))

def main():
	args = sys.argv[1:]
	if len(args) != 3:
		usage()
		exit(1)

	crfilename	= args[0]
	gtffilename	= args[1]
	outfilename	= args[2]

	temp, gtf_dict_beg, gtf_dict_end = load_gtf(gtffilename)
	annotate(crfilename, gtf_dict_beg, gtf_dict_end, outfilename)

if __name__ == '__main__':
	main()
