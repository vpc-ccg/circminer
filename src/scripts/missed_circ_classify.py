import re
import sys
import subprocess

def get_len(cigar):
	mlen = 0
	mysum = 0
	ind = 0
	lens = re.split("M|D|N|=|X|I|S|H|P", cigar)
	if lens[-1] == '':
		del lens[-1]

	#print cigar
	#print lens
	exon_starts_rel = [0]
	exon_mlen = [ ]
	mlen = 0
	for x in lens:
		ind += len(x)

		if cigar[ind] in 'MDN=X':
			mysum += int(x)
			if cigar[ind] == 'N':
				exon_starts_rel.append(mysum)
				exon_mlen.append(mlen)
				mlen = 0
			else:
				mlen += int(x)

		ind += 1
	
	exon_mlen.append(mlen)

	return mysum, exon_starts_rel, exon_mlen

samf = sys.argv[1]
gtff = sys.argv[2]
reasonf = sys.argv[3]

reason_file = open(reasonf, 'w')

def find_fault(chrom, pos, exon_starts_rel, exon_mlen):
	fault = []
	exon_number = [ {}, {}, {}, {} ]
	for i in range(len(exon_starts_rel)):
		rpos = pos + exon_starts_rel[i]
		mlen = exon_mlen[i]
		cmd = '$3=="exon" && $1 == "{}" && $5 >= {} && $4 <= {} {{print $0}}'.format(chrom, rpos, rpos+mlen-1)
		print cmd
		
		ret = subprocess.check_output(['awk', cmd, gtff])

		retl = re.split('\n', ret)
		retll = []
		for r in retl:
			if r == '':
				continue
			rr = re.split('\t', r)
			retll.append(rr)
		
		print len(retll)
		print retll
		
		if len(retll) == 0:
			fault.append('intron')
		
		fully_in = False
		same_start = False
		same_end = False
		for r in retll:
			spos = int(r[3])
			epos = int(r[4])
			rest = r[8].split(' ')
			#print rest
			
			trans_id = rest[rest.index('transcript_id') + 1]
			ex_num = int(rest[rest.index('exon_number') + 1][1:-2])
			
			if spos <= rpos and epos >= rpos+mlen-1:
				fully_in = True
			if spos == rpos and epos >= rpos+mlen-1:
				if i > 0:
					exon_number[i][trans_id] = ex_num
				same_start = True
			if spos <= rpos and epos == rpos+mlen-1:
				if i == 0:
					exon_number[i][trans_id] = ex_num
				same_end = True
		
		if (i == 0 and same_end) or (i == len(exon_starts_rel)-1 and same_start) or (i > 0 and i < len(exon_starts_rel)-1 and same_start and same_end):
			fault.append('boundryok')
		elif fully_in:
			fault.append('middleexon')
		else:
			fault.append('retention')


	print fault
	final_cat = '?'
	if len(exon_starts_rel) == 1:
		if fault[0] == 'intron':
			final_cat = 'IntronicMate'
		elif fault[0] == 'retention':
			final_cat = 'IntronRetention'
		elif fault[0] in ['boundryok', 'middleexon']:
			final_cat = 'Normal'
	
	else:
		if 'intron' in fault:
			final_cat = 'JunctionToIntron'
		elif 'middleexon' in fault:
			final_cat = 'JunctionToMiddleOfExon'
		elif 'retention' in fault:
			final_cat = 'IntronRetention'
		else:
			final_cat = 'Normal'

	if final_cat == 'Normal' and len(fault) > 1:
		print exon_number
		consec = False
		for trans_id in exon_number[0].keys():
			consec_minor = True
			pre_ex_num = exon_number[0][trans_id]
			for i in range(1, len(fault)):
				x = exon_number[i]
				if not trans_id in x.keys():
					consec_minor = False 
					break
				if abs(x[trans_id] - pre_ex_num) != 1:
					consec_minor = False
					break
				pre_ex_num = x[trans_id]

			if consec_minor:
				consec = True
				break

		if not consec:
			final_cat = 'ExonSkipping'

	return final_cat

with open(samf) as f:
	sam_rec = []
	for l in f:
		if len(sam_rec) == 3:
			sam_rec = []

		l = l.strip().split()
		sam_rec.append(l)
		if len(sam_rec) < 3:
			continue
	
		pos_in_pair = [(int(x[1]) >> 6) % 2 for x in sam_rec]
		if pos_in_pair.count(0) == 1:
			k = 0
		else:
			k = 1
			
		ind = pos_in_pair.index(k)
		l = sam_rec[ind]	# the mate that is not split

		rname = l[0]
		chrom = l[2]
		pos = int(l[3])
		cigar = l[5]

		rest = l[11:]
		md = [x[5:] for x in rest if x[:2] == 'MD'][0]
		nm = int([x[5:] for x in rest if x[:2] == 'NM'][0])
		
		clen, exon_starts_rel, exon_mlen = get_len(cigar)

		print '{}:{}-{}, cigar: {}, md: {}, nm:{}'.format(chrom, pos, pos + clen-1, cigar, md, nm)
		print exon_starts_rel
		print exon_mlen

		final_cat = find_fault(chrom, pos, exon_starts_rel, exon_mlen)
		
		if final_cat == 'Normal' and nm > 4:
			final_cat = 'HighError_{}'.format(nm)
		elif final_cat == 'Normal' and sum(exon_mlen) < (76 - 8):
			final_cat = 'HighSoftClip_{}'.format(76 - sum(exon_mlen))

		print '{}\t{}'.format(rname, final_cat)
		print >>reason_file, '{}\t{}'.format(rname, final_cat)
