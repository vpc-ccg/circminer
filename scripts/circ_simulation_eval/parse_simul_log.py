import sys
import re

def parse_old_simul(simul_log_fname, out):
	ch = '-'
	st = 0
	en = 0
	sup = 0
	split_reads = 0

	with open(simul_log_fname) as f:
		for l in f:
			if l[:3] == 'iso':
				continue

			if l[:2] == '**':
				split_reads += 1
				continue

			if l[0] == '!':
				continue

			ll = re.split('\t| |:|\||\n', l)

			if ll[0] == '>':
				sup = ll[1]
			
			else:
				if ch != '-':
					print >>out, '{}\t{}\t{}\t{}\t{}'.format(ch, st, en, split_reads, sup)
				ch = ll[0]
				st = ll[5]
				en = ll[6]
				sup = 0
				split_reads = 0
	if ch != '-':
		print >>out, '{}\t{}\t{}\t{}\t{}'.format(ch, st, en, split_reads, sup)

def parse_new_simul(simul_log_fname):
	event_reads = {}
	reads = {}

	rid = -1
	mate = 0
	mate1_spos = -1
	mate2_spos = -1
	split_mates = 0
	split_reads = 0
	with open(simul_log_fname) as sf:
		for l in sf:
			ll = l.strip().split()
			if ll[0] == '>':
				if rid != -1:	# not the first rid in this event
					if split_mates > 2:
						print >>sys.stderr, 'Error: {} split mates!'.format(split_mates)
					if split_mates > 0:
						split_reads += 1
					reads[rid] = (split_mates, mate, mate1_spos, mate2_spos)
				
				rid = ll[2]
				mate = 0
				mate1_spos = -1
				mate2_spos = -1
				split_mates = 0
			
			elif ll[0] == '**':
				split_mates += 1
				mate = ll[1]
				mate1_spos = ll[2]
				mate2_spos = ll[3]

			elif ll[0] == 'summary:':
				if split_mates > 2:
					print >>sys.stderr, 'Error: {} split mates!'.format(split_mates)
				if split_mates > 0:
					split_reads += 1
				reads[rid] = (split_mates, mate, mate1_spos, mate2_spos)

				gid = ll[1]
				tid = ll[2]
				circRNA = ll[3][8:]
				spliced_len = ll[4][14:]
				support = int(ll[5][23:])
				exon_index = ll[6][10:].split('-')
				if len(exon_index) != 2:
					print >>sys.stderr, 'Error: Invalid exon index!!!'
					start_exon = -1
					end_exon = -1
				else:
					start_exon = int(exon_index[0])
					end_exon = int(exon_index[1])
				type = ll[8:]

				event_reads[circRNA] = (gid, tid, spliced_len, split_reads, start_exon, end_exon, type, reads)
				# print gid, tid, circRNA, spliced_len, split_reads, start_exon, end_exon, type
				# print reads
				if support != split_reads:
					print >>sys.stderr, 'Error: support does not match!', support, split_reads

				# reset variables
				reads = {}
				rid = -1
				split_reads = 0
	return event_reads

def search_circRNA(event_reads, circ_fname, out_fname, offset):
	outf = open(out_fname, 'w')

	with open(circ_fname) as cf:
		for l in cf:
			ll = l.strip().split()
			if len(ll[0]) < 3 or ll[0][:3] != 'chr':
				ch = 'chr' + ll[0]
			else:
				ch = ll[0]
			st = int(ll[1]) + offset
			en = ll[2]
			cr = '{}:{}-{}'.format(ch, st, en)
			print cr

			if cr in event_reads.keys():
				(gid, tid, spliced_len, split_reads, start_exon, end_exon, type, reads) = event_reads[cr]
				for rid in reads.keys():
					(split_mates, mate, mate1_spos, mate2_spos) = reads[rid]
					if split_mates > 0:
						print >>outf, 'simulate:{}'.format(rid)
			else:
				print "not found"

def search_rids(event_reads, rid_fname):
	with open(rid_fname) as rf:
		for l in rf:
			l = l.strip()
			if l[:8] == 'simulate':
				l = l[9:]
			for ev in event_reads.keys():
				if rid in event_reads[ev].reads.keys():
					r = event_reads[ev].reads[rid]

def main():
	simul_log_fname = sys.argv[1]
	circ_fname = sys.argv[2]
	out_fname = sys.argv[3]

	event_reads = parse_new_simul(simul_log_fname)
	print >>sys.stderr, len(event_reads)
	search_circRNA(event_reads, circ_fname, out_fname, 1)

if __name__ == '__main__':
	main()