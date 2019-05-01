import sys

def convert_to_list(rid_dic):
	vmax = 0
	for r in rid_dic.keys():
		if r.isdigit():
			vmax = max(vmax, int(r))
		else:
			vmax = max(vmax, int(r.split(':')[1]))
	rid_list = [ None ] * (vmax+1)

	for r in rid_dic.keys():
		rid = 0
		if r.isdigit():
			rid = int(r)
		else:
			rid = int(r.split(':')[1])
		rid_list[rid] = rid_dic[r]

	return rid_list


def load_slog(slog):
	rid2ev = {}
	ev_rids = {}
	rid2sides = {}
	rid = ''
	rids = []
	with open(slog) as sf:
		for l in sf:
			ll = l.strip().split()
			if ll[0] == '>':
				rid = ll[2]
			# elif ll[0] in ['**', '+']:
			elif ll[0] == '+':
				side1 = ll[1]
				side2 = ll[2] if len(ll) > 2 else '-'
				rids.append((rid, side1, side2))

			elif ll[0] == 'summary:':
				ev = ll[3][8:]
				just_rids = [x[0] for x in rids]
				just_rids.sort()
				# rids = sorted(rids, key=lambda tup: tup[0])
				ev_rids[ev] = just_rids
				for (r, s1, s2) in rids:
					rid2ev[r] = ev
					rid2sides[r] = s1, s2
				rids = []

	return rid2sides, rid2ev, ev_rids

def load_report(rep_fname):
	rid2rep = {}
	rep_rids = {}
	with open(rep_fname) as rf:
		for l in rf:
			ll = l.strip().split()
			# circ = 'chr{}:{}-{}'.format(ll[0], ll[1], ll[2])
			circ = '{}:{}-{}'.format(ll[0], ll[1], ll[2])
			rids = ll[5].split(',')
			con_rids = []
			for r in rids:
				conr = r.split(':')[1]
				rid2rep[conr] = circ
				con_rids.append(conr)
			con_rids.sort()
			rep_rids[circ] = con_rids
	
	return rid2rep, rep_rids

def read_categorize(rid2ev, rid2rep, ev_rids, reado):
	typ = ''
	rev = ''
	rrep = ''
	for r in sorted(rid2ev.keys()):
		rev = rid2ev[r]
		if r not in rid2rep.keys():
			typ = 'not_found'
			rrep = 'x'
		elif rid2ev[r] == rid2rep[r]:
			typ = 'correct'
			rrep = rev
		elif rid2rep[r] in ev_rids.keys():
			typ = 'incorrect_existing_event'
			rrep = rid2rep[r]
		else:
			typ = 'incorrect_nonexisting_event'
			rrep = rid2rep[r]

		print >>reado, '{}\t{}\t{}\t{}'.format(r, typ, rev, rrep)

def read_categorize_list(rid2sides, rid2ev, rid2rep, ev_rids, reado):
	typ = ''
	rev = ''
	rrep = ''
	for r in range(len(rid2ev)):
		if rid2ev[r] == None:
			continue
		rev = rid2ev[r]
		if r >= len(rid2rep) or rid2rep[r] == None:
			typ = 'not_found'
			rrep = 'x'
		elif rid2ev[r] == rid2rep[r]:
			typ = 'correct'
			rrep = rev
		elif rid2rep[r] in ev_rids.keys():
			typ = 'incorrect_existing_event'
			rrep = rid2rep[r]
		else:
			typ = 'incorrect_nonexisting_event'
			rrep = rid2rep[r]

		(s1, s2) = rid2sides[r]
		print >>reado, '{}\t{}\t{}\t{}\t{}\t{}'.format(r, typ, rev, rrep, s1, s2)

def main():
	log_fname = sys.argv[1]
	rep_fname = sys.argv[2]
	out_fname = sys.argv[3]

	outf_reads = open(out_fname+'.reads', 'w')

	rid2sides, rid2ev, ev_rids = load_slog(log_fname)
	rid2rep, rep_rids = load_report(rep_fname)

	rid2sides_l = convert_to_list(rid2sides)
	rid2ev_l = convert_to_list(rid2ev)
	rid2rep_l = convert_to_list(rid2rep)

	read_categorize_list(rid2sides_l, rid2ev_l, rid2rep_l, ev_rids, outf_reads)

	#print rid2ev
	#print ev_rids

if __name__ == '__main__':
		main()

