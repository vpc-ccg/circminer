import sys

gt_fname = sys.argv[1]
my_fname = sys.argv[2]

bp_res = 10

gt_list = []
gt_mark = []

def search_circ(gt_list, ch, sp, ep):
	i = -1
	for ev in gt_list:
		i += 1
		(ech, esp, eep, esup, etyp) = ev
		if ch == ech and abs(esp - sp) <= bp_res and abs(eep - ep) <= bp_res and gt_mark[i] == 0:
			gt_mark[i] += 1
			l = [ech, str(esp), str(eep), str(esup), etyp]
			print '\t'.join(l)
			return
	
	print 'FP'

with open(gt_fname) as gf:
	for l in gf:
		ll = l.strip().split()
		ch = ll[0]
		sp = int(ll[1])
		ep = int(ll[2])
		sup = int(ll[3])
		typ = ' '.join(ll[4:])
		gt_list.append((ch, sp, ep, sup, typ))
	
gt_mark = [ 0 ] * len(gt_list)

with open(my_fname) as mf:
	for l in mf:
		ll = l.strip().split()
		ch = ll[0]
		sp = int(ll[1])
		ep = int(ll[2])
		sup = int(ll[3])
		typ = ll[4]

		sys.stdout.write(l.strip() + '\t')
		search_circ(gt_list, ch, sp, ep)
