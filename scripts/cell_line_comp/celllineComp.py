import sys

read_count = [ 322475040, 147261832, 206362733, 199922486 ]

def load_file(fname):
	d = {}
	with open(fname) as fin:
		for l in fin:
			ll = l.strip().split()
			k = '_'.join(ll[:3])
			sup = int(ll[3])
			d[k] = sup

	return d

def sort_dict(d):
	ds = {k: v for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True)}
	return ds

def print_top_x(pred, postd, pre_rc, post_rc, x=100, enr_ratio=5.0):
	preds = sort_dict(pred)
	combined = {}
	for k in preds.keys():
		pre_sup = preds[k]
		post_sup = 0
		if k in postd.keys():
			post_sup = postd[k]
		#---
		combined[k] = (pre_sup, -1*post_sup)

	combined = sort_dict(combined)
	i = 0
	for k in combined.keys():
		(pre_sup, post_sup) = combined[k]
		#---
		post_sup = -1 * post_sup

		not_dep = 'Y' if ((post_sup * pre_rc) / (pre_sup * post_rc)) >= 1.0 else 'N'
		enr = 'Y' if ((post_sup * pre_rc) / (pre_sup * post_rc)) >= enr_ratio else 'N'

		event = k.split('_')
		print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(event[0], event[1], event[2], pre_sup, post_sup, not_dep, enr))

		i += 1
		if i >= x:
			break

def main():
	args = sys.argv[1:]
	cell = args[0]
	pref = args[1]
	postf = args[2]

	if cell not in [ 'HeLa', 'Hs68' ]:
		print('Error: cell line not found')
		exit(1)

	#print(pref, '---->', postf)

	pre_dict = load_file(pref)
	post_dict = load_file(postf)

	rc_ind = [0, 1] if cell == 'HeLa' else [2, 3]

	print_top_x(pre_dict, post_dict, pre_rc=read_count[rc_ind[0]], post_rc=read_count[rc_ind[1]], x=100)

if __name__ == '__main__':
	main()
