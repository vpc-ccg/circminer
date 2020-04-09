import sys
import re

tools = [ 'CE', 'CE2', 'CircMarker', 'CircMiner', 'circRNA_finder', 'CIRI', 'DCC', 'find_circ', 'KNIFE', 'PTESFinder', 'SEGEMEHL' , 'NCLscan']

celllines = [ 'HeLa-', 'HeLa+', 'Hs68-', 'Hs68+' ]

read_count = [ 322475040, 147261832, 206362733, 199922486 ]
	
basepath = ''
outputpath = ''

def load_toolname(fname):
	tdict = {}
	with open(fname) as f:
		for l in f:
			ll = l.strip().split()
			tname = ll[0]
			pathdir = ll[1:]
			tdict[tname] = pathdir

	return tdict

def extract_events(fname, outname, col, drop_head=False):
	with open(outname, 'w') as fo:
		with open(fname) as f:
			first = True
			for l in f:
				if first and drop_head:
					continue
				first = False

				ll = l.strip().split()
				fo.write('{}\t{}\t{}\t{}\n'.format(ll[col[0]], ll[col[1]], ll[col[2]], ll[col[3]]))

### CircMiner
def extract_events_circminer(fname, outname, col, drop_head=False):
	with open(outname, 'w') as fo:
		with open(fname) as f:
			first = True
			for l in f:
				if first and drop_head:
					first = False
					continue

				ll = l.strip().split()
				if ll[4] != 'STC' or ll[7] != 'Pass':
					continue

				fo.write('{}\t{}\t{}\t{}\n'.format(ll[col[0]], ll[col[1]], ll[col[2]], ll[col[3]]))

### CIRCexplorer, CIRI, PTESFinder, CircMarker
def extract_events_ce(fname, outname, col, start_shift=0, drop_head=False):
	with open(outname, 'w') as fo:
		with open(fname) as f:
			first = True
			for l in f:
				if first and drop_head:
					first = False
					continue

				ll = l.strip().split()

				fo.write('{}\t{}\t{}\t{}\n'.format(ll[col[0]], int(ll[col[1]])+start_shift, ll[col[2]], ll[col[3]]))

### KNIFE
def extract_events_knife(fname, outname, col, drop_head=False):
	with open(outname, 'w') as fo:
		with open(fname) as f:
			first = True
			for l in f:
				if first and drop_head:
					first = False
					continue

				ll = l.strip().split()

				if float(ll[2]) < 0.0 and float(ll[4]) < 0.0:
					continue

				ll0 = re.split(':|\|', ll[0])
				ch = ll0[0]
				st = ll0[2] if ll0[6] == '-' else ll0[4]
				en = ll0[4] if ll0[6] == '-' else ll0[2]

				fo.write('{}\t{}\t{}\t{}\n'.format(ch, st, en, ll[col[3]]))

### NCLscan
def extract_events_ncls(fname, outname, col, drop_head=False):
	with open(outname, 'w') as fo:
		with open(fname) as f:
			first = True
			for l in f:
				if first and drop_head:
					first = False
					continue

				ll = l.strip().split()

				ch = ll[col[0]]
				st = ll[col[1]] if ll[col[4]] == '-' else ll[col[2]]
				en = ll[col[2]] if ll[col[4]] == '-' else ll[col[1]]

				fo.write('{}\t{}\t{}\t{}\n'.format(ch, st, en, ll[col[3]]))

def extract_each_cellline(toolname, pathlist):
	for i in range(0, 4):
		print(toolname)
		if toolname == 'CircMiner':
			extract_events_circminer(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 3])

		if toolname == 'KNIFE':
			extract_events_knife(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 2, 4, 5], drop_head=True)

		if toolname == 'CE' or toolname == 'CE2':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 12], start_shift=1)

		if toolname == 'PTESFinder':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 4], start_shift=1, drop_head=True)

		if toolname == 'CircMarker':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 3], start_shift=0)

		if toolname == 'DCC':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 3], start_shift=0, drop_head=True)

		if toolname == 'find_circ':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 4], start_shift=1)

		if toolname == 'circRNA_finder':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 2, 4], start_shift=1)

		if toolname == 'SEGEMEHL':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [1, 2, 3, 0], start_shift=0)

		if toolname == 'CIRI':
			extract_events_ce(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [1, 2, 3, 16], start_shift=0, drop_head=True)

		if toolname == 'NCLscan':
			extract_events_ncls(basepath + '/' + toolname + '/' + pathlist[i], outputpath + '/' + toolname + '.' + celllines[i], [0, 1, 4, 10, 2])

def extract_events_all(tdict):
	for t in tools:
		extract_each_cellline(t, tdict[t])

###################################################################################################
def load_event(fname):
	evl = []
	with open(fname) as fin:
		for l in fin:
			ll = l.strip().split()
			evl.append(ll)

	return evl

def load_file(fname, sup_th=1, head=None):
	d = {}
	with open(fname) as fin:
		i = 0
		for l in fin:
			i += 1
			ll = l.strip().split()
			k = '_'.join(ll[:3])
			sup = int(ll[3])

			if sup >= sup_th:
				d[k] = sup

			if head != None and i >= head:
				break

	return d


def make_table(cellline, ev_fname, top_only=True):
	if cellline not in [ 'HeLa', 'Hs68' ]:
		print('Error: Invalid cellline!')
		return

	cell_ind = [0, 1] if cellline == 'HeLa' else [2, 3]
	pre_rc = read_count[cell_ind[0]]
	post_rc = read_count[cell_ind[1]]

	print('chr\t{:8s}\t{:8s}\t'.format('start', 'end'), end = '')
	all_pre = {}
	all_post = {}
	for t in tools:
		if top_only:
			all_pre[t]  = load_file(topdir + '/' + t + '.top100')
		else:
			all_pre[t]  = load_file(tempdir + '/' + t + '.' + celllines[cell_ind[0]])
		all_post[t] = load_file(tempdir + '/' + t + '.' + celllines[cell_ind[1]])
		print('{:9s}\t'.format(t), end='')

	print('')

	evl = load_event(ev_fname)
	for ev in evl:
		print('{}\t{:8s}\t{:8s}\t'.format(ev[0], ev[1], ev[2]), end = '')
		k = '_'.join(ev[:3])

		for t in tools:
			pre_sup = 0
			if k in all_pre[t].keys():
				pre_sup = all_pre[t][k]
			post_sup = 0
			if k in all_post[t].keys():
				post_sup = all_post[t][k]

			if pre_sup == 0:
				ratio = -1
			else:
				ratio = post_sup * pre_rc / (pre_sup * post_rc) if pre_sup > 0 else 0

			res = '-'	# depleted
			if ratio >= 5.0:
				res = 'Y'	# enriched
			elif ratio >= 1.0:
				res = 'X'	# not depleted and not enriched

			#print('{:s}({:.2f})({:3d}-{:3d})\t'.format(res, ratio, pre_sup, post_sup), end='')
			print('{:2.2f}{:4s}\t'.format(ratio, ' '), end='')
		print('')

################################################################################################################
def comm(evs1, evs2):
	comm = 0
	for k in evs1.keys():
		if k in evs2.keys():
			comm += 1
	return comm

def make_table_comm(x=100):
	topx = {}
	print('\t\t', end='')
	for t in tools:
		print('\t{:9s}'.format(t), end='')
		topx[t] = load_file(topdir + '/' + t + '.top' + str(x))
	print('')

	for t1 in tools:
		print('{:9s}'.format(t1), end='')
		for t2 in tools:
			c = comm(topx[t1], topx[t2])
			print('\t{:9d}'.format(c), end='')
		print('')

##########################################################################################################################
def calc_ratio(pred, postd, pre_rc, post_rc):
	ratiod = {}
	for k in pred.keys():
		pre_sup = pred[k]
		post_sup = 0
		if k in postd.keys():
			post_sup = postd[k]

		if pre_sup == 0:
			ratio = -1
		else:
			ratio = post_sup * pre_rc / (pre_sup * post_rc) if pre_sup > 0 else 0

		ratiod[k] = ratio

	return ratiod

def make_paper_table(cellline, not_dep_th=1.0, enr_th=5.0):
	if cellline not in [ 'HeLa', 'Hs68' ]:
		print('Error: Invalid cellline!')
		return

	cell_ind = [0, 1] if cellline == 'HeLa' else [2, 3]
	pre_rc = read_count[cell_ind[0]]
	post_rc = read_count[cell_ind[1]]

	print('{:13s}\t{:7s}\t{:7s}\t{}\t{}\t{}\t{}\t{}\t{}'.format(cellline, 'RNaseR-', 'RNaseR+', 'Not depleted', 'Percentage %', 'Top 10 not depleted', 'Top 10 enriched', 'Top 100 not depleted', 'Top 100 enriched'))

	all_pre = {}
	all_post = {}
	top10_pre = {}
	top100_pre = {}
	for t in tools:
		all_pre[t]  = load_file(tempdir + '/' + t + '.' + celllines[cell_ind[0]], sup_th=2)
		all_post[t] = load_file(tempdir + '/' + t + '.' + celllines[cell_ind[1]], sup_th=2)

		ratiod = calc_ratio(all_pre[t], all_post[t], pre_rc, post_rc)
		total_not_dep = 0
		for k in ratiod.keys():
			if ratiod[k] >= not_dep_th:
				total_not_dep += 1

		top100_pre[t] = load_file(topdir + '/' + t + '.top100')
		top100_ratiod = calc_ratio(top100_pre[t], all_post[t], pre_rc, post_rc)
		top100_not_dep = 0
		top100_enr = 0
		for k in top100_ratiod.keys():
			if top100_ratiod[k] >= not_dep_th:
				top100_not_dep += 1
			if top100_ratiod[k] >= enr_th:
				top100_enr += 1

		top10_pre[t] = load_file(topdir + '/' + t + '.top100', head=10)
		top10_ratiod = calc_ratio(top10_pre[t], all_post[t], pre_rc, post_rc)
		top10_not_dep = 0
		top10_enr = 0
		for k in top10_ratiod.keys():
			if top10_ratiod[k] >= not_dep_th:
				top10_not_dep += 1
			if top10_ratiod[k] >= enr_th:
				top10_enr += 1

		print('{:13s}\t{:7d}\t{:7d}\t{}\t{:.2f}\t{}\t{}\t{}\t{}'.format(t, len(all_pre[t]), len(all_post[t]), total_not_dep, total_not_dep*100.0 / len(all_pre[t]), top10_not_dep, top10_enr, top100_not_dep, top100_enr ))

def main():
	args = sys.argv[1:]
	print(args)

	global tempdir
	global topdir
	
	mode = args[0]
	
	if mode == 'init':
		fname = args[1]
		global basepath 
		basepath = args[2]
		global outputpath
		outputpath = args[3]

		tdict = load_toolname(fname)
		extract_events_all(tdict)
	
	elif mode == 'event':
		cellline = args[1]
		ev_fname = args[2]
		tempdir = args[3]
		topdir = args[4]
		make_table(cellline, ev_fname, top_only=True)

	elif mode == 'comm':
		topdir = args[1]
		make_table_comm(100)
	
	elif mode == 'paper_table':
		cellline = args[1]
		tempdir = args[2]
		topdir = args[3]
		not_dep_ratio = float(args[4])
		enr_ratio = float(args[5])

		make_paper_table(cellline, not_dep_th=not_dep_ratio, enr_th=enr_ratio)

	else:
		print('Error: mode {} not recognized'.format(mode))

if __name__ == '__main__':
	main()
