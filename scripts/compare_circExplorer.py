import sys

my_report = sys.argv[1]
ce_report = sys.argv[2]
out_name = sys.argv[3]

my_calls = {}
ce_calls = {}

def load_my_file(my_report):
	with open(my_report) as mf:
		for l in mf:
			ll = l.strip().split()
			chrom = ll[0]
			spos = int(ll[1]) - 1
			epos = int(ll[2])
			support = int(ll[3])

			my_calls[(chrom, spos, epos)] = support

def load_ce_file(ce_report):
	with open(ce_report) as ce:
		for l in ce:
			ll = l.strip().split()
			chrom = ll[0]
			spos = int(ll[1])
			epos = int(ll[2])
			lll = ll[3].split('/')
			support = int(lll[1])
			#support = int(ll[3])

			ce_calls[(chrom, spos, epos)] = support

# c1 - c2
def compare(c1, c2, outf):
	for k in c1.keys():
		if k not in c2.keys():
			(c, s, e) = k
			print >>outf, '{}\t{}\t{}\t{}'.format(c, s, e, c1[k])

def union(c1, c2, outf):
	for k in c1.keys():
		if k in c2.keys():
			(c, s, e) = k
			print >>outf, '{}\t{}\t{}\t{}\t{}'.format(c, s, e, c1[k], c2[k])

load_my_file(my_report)
load_ce_file(ce_report)

us_adv = open('{}.us.adv'.format(out_name), 'w')
compare(my_calls, ce_calls, us_adv)
us_adv.close()

ce_adv = open('{}.ce.adv'.format(out_name), 'w')
compare(ce_calls, my_calls, ce_adv)
ce_adv.close()

both = open('{}.both'.format(out_name), 'w')
union(my_calls, ce_calls, both)
both.close()
