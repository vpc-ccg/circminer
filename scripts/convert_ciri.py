import sys
import re

def clean_list(l):
	while '' in l: l.remove('')
	if 'n/a' in l: l.remove('n/a')

def extract_bsj_reads(bsj_col):
	res = re.sub("<.*?::", "", bsj_col)
	res = re.sub("\(.*?\)", "$", res)
	res = re.sub(">", "$", res)
	res = res.split("$")
	clean_list(res)

	return res

def extract_ro_reads(ro_col):
	res = re.sub("##.*?&&", "$", ro_col)
	res = res.split("$")
	
	clean_list(res)

	return res

def add_count(fname):
	with open(fname) as fin:
		for l in fin:
			ll = l.strip().split()
			if ll[0] == 'BSJ':
				print(l.strip())
				continue

			bsj_col = ll[7]
			ro_col = ll[8]

			bsj_cnt = bsj_col.count('(') / 2
			ro_cnt = ro_col.count('&&')

			# print(bsj_col)
			# print(res)

			bsj_reads = extract_bsj_reads(bsj_col)
			ro_reads = extract_ro_reads(ro_col)

			all_reads = bsj_reads + ro_reads
			all_uniqe_reads = set(all_reads)

			print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(l.strip(), bsj_cnt, len(bsj_reads), ro_cnt, len(ro_reads), len(all_reads), len(all_uniqe_reads)))

def main():
	if len(sys.argv) != 2:
		print("Usage: {} ciri_final".format(sys.argv[0]))
		exit(1)

	fname = sys.argv[1]
	add_count(fname)

if __name__ == '__main__':
	main()
