import sys
from collections import defaultdict
from itertools import chain

def rec_dict():
	return defaultdict(rec_dict)

def load_crfile(crfile, minsup=1):
	crdic = defaultdict(int)

	with open(crfile) as crf:
		for l in crf:
			ltok = l.strip().split()

			ch  = ltok[0]
			beg = int(ltok[1])
			end = int(ltok[2])
			sup = int(ltok[3])
			typ = ltok[4]
			reads = ltok[5]
			
			if sup >= minsup:
				crdic[(ch, beg, end)] = sup

	return crdic

def write_agg(sample_cnt, agg_dict, outfilename):
	with open(outfilename, 'w') as fout:
		fout.write('chr\tstart\tend')
		for i in range(sample_cnt):
			fout.write('\tS{}_sup'.format(i+1))
		fout.write('\n')

		for k in agg_dict.keys():
			(ch, beg, end) = k
			fout.write('{}\t{}\t{}'.format(ch, beg, end))
			for sup in agg_dict[k]:
				fout.write('\t{}'.format(sup))
			fout.write('\n')

def aggregate(crdics, minsup=2):
	all_keys = chain.from_iterable(crdics[i].keys() for i in range(len(crdics)))
	all_keys = list(set(all_keys))

	print('#circRNAs: {}'.format(len(all_keys)))

	agg_dict = defaultdict(list)
	for k in all_keys:
		(ch, beg, end) = k
		vals = [ dic[k] if k in dic.keys() else 0 for dic in crdics ]

		if all(val < minsup for val in vals):
			continue

		agg_dict[k] = vals
		
	return agg_dict

def usage():
	print('\nUsage:\npython {} circ_report_file1 circ_report_file2 circ_report_file3 ... output_file'.format(sys.argv[0]))

def main():
	args = sys.argv[1:]
	if len(args) < 1:
		usage()
		exit(1)
	
	outfilename = args[-1]
	crfilenames = args[:-1]
	crdics = []
	for crfile in crfilenames:
		crdics.append(load_crfile(crfile))

	agg_dict = aggregate(crdics)
	write_agg(len(crdics), agg_dict, outfilename)

if __name__ == '__main__':
	main()
