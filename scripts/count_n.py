import sys

seq_fname = sys.argv[1]

high_th = 50

with open(seq_fname) as sf:
	rname = '?'
	seq = ''
	highn = 0
	for l in sf:
		if l[0] == '>':
			if rname != l[1:]:
				if rname != '?':
					print '{}\t{}\t{}\t{}'.format(rname.strip(), len(seq), seq.count('N'), 100.0 * seq.count('N') / len(seq))
					if 100.0 * seq.count('N') / len(seq) >= high_th:
						highn += 1
				rname = l[1:]
				seq = ''
			continue
		seq += l.strip()
	print '{}\t{}\t{}\t{}'.format(rname.strip(), len(seq), seq.count('N'), 100.0 * seq.count('N') / len(seq))
	if 100.0 * seq.count('N') / len(seq) >= high_th:
		highn += 1
	print '>= {}% N:\t{}'.format(high_th, highn)
