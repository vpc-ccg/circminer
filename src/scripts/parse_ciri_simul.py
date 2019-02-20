import sys
import re

ciri = sys.argv[1]
my_out = sys.argv[2]

out = open(my_out, 'w')
ch = '-'
st = 0
en = 0
sup = 0

with open(ciri) as f:
	for l in f:
		if l[:3] == 'iso':
			continue

		if l[:2] == '**':
			continue

		if l[0] == '!':
			continue

		ll = re.split('\t| |:|\||\n', l)

		if ll[0] == '>':
			sup = ll[1]
		
		else:
			if ch != '-':
				print >>out, '{}\t{}\t{}\t{}'.format(ch, st, en, sup)
			ch = ll[0]
			st = ll[5]
			en = ll[6]
			sup = 0
if ch != '-':
	print >>out, '{}\t{}\t{}\t{}'.format(ch, st, en, sup)
