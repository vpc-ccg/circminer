import sys

fq_file = sys.argv[1]
read_name_file = sys.argv[2]
save_file = sys.argv[3]

read_names = []
prefix = ''

with open(read_name_file) as rnf:
	for rn in rnf:
		rn = rn.strip().split('.')
		prefix = rn[0]
		read_names.append(int(rn[1]))
		#read_names.append(int(rn))

sorted_read_names = sorted(read_names)
print sorted_read_names

with open(save_file, 'w') as sf:
	with open(fq_file, 'r') as fqf:
		ind = 0
		j = 0
		for l in fqf:
			j += 1
			if j == 1:
				rname = l
				if '{}.{}'.format(prefix, sorted_read_names[ind]) in l:
					ind += 1
					sel = True
				else:
					sel = False
			elif j == 2:
				seq = l
			elif j == 3:
				comment = l
			elif j == 4:
				qual = l
				if sel:
					sf.write('{}{}{}{}'.format(rname, seq, comment, qual))
				if (ind >= len(sorted_read_names)):
					break
				j = 0
