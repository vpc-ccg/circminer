import sys

gtf_fname = sys.argv[1]

genes = {}
enums = []

max_enum = 0

with open(gtf_fname) as gf:
	for l in gf:
		if l[0] == '#':
			continue
		ll = l.strip().split()
		if ll[2] == "gene":
			enums.append(max_enum)
			max_enum = 0
		if ll[2] != "exon":
			continue
		exon_num = int(ll[13][1:-2])
		max_enum = max(max_enum, exon_num)

enums.append(max_enum)

print enums
print len(enums)
print sum(enums)
