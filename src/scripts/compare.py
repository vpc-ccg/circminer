import sys

circ_exp_fname = sys.argv[1]
my_fname = sys.argv[2]

def get_sorted_rnames(fname):
	res = [ ]
	with open(fname) as f:
		for l in f:
			l = l.strip()
			res.append(int(l))
	res.sort()
	return res

circ_exp_rname = get_sorted_rnames(circ_exp_fname)
my_rname = get_sorted_rnames(my_fname)

print len(circ_exp_rname)
print len(my_rname)

i = 0
j = 0

comman = 0
exp = 0
mine = 0
fpf = open('fp.rname', 'w')

while i < len(circ_exp_rname) and j < len(my_rname):
	if circ_exp_rname[i] == my_rname[j]:
		comman += 1
		i += 1
		j += 1
	elif circ_exp_rname[i] < my_rname[j]:
		exp += 1
		i += 1
	else:
		mine += 1
		print >>fpf, my_rname[j]
		j += 1

if i < len(circ_exp_rname):
	exp += len(circ_exp_rname) - i

while j < len(my_rname):
	mine += 1
	print >>fpf, my_rname[j]
	j += 1


print "Comman:", comman
print "circExplorer:", exp
print "mine:", mine
