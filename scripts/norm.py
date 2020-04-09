import sys

cnt = {}
with open(sys.argv[1]) as f:
	for l in f:
		ll = l.strip().split()
		cnt[ll[0]] = int(ll[1])


with open(sys.argv[2]) as f:
	for l2 in f:
		l = l2.strip()
		ll = l.split()

		# header
		if ll[0] == 'chr':
			mstr = '_'.join(ll[:3]) + '\t' + '\t'.join(ll[3:])
			print(mstr)
			samp = []
			for i in range(3, len(ll)):
				samp.append(ll[i])
			continue

		mstr = ''
		#for v in ll[:3]:
		#	mstr += "{}_".format(v)
		mstr += '_'.join(ll[:3]) 


		i = 0
		for v in ll[3:]:
			val = 1000000.0 * int(v) / cnt[samp[i]] 
			if val > 0.5:
				val = 0.5
			i += 1
			mstr += "\t{}".format(val)

		print(mstr)

