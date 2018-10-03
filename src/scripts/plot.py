import matplotlib.pyplot as plt
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

tlen1 = []
cnt1 = []

tlen2 = []
cnt2 = []

with open(f1) as f:
	for l in f:
		ll = l.strip().split()
		tlen1.append(ll[1])
		cnt1.append(ll[0])

with open(f2) as f:
	for l in f:
		ll = l.strip().split()
		tlen2.append(ll[1])
		cnt2.append(ll[0])

plt.plot(tlen1, cnt1, color='k')
plt.plot(tlen2, cnt2, color='g')

plt.save('tlen.png')
