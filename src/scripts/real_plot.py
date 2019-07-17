import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from pylab import *
import seaborn as sns

def prepare_for_simul():
	ax = plt.figure()
	sns.set()
	sns.set_context("paper")

	#sns.set_palette("Set2")
	#sns.set_palette(sns.hls_palette(7, h=.4), desat=1.0)
	#sns.set_palette(sns.color_palette("bright", n_colors=7))
	col_list = ["blue", "red", "magenta", "orange", "turquoise", "green"]
	sns.set_palette(sns.xkcd_palette(col_list))
	tool_order = ['circMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'SEGEMEHL (C)', 'CircMarker']

	return ax, tool_order, col_list

ax, tool_order, col_list = prepare_for_simul()
x=  [12.00842862, 12.56843119, 11.5112585, 16.04580231, 12.49909844, 11.88950396]
y=  [0.6791, 0.6304, 0.6817, 0.1626, 0.5457, 0.6921]#, 0.9650, 0.8670, 0.69, 0.944, 0.985]
label=["KNIFE", "CIRI", "CIRCExplorer", "SEGEMEHL", "CircMarker", "CircMiner"]
#label = tool_order
#s=[190, 176, 206, 100, 148, 330]#173, 279, 55, 67, 300]
ss=[76, 73, 79, 55, 67, 100]
s = [x**2/20 for x in ss]
n=[0.6, 0.5, 0.5, 0.3, 0.5, 0.6]

fig, ax = plt.subplots()
for i in range( len(x)):
	ax.scatter( x[i],y[i], label=label[i], s=s[i], norm=n[i], alpha=0.6)
plt.ylim(0.15, 0.75)
plt.xlim( 11, 17 )#, 68000)
plt.title('Comparison of RNaseR- and RNaseR+ for the same patient')#Simulation datasets of 1000 transcripts')
plt.xlabel('Number of Detected circRNAs in RNaseR- (log 2)')
plt.ylabel('Ratio of not depleted circRNAs in RNaseR+')
lgnd = ax.legend()

# Set same size for points in legend
for handle in lgnd.legendHandles:
	handle.set_sizes([50.0])
ax.get_legend().remove()

# Create the second legend and add the artist manually.
from matplotlib.legend import Legend
leg = Legend(ax, lgnd.legendHandles, label, loc='lower right', frameon=False)


# Here we create a legend:
# we'll plot empty lists with the desired size and label
for area in [100, 300, 900]:
    plt.scatter([], [], c='k', alpha=0.3, s=area,
                label=str(area))
plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='')

ax.add_artist(leg)

ax.grid(True)

plt.savefig('real.Hs68.2leg.2.png')
