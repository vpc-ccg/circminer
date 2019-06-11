import sys
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import math
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad


def tlen_plot(args):
	f1 = args[0]
	f2 = args[1]
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

def prepare_for_simul():
	ax = plt.figure()
	sns.set()
	sns.set_context("paper")

	# sns.set_palette("Set2")
	# sns.set_palette(sns.hls_palette(7, h=.4), desat=1.0)
	# sns.set_palette(sns.color_palette("bright", n_colors=7))
	col_list = ["green", "magenta", "blue", "red", "orange", "turquoise"]
	sns.set_palette(sns.xkcd_palette(col_list))
	tool_order = ['circMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'SEGEMEHL (C)', 'CircMarker']

	return ax, tool_order, col_list

def circarrowdraw(x0, y0, kw, radius=1, aspect=1, direction=270, closingangle=-330, rotate_head = 0.0, *args):
	"""
	Circular arrow drawing. x0 and y0 are the anchor points.
	direction gives the angle of the circle center relative to the anchor
	in degrees. closingangle indicates how much of the circle is drawn
	in degrees with positive being counterclockwise and negative being
	clockwise. aspect is important to make the aspect of the arrow 
	fit the current figure. rotate_head is used to rotate the arrow head
	by increasing the y value of the arrow's tail coordinate.
	"""

	# Center of circle
	xc = x0 + radius * np.cos(direction * np.pi / 180)
	yc = y0 + aspect * radius * np.sin(direction * np.pi / 180)	
	# Draw circle
	if closingangle < 0:
		step = -1
	else:
		step = 1
	x = [xc + radius * np.cos((ang + 180 + direction) * np.pi / 180)
		for ang in np.arange(0, closingangle, step)]
	y = [yc + aspect * radius * np.sin((ang + 180 + direction) * np.pi / 180)
		for ang in np.arange(0, closingangle, step)]
	plt.plot(x, y, *args, color=kw['color'])

	# Draw arrow head
	arc_arrow_head = matplotlib.patches.FancyArrowPatch((x[-1], y[-1] + rotate_head),
											(x[0], y[0]), zorder = 10, **kw)
	plt.gca().add_patch(arc_arrow_head)

def positive_mixed_plot(df, total_circ_cnt, pname, dist_ratio=100, rot_head=-0.0003):
	mydf = df[:6]
	mydf2 = df[6:12]
	print(mydf)
	print('-------')
	print(mydf2)

	keys = list(mydf.columns.values)
	mydf['Recall'] = mydf['TP'] / total_circ_cnt
	mydf['Precision'] = mydf['TP'] / mydf['#Detected']

	keys = list(mydf2.columns.values)
	mydf2['Recall'] = mydf2['TP'] / total_circ_cnt
	mydf2['Precision'] = mydf2['TP'] / mydf2['#Detected']

	print(mydf)
	print('------')
	print(mydf2)

	ax, tool_order, col_list = prepare_for_simul()

	# plot points
	light_col_list = [ "light " + x for x in col_list]
	sns.set_palette(sns.xkcd_palette(light_col_list))
	ax = sns.scatterplot(x="Recall", y="Precision", style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
							s=120, alpha=1.0, legend=False, data=mydf2)

	sns.set_palette(sns.xkcd_palette(col_list))
	ax = sns.scatterplot(x="Recall", y="Precision", style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
							s=50, alpha=1.0, data=mydf)

	# plotting arrows
	ax2 = plt.gca()
	style = "Simple,tail_width=0.5,head_width=4,head_length=8"
	style2 = "Simple,tail_width=0.0,head_width=0.01,head_length=0.02"
	kw = dict(arrowstyle=style, color="k")
	kw2 = dict(arrowstyle=style2, color="k")
	for t in tool_order:
		x = [ mydf.loc[mydf['Tools']==t]["Recall"].item(), mydf2.loc[mydf2['Tools']==t]["Recall"].item() ]
		y = [ mydf.loc[mydf['Tools']==t]["Precision"].item(), mydf2.loc[mydf2['Tools']==t]["Precision"].item() ]
		# plt.plot(x, y)
		px0 = (x[0], y[0])
		px1 = (x[1], y[1])
		dist = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2)
		print(t, dist)
		if dist * dist_ratio < 0.1**4:
			bb = plt.gca().axis()
			asp = (bb[3] - bb[2]) / (bb[1] - bb[0])
			circarrowdraw(x[0], y[0], kw, radius=.01, aspect=asp, direction=-90, closingangle=-345, rotate_head=rot_head)

		elif dist < 0.1**3:
			bb = plt.gca().axis()
			# asp = (bb[3] - bb[2]) / (bb[1] - bb[0])
			# circarrowdraw(x[0], y[0], kw, radius=.01, aspect=asp, direction=-90, closingangle=-345, rotate_head=rot_head)

			# pa = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=6.4", **kw)
			# ax2.add_patch(pa)

			tmp0 = [x[1] - .00001, y[1] + .0000001]

			pa1 = matplotlib.patches.FancyArrowPatch(tmp0, px1, connectionstyle="arc3,rad=0.", **kw)
			pa2 = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=6.4", **kw2)
			plt.gca().add_patch(pa1)
			plt.gca().add_patch(pa2)
		else:
			pa = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=.4", **kw)
			plt.gca().add_patch(pa)


	plt.gca().set_ylim(top=1.003)
	plt.gca().set_xlim(left=0.65)

	

	plt.savefig(pname+'.png')


def background_plot(df, pname):
	print(df)

	ax, tool_order, col_list = prepare_for_simul()

	ax = sns.barplot(x="Tools", y="#Detected", alpha=1.0, order=tool_order, data=df)
	ax.set(xlabel='Method name', ylabel='Number of False Positives')

	plt.savefig(pname+'.png')

def simlu1_plot(args):
	xl = args[0]
	df = pd.read_excel(xl, skiprows=1)
	mydf = df.dropna()[:6]
	mydf2 = df.dropna()[14:20]

	final_df1 = pd.concat([mydf, mydf2])
	positive_mixed_plot(final_df1, 1000, 'pos-mix1-n', dist_ratio=.01, rot_head=-0.0003)

	mydf3 = df.dropna()[7:13]
	mydf4 = df.dropna()[21:27]

	final_df2 = pd.concat([mydf3, mydf4])
	positive_mixed_plot(final_df2, 110128, 'pos-mix2-n', dist_ratio=1, rot_head=-0.0001)

	background_plot(df[18:24], 'background1')
	background_plot(df[27:33], 'background2')

def resource_plot(args):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "real(1_thread)")
	# df = pd.parse("real(1_thread)")
	print(df)

	x = ['CircMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'SEGEMEHL', 'CircMarker']
	y = [10800, 28218, 209823, 111934, 345040, 221867]

	ax, tool_order, col_list = prepare_for_simul()

	# ax = sns.barplot(x="Tools", y="#Detected", alpha=1.0, data=mydf)
	# ax.set(xlabel='Method name', ylabel='Time')

	fig = plt.bar(x, y, color=col_list)
	# plt.legend(fig, ['First','Second','Third'], loc = "upper left", title = "cat")

	plt.ylabel('Time (s)')

	plt.savefig('time.png')

def main():
	mode = sys.argv[1]
	args = sys.argv[2:]
	print(mode)
	print(args)
	if mode == 'tlen':
		tlen_plot(args)
	elif mode == 'simul1':
		simlu1_plot(args)
	elif mode == 'resource':
		resource_plot(args)

if __name__ == '__main__':
		main()
