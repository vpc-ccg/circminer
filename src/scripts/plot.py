import sys
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.interpolate import spline
from scipy.interpolate import interp1d

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

def positive_mixed_plot(df, total_circ_cnt, pname):
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

	# plotting lines
	style="Simple,tail_width=0.5,head_width=4,head_length=8"
	kw = dict(arrowstyle=style, color="k")
	for t in tool_order:
		x = [ mydf.loc[mydf['Tools']==t]["Recall"].item(), mydf2.loc[mydf2['Tools']==t]["Recall"].item() ]
		y = [ mydf.loc[mydf['Tools']==t]["Precision"].item(), mydf2.loc[mydf2['Tools']==t]["Precision"].item() ]
		# plt.plot(x, y)
		px0 = (x[0], y[0])
		px1 = (x[1], y[1])
		pa = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=.5", **kw)
		plt.gca().add_patch(pa)

	# plot points
	light_col_list = [ "light " + x for x in col_list]
	sns.set_palette(sns.xkcd_palette(light_col_list))
	ax = sns.scatterplot(x="Recall", y="Precision", style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
							s=225, alpha=1.0, legend=False, data=mydf2)

	sns.set_palette(sns.xkcd_palette(col_list))
	ax = sns.scatterplot(x="Recall", y="Precision", style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
							s=75, alpha=1.0, data=mydf)

	# ax = sns.scatterplot(x="Recall", y="Precision", style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
	# 						marker="o", color='r', s=250, alpha=1.0, legend=False, data=mydf2)

	
	# for t in tool_order:
	# 	circle = plt.Circle((mydf2.loc[mydf2['Tools']==t]["Recall"], mydf2.loc[mydf2['Tools']==t]["Precision"]), .005, lw=2, color='r', fill=False)
	# 	ax.add_artist(circle)

	# x = [ mydf.loc[mydf['Tools']=='CircMarker']["Recall"].item(), mydf2.loc[mydf2['Tools']=='CircMarker']["Recall"].item() ]
	# y = [ mydf.loc[mydf['Tools']=='CircMarker']["Precision"].item(), mydf2.loc[mydf2['Tools']=='CircMarker']["Precision"].item() ]

	# y_curved = np.linspace(y[0], y[1], 300)
	# x_curved = interp1d(y, x, axis=0)(y_curved)
	# print(x, y)
	# plt.plot(x, y)
	# plt.plot(x_curved, y_curved)

	plt.savefig(pname+'.png')


def background_plot(df, pname):
	print(df)

	ax, tool_order, col_list = prepare_for_simul()

	ax = sns.barplot(x="Tools", y="#Detected", alpha=1.0, order=tool_order, data=df)
	ax.set(xlabel='Tool name', ylabel='#FP')

	plt.savefig(pname+'.png')

def simlu1_plot(args):
	xl = args[0]
	df = pd.read_excel(xl, skiprows=1)
	mydf = df.dropna()[:6]
	mydf2 = df.dropna()[14:20]

	final_df1 = pd.concat([mydf, mydf2])
	positive_mixed_plot(final_df1, 1000, 'pos-mix1-new')

	mydf3 = df.dropna()[7:13]
	mydf4 = df.dropna()[21:27]

	final_df2 = pd.concat([mydf3, mydf4])
	positive_mixed_plot(final_df2, 110128, 'pos-mix2-new')

	background_plot(df[18:24], 'background1')
	background_plot(df[27:33], 'background2')


def main():
	mode = sys.argv[1]
	args = sys.argv[2:]
	print(mode)
	print(args)
	if mode == 'tlen':
		tlen_plot(args)
	elif mode == 'simul1':
		simlu1_plot(args)

if __name__ == '__main__':
		main()
