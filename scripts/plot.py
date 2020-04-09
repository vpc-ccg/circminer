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

def prepare_for_stack():
	ax = plt.figure()
	sns.set()
	sns.set_context("paper")
	#sns.set_context("talk")

	sns.set_palette("Set2")
	sns.set_style("whitegrid")
	# sns.set_palette(sns.hls_palette(7, h=.4), desat=1.0)
	# sns.set_palette(sns.color_palette("bright", n_colors=7))
	# col_list = ["green", "magenta", "blue", "red", "orange", "turquoise"]
	# sns.set_palette(sns.xkcd_palette(col_list))

	plt.gca().xaxis.grid(False)

	return ax

def prepare_for_simul(typ=1, font_scale=1.6):
	ax = plt.figure()
	#sns.set()
	#sns.set(font_scale=font_scale)
	sns.set_context("paper")
	#sns.set_context("talk")
	sns.set_style("whitegrid")

	# sns.set_palette("Set2")
	# sns.set_palette(sns.hls_palette(7, h=.4), desat=1.0)
	# sns.set_palette(sns.color_palette("bright", n_colors=7))
	col_list = ["green", "magenta", "blue", "red", "turquoise", "brown", "orange"]
	sns.set_palette(sns.xkcd_palette(col_list))
	if typ == 1:
		tool_order = ['CircMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'CircMarker', 'CIRCexplorer2', 'SEGEMEHL (C)' ]
	elif typ == 3:
		tool_order = ['CircMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'CircMarker', 'CIRCexplorer2']
	elif typ == 2:
		tool_order = ['CircMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'CircMarker', 'CIRCexplorer2', 'SEGEMEHL']
	elif typ == 4:
		tool_order = ['CircMiner', 'CIRCexplorer', 'CIRCexplorer2', 'KNIFE', 'CIRI', 'CircMarker', 'SEGEMEHL', 'circRNA_finder', 'DCC', 'find_circ', 'PTESFinder', 'NCLscan']
		col_list = ["green", "magenta", "brown", "blue", "red", "turquoise", "orange", "navy" , "maroon", "pink", "violet", "gold"]
	elif typ == 5:
		tool_order = ['CircMiner', 'CIRCexplorer', 'CIRCexplorer2', 'KNIFE', 'CIRI', 'CircMarker', 'SEGEMEHL (C)', 'circRNA_finder', 'DCC', 'find_circ', 'PTESFinder', 'NCLscan']
		col_list = ["green", "magenta", "brown", "blue", "red", "turquoise", "orange", "navy" , "maroon", "pink", "violet", "gold"]

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
	plt.plot(x, y, *args, color=kw['color'], alpha=kw['alpha'])

	# Draw arrow head
	arc_arrow_head = matplotlib.patches.FancyArrowPatch((x[-1], y[-1] + rotate_head),
											(x[0], y[0]), zorder = 10, **kw)
	plt.gca().add_patch(arc_arrow_head)

def positive_mixed_plot(df, total_circ_cnt, pname, dist_ratio=100, rot_head=-0.0003):
	mydf = df[:12]
	mydf2 = df[12:24]
	print(mydf)
	print('-------')
	print(mydf2)

	keys = list(mydf.columns.values)
	mydf['Recall'] = mydf['TP'] / total_circ_cnt
	mydf['Precision'] = mydf['TP'] / mydf['#Detected']

	keys = list(mydf2.columns.values)
	mydf2['Recall'] = mydf2['TP'] / total_circ_cnt
	mydf2['Precision'] = mydf2['TP'] / mydf2['#Detected']


	#mydf.drop(mydf[mydf.Tools == "SEGEMEHL (C)"].index, inplace=True)
	#mydf2.drop(mydf2[mydf2.Tools == "SEGEMEHL (C)"].index, inplace=True)

	print(mydf)
	print('------')
	print(mydf2)

	ax, tool_order, col_list = prepare_for_simul(5)

	print(tool_order)
	print(col_list)

	#markers = [ 'o', 'x', 's', '+', 'd', '^', 'v', '>', '<', 'p', 'h' ]
	markers = [ "o", "x", "s", "+", "d", "^", "v", ">", "<", "p", "h", "D" ]

	# plot points
	light_col_list = [ "light " + x for x in col_list]
	sns.set_palette(sns.xkcd_palette(light_col_list))
	#ax = sns.scatterplot(x="Recall", y="Precision", marker=markers, style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
	ax = sns.scatterplot(x="Recall", y="Precision", marker="X", style_order=tool_order, hue="Tools", hue_order=tool_order,
							s=100, alpha=1.0, legend=False, data=mydf2)

	sns.set_palette(sns.xkcd_palette(col_list))
	#ax = sns.scatterplot(x="Recall", y="Precision", marker=markers, style="Tools", style_order=tool_order, hue="Tools", hue_order=tool_order,
	ax = sns.scatterplot(x="Recall", y="Precision", markers=markers, style_order=tool_order, hue="Tools", hue_order=tool_order,
							s=20, alpha=1.0, data=mydf)

	# plotting arrows
	ax2 = plt.gca()
	style = "Simple,tail_width=0.5,head_width=4,head_length=8"
	style2 = "Simple,tail_width=0.0,head_width=0.01,head_length=0.02"
	arrow_alpha = 0.95
	kw = dict(arrowstyle=style, color="gray", alpha=arrow_alpha)
	kw2 = dict(arrowstyle=style2, color="gray", alpha=arrow_alpha)
	for t in tool_order:
		x = [ mydf.loc[mydf['Tools']==t]["Recall"].item(), mydf2.loc[mydf2['Tools']==t]["Recall"].item() ]
		y = [ mydf.loc[mydf['Tools']==t]["Precision"].item(), mydf2.loc[mydf2['Tools']==t]["Precision"].item() ]
		# plt.plot(x, y)
		px0 = (x[0], y[0])
		px1 = (x[1], y[1])
		dist = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2)
		print(t, dist)
		#if dist * dist_ratio < 0.1**4 or t == 'CIRCexplorer2':
		if dist * dist_ratio < 0.1**4:
			bb = plt.gca().axis()
			asp = (bb[3] - bb[2]) / (bb[1] - bb[0])
			radius = .01 if t != 'CIRCexplorer2' else .0075
			circarrowdraw(x[0], y[0], kw, radius=radius, aspect=asp, direction=-90, closingangle=-345, rotate_head=rot_head)

		elif dist * dist_ratio < 0.1**3 or t in ['find_circ']:
			bb = plt.gca().axis()
			# asp = (bb[3] - bb[2]) / (bb[1] - bb[0])
			# circarrowdraw(x[0], y[0], kw, radius=.01, aspect=asp, direction=-90, closingangle=-345, rotate_head=rot_head)

			# pa = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=6.4", **kw)
			# ax2.add_patch(pa)

			tmp0 = [x[1] - .00001, y[1] + .0000001]

			if t in ['CIRCexplorer2', 'DCC'] and total_circ_cnt == 1000:
				pa1 = matplotlib.patches.FancyArrowPatch(tmp0, px1, connectionstyle="arc3,rad=-0.2", **kw)
				pa2 = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=10.4", **kw2)
			elif t in ['CIRCexplorer2', 'DCC'] and total_circ_cnt > 1000:
				tmp0 = [x[0] - .00001, y[0] + .0000001]
				pa1 = matplotlib.patches.FancyArrowPatch(tmp0, px0, connectionstyle="arc3,rad=-0.15", **kw)
				pa2 = matplotlib.patches.FancyArrowPatch(px1, px0, connectionstyle="arc3,rad=10.4", **kw2)
			elif t in ['find_circ']:
				if total_circ_cnt == 1000:
					pa1 = matplotlib.patches.FancyArrowPatch(tmp0, px1, connectionstyle="arc3,rad=-0.3", **kw)
					pa2 = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=.9", **kw2)
				else:
					pa1 = matplotlib.patches.FancyArrowPatch(tmp0, px1, connectionstyle="arc3,rad=0.0", **kw)
					pa2 = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=2.9", **kw2)

			else:
				pa1 = matplotlib.patches.FancyArrowPatch(tmp0, px1, connectionstyle="arc3,rad=0.", **kw)
				pa2 = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=6.4", **kw2)
			plt.gca().add_patch(pa1)
			plt.gca().add_patch(pa2)
		else:
			pa = matplotlib.patches.FancyArrowPatch(px0, px1, connectionstyle="arc3,rad=.4", **kw)
			plt.gca().add_patch(pa)


	plt.gca().set_ylim(top=1.004)
	plt.gca().set_xlim(left=0.54)

	# Create first legend and add artist manually
	handles, labels = ax.get_legend_handles_labels()
	labels[7] = 'SEGEMEHL'
	leg = plt.legend(handles=handles[1:], labels=labels[1:], loc='lower left', bbox_to_anchor=(1.0, 0.0), title='Methods')
	leg._legend_box.align = "left"
	ax.add_artist(leg)

	# Create the second legend and add the artist manually
	# from matplotlib.legend import Legend
	# leg = Legend(ax, lgnd.legendHandles, tool_order, loc='lower right', frameon=False)

	mk = [ "o","o","o","o","o","o", "X", "X", "X"]
	h = [plt.plot([],[], color='k', marker=mk[i], ms=i, ls="")[0] for i in range(5, 9, 3)]
	# lgd = plt.legend(handles=h, labels=['Positive', 'Mixed'], loc='lower left', title="Dataset", bbox_to_anchor=(1.0, 0.4))
	lgd = plt.legend(handles=h, labels=['Positive Control', 'Diluted'], loc='lower left', title="Dataset", bbox_to_anchor=(1.0, 0.7))
	lgd._legend_box.align = "left"


	#plt.savefig(pname+'.svg', bbox_extra_artists=(leg,lgd), bbox_inches='tight', dpi=1200)
	#plt.savefig(pname+'.png', bbox_extra_artists=(leg,lgd), bbox_inches='tight', dpi=800)
	plt.savefig(pname+'.pdf', bbox_extra_artists=(leg,lgd), bbox_inches='tight')


def false_pos_plot(df, pname, col):
	plt.gca().xaxis.grid(False)
	print(df)

	ax, tool_order, col_list = prepare_for_simul()

	ax = sns.barplot(x="Tools", y=col, alpha=1.0, order=tool_order, data=df)
	ax.set(xlabel='Method name', ylabel='Number of False Positives')

	ax.get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
	plt.grid(b=True, which='major', color='gray', linewidth=1.0)
	plt.grid(b=True, which='minor', color='gray', linewidth=0.5)
	plt.gca().xaxis.grid(False)

	plt.savefig(pname+'.svg', dpi=1200)

def simul1_plot(args):
	xl = args[0]
	df = pd.read_excel(xl, skiprows=1)
	mydf = df.dropna()[:12]
	mydf2 = df.dropna()[26:38]

	final_df1 = pd.concat([mydf, mydf2])
	positive_mixed_plot(final_df1, 1000, 'pos-mix-1k', dist_ratio=.1, rot_head=-0.0003)

	mydf3 = df.dropna()[13:25]
	mydf4 = df.dropna()[39:51]

	final_df2 = pd.concat([mydf3, mydf4])
	positive_mixed_plot(final_df2, 110128, 'pos-mix-110k', dist_ratio=1.0, rot_head=-0.0001)

	mydf2['FP'] = mydf2['#Detected'] - mydf2['TP']
	mydf4['FP'] = mydf4['#Detected'] - mydf4['TP']

	#false_pos_plot(mydf2, 'fp_rate1k', 'FP')
	#false_pos_plot(mydf4, 'fp_rate110k', 'FP')

	#false_pos_plot(df[18:24], 'background1k', '#Detected')
	#false_pos_plot(df[27:33], 'background110k', '#Detected')

def confidence_plot(df,samp):
	df['size'] = (1.3 ** (df['Top 100 not depleted'] * df['Top 100 enriched'] / 100.0)).astype('float64')
	df['RNaseR-_log'] = np.log2(df['RNaseR-'].astype('float64'))
	print(df)

	ax, tool_order, col_list = prepare_for_simul(2)
	# font setting
	plt.rcParams['font.family'] = "Liberation Serif"
	plt.rcParams['mathtext.fontset'] = 'dejavuserif'
	plt.rcParams['font.size'] = 16
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['legend.fontsize'] = 9
	plt.rcParams['xtick.labelsize'] = 14
	plt.rcParams['ytick.labelsize'] = 14

	font = {'family': 'serif',
			'size': 12,
			'weight': 'bold',
        # 'color':  'darkred',
        }


	#ax = sns.scatterplot(x="RNaseR-", y="Percentage %", hue="Tool", size="Top 100 not depleted", data=df, sizes=(20,120), linewidth=0)
	ax = sns.scatterplot(x="RNaseR-_log", y="Percentage %", hue="Tool", hue_order=tool_order, size="size", sizes=(60, 200), alpha=0.6, data=df)

	#plt.gca().set_ylim(bottom=95.0)

	handles, labels = ax.get_legend_handles_labels()
	print(handles)
	del(handles[-1])
	del(labels[-1])
	labels[0] = 'Methods'
	labels[7] = 'Concordane of\nhigh-confidence calls'
	labels[8] = 'Low'
	labels[9] = 'Medium'
	labels[10] = 'High'
	print(labels)
	leg = plt.legend(handles=handles[:], labels=labels[:], loc='upper right', bbox_to_anchor=(1.0, 1.0))
	leg._legend_box.align = "right"
	ax.add_artist(leg)

	plt.title("Comparison of calls on RNaseR- and RNaseR+ datasets", fontdict=font)
	ax.set_ylabel("Ratio of self-consistent circRNA in RNaseR+")
	ax.set_xlabel("Number of detected circRNAs in RNaseR-")

	
	#plt.savefig('{}_perf.pdf'.format(samp), bbox_inches='tight')
	plt.savefig('{}_perf.png'.format(samp), bbox_inches='tight', dpi=900)

def real1_plot(args):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "real(results)")
	#print(df)
	
	df_hs68 = df[21:27]
	#confidence_plot(df_hs68, 'Hs68')
	#top_calls_plot(df_hs68, 'HS68')

	df_hela = df[29:35]
	#confidence_plot(df_hela, 'HeLa')

	hs68_cpu, hs68_mem = prepare_resource_real_data(args, 'hs68')
	hela_cpu, hela_mem = prepare_resource_real_data(args, 'hela')

	#df_list = [ hs68_cpu, hs68_mem, hela_cpu, hela_mem ]
	#resource_real_plot(df_list)
	
	hs68_cpu['cell line'] = 'Hs68'
	hs68_mem['cell line'] = 'Hs68'

	hela_cpu['cell line'] = 'HeLa'
	hela_mem['cell line'] = 'HeLa'

	df_cpu = pd.concat([hs68_cpu, hela_cpu])
	df_mem = pd.concat([hs68_mem, hela_mem])

	resource_real_plot2(df_cpu, 'time')
	resource_real_plot2(df_mem, 'mem')

def top_calls_plot(df, sample):
	df["Top 100 not enriched not depleted"] = df["Top 100 not depleted"] - df["Top 100 enriched"]
	print(df)	

	ax, tool_order, col_list = prepare_for_simul(2)
	# Create correct bars
	#ax1 = plt.bar(r, meandf['Correct'], edgecolor='gray', width=bar_width, label="Correct mapping", yerr=stddf['AllCorrectPerc'])
	# Create incorrect bars
	#ax2 = plt.bar(r, meandf['Incorrect'], bottom=meandf['Correct'], edgecolor='gray', width=bar_width, label="Incorrect mapping", yerr=stddf['AllIncorrectPerc'])
	#ax = sns.barplot(x = 'Tool', y = 'Top 100 enriched', alpha=1.0, data=df, order=tool_order)
	
	ax = sns.barplot(x = 'Tool', y = 'Top 100 not depleted', alpha=1.0, data=df, order=tool_order, color='b')
	ax = sns.barplot(x = 'Tool', y = 'Top 100 enriched', alpha=1.0, data=df, order=tool_order, color='g')

	sns.despine(left=True, right=True, top=True)

	matplotlib.pyplot.sca(ax)
	plt.xticks(rotation=90)

	#ax.get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
	#plt.grid(b=True, which='major', color='gray', linewidth=1.0)
	#plt.grid(b=True, which='minor', color='gray', linewidth=0.5)
	#plt.gca().xaxis.grid(False)
	
	ax.set(xlabel='Methods', ylabel='Number of calls')
	ax.set_ylim(top=100)

	matplotlib.style.use('seaborn')
	name_to_color = {
					'Enriched':   'green',
					'Not depleted':   'blue',
					}

	patches = [matplotlib.patches.Patch(color=v, label=k) for k,v in name_to_color.items()]
	matplotlib.pyplot.legend(handles=patches,loc = "upper right", bbox_to_anchor=(1,1) )
			

	## There is no labels, need to define the labels
	#legend_labels  = ['Enriched', 'Not depleted']

	## Create the legend patches
	#legend_patches = [matplotlib.patches.Patch(color=C, label=L) for
    #              C, L in zip([item.get_facecolor() for item in Boxes],
    #                          legend_labels)]

	## Plot the legend
	#plt.legend(handles=legend_patches)

	plt.savefig('{}_top100_enriched_notdepleted.png'.format(sample), bbox_inches='tight', dpi=600)

def prepare_resource_real_data(args, sample):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "real(1_thread)")
	#df = df.iloc[:,:-3]
	print(df)

	if sample == 'hela':
		df_cput = df.iloc[[5,6,7], 2:].T
	else:
		df_cput = df.iloc[[5,8,9], 2:].T
	df_cput.columns = ['Tool', 'R-', 'R+']
	df_cput = df_cput.reset_index().replace({'index': {'STAR': 'CIRCexplorer'}}).groupby('index', sort=False).sum()

	df_cput.loc[df_cput["Tool"] == "STARCIRCexplorer", "Tool"] = "CIRCexplorer" 
	df_cput.loc[df_cput["Tool"] == "v 0.2", "Tool"] = "CircMiner" 

	#df_cput['RNaseR-'] = df_cput['R-'] / 3600.0
	#df_cput['RNaseR+'] = df_cput['R+'] / 3600.0
	df_cput['RNaseR-'] = df_cput['R-'] / 60.0
	df_cput['RNaseR+'] = df_cput['R+'] / 60.0

	df_cput = pd.melt(df_cput, id_vars=['Tool', 'R-', 'R+'], value_vars=['RNaseR-', 'RNaseR+'], var_name='dataset', value_name='time')
	print(df_cput)


	if sample == 'hela':
		df_mem = df.iloc[[11,12,13], 2:].T
	else:
		df_mem = df.iloc[[11,14,15], 2:].T
	df_mem.columns = ['Tool', 'R-', 'R+']
	df_mem = df_mem.reset_index().replace({'index': {'STAR': 'CIRCexplorer'}}).groupby('index', sort=False).max()

	df_mem.loc[df_mem["Tool"] == "STAR", "Tool"] = "CIRCexplorer" 
	df_mem.loc[df_mem["Tool"] == "v 0.2", "Tool"] = "CircMiner" 

	df_mem['RNaseR-'] = df_mem['R-'] / (1024.0 * 1024.0)
	df_mem['RNaseR+'] = df_mem['R+'] / (1024.0 * 1024.0)

	df_mem = pd.melt(df_mem, id_vars=['Tool', 'R-', 'R+'], value_vars=['RNaseR-', 'RNaseR+'], var_name='dataset', value_name='mem')
	print(df_mem)

	return df_cput, df_mem

def resource_real_plot2(df, typ):
	ax, tool_order, col_list = prepare_for_simul(4)
	x = tool_order
	
	#ax = sns.barplot(x = 'Tool', y = typ, hue = 'dataset', alpha=1.0, data=df, order=tool_order)
	ax = sns.catplot(x = 'Tool', y = typ, hue = 'dataset', col='cell line', alpha=1.0, data=df, order=tool_order, kind='bar', legend=True, sharey=False)
	
	if typ == 'time':
		ax.set(ylabel='Time (min)')
	else:
		ax.set(ylabel='Memory usage (GB)')

	sns.despine(left=True, right=True, top=True)

	#matplotlib.pyplot.sca(ax)
	for i in range(2):
		xlabels = ax.axes[0,i].get_xticklabels()
		ax.axes[0,i].set_xticklabels(xlabels, rotation=90)
		ax.axes[0,i].get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
		ax.axes[0,i].grid(b=True, which='major', color='gray', linewidth=1.0)
		ax.axes[0,i].grid(b=True, which='minor', color='gray', linewidth=0.5)
		ax.axes[0,i].xaxis.grid(False)


	ax.set(xlabel='Methods')

	#if sample == 'hela' and t == 'mem':
	#	ax.fig.get_axes()[0].legend(loc='lower left', bbox_to_anchor=(1,0.5), title='dataset')

	plt.savefig('{}_combined.pdf'.format(typ), bbox_inches='tight')
		
def resource_real_plot(df_list):
	ax, tool_order, col_list = prepare_for_simul(2)

	#fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, squeeze=False, figsize=(10,5))
	fig, axes = plt.subplots(1, 4)

	res_type = [ 'time', 'mem', 'time', 'mem' ]

	for t, df, ax_ind in zip(res_type, df_list, range(len(res_type))):
		#ax, tool_order, col_list = prepare_for_simul(4)
		x = tool_order
		
		ax = sns.catplot(x = 'Tool', y = t, hue = 'dataset', alpha=1.0, data=df, order=tool_order, kind='bar', legend=False, ax=axes[ax_ind])
		
		#ax = axes[ax_ind]


		#if t == 'time':
		#	ax.set(ylabel='Time (min)')
		#else:
		#	ax.set(ylabel='Memory usage (GB)')

		#sns.despine(left=True, right=True, top=True)

		##matplotlib.pyplot.sca(ax)
		#plt.xticks(rotation=90)

		##axes[0, 0].axes[0,0].get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
		##ax.axes[0, 0].get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
		##plt.grid(b=True, which='major', color='gray', linewidth=1.0)
		##plt.grid(b=True, which='minor', color='gray', linewidth=0.5)
		##plt.gca().xaxis.grid(False)
		#
		#ax.set(xlabel='Methods')

		##if sample == 'hela' and t == 'mem':
		##	ax.fig.get_axes()[0].legend(loc='lower left', bbox_to_anchor=(1,0.5), title='dataset')

		##plt.savefig('{}_{}_combined.pdf'.format(sample, t), bbox_inches='tight')
		
	plt.savefig('real_time_mem.pdf', bbox_inches='tight')


def resource_plot_single(args):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "simul(1_thread)")
	# df = pd.parse("real(1_thread)")
	print(df)


	df_cput = df.iloc[[7, 12,13], 2:].T
	df_cput.columns = ['Tool', '1k', '110k']
	df_cput = df_cput.reset_index().replace({'index': {'STAR': 'CIRCexplorer'}}).groupby('index', sort=False).sum()

	df_cput.loc[df_cput["Tool"] == "STARCIRCexplorer", "Tool"] = "CIRCexplorer" 
	df_cput.loc[df_cput["Tool"] == "v 0.2", "Tool"] = "CircMiner" 

	df_cput['1k-hr'] = df_cput['1k'] / 3600.0
	df_cput['110k-hr'] = df_cput['110k'] / 3600.0

	print(df_cput)

	df_mem = df.iloc[[15, 20,21], 2:].T
	df_mem.columns = ['Tool', '1k', '110k']
	df_mem = df_mem.reset_index().replace({'index': {'STAR': 'CIRCexplorer'}}).groupby('index', sort=False).max()

	df_mem.loc[df_mem["Tool"] == "STAR", "Tool"] = "CIRCexplorer" 
	df_mem.loc[df_mem["Tool"] == "v 0.2", "Tool"] = "CircMiner" 

	df_mem['1k-gb'] = df_mem['1k'] / (1024.0 * 1024.0)
	df_mem['110k-gb'] = df_mem['110k'] / (1024.0 * 1024.0)

	print(df_mem)


	for r in ['hr', 'gb']:
		for s in ['1k', '110k']:
			ax, tool_order, col_list = prepare_for_simul(2)

			#fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, squeeze=False)
			fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, squeeze=False, figsize=(10,5))

			#x = ['CircMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'SEGEMEHL', 'CircMarker']
			x = tool_order

			ax = sns.barplot(x = 'Tool', y = '{}-{}'.format(s, r), alpha=1.0, data=df_cput if (r == 'hr') else df_mem, order=tool_order, ax=axes[0,0])

			#ax.set_yscale('log')

			#fig = plt.bar(x, y, color=col_list)
			# plt.legend(fig, ['First','Second','Third'], loc = "upper left", title = "cat")

			if r == 'hr':
				axes[0,0].set_ylabel('Time (hr)')
			else:
				axes[0,0].set_ylabel('Memory usage (GB)')

			sns.despine(left=True, right=True, top=True)

			for ax in fig.axes:
				matplotlib.pyplot.sca(ax)
				plt.xticks(rotation=90)

				ax.get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
				plt.grid(b=True, which='major', color='gray', linewidth=1.0)
				plt.grid(b=True, which='minor', color='gray', linewidth=0.5)
				plt.gca().xaxis.grid(False)
	
			if s == '110k':
				axes[0,0].set(xlabel='Methods\nDiluted Dataset1 (DC1)')
			else:
				axes[0,0].set(xlabel='Methods\nDiluted Dataset2 (DC2)')

			#plt.savefig('time_diluted_1k.pdf')
			#plt.savefig('time_diluted_110k.pdf')
			#plt.savefig('mem_1k.pdf')
			#plt.savefig('mem_110k.pdf')
			#plt.savefig('dc_time_mem.pdf', bbox_inches='tight')
			plt.savefig('{}_{}.png'.format(s, r), bbox_inches='tight', dpi=500)

def resource_plot(args):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "simul(1_thread)")
	# df = pd.parse("real(1_thread)")
	print(df)


	df_cput = df.iloc[[7, 12,13], 2:].T
	df_cput.columns = ['Tool', '1k', '110k']
	df_cput = df_cput.reset_index().replace({'index': {'STAR': 'CIRCexplorer'}}).groupby('index', sort=False).sum()

	df_cput.loc[df_cput["Tool"] == "STARCIRCexplorer", "Tool"] = "CIRCexplorer" 
	df_cput.loc[df_cput["Tool"] == "v 0.2", "Tool"] = "CircMiner" 

	df_cput['1k-hr'] = df_cput['1k'] / 3600.0
	df_cput['110k-hr'] = df_cput['110k'] / 3600.0

	print(df_cput)

	df_mem = df.iloc[[15, 20,21], 2:].T
	df_mem.columns = ['Tool', '1k', '110k']
	df_mem = df_mem.reset_index().replace({'index': {'STAR': 'CIRCexplorer'}}).groupby('index', sort=False).max()

	df_mem.loc[df_mem["Tool"] == "STAR", "Tool"] = "CIRCexplorer" 
	df_mem.loc[df_mem["Tool"] == "v 0.2", "Tool"] = "CircMiner" 

	df_mem['1k-gb'] = df_mem['1k'] / (1024.0 * 1024.0)
	df_mem['110k-gb'] = df_mem['110k'] / (1024.0 * 1024.0)

	print(df_mem)


	ax, tool_order, col_list = prepare_for_simul(4)

	#fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, squeeze=False)
	fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, squeeze=False, figsize=(10,5))

	#x = ['CircMiner', 'CIRCexplorer', 'KNIFE', 'CIRI', 'SEGEMEHL', 'CircMarker']
	x = tool_order

	ax = sns.barplot(x = 'Tool', y = '110k-hr', alpha=1.0, data=df_cput, order=tool_order, palette=col_list, ax=axes[0,0])
	ax = sns.barplot(x = 'Tool', y = '1k-hr', alpha=1.0, data=df_cput, order=tool_order, palette=col_list, ax=axes[0,1])
	ax = sns.barplot(x = 'Tool', y = '110k-gb', alpha=1.0, data=df_mem, order=tool_order, palette=col_list, ax=axes[1,0])
	ax = sns.barplot(x = 'Tool', y = '1k-gb', alpha=1.0, data=df_mem, order=tool_order, palette=col_list, ax=axes[1,1])

	#ax.set_yscale('log')

	#fig = plt.bar(x, y, color=col_list)
	# plt.legend(fig, ['First','Second','Third'], loc = "upper left", title = "cat")

	axes[0,0].set_ylabel('Time (hr)')
	axes[1,0].set_ylabel('Memory usage (GB)')

	axes[0,1].yaxis.label.set_visible(False)
	axes[1,1].yaxis.label.set_visible(False)
	axes[0,0].xaxis.label.set_visible(False)
	axes[0,1].xaxis.label.set_visible(False)

	sns.despine(left=True, right=True, top=True)

	for ax in fig.axes:
		matplotlib.pyplot.sca(ax)
		plt.xticks(rotation=90)

		ax.get_yaxis().set_minor_locator(matplotlib.ticker.AutoMinorLocator())
		plt.grid(b=True, which='major', color='gray', linewidth=1.0)
		plt.grid(b=True, which='minor', color='gray', linewidth=0.5)
		plt.gca().xaxis.grid(False)
	
	axes[1,0].set(xlabel='Methods\nDiluted Dataset1 (DC1)')
	axes[1,1].set(xlabel='Methods\nDiluted Dataset2 (DC2)')

	#plt.savefig('time_diluted_1k.pdf')
	#plt.savefig('time_diluted_110k.pdf')
	#plt.savefig('mem_1k.pdf')
	#plt.savefig('mem_110k.pdf')
	plt.savefig('dc_time_mem.pdf', bbox_inches='tight')
	#plt.savefig('dc_time_mem.png', bbox_inches='tight', dpi=700)

def seedlim_single_sample(mydf, id, scale=90):
	keys = list(mydf.columns.values)
	print(keys)
	mydf['AllCorrectPerc'] = 1000.0 * (mydf['Correct (Us)'] + mydf['Correct (Us).1']) / mydf['Total read cnt']
	mydf['AllIncorrectPerc'] = 1000.0 * (mydf['Concordant Cnt'] + mydf['Discordant Cnt'] - (mydf['Correct (Us)'] + mydf['Correct (Us).1'])) / mydf['Total read cnt']
	gk = mydf.groupby('Seed Limit')
	print(gk.first())
	meandf = gk.mean()
	stddf = gk.std()
	print('STD')
	print(stddf)

	meandf['Unmapped'] = 100.0 * (meandf['Total read cnt'] - meandf['Concordant Cnt'] - meandf['Discordant Cnt']) / meandf['Total read cnt']
	# meandf['Correct'] = 100.0 * (meandf['Correct (Us)'] + meandf['Correct (Us).1']) / meandf['Total read cnt']
	meandf['Correct'] = 100.0 * (meandf['Correct (Us)'] + meandf['Correct (Us).1']) / meandf['Total read cnt'] - scale
	meandf['Incorrect'] = 100.0 * ((meandf['Concordant Cnt'] + meandf['Discordant Cnt']) - (meandf['Correct (Us)'] + meandf['Correct (Us).1'])) / meandf['Total read cnt']
	print(meandf)

	# plot
	ax = prepare_for_stack()

	# font setting
	plt.rcParams['font.family'] = "Liberation Serif"
	plt.rcParams['mathtext.fontset'] = 'dejavuserif'
	plt.rcParams['font.size'] = 22
	plt.rcParams['axes.labelsize'] = 16
	plt.rcParams['legend.fontsize'] = 14
	plt.rcParams['xtick.labelsize'] = 22
	plt.rcParams['ytick.labelsize'] = 16

	font = {'family': 'serif',
			'size': 14,
			'weight': 'bold',
        # 'color':  'darkred',
        }

	bar_width = 0.55
	# names = ('SL=100','SL=500','SL=1000', 'SL=5000')
	names = ('100','500','1000', '5000')
	r = [0, 1, 2, 3]
	# Create correct bars
	ax1 = plt.bar(r, meandf['Correct'], edgecolor='gray', width=bar_width, label="Correct mapping", yerr=stddf['AllCorrectPerc'])
	# Create incorrect bars
	ax2 = plt.bar(r, meandf['Incorrect'], bottom=meandf['Correct'], edgecolor='gray', width=bar_width, label="Incorrect mapping", yerr=stddf['AllIncorrectPerc'])
	# Create unmapped bars
	ax3 = plt.bar(r, meandf['Unmapped'], bottom=[i+j for i,j in zip(meandf['Correct'], meandf['Incorrect'])], edgecolor='gray', width=bar_width, label="Unmapped")
	
	# add percentage on top
	for r1, r2, r3 in zip(ax1, ax2, ax3):
		h1 = r1.get_height()
		h2 = r2.get_height()
		h3 = r3.get_height()
		plt.text(r1.get_x() + r1.get_width() / 2., h1 / 2., "%.2lf%%" % (h1+scale), ha="center", va="center", color="black", fontsize=font['size'], fontweight="bold")
		plt.text(r2.get_x() + r2.get_width() / 2., h1 + h2 / 2., "%.2lf%%" % h2, ha="center", va="center", color="black", fontsize=font['size'], fontweight="bold")
		plt.text(r3.get_x() + r3.get_width() / 2., h1 + h2 + h3 / 2., "%.2lf%%" % h3, ha="center", va="center", color="black", fontsize=font['size'], fontweight="bold")

	# Custom x, y axis
	plt.xticks(r, names, fontname=plt.rcParams['font.family'], fontsize=font['size'])
	plt.xlabel("Seed limit value", fontdict=font)
	plt.yticks(np.arange(0, 101-scale, 1.0), np.arange(scale, 101, 1.0), fontname=plt.rcParams['font.family'], fontsize=font['size'])
	
	plt.title("{} Mapping Accuracy".format(id), fontdict=font)
	# if id == 'S3':
	plt.legend(bbox_to_anchor=(1.00,1.10), loc="lower right", borderaxespad=0)

	plt.savefig('{}_stack.pdf'.format(id), bbox_inches='tight')

def seedlim_single_sample_wstar(mydf, id, subax, scale=90):
	subax.xaxis.grid(False)

	keys = list(mydf.columns.values)
	print(keys)
	mydf['AllCorrectPerc'] = 1000.0 * (mydf['Correct (Us)'] + mydf['Correct (Us).1']) / mydf['Total read cnt']
	mydf['AllIncorrectPerc'] = 1000.0 * (mydf['Concordant Cnt'] + mydf['Discordant Cnt'] - (mydf['Correct (Us)'] + mydf['Correct (Us).1'])) / mydf['Total read cnt']
	gk = mydf.groupby('Seed Limit')
	print(gk.first())
	meandf = gk.mean()
	stddf = gk.std()
	print('STD')
	print(stddf)

	meandf['Unmapped'] = 100.0 * (meandf['Total read cnt'] - meandf['Concordant Cnt'] - meandf['Discordant Cnt']) / meandf['Total read cnt']
	meandf['Correct'] = 100.0 * (meandf['Correct (Us)'] + meandf['Correct (Us).1']) / meandf['Total read cnt'] - scale
	meandf['Incorrect'] = 100.0 * ((meandf['Concordant Cnt'] + meandf['Discordant Cnt']) - (meandf['Correct (Us)'] + meandf['Correct (Us).1'])) / meandf['Total read cnt']

	meandf['STARUnmapped'] = 100.0 * (meandf['Total read cnt'] - meandf['All STAR Mappings(ed<=8)']) / meandf['Total read cnt']
	meandf['STARCorrect'] = 100.0 * (meandf['All STAR Correct(ed<=8)']) / meandf['Total read cnt'] - scale
	meandf['STARIncorrect'] = 100.0 * (meandf['All STAR Mappings(ed<=8)'] - meandf['All STAR Correct(ed<=8)']) / meandf['Total read cnt']
	print(meandf)

	# font setting
	plt.rcParams['font.family'] = "Liberation Serif"
	plt.rcParams['mathtext.fontset'] = 'dejavuserif'
	plt.rcParams['font.size'] = 22
	plt.rcParams['axes.labelsize'] = 16
	plt.rcParams['legend.fontsize'] = 16
	plt.rcParams['xtick.labelsize'] = 22
	plt.rcParams['ytick.labelsize'] = 16

	font = {'family': 'serif',
			'size': 14,
			'weight': 'bold',
        # 'color':  'darkred',
        }

	bar_width = 0.65
	names = ('SL=100','SL=500','SL=1000', 'SL=5000', 'STAR')
	# names = ('100','500','1000', '5000')
	r = [0, 1, 2, 3, 4]
	# Create correct bars
	corrdf = meandf['Correct'].append(meandf['STARCorrect'])[:5]
	print (corrdf)
	# ax1 = subax.bar(r, corrdf, edgecolor='gray', width=bar_width, label="Correct mapping", yerr=stddf['AllCorrectPerc'])
	ax1 = subax.bar(r, corrdf, edgecolor='gray', width=bar_width, label="Correct mapping")
	# Create incorrect bars
	incorrdf = meandf['Incorrect'].append(meandf['STARIncorrect'])[:5]
	# ax2 = subax.bar(r, incorrdf, bottom=meandf['Correct'], edgecolor='gray', width=bar_width, label="Incorrect mapping", yerr=stddf['AllIncorrectPerc'])
	ax2 = subax.bar(r, incorrdf, bottom=corrdf, edgecolor='gray', width=bar_width, label="Incorrect mapping")
	# Create unmapped bars
	unmappeddf = meandf['Unmapped'].append(meandf['STARUnmapped'])[:5]
	ax3 = subax.bar(r, unmappeddf, bottom=[i+j for i,j in zip(corrdf, incorrdf)], edgecolor='gray', width=bar_width, label="Unmapped")
	
	# add percentage on top
	for r1, r2, r3 in zip(ax1, ax2, ax3):
		h1 = r1.get_height()
		h2 = r2.get_height()
		h3 = r3.get_height()
		subax.text(r1.get_x() + r1.get_width() / 2., h1 / 2., "%.2lf%%" % (h1+scale), ha="center", va="center", color="black", fontsize=font['size'], fontweight="bold")
		subax.text(r2.get_x() + r2.get_width() / 2., h1 + h2 / 2., "%.2lf%%" % h2, ha="center", va="center", color="black", fontsize=font['size'], fontweight="bold")
		subax.text(r3.get_x() + r3.get_width() / 2., h1 + h2 + h3 / 2., "%.2lf%%" % h3, ha="center", va="center", color="black", fontsize=font['size'], fontweight="bold")

	# subax.set_yticklabels(np.arange(0, 101-scale, 1.0), np.arange(scale, 101, 1.0), fontname=plt.rcParams['font.family'], fontsize=font['size'])
	# subax.set_xticklabels(r, names, fontname=plt.rcParams['font.family'], fontsize=font['size'])
	fontdic = {'fontname': plt.rcParams['font.family'], 'fontsize': font['size']}
	subax.set_xticks(r)
	subax.set_xticklabels(names, fontdict=fontdic)
	subax.set_yticks(np.arange(0, 101-scale, 1.0))
	subax.set_yticklabels(np.arange(scale, 101, 1.0), fontdict=fontdic)
	subax.set_title("{}".format(id), fontdict=font)

	return font

def seedlim_plot(args):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "All")

	# setup plot and subplots
	scale = 90
	ax = prepare_for_stack()
	fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(24, 8))

	mydf = df[0:20]
	# print(mydf)
	seedlim_single_sample_wstar(mydf, 'S1', ax[0], scale=scale)

	mydf = df[20:40]
	# print(mydf)
	seedlim_single_sample_wstar(mydf, 'S2', ax[1], scale=scale)

	mydf = df[40:]
	# print(mydf)
	font = seedlim_single_sample_wstar(mydf, 'S3', ax[2], scale=scale)

	# add a big axis, hide frame
	fig.add_subplot(111, frameon=False)
	plt.gca().xaxis.grid(False)
	plt.gca().yaxis.grid(False)
	# hide tick and tick label of the big axis
	plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	plt.xlabel("CircMiner (Seed Limit, SL) / STAR", fontdict=font, labelpad=30)

	plt.title("Mapping Accuracy", fontdict=font, pad=30.0)
	ax[2].legend(bbox_to_anchor=(1.00,1.10), loc="lower right", borderaxespad=0)

	for ax in fig.axes:
		plt.sca(ax)
		plt.xticks(rotation=30)

	plt.savefig('all_stack_wstar.pdf', bbox_inches='tight')

def seedlim_plot2(args):
	xl = pd.ExcelFile(args[0])
	print(xl.sheet_names)
	df = pd.io.excel.ExcelFile.parse(xl, "CM")
	#df["Seed Limit"] = math.sqrt(int(df["SL"]))
	df["Seed Limit"] = df["SL"] 
	
	df.loc[df["Seed Limit"] == 100, "size"] = 1
	df.loc[df["Seed Limit"] == 500, "size"] = 2
	df.loc[df["Seed Limit"] == 1000, "size"] = 4
	df.loc[df["Seed Limit"] == 5000, "size"] = 8

	print(df)

	ax, tool_order, col_list = prepare_for_simul()
	#plt.gca().xaxis.grid(True)
	#plt.gca().yaxis.grid(True)
	
	kw = dict(edgecolor="none")
	
	#ax = sns.lineplot(x="Recall", y="Precision", hue="SL", style="Sample", markers=['*', 'D', 'o'], dashes=False, data=df)#, legend_out = True)
	#ax = sns.lineplot(x="Recall", y="Precision", hue="Sample", style="Sample", markers=['.', '.', '.'], dashes=False, data=df, alpha=0.7, lw=0.6, legend=False)
	ax = sns.lineplot(x="Recall", y="Precision", color="k", style="Sample", markers=[',', ',', ','], dashes=False, data=df, alpha=1.0, lw=0.6, legend=False)
	ax = sns.scatterplot(x="Recall", y="Precision", hue="Sample", size="size", style="Sample", markers=['>', 'D', 'p'], data=df, sizes=(20,120), linewidth=0)

	#for line in range(0,df.shape[0]):
	#	#ax.text(df.Recall[line]+0.2, df.Precision[line], df.SL[line], horizontalalignment='left', size='medium', color='black', weight='semibold')
	#	yval = df.Precision[line] + 0.08 if line % 2 == 0 else df.Precision[line] - 0.1
	#	ax.text(df.Recall[line], yval, df.SL[line], horizontalalignment='left', color='black', fontsize=6)

	#plt.title("Mapping Accuracy", fontdict=font, pad=30.0)
	#ax[2].legend(bbox_to_anchor=(1.00,1.10), loc="lower right", borderaxespad=0)
	
	#plt.gca().set_xlim(right=100.0)
	#plt.gca().set_ylim(top=100.0)
	plt.gca().set_ylim(bottom=95.0)

	handles, labels = ax.get_legend_handles_labels()
	print(handles)
	del(handles[-1])
	del(labels[-1])
	labels[0] = 'Simulated Dataset'
	labels[4] = 'Seed Limit'
	labels[5] = 100
	labels[6] = 500
	labels[7] = 1000
	labels[8] = 5000
	print(labels)
	leg = plt.legend(handles=handles[:], labels=labels[:], loc='upper left', bbox_to_anchor=(0.0, 1.0))
	leg._legend_box.align = "left"
	ax.add_artist(leg)

	ax.set(xlabel='Recall (avg. over 5 replicates)', ylabel='Precision (avg. over 5 replicates)')
	
	#plt.savefig('sl_prec_recall.pdf', bbox_inches='tight')
	plt.savefig('sl_prec_recall.png', bbox_inches='tight', dpi=800)


def main():
	# import  matplotlib.font_manager
	# flist = matplotlib.font_manager.get_fontconfig_fonts()
	# names = [matplotlib.font_manager.FontProperties(fname=fname).get_name() for fname in flist]
	# print(names)

	mode = sys.argv[1]
	args = sys.argv[2:]
	print(mode)
	print(args)
	if mode == 'tlen':
		tlen_plot(args)
	elif mode == 'simul1':
		simul1_plot(args)
	elif mode == 'resource':
		resource_plot(args)
		resource_plot_single(args)
	elif mode == 'seedlim':
		seedlim_plot(args)
	elif mode == 'seedlim2':
		seedlim_plot2(args)
	elif mode == 'real1':
		real1_plot(args)

if __name__ == '__main__':
	main()
