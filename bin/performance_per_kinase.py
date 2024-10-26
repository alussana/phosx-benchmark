import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from math import pi
from scipy.signal import resample
import sys

font = {'size'   : 11.}

plt.rc('font', **font)

rank_phuego={}
rank_seeds={}
rank_rwr={}
rank_pcsf={}
pie_count={}
for comp in [sys.argv[1]]:
	f1=open("enrich_comparison/"+comp+".txt")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")

		if seq[0]=="seeds":
			if seq[3]=="nan":
				rank_seeds[seq[2]]=11
			elif int(seq[3])>=10:
				rank_seeds[seq[2]]=10
			else:
				rank_seeds[seq[2]]=int(seq[3])

		if seq[0]=="pcsf":
			if seq[3]=="nan":
				rank_pcsf[seq[2]]=11
			elif int(seq[3])>=10:
				rank_pcsf[seq[2]]=10
			else:
				rank_pcsf[seq[2]]=int(seq[3])

		if seq[0]=="rwr":
			if seq[3]=="nan":
				rank_rwr[seq[2]]=11
			elif int(seq[3])>=10:
				rank_rwr[seq[2]]=10
			else:
				rank_rwr[seq[2]]=int(seq[3])

		if seq[0]=="net":
			if seq[-2] in pie_count:
				pie_count[seq[-2]]+=1
			else:
				pie_count[seq[-2]]=1

			if seq[3]=="nan":
				rank_phuego[seq[2]]=11
			elif int(seq[3])>=10:
				rank_phuego[seq[2]]=10
			else:
				rank_phuego[seq[2]]=int(seq[3])

		seq=f1.readline()


categories_seeds = list(rank_seeds.keys())
categories_seeds_1 = [*categories_seeds, categories_seeds[0]]

split_cat=[]
for i in categories_seeds_1:
	author=i.split(",")[0]
	rest=i.split(",")[1:]
	split_cat.append(author+",\n"+" ".join(rest))
split_cat[-1]=""


rankings_seeds =list(rank_seeds.values())
rankings_seeds_1 = [*rankings_seeds, rankings_seeds[0]]

categories_pcsf = list(rank_pcsf.keys())
categories_pcsf_1 = [*categories_pcsf, categories_pcsf[0]]
rankings_pcsf =list(rank_pcsf.values())
rankings_pcsf_1 = [*rankings_pcsf, rankings_pcsf[0]]


categories_rwr = list(rank_rwr.keys())
categories_rwr_1 = [*categories_rwr, categories_rwr[0]]
rankings_rwr =list(rank_rwr.values())
rankings_rwr_1 = [*rankings_rwr, rankings_rwr[0]]

categories = list(rank_phuego.keys())
categories = [*categories, categories[0]]
rankings =list(rank_phuego.values())
rankings_1 = [*rankings, rankings[0]]




label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(rankings_1))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 12), subplot_kw=dict(polar=True, zorder=0))


ax.plot(label_loc, rankings_seeds_1, label='Seeds',marker="x",c="#009139",alpha=0.5)
ax.fill(label_loc,rankings_seeds_1, '#009139', alpha=0.1)

ax.plot(label_loc, rankings_rwr_1, label='RWR',marker="D",c="#FB9A29",alpha=0.5)
ax.fill(label_loc,rankings_rwr_1, '#FB9A29', alpha=0.1)

ax.plot(label_loc, rankings_pcsf_1, label='PCSF',marker="v",c="#46353A",alpha=0.5)
ax.fill(label_loc,rankings_pcsf_1, '#46353A', alpha=0.1)

ax.plot(label_loc, rankings_1, label='phuEGO',marker="o",c="#882E72",alpha=0.5)
ax.fill(label_loc,rankings_1, '#882E72', alpha=0.1)
ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1.35))

lines, labels = plt.thetagrids(np.degrees(label_loc), labels=split_cat)

ax.grid(True)

angles = np.linspace(0,2*np.pi,len(ax.get_xticklabels())+1)
angles[np.cos(angles) < 0] = angles[np.cos(angles) < 0] + np.pi
angles = np.rad2deg(angles)
labels = []
for label, angle in zip(ax.get_xticklabels(), angles):
	x,y = label.get_position()
	"""
	if label.get_text().startswith("Pan"):
		lab = ax.text(x,y-0.23, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif len(label.get_text())>37:
		lab = ax.text(x,y-0.34, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 32<len(label.get_text())<=37:
		lab = ax.text(x,y-0.28, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 30<=len(label.get_text())<=32:
		lab = ax.text(x,y-0.24, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 28<=len(label.get_text())<30:
		lab = ax.text(x,y-0.16, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 26<=len(label.get_text())<=27:
		lab = ax.text(x,y-0.16, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	else:
		lab = ax.text(x,y-0.1, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	"""
	if label.get_text().startswith("Pan"):
		lab = ax.text(x,y-0.23, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif len(label.get_text())>37:
		lab = ax.text(x,y-0.28, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 32<len(label.get_text())<=37:
		lab = ax.text(x,y-0.25, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 30<=len(label.get_text())<=32:
		lab = ax.text(x,y-0.22, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 28<=len(label.get_text())<30:
		lab = ax.text(x,y-0.16, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	elif 26<=len(label.get_text())<=27:
		lab = ax.text(x,y-0.12, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	else:
		lab = ax.text(x,y-0.05, label.get_text(), transform=label.get_transform(), ha=label.get_ha(), va=label.get_va())
	lab.set_rotation(angle)
	labels.append(lab)

ax.set_xticklabels([])
ax.set_rlabel_position(100) # yticks position

plt.yticks([1,3,5,7,10,11], ["1","3","5","7",">=10","None"], color="black", size=10)

ax.set_rorigin(-3)
ax2 = fig.add_subplot(111, zorder=1)
sizes = [12,1, 10, 4,4, 12, 4, 5]
sizes=pie_count.values()
size=0.1
patches,text=ax2.pie(sizes, radius=0.36 - size, wedgeprops=dict(width=size, edgecolor='w',linewidth= 1.5),startangle = -2)

ax2.legend()
plt.legend(patches, pie_count.keys(), loc="upper right", bbox_to_anchor=(1.45, 1.35))
ax.tick_params(pad=20)
plt.tight_layout()
plt.savefig("figures/figure1a_"+sys.argv[1]+".pdf",dpi=600)
plt.show()
