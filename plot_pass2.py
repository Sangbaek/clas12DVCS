#!/usr/bin/env python3
"""
Script to reproduce some of plots used in Sangbaek's thesis.
The original plots were produced in the jupyter notebook.
Some of them had problems in that
(1) they are not clearly visible,
(2) they do not have color bar scale for the 2d histogram.
"""
import gc
import matplotlib.pyplot as plt
from copy import copy
cmap = copy(plt.cm.get_cmap("jet"))
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from utils.const import *
from utils.physics import *
from utils.fiducial import *
from matplotlib.colors import LogNorm
import argparse

cmap.set_under('w',0)
cmap.set_bad('w',0)

degree = r"${}^{\circ}$"
GeV = "GeV"
GeV2 = "GeV"+r"${}^{2}$"
GeVc = "GeV/c"
GeVc2 = "(GeV/c)"+r"${}^{2}$"

import matplotlib
# initial settings
pgf_with_latex = {
		"pgf.texsystem": "pdflatex",
		"text.usetex": True,			# use LaTeX to write all text
		"font.family": "sans-serif",		
		"font.sans-serif": "Helvetica",
		"font.size": 25,				# default font size
		"axes.titlepad": 20,			# x and y label size
		"axes.labelsize": 24,			# x and y label size
		"axes.titlesize": 24,		  # subfigure title size, i.e. title size when one figure
		"legend.fontsize": 22,			# legend size
		"xtick.labelsize": 23,			# x axis tick label size
		"ytick.labelsize": 23,			# y axis tick label 
		"figure.titlesize": 25,         # Figure title size, useful when you have multiple plots in one canvas.
		"pgf.preamble": r"\usepackage{xcolor}",     # xcolor for colours
		"figure.autolayout": False
}
matplotlib.rcParams.update(pgf_with_latex)

parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-ch","--chapter", help="chapter", default = "2")
parser.add_argument("-pol","--polarity", help="polarity", default = "inbending")
parser.add_argument("-fom","--figureofmerit", help="figure of merit", default = None)
parser.add_argument("-det","--detector", help="detector", default = None)
args = parser.parse_args()


point = r'$.$'
zero = r'$0$'
one = r"$1$"
two = r"$2$"
three = r"$3$"
four = r"$4$"
five = r"$5$"
six = r"$6$"
seven = r"$7$"
eight = r"$8$"
nine = r"$9$"
times = r"$\times$"
ten = r"$10$"
hundred = r"$10^2$"
thousand = r"$10^3$"
tenthousands = r"$10^4$"
hundredthousands = r'$10^5$'
tenth = r"$10^{-1}$"
hundredth = r"$10^{-2}$"
thousandth = r"$10^{-3}$"

parent_exp = "/volatile/clas12/sangbaek/pass1_test/convPkl_full/inb/exp/"

#epg Exp
epgExpInb = pd.read_pickle(parent_exp + "dvcs_0.pkl")
pi0ExpInb = pd.read_pickle(parent_exp + "pi0_0.pkl")

epgExpInbCDFT = epgExpInb.loc[epgExpInb.config == 3]
pi0ExpInbCDFT = pi0ExpInb.loc[(pi0ExpInb.config == 3)]

epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
pi0ExpInbCD = pi0ExpInb.loc[(pi0ExpInb.config == 2)]

epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]
pi0ExpInbFD = pi0ExpInb.loc[(pi0ExpInb.config == 1)]

parent_exp = "/volatile/clas12/sangbaek/pass2_test/convPkl_full/inb/exp/"

#epg Exp
epgExpInb2 = pd.read_pickle(parent_exp + "dvcs_0.pkl")
pi0ExpInb2 = pd.read_pickle(parent_exp + "pi0_0.pkl")

epgExpInbCDFT2 = epgExpInb2.loc[epgExpInb2.config == 3]
pi0ExpInbCDFT2 = pi0ExpInb2.loc[(pi0ExpInb2.config == 3)]

epgExpInbCD2 = epgExpInb2.loc[epgExpInb2.config == 2]
pi0ExpInbCD2 = pi0ExpInb2.loc[(pi0ExpInb2.config == 2)]

epgExpInbFD2 = epgExpInb2.loc[epgExpInb2.config == 1]
pi0ExpInbFD2 = pi0ExpInb2.loc[(pi0ExpInb2.config == 1)]

parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full_noCorr/inb/exp/"

#epg Exp
epgExpInb3 = pd.read_pickle(parent_exp + "dvcs.pkl")
pi0ExpInb3 = pd.read_pickle(parent_exp + "pi0.pkl")

epgExpInbCDFT3 = epgExpInb3.loc[epgExpInb3.config == 3]
pi0ExpInbCDFT3 = pi0ExpInb3.loc[(pi0ExpInb3.config == 3)]

epgExpInbCD3 = epgExpInb3.loc[epgExpInb3.config == 2]
pi0ExpInbCD3 = pi0ExpInb3.loc[(pi0ExpInb3.config == 2)]

epgExpInbFD3 = epgExpInb3.loc[epgExpInb3.config == 1]
pi0ExpInbFD3 = pi0ExpInb3.loc[(pi0ExpInb3.config == 1)]


varstoplot = ["coneAngle", "MM2_eg", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
title = [r"$\theta_{e'\gamma}$", "MM"+r"${}^2_{e'\gamma}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi_{\vec{L}\vec{\Gamma}}$" , "ME"+r"${}_{e'p'\gamma}$", "MM"+r"${}^{2}_{e'p'\gamma}$", "MM"+r"${}^{2}_{e'p'}$", "MPt"+r"${}_{e'p'\gamma}$"]
unit = [degree, GeV2, degree, degree, GeV, GeV2, GeV2, GeVc]
df1 = epgExpInbCDFT
df2 = epgExpInbCDFT2
df3 = epgExpInbCDFT3
xlb = [10, 0.3, 0, 0, -.3, -0.01, -0.3, 0]
xub = [35, 1.5, 0.75, 6, .3, 0.01, 0.3, 0.1]
yub = [0.15, 2, 3, 0.4, 3, 250, 4, 25]
fig, axs = plt.subplots(2, 4, figsize = (16, 8))
for yind in range(0, 2):
	for xind in range(0, 4):
		ind = 4*yind + xind
		# simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
		# simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
		# simDist = (1-contInbCDFT)*simDist_dvcs + contInbCDFT*simDist_dvpi0
		_, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(df1.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "spring2019 pass1, 10.2 GeV")
		axs[yind, xind].hist(df2.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "spring2019 pass2, 10.2 GeV")
		axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='g', density=True, linewidth=1, label = "fall2018 pass1, 10.6 GeV")
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])
axs[0, 0].set_xlim([10, 30])
axs[0, 0].set_xticks([10, 20, 30])
axs[0, 1].set_xlim([0.5, 1.5])
axs[0, 1].set_xticks([0.5, 1, 1.5])
axs[0, 2].set_xticks([0, 0.2, 0.4, 0.6])
axs[0, 3].set_xticks([0, 2, 4, 6])
axs[1, 0].set_xlim([-0.3, 0.6])
axs[1, 0].set_xticks([-0.3, 0, 0.3, 0.6])
axs[1, 1].set_xlim([-0.01, 0.01])
axs[1, 1].set_xticks([-0.01, 0, 0.01])
axs[1, 2].set_xlim([-0.3, 0.3])
axs[1, 2].set_xticks([-0.3, 0, 0.3])
axs[1, 3].set_xlim([0, 0.08])
axs[1, 3].set_xticks([0, 0.04, 0.08])
# axs[0, 0].set_yticks([0.02, 0.04, 0.06, 0.08, 0.1, 0.12])
# axs[0, 1].set_yticks([0.5, 1, 1.5, 2])
# axs[0, 2].set_yticks([1, 2, 3])
# axs[0, 3].set_yticks([0.2, 0.4])
# axs[1, 0].set_yticks([1, 2, 3])
# axs[1, 1].set_yticks([50, 100, 150, 200, 250])
# axs[1, 2].set_yticks([1,2,3])
# axs[1, 3].set_yticks([5, 10, 15, 20, 25])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/dvcsInbCDFTexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

df1 = epgExpInbCD
df2 = epgExpInbCD2
df3 = epgExpInbCD3
fig, axs = plt.subplots(2, 4, figsize = (16, 8))
for yind in range(0, 2):
	for xind in range(0, 4):
		ind = 4*yind + xind
		# simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
		# simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
		# simDist = (1-contInbCD)*simDist_dvcs + contInbCD*simDist_dvpi0
		_, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(df1.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "spring2019 pass1, 10.2 GeV")
		axs[yind, xind].hist(df2.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "spring2019 pass2, 10.2 GeV")
		axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='g', density=True, linewidth=1, label = "fall2018 pass1, 10.6 GeV")
		# axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])
axs[0, 0].set_xlim([10, 30])
axs[0, 0].set_xticks([10, 20, 30])
axs[0, 1].set_xlim([0, 2])
axs[0, 1].set_xticks([0, 0.5, 1, 1.5, 2])
axs[0, 2].set_xlim([0, 0.6])
axs[0, 2].set_xticks([0, 0.2, 0.4, 0.6])
axs[0, 3].set_xticks([0, 2.5, 5])
axs[1, 0].set_xlim([-0.6, 0.6])
axs[1, 0].set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
axs[1, 1].set_xlim([-0.015, 0.015])
axs[1, 1].set_xticks([-0.015, 0, 0.015])
axs[1, 2].set_xlim([-0.3, 0.3])
axs[1, 2].set_xticks([-0.3, 0, 0.3])
axs[1, 3].set_xlim([0, 0.08])
axs[1, 3].set_xticks([0, 0.04, 0.08])

# axs[0, 0].set_yticks([0.05, 0.1, 0.15])
# axs[0, 1].set_yticks([0.2, 0.4, 0.6, 0.8, 1, 1.2])
# axs[0, 2].set_yticks([1, 2, 3])
# axs[0, 3].set_yticks([0.25, 0.5])
# axs[1, 0].set_yticks([0.5, 1, 1.5])
# axs[1, 1].set_yticks([50, 100, 150, 200])
# axs[1, 2].set_yticks([1,2,3,4])
# axs[1, 3].set_yticks([5, 10, 15, 20])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/dvcsInbCDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

df1 = epgExpInbFD
df2 = epgExpInbFD2
df3 = epgExpInbFD3
fig, axs = plt.subplots(2, 4, figsize = (16, 8))
for yind in range(0, 2):
	for xind in range(0, 4):
		ind = 4*yind + xind
		# simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
		# simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
		# simDist = (1-contInbFD)*simDist_dvcs + contInbFD*simDist_dvpi0
		_, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(df1.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "spring2019 pass1, 10.2 GeV")
		axs[yind, xind].hist(df2.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "spring2019 pass2, 10.2 GeV")
		axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='g', density=True, linewidth=1, label = "fall2018 pass1, 10.6 GeV")
		# axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])
axs[0, 0].set_xticks([25, 30, 35, 40, 45])
axs[0, 0].set_xlim([25, 45])
axs[0, 1].set_xticks([0, 0.5, 1, 1.5, 2])
axs[0, 2].set_xticks([0, 0.5, 1, 1.5])
axs[0, 3].set_xticks([0, 2, 4, 6, 8])
axs[1, 0].set_xlim([-0.8, 0.8])
axs[1, 0].set_xticks([-0.8, -0.4, 0, 0.4, 0.8])
axs[1, 1].set_xlim([-0.015, 0.015])
axs[1, 1].set_xticks([-0.015, 0, 0.015])
axs[1, 2].set_xlim([-0.2, 0.2])
axs[1, 2].set_xticks([-0.2, 0, 0.2])
axs[1, 3].set_xlim([0, 0.3])
axs[1, 3].set_xticks([0, 0.1, 0.2, 0.3])
# axs[0, 0].set_yticks([0.02, 0.04, 0.06, 0.08, 0.1, 0.12])
# axs[0, 1].set_yticks([0.2, 0.4, 0.6, 0.8, 1, 1.2])
# axs[0, 2].set_yticks([0.5, 1, 1.5])
# axs[0, 3].set_yticks([0.2, 0.4])
# axs[1, 0].set_yticks([0.5, 1, 1.5])
# axs[1, 1].set_yticks([50, 100, 150])
# axs[1, 2].set_yticks([2, 4, 6])
# axs[1, 3].set_yticks([2, 4, 6, 8])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/dvcsInbFDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (10, 6))

h = axs.hist2d(epgExpInbCDFT2.Ptheta, epgExpInbCDFT2.Pp, bins = [np.linspace(40, 80, 101), np.linspace(0, 1.5, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
cbar = plt.colorbar(h[3], ax = axs)
axs.set_ylim([0, 1.5])
axs.set_xlim([40, 80])
axs.set_ylabel(r"$p_{p'}$" + " [" + GeVc + "]")
axs.set_xlabel(r"$\theta_{p'}$" + " [" + degree + "]")
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pass2_CDFT_Pp_Ptheta.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (10, 6))

h = axs.hist2d(epgExpInbCD2.Ptheta, epgExpInbCD2.Pp, bins = [np.linspace(40, 80, 101), np.linspace(0, 1.5, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
cbar = plt.colorbar(h[3], ax = axs)
axs.set_ylim([0, 1.5])
axs.set_xlim([40, 80])
axs.set_ylabel(r"$p_{p'}$" + " [" + GeVc + "]")
axs.set_xlabel(r"$\theta_{p'}$" + " [" + degree + "]")
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pass2_CDFD_Pp_Ptheta.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (10, 6))

h = axs.hist2d(epgExpInb2.loc[epgExpInb2.Psector>7].Pp, epgExpInb2.loc[epgExpInb2.Psector>7].MM2_ep, bins = [np.linspace(0, 1.5, 101), np.linspace(-.3, .3, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
cbar = plt.colorbar(h[3], ax = axs)
axs.set_xlim([0, 1.5])
axs.set_ylim([-.3, .3])
axs.set_xlabel(r"$p_{p'}$" + " [" + GeVc + "]")
axs.set_ylabel("MM"+r"${}^{2}_{e'p'}$" +" [" + GeV2 + "]")
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pass2_CDprotons_Pp_MM2ep.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (10, 6))

h = axs.hist2d(epgExpInb.loc[epgExpInb.Psector>7].Pp, epgExpInb.loc[epgExpInb.Psector>7].MM2_ep, bins = [np.linspace(0, 1.5, 101), np.linspace(-.3, .3, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
cbar = plt.colorbar(h[3], ax = axs)
axs.set_xlim([0, 1.5])
axs.set_ylim([-.3, .3])
axs.set_xlabel(r"$p_{p'}$" + " [" + GeVc + "]")
axs.set_ylabel("MM"+r"${}^{2}_{e'p'}$" +" [" + GeV2 + "]")
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pass1_CDprotons_Pp_MM2ep.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (10, 6))

h = axs.hist2d(epgExpInb2.loc[epgExpInb2.Psector>7].Ptheta, epgExpInb2.loc[epgExpInb2.Psector>7].MM2_ep, bins = [np.linspace(40, 80, 101), np.linspace(-.3, .3, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
cbar = plt.colorbar(h[3], ax = axs)
axs.set_xlim([40, 80])
axs.set_ylim([-.3, .3])
axs.set_xlabel(r"$\theta_{p'}$" + " [" + degree + "]")
axs.set_ylabel("MM"+r"${}^{2}_{e'p'}$" +" [" + GeV2 + "]")
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pass2_CDprotons_Ptheta_MM2ep.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (10, 6))

h = axs.hist2d(epgExpInb.loc[epgExpInb.Psector>7].Ptheta, epgExpInb.loc[epgExpInb.Psector>7].MM2_ep, bins = [np.linspace(40, 80, 101), np.linspace(-.3, .3, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
cbar = plt.colorbar(h[3], ax = axs)
axs.set_xlim([40, 80])
axs.set_ylim([-.3, .3])
axs.set_xlabel(r"$\theta_{p'}$" + " [" + degree + "]")
axs.set_ylabel("MM"+r"${}^{2}_{e'p'}$" +" [" + GeV2 + "]")

plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pass1_CDprotons_Ptheta_MM2ep.pdf")
plt.clf()


varstoplot = ["Mpi0", "MM2_egg", "reconPi", "coplanarity", "ME_epgg", "MM2_epgg", "MM2_ep", "MPt"]
title = [r"$IM_{\pi^0}$", "MM"+r"${}^2_{e'\pi^0}$", r"$\theta_{\pi^0_{det.}\pi^0_{rec.}}$", r"$\Delta\phi_{e'p'\pi^0}$" , "ME"+r"${}_{e'p'\pi^0}$", "MM"+r"${}^{2}_{e'p'\pi^0}$", "MM"+r"${}^{2}_{e'p'}$", "MPt"+r"${}_{e'p'\pi^0}$"]
unit = [GeV, GeV2, degree, degree, GeV, GeV2, GeV2, GeVc]

df4 = pi0ExpInbCDFT
df5 = pi0ExpInbCDFT2
df6 = pi0ExpInbCDFT3

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
for yind in range(0, 2):
	for xind in range(0, 4):
		ind = 4*yind + xind
		# simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
		_, bins = np.histogram(df4.loc[:, varstoplot[ind]], 100, density = True)
		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "spring2019 pass1, 10.2 GeV")
		axs[yind, xind].hist(df5[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "spring2019 pass2, 10.2 GeV")
		axs[yind, xind].hist(df6[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='g', density=True, linewidth=1, label = "fall2018 pass1, 10.6 GeV")
		# axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])
axs[0, 0].set_xticks([0.1, 0.13, 0.16])
axs[0, 0].set_xlim([0.1, 0.16])
axs[0, 1].set_xticks([0, 0.4, 0.8, 1.2, 1.6])
axs[0, 2].set_xticks([0, 0.25, 0.5 ,0.75, 1])
axs[0, 3].set_xticks([0, 2, 4, 6, 8])
axs[1, 0].set_xlim([-0.4, 0.4])
axs[1, 0].set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
axs[1, 1].set_xlim([-0.02, 0.02])
axs[1, 1].set_xticks([-0.02, 0, 0.02])
axs[1, 2].set_xlim([-0.5, 0.5])
axs[1, 2].set_xticks([-0.5, 0, 0.5])
axs[1, 3].set_xlim([0, 0.15])
axs[1, 3].set_xticks([0, 0.05, 0.1, 0.15])
# axs[0, 0].set_yticks([50, 100])
# axs[0, 1].set_yticks([0.5, 1, 1.5, 2])
# axs[0, 2].set_yticks([.5, 1, 1.5, 2, 2.5])
# axs[0, 3].set_yticks([0.25, 0.5])
# axs[1, 0].set_yticks([1, 2, 3])
# axs[1, 1].set_yticks([50, 100, 150, 200])
# axs[1, 2].set_yticks([2, 4])
# axs[1, 3].set_yticks([5, 10, 15, 20, 25])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pi0InbCDFTexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()


df4 = pi0ExpInbCD
df5 = pi0ExpInbCD2
df6 = pi0ExpInbCD3

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
for yind in range(0, 2):
	for xind in range(0, 4):
		ind = 4*yind + xind
		# simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
		_, bins = np.histogram(df4.loc[:, varstoplot[ind]], 100, density = True)
		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "spring2019 pass1, 10.2 GeV")
		axs[yind, xind].hist(df5[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "spring2019 pass2, 10.2 GeV")
		axs[yind, xind].hist(df6[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='g', density=True, linewidth=1, label = "fall2018 pass1, 10.6 GeV")
		# axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])
axs[0, 0].set_xticks([0.1, 0.13, 0.16])
axs[0, 0].set_xlim([0.1, 0.16])
axs[0, 1].set_xticks([0, 0.5, 1, 1.5, 2])
axs[0, 2].set_xticks([0, 0.5, 1, 1.5])
axs[0, 3].set_xticks([0, 5, 10])
axs[1, 0].set_xlim([-0.6, 0.6])
axs[1, 0].set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
axs[1, 1].set_xlim([-0.02, 0.02])
axs[1, 1].set_xticks([-0.02, 0, 0.02])
axs[1, 2].set_xlim([-0.3, 0.3])
axs[1, 2].set_xticks([-0.3, 0, 0.3])
axs[1, 3].set_xlim([0, 0.2])
axs[1, 3].set_xticks([0, 0.05, 0.1, 0.15, 0.2])
# axs[0, 0].set_yticks([20, 40])
# axs[0, 1].set_yticks([0.2, 0.4, 0.6, 0.8, 1, 1.2])
# axs[0, 2].set_yticks([0.5, 1, 1.5])
# axs[0, 3].set_yticks([0.2, 0.4])
# axs[1, 0].set_yticks([0.4, 0.8, 1.2, 1.6])
# axs[1, 1].set_yticks([50, 100, 150])
# axs[1, 2].set_yticks([2.5, 5])
# axs[1, 3].set_yticks([5, 10, 15])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pi0InbCDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()


df4 = pi0ExpInbFD
df5 = pi0ExpInbFD2
df6 = pi0ExpInbFD3

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
for yind in range(0, 2):
	for xind in range(0, 4):
		ind = 4*yind + xind
		# simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
		_, bins = np.histogram(df4.loc[:, varstoplot[ind]], 100, density = True)
		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "spring2019 pass1, 10.2 GeV")
		axs[yind, xind].hist(df5[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "spring2019 pass2, 10.2 GeV")
		axs[yind, xind].hist(df6[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='g', density=True, linewidth=1, label = "fall2018 pass1, 10.6 GeV")
		# axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])
axs[0, 0].set_xticks([0.1, 0.13, 0.16])
axs[0, 0].set_xlim([0.1, 0.16])
axs[0, 1].set_xticks([0, 0.5, 1, 1.5, 2])
axs[0, 2].set_xticks([0, 0.5, 1, 1.5])
axs[0, 3].set_xticks([0, 5, 10])
axs[1, 0].set_xlim([-0.6, 0.6])
axs[1, 0].set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
axs[1, 1].set_xlim([-0.02, 0.02])
axs[1, 1].set_xticks([-0.02, 0, 0.02])
axs[1, 2].set_xlim([-0.3, 0.3])
axs[1, 2].set_xticks([-0.3, 0, 0.3])
axs[1, 3].set_xlim([0, 0.2])
axs[1, 3].set_xticks([0, 0.05, 0.1, 0.15, 0.2])
# axs[0, 0].set_yticks([25, 50])
# axs[0, 1].set_yticks([0.5, 1, 1.5])
# axs[0, 2].set_yticks([0.5, 1, 1.5])
# axs[0, 3].set_yticks([0.2, 0.4])
# axs[1, 0].set_yticks([0.5, 1, 1.5, 2])
# axs[1, 1].set_yticks([50, 100, 150, 200])
# axs[1, 2].set_yticks([2.5, 5])
# axs[1, 3].set_yticks([2,4,6,8,10,12])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.tight_layout()
plt.savefig("/volatile/clas12/sangbaek/pass2_test/plots/spring2019/pi0InbFDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

