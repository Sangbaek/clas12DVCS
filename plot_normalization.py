#!/usr/bin/env python3
"""
Script to study the normalization issues.
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

# parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# parser.add_argument("-ch","--chapter", help="chapter", default = "2")
# args = parser.parse_args()

# Electron reconstruction
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/"
parent_MC_BH = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bh/"
parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_1g/"
parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_2g/"
parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/exp/"

epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
bhSimInb = pd.read_pickle(parent_MC_BH + "4892.pkl")
dvcsSimInb = pd.read_pickle(parent_MC + "4893.pkl")
bkgSimInb = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
pi0ExpInb = pd.read_pickle(parent_exp + "pi0.pkl")
pi0SimInb = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")

epgExpInbCDFT = epgExpInb.loc[epgExpInb.config >= 3]
bhSimInbCDFT = bhSimInb.loc[bhSimInb.config >= 3]
dvcsSimInbCDFT = dvcsSimInb.loc[dvcsSimInb.config >= 3]
bkgSimInbCDFT = bkgSimInb.loc[bkgSimInb.config >= 3]
pi0ExpInbCDFT = pi0ExpInb.loc[(pi0ExpInb.config >= 3)]
pi0SimInbCDFT = pi0SimInb.loc[(pi0SimInb.config >= 3)]

epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
bhSimInbCD = dvbhmInb.loc[dvbhmInb.config == 2]
dvcsSimInbCD = dvcsSimInb.loc[dvcsSimInb.config == 2]
bkgSimInbCD = bkgSimInb.loc[bkgSimInb.config == 2]
pi0ExpInbCD = pi0ExpInb.loc[(pi0ExpInb.config == 2)]
pi0SimInbCD = pi0SimInb.loc[(pi0SimInb.config == 2)]

epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]
bhSimInbFD = dvbhmInb.loc[dvbhmInb.config == 1]
dvcsSimInbFD = dvcsSimInb.loc[dvcsSimInb.config == 1]
bkgSimInbFD = bkgSimInb.loc[bkgSimInb.config == 1]
pi0ExpInbFD = pi0ExpInb.loc[(pi0ExpInb.config == 1)]
pi0SimInbFD = pi0SimInb.loc[(pi0SimInb.config == 1)]

contInbCDFT = len(pi0ExpInbCDFT)*len(bkgSimInbCDFT)/len(pi0SimInbCDFT)/len(epgExpInbCDFT)
contInbCD = (len(pi0ExpInbCD.loc[(pi0ExpInbCD.phi1<30)|(pi0ExpInbCD.phi1>330)]) *
			len(bkgSimInbCD.loc[(bkgSimInbCD.phi1<30)|(bkgSimInbCD.phi1>330)]) /
			len(pi0SimInbCD.loc[(pi0SimInbCD.phi1<30)|(pi0SimInbCD.phi1>330)]) /
			len(epgExpInbCD.loc[(epgExpInbCD.phi1<30)|(epgExpInbCD.phi1>330)]))
contInbFD = len(pi0ExpInbFD)*len(bkgSimInbFD)/len(pi0SimInbFD)/len(epgExpInbFD)

varstoplot = ["xB", "Q2", "t1", "Ep", "Etheta", "Ephi", "Pp", "Ptheta", "Pphi", "Gp", "Gtheta", "Gphi"]
title = [r"$x_B$", r"$Q^2$", r"$|t|$", r"$p_{e'}$", r"$\theta_{e'}", r"$\phi_{e'}$", r"$p_{p'}$", r"$\theta_{p'}", r"$\phi_{p'}$", r"$p_{\gamma}$", r"$\theta_{\gamma}", r"$\phi_{\gamma}$"]
binstoplot = [np.linspace(0.05, 0.7, 101), np.linspace(1, 7, 101), np.linspace(0, 1, 101), 100, 100, 100, 100, 100, 100, 100, 100, 100]
unit = ["", GeVc2, GeV2, GeVc, degree, degree, GeVc, degree, degree, GeVc, degree, degree]

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
for yind in range(0, 4):
	for xind in range(0, 3):
		ind = 4*yind + xind
		simDist_dvcs, bins = np.histogram(dvcsSimInbCDFT.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
		simDist_dvpi0, _ = np.histogram(bkgSimInbCDFT.loc[:, varstoplot[ind]], bins, density = True)
		simDist = (1-contInbCDFT)*simDist_dvcs + contInbCDFT*simDist_dvpi0

		simDist_bh, _ = np.histogram(bhSimInbCDFT.loc[:, varstoplot[ind]], bins, density = True)
		simDist2 = (1-contInbCDFT)*simDist_bh + contInbCDFT*simDist_dvpi0

		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(epgExpInbCDFT.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
		axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
		axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])

plt.savefig("plots/normalization/CDFT_particle_kine.pdf")
plt.clf()

fig, axs = plt.subplots(3, 1, figsize = (12, 6))
for yind in range(0, 4):
	for xind in range(0, 3):
		ind = 4*yind + xind
		simDist_dvcs, bins = np.histogram(dvcsSimInbCD.loc[(dvcsSimInbCD.phi1<30)|(dvcsSimInbCD.phi1>30), varstoplot[ind]], binstoplot[ind], density = True)
		simDist_dvpi0, _ = np.histogram(bkgSimInbCD.loc[(bkgSimInbCD.phi1<30)|(bkgSimInbCD.phi1>30), varstoplot[ind]], bins, density = True)
		simDist = (1-contInbCD)*simDist_dvcs + contInbCD*simDist_dvpi0

		simDist_bh, _ = np.histogram(bhSimInbCD.loc[(bhSimInbCD.phi1<30)|(bhSimInbCD.phi1>30), varstoplot[ind]], bins, density = True)
		simDist2 = (1-contInbCD)*simDist_bh + contInbCD*simDist_dvpi0

		bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		axs[yind, xind].hist(epgExpInbCD.loc[(epgExpInbCD.phi1<30)|(epgExpInbCD.phi1>30),varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
		axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
		axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
		axs[yind, xind].set_title(title[ind])
		if (unit[ind]):
			axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
		else:
			axs[yind, xind].set_xlabel(title[ind])
		axs[yind, xind].set_xlim([bins[0], bins[-1]])

plt.savefig("plots/normalization/CD_particle_kine.pdf")
plt.clf()
