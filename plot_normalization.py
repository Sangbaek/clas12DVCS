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
        "axes.titlesize": 24,			# subfigure title size, i.e. title size when one figure
        "legend.fontsize": 22,			# legend size
        "xtick.labelsize": 23,			# x axis tick label size
        "ytick.labelsize": 23,			# y axis tick label 
        "figure.titlesize": 25,			# Figure title size, useful when you have multiple plots in one canvas.
        "pgf.preamble": r"\usepackage{xcolor}",     # xcolor for colours
        "figure.autolayout": False
}
matplotlib.rcParams.update(pgf_with_latex)

# parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# parser.add_argument("-ch","--chapter", help="chapter", default = "2")
# args = parser.parse_args()

# Read the data files
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

epgExpInbCR = epgExpInb.loc[epgExpInb.config == 4]
bhSimInbCR = bhSimInb.loc[bhSimInb.config == 4]
dvcsSimInbCR = dvcsSimInb.loc[dvcsSimInb.config == 4]
bkgSimInbCR = bkgSimInb.loc[bkgSimInb.config == 4]

epgExpInbCDFT = epgExpInb.loc[epgExpInb.config >= 3]
bhSimInbCDFT = bhSimInb.loc[bhSimInb.config >= 3]
dvcsSimInbCDFT = dvcsSimInb.loc[dvcsSimInb.config >= 3]
bkgSimInbCDFT = bkgSimInb.loc[bkgSimInb.config >= 3]
pi0ExpInbCDFT = pi0ExpInb.loc[(pi0ExpInb.config >= 3)]
pi0SimInbCDFT = pi0SimInb.loc[(pi0SimInb.config >= 3)]

epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
bhSimInbCD = bhSimInb.loc[bhSimInb.config == 2]
dvcsSimInbCD = dvcsSimInb.loc[dvcsSimInb.config == 2]
bkgSimInbCD = bkgSimInb.loc[bkgSimInb.config == 2]
pi0ExpInbCD = pi0ExpInb.loc[(pi0ExpInb.config == 2)]
pi0SimInbCD = pi0SimInb.loc[(pi0SimInb.config == 2)]

epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]
bhSimInbFD = bhSimInb.loc[bhSimInb.config == 1]
dvcsSimInbFD = dvcsSimInb.loc[dvcsSimInb.config == 1]
bkgSimInbFD = bkgSimInb.loc[bkgSimInb.config == 1]
pi0ExpInbFD = pi0ExpInb.loc[(pi0ExpInb.config == 1)]
pi0SimInbFD = pi0SimInb.loc[(pi0SimInb.config == 1)]

contInbCR = 0
contInbCDFT = len(pi0ExpInbCDFT)*len(bkgSimInbCDFT)/len(pi0SimInbCDFT)/len(epgExpInbCDFT)
contInbCD = (len(pi0ExpInbCD.loc[(pi0ExpInbCD.phi1<30)|(pi0ExpInbCD.phi1>330)]) *
            len(bkgSimInbCD.loc[(bkgSimInbCD.phi1<30)|(bkgSimInbCD.phi1>330)]) /
            len(pi0SimInbCD.loc[(pi0SimInbCD.phi1<30)|(pi0SimInbCD.phi1>330)]) /
            len(epgExpInbCD.loc[(epgExpInbCD.phi1<30)|(epgExpInbCD.phi1>330)]))
contInbFD = (len(pi0ExpInbFD.loc[(pi0ExpInbFD.phi1<30)|(pi0ExpInbFD.phi1>330)]) *
            len(bkgSimInbFD.loc[(bkgSimInbFD.phi1<30)|(bkgSimInbFD.phi1>330)]) /
            len(pi0SimInbFD.loc[(pi0SimInbFD.phi1<30)|(pi0SimInbFD.phi1>330)]) /
            len(epgExpInbFD.loc[(epgExpInbFD.phi1<30)|(epgExpInbFD.phi1>330)]))

varstoplot = ["xB", "Q2", "t1", "Ep", "Etheta", "Ephi", "Pp", "Ptheta", "Pphi", "Gp", "Gtheta", "Gphi"]
title = [r"$x_B$", r"$Q^2$", r"$|t|$", r"$p_{e'}$", r"$\theta_{e'}$", r"$\phi_{e'}$", r"$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", r"$p_{\gamma}$", r"$\theta_{\gamma}$", r"$\phi_{\gamma}$"]
binstoplot = [np.linspace(0.05, 0.7, 101), np.linspace(1, 7, 101), np.linspace(0, 1, 101), 100, 100, 100, 100, 100, 100, 100, 100, 100]
unit = ["", GeVc2, GeV2, GeVc, degree, degree, GeVc, degree, degree, GeVc, degree, degree]

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimInbCR.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCR)*simDist_dvcs + contInbCR*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCR)*simDist_bh + contInbCR*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCR.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CR, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CR_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
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
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFT, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CDFT_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
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
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFD, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CD_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimInbFD.loc[(dvcsSimInbFD.phi1<30)|(dvcsSimInbFD.phi1>30), varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbFD.loc[(bkgSimInbFD.phi1<30)|(bkgSimInbFD.phi1>30), varstoplot[ind]], bins, density = True)
        simDist = (1-contInbFD)*simDist_dvcs + contInbFD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimInbFD.loc[(bhSimInbFD.phi1<30)|(bhSimInbFD.phi1>30), varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbFD)*simDist_bh + contInbFD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbFD.loc[(epgExpInbFD.phi1<30)|(epgExpInbFD.phi1>30),varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'FDFD, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/FD_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

# Read the data files
parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/dvcs/"
parent_MC_BH = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/bh/"
parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/bkg_1g/"
parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/bkg_2g/"
parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/exp/"

epgExpOutb = pd.read_pickle(parent_exp + "dvcs.pkl")
bhSimOutb = pd.concat([pd.read_pickle(parent_MC_BH + "4902.pkl")
    , pd.read_pickle(parent_MC_BH + "4903.pkl"), pd.read_pickle(parent_MC_BH + "4905.pkl")])
dvcsSimOutb = pd.concat([pd.read_pickle(parent_MC + "4907.pkl")
    , pd.read_pickle(parent_MC + "4909.pkl"), pd.read_pickle(parent_MC + "4912.pkl")])
bkgSimOutb = pd.read_pickle(parent_MC_bkg1g + "4243.pkl")
pi0ExpOutb = pd.read_pickle(parent_exp + "pi0.pkl")
pi0SimOutb = pd.read_pickle(parent_MC_bkg2g + "4243.pkl")

epgExpOutbCR = epgExpOutb.loc[epgExpOutb.config == 4]
bhSimOutbCR = bhSimOutb.loc[bhSimOutb.config == 4]
dvcsSimOutbCR = dvcsSimOutb.loc[dvcsSimOutb.config == 4]
bkgSimOutbCR = bkgSimOutb.loc[bkgSimOutb.config == 4]

epgExpOutbCDFT = epgExpOutb.loc[epgExpOutb.config == 3]
bhSimOutbCDFT = bhSimOutb.loc[bhSimOutb.config == 3]
dvcsSimOutbCDFT = dvcsSimOutb.loc[dvcsSimOutb.config == 3]
bkgSimOutbCDFT = bkgSimOutb.loc[bkgSimOutb.config == 3]
pi0ExpOutbCDFT = pi0ExpOutb.loc[(pi0ExpOutb.config == 3)]
pi0SimOutbCDFT = pi0SimOutb.loc[(pi0SimOutb.config == 3)]

epgExpOutbCD = epgExpOutb.loc[epgExpOutb.config == 2]
bhSimOutbCD = bhSimOutb.loc[bhSimOutb.config == 2]
dvcsSimOutbCD = dvcsSimOutb.loc[dvcsSimOutb.config == 2]
bkgSimOutbCD = bkgSimOutb.loc[bkgSimOutb.config == 2]
pi0ExpOutbCD = pi0ExpOutb.loc[(pi0ExpOutb.config == 2)]
pi0SimOutbCD = pi0SimOutb.loc[(pi0SimOutb.config == 2)]

epgExpOutbFD = epgExpOutb.loc[epgExpOutb.config == 1]
bhSimOutbFD = bhSimOutb.loc[bhSimOutb.config == 1]
dvcsSimOutbFD = dvcsSimOutb.loc[dvcsSimOutb.config == 1]
bkgSimOutbFD = bkgSimOutb.loc[bkgSimOutb.config == 1]
pi0ExpOutbFD = pi0ExpOutb.loc[(pi0ExpOutb.config == 1)]
pi0SimOutbFD = pi0SimOutb.loc[(pi0SimOutb.config == 1)]

contOutbCR = 0
contOutbCDFT = len(pi0ExpOutbCDFT)*len(bkgSimOutbCDFT)/len(pi0SimOutbCDFT)/len(epgExpOutbCDFT)
contOutbCD = (len(pi0ExpOutbCD.loc[(pi0ExpOutbCD.phi1<30)|(pi0ExpOutbCD.phi1>330)]) *
            len(bkgSimOutbCD.loc[(bkgSimOutbCD.phi1<30)|(bkgSimOutbCD.phi1>330)]) /
            len(pi0SimOutbCD.loc[(pi0SimOutbCD.phi1<30)|(pi0SimOutbCD.phi1>330)]) /
            len(epgExpOutbCD.loc[(epgExpOutbCD.phi1<30)|(epgExpOutbCD.phi1>330)]))
contOutbFD = (len(pi0ExpOutbFD.loc[(pi0ExpOutbFD.phi1<30)|(pi0ExpOutbFD.phi1>330)]) *
            len(bkgSimOutbFD.loc[(bkgSimOutbFD.phi1<30)|(bkgSimOutbFD.phi1>330)]) /
            len(pi0SimOutbFD.loc[(pi0SimOutbFD.phi1<30)|(pi0SimOutbFD.phi1>330)]) /
            len(epgExpOutbFD.loc[(epgExpOutbFD.phi1<30)|(epgExpOutbFD.phi1>330)]))

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCR.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCR)*simDist_dvcs + contOutbCR*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCR)*simDist_bh + contOutbCR*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCR.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CR, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CR_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCDFT.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCDFT)*simDist_dvcs + contOutbCDFT*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCDFT)*simDist_bh + contOutbCDFT*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCDFT.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFT, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CDFT_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCD.loc[(dvcsSimOutbCD.phi1<30)|(dvcsSimOutbCD.phi1>30), varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCD.loc[(bkgSimOutbCD.phi1<30)|(bkgSimOutbCD.phi1>30), varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCD)*simDist_dvcs + contOutbCD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimOutbCD.loc[(bhSimOutbCD.phi1<30)|(bhSimOutbCD.phi1>30), varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCD)*simDist_bh + contOutbCD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCD.loc[(epgExpOutbCD.phi1<30)|(epgExpOutbCD.phi1>30),varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFD, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CD_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 3, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
for yind in range(0, 4):
    for xind in range(0, 3):
        ind = 3*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbFD.loc[(dvcsSimOutbFD.phi1<30)|(dvcsSimOutbFD.phi1>30), varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbFD.loc[(bkgSimOutbFD.phi1<30)|(bkgSimOutbFD.phi1>30), varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbFD)*simDist_dvcs + contOutbFD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimOutbFD.loc[(bhSimOutbFD.phi1<30)|(bhSimOutbFD.phi1>30), varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbFD)*simDist_bh + contOutbFD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbFD.loc[(epgExpOutbFD.phi1<30)|(epgExpOutbFD.phi1>30),varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'FDFD, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/FD_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()
