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

epgExpInb  = epgExpInb.loc [(epgExpInb.phi1<30)   | (epgExpInb.phi1>330)]
bhSimInb   = bhSimInb.loc  [(bhSimInb.phi1<30)    | (bhSimInb.phi1 >330)]
dvcsSimInb = dvcsSimInb.loc[(dvcsSimInb.phi1<30)  | (dvcsSimInb.phi1>330)]
bkgSimInb  = bkgSimInb.loc [(bkgSimInb.phi1<30)   | (bkgSimInb.phi1>330)]
pi0ExpInb  = pi0ExpInb.loc [(pi0ExpInb.phi1<30)   | (pi0ExpInb.phi1>330)]
pi0SimInb  = pi0SimInb.loc [(pi0SimInb.phi1<30)   | (pi0SimInb.phi1>330)]

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
contInbCD   = len(pi0ExpInbCD)*len(bkgSimInbCD)/len(pi0SimInbCD)/len(epgExpInbCD)
contInbFD   = len(pi0ExpInbFD)*len(bkgSimInbFD)/len(pi0SimInbFD)/len(epgExpInbFD)

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

epgExpOutb  = epgExpOutb.loc [(epgExpOutb.phi1<30)   | (epgExpOutb.phi1>330)]
bhSimOutb   = bhSimOutb.loc  [(bhSimOutb.phi1<30)    | (bhSimOutb.phi1 >330)]
dvcsSimOutb = dvcsSimOutb.loc[(dvcsSimOutb.phi1<30)  | (dvcsSimOutb.phi1>330)]
bkgSimOutb  = bkgSimOutb.loc [(bkgSimOutb.phi1<30)   | (bkgSimOutb.phi1>330)]
pi0ExpOutb  = pi0ExpOutb.loc [(pi0ExpOutb.phi1<30)   | (pi0ExpOutb.phi1>330)]
pi0SimOutb  = pi0SimOutb.loc [(pi0SimOutb.phi1<30)   | (pi0SimOutb.phi1>330)]

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
contOutbCD = len(pi0ExpOutbCD)*len(bkgSimOutbCD)/len(pi0SimOutbCD)/len(epgExpOutbCD)
contOutbFD = len(pi0ExpOutbFD)*len(bkgSimOutbFD)/len(pi0SimOutbFD)/len(epgExpOutbFD)

#particle kinematics
varstoplot = ["xB", "Q2", "t1", "phi1", "Ep", "Etheta", "Ephi", "Pp", "Ptheta", "Pphi", "Gp", "Gtheta", "Gphi"]
title = [r"$x_B$", r"$Q^2$", r"$|t|$", r"$\phi_{H\Gamma}$", r"$p_{e'}$", r"$\theta_{e'}$", r"$\phi_{e'}$", r"$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", r"$p_{\gamma}$", r"$\theta_{\gamma}$", r"$\phi_{\gamma}$"]
binstoplot_CR = [np.linspace(0.05, 0.7, 101), np.linspace(1, 2.5, 101), np.linspace(0.6, 1.4, 101), 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
binstoplot = [np.linspace(0.05, 0.7, 101), np.linspace(1, 7, 101), np.linspace(0, 1, 101), 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
unit = ["", GeVc2, GeV2, degree, GeVc, degree, degree, GeVc, degree, degree, GeVc, degree, degree]

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimInbCR.loc[:, varstoplot[ind]], binstoplot_CR[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCR)*simDist_dvcs + contInbCR*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCR)*simDist_bh + contInbCR*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCR.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CR, Inb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/CR_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimInbCDFT.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCDFT)*simDist_dvcs + contInbCDFT*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCDFT)*simDist_bh + contInbCDFT*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCDFT.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFT, Inb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/CDFT_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimInbCD.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCD)*simDist_dvcs + contInbCD*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCD)*simDist_bh + contInbCD*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFD, Inb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/CD_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimInbFD.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbFD)*simDist_dvcs + contInbFD*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbFD)*simDist_bh + contInbFD*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbFD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'FDFD, Inb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/FD_Inb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCR.loc[:, varstoplot[ind]], binstoplot_CR[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCR)*simDist_dvcs + contOutbCR*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCR)*simDist_bh + contOutbCR*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCR.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CR, Outb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/CR_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCDFT.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCDFT)*simDist_dvcs + contOutbCDFT*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCDFT)*simDist_bh + contOutbCDFT*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCDFT.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFT, Outb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/CDFT_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCD.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCD)*simDist_dvcs + contOutbCD*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCD)*simDist_bh + contOutbCD*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFD, Outb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/CD_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(4, 4, figsize = (18, 30))
fig.subplots_adjust(wspace = .3, hspace = .3)
ind = -1
for yind in range(0, 4):
    for xind in range(0, 4):
        if (xind==3) and (yind>0):
            axs[yind, xind].yaxis.set_visible(False)
            axs[yind, xind].xaxis.set_visible(False)
            axs[yind, xind].axis('off')
            continue
        ind += 1
        simDist_dvcs, bins = np.histogram(dvcsSimOutbFD.loc[:, varstoplot[ind]], binstoplot[ind], density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbFD)*simDist_dvcs + contOutbFD*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbFD)*simDist_bh + contOutbFD*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbFD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'FDFD, Outb.', title_fontsize = 30, bbox_to_anchor = (0.7, 0.5))
plt.savefig("plots/normalization/FD_Outb_particle_kine.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

#exclusivity
varstoplot = dvcsvars
title = dvcstitles
unit = dvcsunits

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimInbCR.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCR)*simDist_dvcs + contInbCR*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCR)*simDist_bh + contInbCR*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCR.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CR, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CR_Inb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimInbCDFT.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCDFT)*simDist_dvcs + contInbCDFT*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimInbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCDFT)*simDist_bh + contInbCDFT*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCDFT.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFT, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CDFT_Inb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimInbCD.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCD)*simDist_dvcs + contInbCD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimInbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbCD)*simDist_bh + contInbCD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbCD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFD, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CD_Inb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimInbFD.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimInbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbFD)*simDist_dvcs + contInbFD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimInbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contInbFD)*simDist_bh + contInbFD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpInbFD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'FDFD, Inb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/FD_Inb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCR.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCR)*simDist_dvcs + contOutbCR*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCR.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCR)*simDist_bh + contOutbCR*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCR.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CR, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CR_Outb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCDFT.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCDFT)*simDist_dvcs + contOutbCDFT*simDist_dvpi0
        simDist_bh, _ = np.histogram(bhSimOutbCDFT.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCDFT)*simDist_bh + contOutbCDFT*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCDFT.loc[:, varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFT, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CDFT_Outb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbCD.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbCD)*simDist_dvcs + contOutbCD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimOutbCD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbCD)*simDist_bh + contOutbCD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbCD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'CDFD, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/CD_Outb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

fig, axs = plt.subplots(2, 4, figsize = (16, 8))
fig.subplots_adjust(wspace = .3, hspace = .6)
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        simDist_dvcs, bins = np.histogram(dvcsSimOutbFD.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(bkgSimOutbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contOutbFD)*simDist_dvcs + contOutbFD*simDist_dvpi0

        simDist_bh, _ = np.histogram(bhSimOutbFD.loc[:, varstoplot[ind]], bins, density = True)
        simDist2 = (1-contOutbFD)*simDist_bh + contOutbFD*simDist_dvpi0

        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].hist(epgExpOutbFD.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = "Experimental Data")
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = 'Simulation based on VGG')
        axs[yind, xind].hist(bins[:-1], bins, weights = simDist2, histtype = 'step', color='g', linewidth=1, label = 'Simulation based on pure BH')
        axs[yind, xind].set_title(title[ind])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
        axs[yind, xind].set_xlim([bins[0], bins[-1]])
handles, labels = axs[0, 0].get_legend_handles_labels()
lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title = 'FDFD, Outb.', title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
plt.savefig("plots/normalization/FD_Outb_excl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
plt.clf()

print("\& Sector \& Polarity \& CR \& (CD, FT) \& (CD, FD), \& (FD, FD)" + "\\\\")
print("Contamination Ratio        \& All  \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(contInbCR, contInbCDFT, contInbCD, contInbFD) + "\\\\")
print("Experimental Data          \&  1   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpInbCR.Esector==1)/len(epgExpInbCR), sum(epgExpInbCDFT.Esector==1)/len(epgExpInbCDFT)
    , sum(epgExpInbCD.Esector==1)/len(epgExpInbCD), sum(epgExpInbFD.Esector==1)/len(epgExpInbFD)) + "\\\\")
print("Experimental Data          \&  2   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpInbCR.Esector==2)/len(epgExpInbCR), sum(epgExpInbCDFT.Esector==2)/len(epgExpInbCDFT)
    , sum(epgExpInbCD.Esector==2)/len(epgExpInbCD), sum(epgExpInbFD.Esector==2)/len(epgExpInbFD)) + "\\\\")
print("Experimental Data          \&  3   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpInbCR.Esector==3)/len(epgExpInbCR), sum(epgExpInbCDFT.Esector==3)/len(epgExpInbCDFT)
    , sum(epgExpInbCD.Esector==3)/len(epgExpInbCD), sum(epgExpInbFD.Esector==3)/len(epgExpInbFD)) + "\\\\")
print("Experimental Data          \&  4   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpInbCR.Esector==4)/len(epgExpInbCR), sum(epgExpInbCDFT.Esector==4)/len(epgExpInbCDFT)
    , sum(epgExpInbCD.Esector==4)/len(epgExpInbCD), sum(epgExpInbFD.Esector==4)/len(epgExpInbFD)) + "\\\\")
print("Experimental Data          \&  5   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpInbCR.Esector==5)/len(epgExpInbCR), sum(epgExpInbCDFT.Esector==5)/len(epgExpInbCDFT)
    , sum(epgExpInbCD.Esector==5)/len(epgExpInbCD), sum(epgExpInbFD.Esector==5)/len(epgExpInbFD)) + "\\\\")
print("Experimental Data          \&  6   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpInbCR.Esector==6)/len(epgExpInbCR), sum(epgExpInbCDFT.Esector==6)/len(epgExpInbCDFT)
    , sum(epgExpInbCD.Esector==6)/len(epgExpInbCD), sum(epgExpInbFD.Esector==6)/len(epgExpInbFD)) + "\\\\")

print("BH Simulation          \&  1   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimInbCR.Esector==1)/len(bhSimInbCR), sum(bhSimInbCDFT.Esector==1)/len(bhSimInbCDFT)
    , sum(bhSimInbCD.Esector==1)/len(bhSimInbCD), sum(bhSimInbFD.Esector==1)/len(bhSimInbFD))  + "\\\\")
print("BH Simulation          \&  2   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimInbCR.Esector==2)/len(bhSimInbCR), sum(bhSimInbCDFT.Esector==2)/len(bhSimInbCDFT)
    , sum(bhSimInbCD.Esector==2)/len(bhSimInbCD), sum(bhSimInbFD.Esector==2)/len(bhSimInbFD))  + "\\\\")
print("BH Simulation          \&  3   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimInbCR.Esector==3)/len(bhSimInbCR), sum(bhSimInbCDFT.Esector==3)/len(bhSimInbCDFT)
    , sum(bhSimInbCD.Esector==3)/len(bhSimInbCD), sum(bhSimInbFD.Esector==3)/len(bhSimInbFD))  + "\\\\")
print("BH Simulation          \&  4   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimInbCR.Esector==4)/len(bhSimInbCR), sum(bhSimInbCDFT.Esector==4)/len(bhSimInbCDFT)
    , sum(bhSimInbCD.Esector==4)/len(bhSimInbCD), sum(bhSimInbFD.Esector==4)/len(bhSimInbFD))  + "\\\\")
print("BH Simulation          \&  5   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimInbCR.Esector==5)/len(bhSimInbCR), sum(bhSimInbCDFT.Esector==5)/len(bhSimInbCDFT)
    , sum(bhSimInbCD.Esector==5)/len(bhSimInbCD), sum(bhSimInbFD.Esector==5)/len(bhSimInbFD))  + "\\\\")
print("BH Simulation          \&  6   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimInbCR.Esector==6)/len(bhSimInbCR), sum(bhSimInbCDFT.Esector==6)/len(bhSimInbCDFT)
    , sum(bhSimInbCD.Esector==6)/len(bhSimInbCD), sum(bhSimInbFD.Esector==6)/len(bhSimInbFD))  + "\\\\")

print("VGG Simulation          \&  1   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimInbCR.Esector==1)/len(dvcsSimInbCR), sum(dvcsSimInbCDFT.Esector==1)/len(dvcsSimInbCDFT)
    , sum(dvcsSimInbCD.Esector==1)/len(dvcsSimInbCD), sum(dvcsSimInbFD.Esector==1)/len(dvcsSimInbFD)) + "\\\\")
print("VGG Simulation          \&  2   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimInbCR.Esector==2)/len(dvcsSimInbCR), sum(dvcsSimInbCDFT.Esector==2)/len(dvcsSimInbCDFT)
    , sum(dvcsSimInbCD.Esector==2)/len(dvcsSimInbCD), sum(dvcsSimInbFD.Esector==2)/len(dvcsSimInbFD)) + "\\\\")
print("VGG Simulation          \&  3   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimInbCR.Esector==3)/len(dvcsSimInbCR), sum(dvcsSimInbCDFT.Esector==3)/len(dvcsSimInbCDFT)
    , sum(dvcsSimInbCD.Esector==3)/len(dvcsSimInbCD), sum(dvcsSimInbFD.Esector==3)/len(dvcsSimInbFD)) + "\\\\")
print("VGG Simulation          \&  4   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimInbCR.Esector==4)/len(dvcsSimInbCR), sum(dvcsSimInbCDFT.Esector==4)/len(dvcsSimInbCDFT)
    , sum(dvcsSimInbCD.Esector==4)/len(dvcsSimInbCD), sum(dvcsSimInbFD.Esector==4)/len(dvcsSimInbFD)) + "\\\\")
print("VGG Simulation          \&  5   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimInbCR.Esector==5)/len(dvcsSimInbCR), sum(dvcsSimInbCDFT.Esector==5)/len(dvcsSimInbCDFT)
    , sum(dvcsSimInbCD.Esector==5)/len(dvcsSimInbCD), sum(dvcsSimInbFD.Esector==5)/len(dvcsSimInbFD)) + "\\\\")
print("VGG Simulation          \&  6   \& Inb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimInbCR.Esector==6)/len(dvcsSimInbCR), sum(dvcsSimInbCDFT.Esector==6)/len(dvcsSimInbCDFT)
    , sum(dvcsSimInbCD.Esector==6)/len(dvcsSimInbCD), sum(dvcsSimInbFD.Esector==6)/len(dvcsSimInbFD)) + "\\\\")

print("Contamination Ratio        \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(contOutbCR, contOutbCDFT, contOutbCD, contOutbFD) + "\\\\")
print("Experimental Data          \&  1   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpOutbCR.Esector==1)/len(epgExpOutbCR), sum(epgExpOutbCDFT.Esector==1)/len(epgExpOutbCDFT)
    , sum(epgExpOutbCD.Esector==1)/len(epgExpOutbCD), sum(epgExpOutbFD.Esector==1)/len(epgExpOutbFD)) + "\\\\")
print("Experimental Data          \&  2   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpOutbCR.Esector==2)/len(epgExpOutbCR), sum(epgExpOutbCDFT.Esector==2)/len(epgExpOutbCDFT)
    , sum(epgExpOutbCD.Esector==2)/len(epgExpOutbCD), sum(epgExpOutbFD.Esector==2)/len(epgExpOutbFD)) + "\\\\")
print("Experimental Data          \&  3   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpOutbCR.Esector==3)/len(epgExpOutbCR), sum(epgExpOutbCDFT.Esector==3)/len(epgExpOutbCDFT)
    , sum(epgExpOutbCD.Esector==3)/len(epgExpOutbCD), sum(epgExpOutbFD.Esector==3)/len(epgExpOutbFD)) + "\\\\")
print("Experimental Data          \&  4   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpOutbCR.Esector==4)/len(epgExpOutbCR), sum(epgExpOutbCDFT.Esector==4)/len(epgExpOutbCDFT)
    , sum(epgExpOutbCD.Esector==4)/len(epgExpOutbCD), sum(epgExpOutbFD.Esector==4)/len(epgExpOutbFD)) + "\\\\")
print("Experimental Data          \&  5   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpOutbCR.Esector==5)/len(epgExpOutbCR), sum(epgExpOutbCDFT.Esector==5)/len(epgExpOutbCDFT)
    , sum(epgExpOutbCD.Esector==5)/len(epgExpOutbCD), sum(epgExpOutbFD.Esector==5)/len(epgExpOutbFD)) + "\\\\")
print("Experimental Data          \&  6   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(epgExpOutbCR.Esector==6)/len(epgExpOutbCR), sum(epgExpOutbCDFT.Esector==6)/len(epgExpOutbCDFT)
    , sum(epgExpOutbCD.Esector==6)/len(epgExpOutbCD), sum(epgExpOutbFD.Esector==6)/len(epgExpOutbFD)) + "\\\\")

print("BH Simulation          \&  1   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimOutbCR.Esector==1)/len(bhSimOutbCR), sum(bhSimOutbCDFT.Esector==1)/len(bhSimOutbCDFT)
    , sum(bhSimOutbCD.Esector==1)/len(bhSimOutbCD), sum(bhSimOutbFD.Esector==1)/len(bhSimOutbFD)) + "\\\\")
print("BH Simulation          \&  2   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimOutbCR.Esector==2)/len(bhSimOutbCR), sum(bhSimOutbCDFT.Esector==2)/len(bhSimOutbCDFT)
    , sum(bhSimOutbCD.Esector==2)/len(bhSimOutbCD), sum(bhSimOutbFD.Esector==2)/len(bhSimOutbFD)) + "\\\\")
print("BH Simulation          \&  3   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimOutbCR.Esector==3)/len(bhSimOutbCR), sum(bhSimOutbCDFT.Esector==3)/len(bhSimOutbCDFT)
    , sum(bhSimOutbCD.Esector==3)/len(bhSimOutbCD), sum(bhSimOutbFD.Esector==3)/len(bhSimOutbFD)) + "\\\\")
print("BH Simulation          \&  4   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimOutbCR.Esector==4)/len(bhSimOutbCR), sum(bhSimOutbCDFT.Esector==4)/len(bhSimOutbCDFT)
    , sum(bhSimOutbCD.Esector==4)/len(bhSimOutbCD), sum(bhSimOutbFD.Esector==4)/len(bhSimOutbFD)) + "\\\\")
print("BH Simulation          \&  5   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimOutbCR.Esector==5)/len(bhSimOutbCR), sum(bhSimOutbCDFT.Esector==5)/len(bhSimOutbCDFT)
    , sum(bhSimOutbCD.Esector==5)/len(bhSimOutbCD), sum(bhSimOutbFD.Esector==5)/len(bhSimOutbFD)) + "\\\\")
print("BH Simulation          \&  6   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(bhSimOutbCR.Esector==6)/len(bhSimOutbCR), sum(bhSimOutbCDFT.Esector==6)/len(bhSimOutbCDFT)
    , sum(bhSimOutbCD.Esector==6)/len(bhSimOutbCD), sum(bhSimOutbFD.Esector==6)/len(bhSimOutbFD)) + "\\\\")

print("VGG Simulation         \&  1   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimOutbCR.Esector==1)/len(dvcsSimOutbCR), sum(dvcsSimOutbCDFT.Esector==1)/len(dvcsSimOutbCDFT)
    , sum(dvcsSimOutbCD.Esector==1)/len(dvcsSimOutbCD), sum(dvcsSimOutbFD.Esector==1)/len(dvcsSimOutbFD)) + "\\\\")
print("VGG Simulation         \&  2   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimOutbCR.Esector==2)/len(dvcsSimOutbCR), sum(dvcsSimOutbCDFT.Esector==2)/len(dvcsSimOutbCDFT)
    , sum(dvcsSimOutbCD.Esector==2)/len(dvcsSimOutbCD), sum(dvcsSimOutbFD.Esector==2)/len(dvcsSimOutbFD)) + "\\\\")
print("VGG Simulation         \&  3   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimOutbCR.Esector==3)/len(dvcsSimOutbCR), sum(dvcsSimOutbCDFT.Esector==3)/len(dvcsSimOutbCDFT)
    , sum(dvcsSimOutbCD.Esector==3)/len(dvcsSimOutbCD), sum(dvcsSimOutbFD.Esector==3)/len(dvcsSimOutbFD)) + "\\\\")
print("VGG Simulation         \&  4   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimOutbCR.Esector==4)/len(dvcsSimOutbCR), sum(dvcsSimOutbCDFT.Esector==4)/len(dvcsSimOutbCDFT)
    , sum(dvcsSimOutbCD.Esector==4)/len(dvcsSimOutbCD), sum(dvcsSimOutbFD.Esector==4)/len(dvcsSimOutbFD)) + "\\\\")
print("VGG Simulation         \&  5   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimOutbCR.Esector==5)/len(dvcsSimOutbCR), sum(dvcsSimOutbCDFT.Esector==5)/len(dvcsSimOutbCDFT)
    , sum(dvcsSimOutbCD.Esector==5)/len(dvcsSimOutbCD), sum(dvcsSimOutbFD.Esector==5)/len(dvcsSimOutbFD)) + "\\\\")
print("VGG Simulation         \&  6   \& Outb. \& {:.2e} \& {:.2e} \& {:.2e}\& {:.2e}".format(sum(dvcsSimOutbCR.Esector==6)/len(dvcsSimOutbCR), sum(dvcsSimOutbCDFT.Esector==6)/len(dvcsSimOutbCDFT)
    , sum(dvcsSimOutbCD.Esector==6)/len(dvcsSimOutbCD), sum(dvcsSimOutbFD.Esector==6)/len(dvcsSimOutbFD)) + "\\\\")


# Q2 [1.000, 1.200, 1.456, 1.912, 2.510, 3.295, 4.326, 5.761, 7.000]
# xB [0.062, 0.090, 0.118, 0.155, 0.204, 0.268, 0.357, 0.446, 0.581]
# t  [0.11, 0.15, 0.25, 0.4, 0.6, 0.8, 1.0]

# Find the bin that is shared by (CD, FT) and (CD, FD): (3, 3, 0), (3, 3, 1), (4, 4, 1), (4, 4, 2)

# (3, 3, 0)
# 0.155 < xB < 0.204
# 1.912 < Q2 < 2.510
# 0.11  < |t| < 0.15

# (3, 3, 1)
# 0.155 < xB < 0.204
# 1.912 < Q2 < 2.510
# 0.15  < |t| < 0.25

# (4, 4, 1)
# 0.204 < xB < 0.268
# 2.510 < Q2 < 3.295
# 0.15  < |t| < 0.25

# (4, 4, 2)
# 0.204 < xB < 0.268
# 2.510 < Q2 < 3.295
# 0.25  < |t| < 0.4

# def pickOneBinOneSector(df, xBbin, Q2bin, tbin, sector = 0):
#     xBcond = (df.xB>=newxBbins2[xBbin]) & (df.xB<newxBbins2[xBbin+1])
#     Q2cond = (df.Q2>=newQ2bins2[Q2bin]) & (df.Q2<newQ2bins2[Q2bin+1])
#     tcond  = (df.t1>=newtbins[tbin])    & (df.t1<newtbins[tbin+1])
#     print(sum(xBcond), sum(Q2cond), sum(tcond))
#     cond   = xBcond & Q2cond & tcond
#     if sector > 0:
#         sectorcond = df.Esector == sector
#         cond   = cond & sectorcond
#     return df.loc[cond, :]

# plt.hist(pickOneBinOneSector(dvcsSimOutbCR, 3, 3, 5).Etheta, bins = 100)

# dist, bins = np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==1, "Etheta"], bins = 100)
# plt.hist(bins[:-1], bins, weights = dist/np.diff(bins)/len(dvcsSimOutbCDFT), histtype = 'step')
# dist, _    = np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==2, "Etheta"], bins = bins)
# plt.hist(bins[:-1], bins, weights = dist/np.diff(bins)/len(dvcsSimOutbCDFT), histtype = 'step')
# plt.hist2d(bhSimOutbCDFT.t1, bhSimOutbCDFT.Pphi, bins = 100, norm = LogNorm())
# plt.hist2d(epgExpOutbCDFT.t1, epgExpOutbCDFT.Pphi, bins = 100, norm = LogNorm())
# plt.show()


fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.Etheta>10, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .Etheta>10, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulaton based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$\theta_{e'}>10$" +" "+ degree, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_highEtheta_Ephi.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.Etheta<10, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .Etheta<10, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$\theta_{e'}<10$" +" "+ degree, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_lowEtheta_Ephi.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.Etheta>10, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .Etheta>10, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulaton based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$\theta_{e'}>10$" +" "+ degree, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_highEtheta_Pphi.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.Etheta<10, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .Etheta<10, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$\theta_{e'}<10$" +" "+ degree, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_lowEtheta_Pphi.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.t1>0.6, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .t1>0.6, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulaton based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$p_{p'}>0.6$" +" "+ GeVc, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_hight_Ephi.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.t1<0.6, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .t1<0.6, "Ephi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$p_{p'}<0.6$" +" "+ GeVc, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_lowt_Ephi.pdf")
plt.clf()


fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.t1>0.6, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .t1>0.6, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulaton based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$p_{p'}>0.6$" +" "+ GeVc, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_hight_Pphi.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[epgExpOutbCD.t1<0.6, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[bhSimOutbCD .t1<0.6, "Pphi"].hist(ax = ax, bins = np.linspace(-180, 180, 91), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
ax.set_xlim([-180, 180])
ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$\phi_{p'}$"+ " ["+ degree + "]")
plt.legend(loc='lower left', title = r"$p_{p'}<0.6$" +" "+ GeVc, bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_lowt_Pphi.pdf")
plt.clf()

print("Experimental Data &  Inb.  & {:.2e} & {:.2e} & {:.2e}\\\\".format(len(epgExpInbCDFT.loc[epgExpInbCDFT.t1<0.6, :]), len(epgExpInbCD.loc[epgExpInbCD.t1<0.6, :]), len(epgExpInbFD.loc[epgExpInbFD.t1<0.6, :])))
print("BH Simulation &  Inb.  & {:.2e} & {:.2e} & {:.2e}\\\\".format(len(bhSimInbCDFT.loc[bhSimInbCDFT.t1<0.6, :]), len(bhSimInbCD.loc[bhSimInbCD.t1<0.6, :]), len(bhSimInbFD.loc[bhSimInbFD.t1<0.6, :])))
print("Experimental Data &  Outb. & {:.2e} & {:.2e} & {:.2e}\\\\".format(len(epgExpOutbCDFT.loc[epgExpOutbCDFT.t1<0.6, :]), len(epgExpOutbCD.loc[epgExpOutbCD.t1<0.6, :]), len(epgExpOutbFD.loc[epgExpOutbFD.t1<0.6, :])))
print("BH Simulation &  Outb. & {:.2e} & {:.2e} & {:.2e}\\\\".format(len(bhSimOutbCDFT.loc[bhSimOutbCDFT.t1<0.6, :]), len(bhSimOutbCD.loc[bhSimOutbCD.t1<0.6, :]), len(bhSimOutbFD.loc[bhSimOutbFD.t1<0.6, :])))

print("Experimental Data &  Inb.  & {:.2e} & {:.2e}\\\\".format(len(epgExpInbCD)/len(epgExpInbCDFT), len(epgExpInbFD)/len(epgExpInbCDFT)))
print("BH Simulation &  Inb.  & {:.2e} & {:.2e}\\\\".format(len(bhSimInbCD)/len(bhSimInbCDFT), len(bhSimInbFD)/len(bhSimInbCDFT)))
print("Experimental Data &  Outb. & {:.2e} & {:.2e}\\\\".format(len(epgExpOutbCD)/len(epgExpOutbCDFT), len(epgExpOutbFD)/len(epgExpOutbCDFT)))
print("BH Simulation &  Outb. & {:.2e} & {:.2e}\\\\".format(len(bhSimOutbCD)/len(bhSimOutbCDFT), len(bhSimOutbFD)/len(bhSimOutbCDFT)))

fig, ax = plt.subplots(1,1, figsize=(8, 5))
ax.hist2d(epgExpOutbCD.t1, epgExpOutbCD.t1 - epgExpOutbCD.t1Orig, bins = [np.linspace(0, 1, 100), 100], cmap = cmap, cmin = 1)
ax.set_xlabel(r"$|t|_{after}$"+ " ["+ GeV2 + "] (after correction)")
ax.set_ylabel(r"$|t|_{after}-|t|_{before}$"+ " ["+ GeV2 + "]")
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_t1diff.pdf")
plt.clf()

bhSimOutbCD.loc[:, "PpxOrig"] = bhSimOutbCD.loc[:, "PpOrig"]*np.sin(np.radians(bhSimOutbCD.loc[:, "PthetaOrig"]))*np.cos(np.radians(bhSimOutbCD.loc[:, "PphiOrig"]))
bhSimOutbCD.loc[:, "PpyOrig"] = bhSimOutbCD.loc[:, "PpOrig"]*np.sin(np.radians(bhSimOutbCD.loc[:, "PthetaOrig"]))*np.sin(np.radians(bhSimOutbCD.loc[:, "PphiOrig"]))
bhSimOutbCD.loc[:, "PpzOrig"] = bhSimOutbCD.loc[:, "PpOrig"]*np.cos(np.radians(bhSimOutbCD.loc[:, "PthetaOrig"]))

VmissG = [-bhSimOutbCD["Epx"] - bhSimOutbCD["PpxOrig"], -bhSimOutbCD["Epy"] - bhSimOutbCD["PpyOrig"],
          self.pbeam - bhSimOutbCD["Epz"] - bhSimOutbCD["PpzOrig"]]
bhSimOutbCD.loc[:,'MM2_epOrig'] = (-M - self.ebeam + bhSimOutbCD["Ee"] + bhSimOutbCD["PeOrig"])**2 - mag2(VmissG)

epgExpOutbCD.loc[:, "PpxOrig"] = epgExpOutbCD.loc[:, "PpOrig"]*np.sin(np.radians(epgExpOutbCD.loc[:, "PthetaOrig"]))*np.cos(np.radians(epgExpOutbCD.loc[:, "PphiOrig"]))
epgExpOutbCD.loc[:, "PpyOrig"] = epgExpOutbCD.loc[:, "PpOrig"]*np.sin(np.radians(epgExpOutbCD.loc[:, "PthetaOrig"]))*np.sin(np.radians(epgExpOutbCD.loc[:, "PphiOrig"]))
epgExpOutbCD.loc[:, "PpzOrig"] = epgExpOutbCD.loc[:, "PpOrig"]*np.cos(np.radians(epgExpOutbCD.loc[:, "PthetaOrig"]))

VmissG = [-epgExpOutbCD["Epx"] - epgExpOutbCD["PpxOrig"], -epgExpOutbCD["Epy"] - epgExpOutbCD["PpyOrig"],
          self.pbeam - epgExpOutbCD["Epz"] - epgExpOutbCD["PpzOrig"]]
epgExpOutbCD.loc[:,'MM2_epOrig'] = (-M - self.ebeam + epgExpOutbCD["Ee"] + epgExpOutbCD["PeOrig"])**2 - mag2(VmissG)


fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[:, "t1Orig"].hist(ax = ax, bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[:, "t1Orig"].hist(ax = ax, bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'g', label = 'Simulaton based on pure BH')
ax.set_xlim([0, 1])
# ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$|t|$"+ " ["+ GeV2 + "]")
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_t1Orig.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[:, "t1Orig"].hist(ax = ax, bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[:, "t1Orig"].hist(ax = ax, bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
ax.set_xlim([0, 1])
# ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$|t|$"+ " ["+ GeV2 + "]")
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_t1Orig.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[:, "MM2ep"].hist(ax = ax, bins = np.linspace(-0.3, 0.3, 51), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[:, "MM2ep"].hist(ax = ax, bins = np.linspace(-0.3, 0.3, 51), density = True, histtype = 'step', color = 'g', label = 'Simulaton based on pure BH')
ax.set_xlim([-0.3, 0.3])
# ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$MM^2_{e'p'}$"+ " ["+ GeV2 + "]")
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_MM2ep.pdf")
plt.clf()

fig, ax = plt.subplots(1,1, figsize=(8, 5))
epgExpOutbCD.loc[:, "MM2ep"].hist(ax = ax, bins = np.linspace(-0.3, 0.3, 51), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
bhSimOutbCD .loc[:, "MM2ep"].hist(ax = ax, bins = np.linspace(-0.3, 0.3, 51), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
ax.set_xlim([-0.3, 0.3])
# ax.set_ylim([0, 0.01])
ax.set_xlabel(r"$MM^2_{e'p'}$"+ " ["+ GeV2 + "]")
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.65), framealpha = 0.5)
plt.tight_layout()
plt.savefig("plots/normalization/CDFD_Outb_MM2ep.pdf")
plt.clf()



# bhSimOutbCDFT .loc[:, "t1Orig"].hist(bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
# epgExpOutbCDFT.loc[:, "t1Orig"].hist(bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
# plt.show()

# bhSimOutbCDFT .loc[:, "t1Orig"].hist(bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'b', label = 'Experimental Data')
# epgExpOutbCDFT.loc[:, "t1Orig"].hist(bins = np.linspace(0, 1, 51), density = True, histtype = 'step', color = 'g', label = 'Simulation based on pure BH')
# plt.show()


# bhSimOutbCD.loc[(bhSimOutbCD.Pphi<-50) & (bhSimOutbCD.Pphi>-100), "t1"].hist(bins = np.linspace(0, 1.5, 100), density = True, histtype = 'step')
# epgExpOutbCD.loc[(epgExpOutbCD.Pphi<-50) & (epgExpOutbCD.Pphi>-100), "t1"].hist(bins = np.linspace(0, 1.5, 100), density = True, histtype = 'step')
# plt.show()

# bhSimOutbCD.loc[(bhSimOutbCD.Pphi<-100) & (bhSimOutbCD.Pphi>-180), "t1"].hist(bins = np.linspace(0, 1.5, 100), density = True, histtype = 'step')
# epgExpOutbCD.loc[(epgExpOutbCD.Pphi<-100) & (epgExpOutbCD.Pphi>-180), "t1"].hist(bins = np.linspace(0, 1.5, 100), density = True, histtype = 'step')
# plt.show()

# bhSimOutbCD.loc[(bhSimOutbCD.Pphi>-50), "t1"].hist(bins = np.linspace(0, 1.5, 100), density = True, histtype = 'step')
# epgExpOutbCD.loc[(epgExpOutbCD.Pphi>-50), "t1"].hist(bins = np.linspace(0, 1.5, 100), density = True, histtype = 'step')
# plt.show()

# var = "Pphi"
# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==1, var], density = True, bins = np.linspace(-180, 180, 180))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==1, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==1, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==1, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==2, var], density = True, bins = np.linspace(0, 45, 46))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==2, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==2, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==2, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==2, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==2, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==2, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==2, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==3, var], density = True, bins = np.linspace(0, 45, 46))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==3, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==3, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==3, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==3, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==3, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==3, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==3, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==4, var], density = True, bins = np.linspace(0, 45, 46))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==4, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==4, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==4, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==4, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==4, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==4, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==4, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==5, var], density = True, bins = np.linspace(0, 45, 46))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==5, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==5, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==5, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==5, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==5, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==5, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==5, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==6, var], density = True, bins = np.linspace(0, 45, 46))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==6, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==6, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==6, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==6, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==6, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==6, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==6, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# distSim, bins    = np.histogram(bhSimOutbCR.loc[bhSimOutbCR.Esector==1, var], density = True, bins = np.linspace(0, 45, 46))
# distExp, _    = np.histogram(epgExpOutbCR.loc[epgExpOutbCR.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step', color = 'k')

# distSim, bins    = np.histogram(bhSimOutbCDFT.loc[bhSimOutbCDFT.Esector==1, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCDFT.loc[epgExpOutbCDFT.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbCD.loc[bhSimOutbCD.Esector==1, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbCD.loc[epgExpOutbCD.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')

# distSim, bins    = np.histogram(bhSimOutbFD.loc[bhSimOutbFD.Esector==1, var], density = True, bins = bins)
# distExp, _    = np.histogram(epgExpOutbFD.loc[epgExpOutbFD.Esector==1, var], density = True, bins = bins)
# plt.hist(bins[:-1], bins, weights = divideHist(distExp, distSim), histtype = 'step')
# plt.yscale('log')
# plt.show()

# dist, _    = np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==4, "Etheta"], bins = bins)
# plt.hist(bins[:-1], bins, weights = dist/np.diff(bins)/len(dvcsSimOutbCDFT), histtype = 'step')
# dist, _    = np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==5, "Etheta"], bins = bins)
# plt.hist(bins[:-1], bins, weights = dist/np.diff(bins)/len(dvcsSimOutbCDFT), histtype = 'step')
# dist, _    = np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==6, "Etheta"], bins = bins)
# plt.hist(bins[:-1], bins, weights = dist/np.diff(bins)/len(dvcsSimOutbCDFT), histtype = 'step')


# np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==2, "Etheta"], bins = 100, histtype = 'step')
# np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==3, "Etheta"], bins = 100, histtype = 'step')
# np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==4, "Etheta"], bins = 100, histtype = 'step')
# np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==5, "Etheta"], bins = 100, histtype = 'step')
# np.histogram(dvcsSimOutbCDFT.loc[dvcsSimOutbCDFT.Esector==6, "Etheta"], bins = 100, histtype = 'step')