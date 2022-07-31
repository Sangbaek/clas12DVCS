#!/usr/bin/env python3
"""
Script to reproduce some of plots used in Sangbaek's thesis.
The original plots were produced in the jupyter notebook.
Some of them has problems in that
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
		"text.usetex": True,            # use LaTeX to write all text
		"font.family": "sans-serif",         
		"font.sans-serif": "Helvetica",
		"font.size": 25,				# default font size
		"axes.titlepad": 20,			# x and y label size
		"axes.labelsize": 24,			# x and y label size
		"axes.titlesize": 24,           # subfigure title size, i.e. title size when one figure
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
args = parser.parse_args()

if args.chapter:
	chapter = int(args.chapter)
else:
	print("no chapter given.")
	exit()

if args.polarity in ["inbending", "outbending"]:
	polarity = args.polarity
else:
	print("The polarity must be within inb and outb.")
	exit()
#chapter 2
if chapter == 2:
	'''
	The chapter 2 is about the fiducial cuts.
	data set 1: without fiducial cut. no momenutm correction cut.
	data set 2: with fiducial cut. no momenutm correction cut.
	'''

	ecal_e_sampl_mu = [[  0.2531,  0.2550,  0.2514,  0.2494,  0.2528,  0.2521 ],[-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703  ],[  4.939,  5.350,  5.102,  6.440,  6.149,  4.957  ]]
	ecal_e_sampl_sigm = [[ 2.726e-3,  4.157e-3,  5.222e-3,  5.398e-3,  8.453e-3,  6.533e-3 ],[1.062,  0.859,  0.5564,  0.6576,  0.3242,  0.4423],[-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274]]

	ecal_e_sampl_mu_mc = [[0.248605, 0.248605, 0.248605, 0.248605, 0.248605, 0.248605 ],
	[-0.844221, -0.844221, -0.844221, -0.844221, -0.844221, -0.844221  ],
	[ 4.87777, 4.87777, 4.87777, 4.87777, 4.87777, 4.87777  ]]
	ecal_e_sampl_sigm_mc = [[7.41575e-3, 7.41575e-3, 7.41575e-3, 7.41575e-3, 7.41575e-3, 7.41575e-3 ],
	[ 0.215861, 0.215861, 0.215861, 0.215861, 0.215861, 0.215861   ],
	[-0.319801, -0.319801, -0.319801, -0.319801, -0.319801, -0.319801]]

	if polarity == "inbending":
		# before the fiducial cuts!
		dvcsSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/dvcs/4893.pkl")
		expSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/exp/dvcs.pkl")

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))

		axs.hist(expSample.Enphe, bins = np.linspace(0, 50, 51), histtype = 'step', color = 'k', density = True, label = 'experimental')
		axs.hist(dvcsSample.Enphe, bins = np.linspace(0, 50, 51), histtype = 'step', color = 'r', density = True, label = 'simulation')
		axs.axvline(2, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlim([0, 50])
		axs.set_xticks([2, 10, 20, 30, 40, 50])
		axs.set_xticklabels([2, 10, 20, 30, 40, 50])
		axs.set_yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1])
		axs.set_yticklabels([0, 0.02, 0.04, 0.06, 0.08, 0.1])
		axs.set_xlabel("Number of Photoelectrons in HTCC (" + r"$n_{phe.})$")
		plt.tight_layout()
		plt.savefig("plots/ch2/Enphe.pdf")
		plt.clf()

		fig, axs = plt.subplots(2, 3, figsize = (15,10))
		partp = np.linspace(2, 8.6, 101)
		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*(xind) + yind + 1
				h = axs[xind, yind].hist2d(expSample.loc[expSample.Esector==Esector, "Ep"],  expSample.loc[expSample.Esector==Esector, "ESamplFrac"], cmin =1, bins = [np.linspace(2, 8.6, 100), np.linspace(0, 0.35, 100)], cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
				ticks = [1, 10, 100, 1000]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticks)
				mean = ecal_e_sampl_mu[0][Esector-1] + ecal_e_sampl_mu[1][Esector-1]/1000*pow(partp-ecal_e_sampl_mu[2][Esector-1],2);
				sigma = ecal_e_sampl_sigm[0][Esector-1] + ecal_e_sampl_sigm[1][Esector-1]/(10*(partp-ecal_e_sampl_sigm[2][Esector-1]));
				axs[xind, yind].plot(partp, mean+3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].plot(partp, mean-3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].set_title(r"$e'$" + " Samp. Frac. Sector {}".format(Esector))
				axs[xind, yind].set_xlim([2, 8.6])
				axs[xind, yind].set_ylim([0, 0.37])
				#         axs[xind, yind].set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
				#         axs[xind, yind].set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_xticks([2, 4, 6, 8, 8.6])
				axs[xind, yind].set_xticklabels([2, 4, 6, "", 8.6])
				axs[xind, yind].set_xlabel(r"$p_{e'}$" + " ["+GeVc+"]")
				axs[xind, yind].set_ylabel(r"$E_{dep.}/p_{e'}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutSampling.pdf")
		plt.clf()

		fig, axs = plt.subplots(2, 3, figsize = (15,10))
		partp = np.linspace(2, 8.6, 101)
		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*xind + yind  +1
				h = axs[xind, yind].hist2d(dvcsSample.loc[dvcsSample.Esector==Esector, "Ep"],  dvcsSample.loc[dvcsSample.Esector==Esector, "ESamplFrac"], cmin =1, bins = [np.linspace(2, 8.6, 100), np.linspace(0, 0.35, 100)], cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticks)
				mean = ecal_e_sampl_mu_mc[0][Esector-1] + ecal_e_sampl_mu_mc[1][Esector-1]/1000*pow(partp-ecal_e_sampl_mu_mc[2][Esector-1],2);
				sigma = ecal_e_sampl_sigm_mc[0][Esector-1] + ecal_e_sampl_sigm_mc[1][Esector-1]/(10*(partp-ecal_e_sampl_sigm_mc[2][Esector-1]));
				axs[xind, yind].plot(partp, mean+3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].plot(partp, mean-3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].set_title(r"$e'$" + " Samp. Frac. Sector {}".format(Esector))
				axs[xind, yind].set_xlim([2, 8.6])
				axs[xind, yind].set_ylim([0, 0.37])
				axs[xind, yind].set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_xticks([2, 4, 6, 8, 8.6])
				axs[xind, yind].set_xticklabels([2, 4, 6, "", 8.6])
				axs[xind, yind].set_xlabel(r"$p_{e'}$" + " ["+GeVc+"]")
				axs[xind, yind].set_ylabel(r"$E_{dep.}/p_{e'}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutSamplingMC.pdf")
		plt.clf()

		fig, axs = plt.subplots(2, 3, figsize = (15,10))
		partp = np.linspace(2, 8.6, 101)

		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*(xind) + yind + 1
				h  = axs[xind, yind].hist2d(expSample.loc[expSample.Esector==Esector, "Eedep1"],  expSample.loc[expSample.Esector==Esector, "Eedep2"]+expSample.loc[expSample.Esector==Esector, "Eedep3"], cmin =1, bins = np.linspace(0, 0.8, 100), cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 200))
				ticks = [1, 10, 100, 200]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticks)
				axs[xind, yind].set_title(r"$e',~E_{dep.}$" + " Sector {}".format(Esector))
				axs[xind, yind].set_xlim([0.03, 0.8])
				axs[xind, yind].set_ylim([0.03, 0.8])
				axs[xind, yind].axvline(0.07, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].set_yticks([0, 0.2, 0.4, 0.6, 0.8])
				axs[xind, yind].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8])
				axs[xind, yind].set_xticks([0, 0.07, 0.2, 0.4, 0.6, 0.8])
				axs[xind, yind].set_xticklabels([0, "", 0.2, 0.4, 0.6, 0.8])
				axs[xind, yind].set_xlabel(r"$E_{dep.,~\mathrm{PCAL}}$" + " ["+GeV+"]")
				axs[xind, yind].set_ylabel(r"$E_{dep.,~\mathrm{ECAL}}$" + " ["+GeV+"]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutEdep.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))

		axs.hist(expSample.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'k', density = True, label = 'experiment')
		axs.hist(dvcsSample.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'r', density = True, label = 'simulation')
		axs.axvline(-13, color = 'k', linestyle = '--', linewidth = 5)
		axs.axvline(12, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlim([-20, 15])
		axs.set_xticks([-20, -15, -13, -10, -5, 0, 5, 10, 12, 15])
		axs.set_xticklabels([-20, "", -13, -10, -5, 0, 5, "",12, 15])
		axs.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
		axs.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2])
		# plt.hist(dvcsSample.Enphe, bins = np.linspace(0, 50, 51), density = True, histtype = 'step')
		plt.legend(loc='upper right', bbox_to_anchor = (1.05, 0.95), title = 'Inbending', framealpha = 1)
		axs.set_xlabel("$vz_{e'}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutVz_inb.pdf")
		plt.clf()

		fig, axs = plt.subplots(2, 3, figsize = (15,10))
		x = np.linspace(0, 0.3, 101)
		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*(xind) + yind + 1
				h = axs[xind, yind].hist2d(expSample.loc[expSample.Esector==Esector, "Eedep2"]/expSample.loc[expSample.Esector==Esector, "Ep"],  expSample.loc[expSample.Esector==Esector, "Eedep1"]/expSample.loc[expSample.Esector==Esector, "Ep"], cmin =1, bins = np.linspace(0, 0.3, 101), cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
				ticks = [1, 10, 100, 1000]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticks)
				axs[xind, yind].set_title(r"$e',~E_{dep.}$" + " Sector {}".format(Esector))
				axs[xind, yind].set_xlim([0.0, 0.3])
				axs[xind, yind].set_ylim([0.0, 0.3])
				axs[xind, yind].arrow(0.2, 0.2, 0.01, 0.01, linewidth = 5)
				axs[xind, yind].plot(x, 0.2-x, linestyle = '--', color = 'k', linewidth = 5)
				axs[xind, yind].set_yticks([0, 0.1, 0.2, 0.3])
				axs[xind, yind].set_yticklabels([0, 0.1, 0.2, 0.3])
				axs[xind, yind].set_xticks([0, 0.1, 0.2, 0.3])
				axs[xind, yind].set_xticklabels(["", 0.1, 0.2, 0.3])
				axs[xind, yind].set_xlabel(r"$E_{dep.,~\mathrm{ECAL-inner}}/p_{e'}$")
				axs[xind, yind].set_ylabel(r"$E_{dep.,~\mathrm{PCAL}}/p_{e'}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutAntiPion.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EDc3Hitx, expSample.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
		ticks = [1, 10, 100, 1000, 5000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title(r"$e'$"+" DC Outmost Layer Hits, Pre-fiducial, Inbending")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidDC.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EcalV1, expSample.ESamplFrac, bins = [np.linspace(0, 30, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 500))
		ticks = [1, 10, 100, 500]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.axvline(9, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlabel(r"$l_{V}$"+ " [cm]")
		axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
		axs.set_xticks([0, 5, 9, 10, 15, 20, 25, 30])
		axs.set_xticklabels([0, 5, 9, "", 15, 20, 25, 30])
		axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
		axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidV.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EcalW1, expSample.ESamplFrac, bins = [np.linspace(0, 30, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 500))
		ticks = [1, 10, 100, 500]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.axvline(9, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlabel(r"$l_{W}$"+ " [cm]")
		axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
		axs.set_xticks([0, 5, 9, 10, 15, 20, 25, 30])
		axs.set_xticklabels([0, 5, 9, "", 15, 20, 25, 30])
		axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
		axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidW.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config>1].Pp, expSample.loc[expSample.config>1].Pchi2pid, bins = [np.linspace(0, 1.8, 101), np.linspace(-6, 6, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 500))
		ticks = [1, 10, 100, 500]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xticks([0, 0.5, 1, 1.5, 2])
		axs.set_xticklabels([0, 0.5, 1, 1.5, 2])
		axs.set_yticks([-6, -4, -2 ,0, 2, 4, 6])
		axs.set_yticklabels([-6, -4, -2 ,0, 2, 4, 6])
		axs.set_xlabel(r"$p_{p'}$"+ " ["+GeVc+"]")
		axs.set_ylabel(r"$\chi$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidchi.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.PDc3Hitx, expSample.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticks = [1, 10, 100, 300]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_title(r"$p'$"+" DC Outmost Layer Hits, Pre-fiducial")
		# axs.set_xticks([-400, -200, 0, 200, 400])
		# axs.set_yticks([-400, -200, 0, 200, 400])
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidDC.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		expSample.loc[(expSample.config>1), "Ptheta"].hist(bins = np.linspace(30, 80, 101), histtype = 'step', density = True, ax = axs, color = 'k', label = 'experiment')
		dvcsSample.loc[(dvcsSample.config>1), "Ptheta"].hist(bins = np.linspace(30, 80, 101), histtype = 'step', density = True, ax = axs, color = 'r', label = 'simulation')
		plt.axvline(64.23, color = 'k', linestyle = '--', linewidth = 2)
		plt.legend(loc='upper left', bbox_to_anchor = (0.0, 0.95), framealpha = 1)
		plt.xlabel(r"$\theta_{p'}$" + " ["+degree+"]")
		axs.set_xticks([30, 40, 50, 60, 64.23, 70, 80])
		axs.set_xticklabels([30, 40, 50, "", 64.23, 70, 80])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		expSample.loc[(expSample.config>1), "PCvt12theta"].hist(bins = np.linspace(40, 90, 101), histtype = 'step', density = True, ax = axs, color = 'k', label = 'experiment')
		dvcsSample.loc[(dvcsSample.config>1), "PCvt12theta"].hist(bins = np.linspace(40, 90, 101), histtype = 'step', density = True, ax = axs, color = 'r', label = 'simulation')
		plt.axvline(44.5, color = 'k', linestyle = '--', linewidth = 2)
		plt.legend(loc='upper right', bbox_to_anchor = (1.05, 0.95), framealpha = 1)
		plt.xlabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
		axs.set_xticks([40, 44.5, 50, 60, 70, 80, 90])
		axs.set_xticklabels(["", 44.5, 50, 60, 70, 80, 90])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd_cvt.pdf")
		plt.clf()

		def linearfit(args, x):
			x = np.array(x)
			a, b = args
			return a + b * x
		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config>1, "Ptheta"], expSample.loc[expSample.config>1, "PCvt12theta"], bins = [np.linspace(30,70, 101), np.linspace(40, 80, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
		ticks = [1, 10, 100, 1000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		plt.ylabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
		plt.xlabel(r"$\theta_{~p'}$" + " ["+degree+"]")
		axs.set_xticks([30, 35,40, 45,50,55, 60, 65, 70])
		axs.set_xticklabels([30, 35,40, 45,50,55, 60, 65, 70])
		axs.set_yticks([40, 50, 60, 70, 80])
		axs.set_yticklabels([40, 50, 60, 70, 80])
		x1 = np.linspace((2.924+44.5)/1.274, 64.23, 101)
		axs.plot(x1, linearfit([-2.924, 1.274], x1), color = 'k', linewidth = 4, linestyle = '--')
		x1 = np.linspace((3.523+44.5)/(1.046), 64.23, 101)
		axs.plot(x1, linearfit([-3.523, 1.046], x1), color = 'k', linewidth = 4, linestyle = '--')
		x1 = np.linspace((2.924+44.5)/1.274, (3.523+44.5)/(1.046), 101)
		axs.plot(x1, 44.5 + 0*x1, color = 'k', linewidth = 4, linestyle = '--')
		y1 = np.linspace(-2.924+1.274*64.23, -3.523+1.046*64.23, 101)
		axs.plot(y1*0 + 64.23, y1, color = 'k', linewidth = 4, linestyle = '--')
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd2.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config>1, "PCvt12phi"], expSample.loc[expSample.config>1, "PCvt12theta"], bins = [np.linspace(-180,180, 361), np.linspace(40, 80, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 100))
		ticks = [1, 10, 100]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		# plt.axvline(62)
		plt.ylabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
		plt.xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
		axs.set_xticks([-180, -90, 0, 90, 180])
		axs.set_xticklabels([-180, -90, 0, 90, 180])
		axs.set_yticks([40, 50, 60, 70, 80])
		axs.set_yticklabels([40, 50, 60, 70, 80])
		plt.axvline(-95,  color = 'k', linestyle = '--', linewidth = 4)
		plt.axvline(-80,  color = 'k', linestyle = '--', linewidth = 4)
		plt.axvline(25,  color = 'k', linestyle = '--', linewidth = 4)
		plt.axvline(40,  color = 'k', linestyle = '--', linewidth = 4)
		plt.axvline(143,  color = 'k', linestyle = '--', linewidth = 4)
		plt.axvline(158,  color = 'k', linestyle = '--', linewidth = 4)
		# plt.axhline(44.5,  color = 'k', linestyle = '--', linewidth = 4)
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd3.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		dvcsPthetahist, bins = np.histogram(dvcsSample.loc[(dvcsSample.config>1), "Ptheta"], bins = np.linspace(60, 70, 101))
		expPthetahist, bins = np.histogram(expSample.loc[(expSample.config>1), "Ptheta"], bins = np.linspace(60, 70, 101))
		bincenters = (bins[:-1] + bins[1:])/2
		axs.errorbar(bincenters, expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist)), xerr = 0.03, yerr =expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))*np.sqrt(1/dvcsPthetahist + 1/expPthetahist), linestyle = '', color = 'k')
		axs.axvline(64.23, color = 'k', linestyle = '--')
		fom = expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))
		for j in range(60, 100):
			moving_average = np.sum(fom[0:j])/j
			current = fom[j]
			if current > moving_average*1.05:
				print(j, bincenters[j], moving_average, fom[j])
				break
		axs.set_xlabel(r"$\theta_{p'}$"+ " ["+degree+"]")
		axs.set_ylabel(r"$n_{exp.}/n_{sim.}$")
		axs.set_xlim([60, 70])
		axs.set_xticks([60, 62, 64.23, 66, 68, 70])
		axs.set_xticklabels([60, 62, 64.23, 66, 68, 70])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd4.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		bins = np.linspace(40, 80, 81)
		dvcsPthetahist, bins = np.histogram(dvcsSample.loc[(dvcsSample.config>1), "PCvt12theta"], bins = bins)
		expPthetahist, bins = np.histogram(expSample.loc[(expSample.config>1), "PCvt12theta"], bins = bins)
		bincenters = (bins[:-1] + bins[1:])/2
		axs.errorbar(bincenters, expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist)), xerr = 0.25, yerr = expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))*np.sqrt(1/dvcsPthetahist + 1/expPthetahist), linestyle = '', color = 'k')
		for j in range(79, 0, -1):
			moving_average = np.sum((expPthetahist/ dvcsPthetahist)[j:])/(len(bins)-1-j)
			current = (expPthetahist/ dvcsPthetahist)[j]
			if current > moving_average*1.05:
				print(j, bincenters[j], moving_average,(expPthetahist/ dvcsPthetahist)[j])
				break
		axs.axvline(44.5, color = 'k', linestyle = '--')
		axs.set_xlim([40, 80])
		axs.set_xticks([40, 44.5, 50, 60, 70, 80])
		axs.set_xticklabels([40, 44.5, 50, 60, 70, 80])
		axs.set_xlabel(r"$\theta_{\mathrm{CVT}}$"+ " ["+degree+"]")
		axs.set_ylabel(r"$n_{exp.}/n_{sim.}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd5.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		expSample.loc[(expSample.config>1), "PCvt12phi"].hist(bins = np.linspace(-180, 180, 181), histtype = 'step', density = True, ax = axs, color = 'k', label = 'experiment')
		dvcsSample.loc[(dvcsSample.config>1), "PCvt12phi"].hist(bins = np.linspace(-180, 180, 181), histtype = 'step', density = True, ax = axs, color = 'r', label = 'simulation')
		plt.xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
		axs.set_xticks([-180, -90, 0, 90, 180])
		axs.set_xticklabels([-180, -90, 0, 90, 180])
		plt.legend(loc='upper right', bbox_to_anchor = (1.2, 0.95), framealpha = 1)
		axs.set_xlim([-180, 360])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd6.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		dvcsPthetahist, bins = np.histogram(dvcsSample.loc[(dvcsSample.config>1), "PCvt12phi"], bins = np.linspace(-180, 180, 181))
		expPthetahist, bins = np.histogram(expSample.loc[(expSample.config>1), "PCvt12phi"], bins = np.linspace(-180, 180, 181))
		bincenters = (bins[:-1] + bins[1:])/2
		expPthetahist = np.where(expPthetahist>0, expPthetahist, np.inf)
		dvcsPthetahist = np.where(dvcsPthetahist>0, dvcsPthetahist, np.nan)
		ratio = np.where(expPthetahist/ dvcsPthetahist>0, expPthetahist/ dvcsPthetahist, 0)
		ratiounc = np.where(ratio>0, ratio*np.sqrt(1/dvcsPthetahist + 1/expPthetahist), 0)
		ratio = ratio/np.sum(expPthetahist)*np.sum(dvcsPthetahist)
		ratiounc = ratiounc/np.sum(expPthetahist)*np.sum(dvcsPthetahist)
		axs.errorbar(bincenters, ratio, xerr = 1, yerr = ratiounc, linestyle = '', color = 'k')
		axs.set_ylim([0, 3])
		axs.set_xlim([-180, 180])
		axs.set_xticks([-180, -90, 0, 90, 180])
		axs.set_xticklabels([-180, -90, 0, 90, 180])
		axs.axvline(-95,  color ='k', linestyle = '--')
		axs.axvline(-80,  color ='k', linestyle = '--')
		axs.axvline(25,  color ='k', linestyle = '--')
		axs.axvline(40,  color ='k', linestyle = '--')
		axs.axvline(143,  color ='k', linestyle = '--')
		axs.axvline(158,  color ='k', linestyle = '--')
		axs.set_xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
		axs.set_ylabel(r"$n_{exp.}/n_{sim.}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd7.pdf")
		plt.clf()

		#after the fiducial cuts
		dvcsSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/inb/dvcs/4893.pkl")
		expSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/inb/exp/dvcs.pkl")

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EDc3Hitx, expSample.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
		ticks = [1, 10, 100, 1000, 5000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title(r"$e'$"+" DC Outmost Layer Hits, Post-fiducial, Inbending")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_efidDC.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.PDc3Hitx, expSample.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 200))
		ticks = [1, 10, 100, 200]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_title(r"$p'$"+" DC Outmost Layer Hits, Post-fiducial")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_pfidDC.pdf")
		plt.clf()

	if polarity == "outbending":
		dvcsSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/dvcs/4907.pkl")
		expSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/exp/dvcs.pkl")


		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		axs.hist(expSampleOutb.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'k', density = True, label = 'experiment')
		axs.hist(dvcsSampleOutb.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'r', density = True, label = 'simulation')
		axs.axvline(-18, color = 'k', linestyle = '--', linewidth = 5)
		axs.axvline(10, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlim([-20, 15])
		axs.set_xticks([-20, -18, -15, -10, -5, 0, 5, 10, 15])
		axs.set_xticklabels(["", -18,  -15, -10, -5, 0, 5, 10, 15])
		axs.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
		axs.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2])
		# plt.hist(dvcsSample.Enphe, bins = np.linspace(0, 50, 51), density = True, histtype = 'step')
		plt.legend(loc='upper right', bbox_to_anchor = (1.05, 0.95), title = 'Outbending', framealpha = 1)
		# axs.set_xlabel("Number of Photoelectrons in HTCC (" + r"$n_{phe.})$")
		axs.set_xlabel("$vz_{e'}$"+ " [cm]")
		plt.tight_layout()	
		plt.savefig("plots/ch2/precutVz_outb.pdf")
		plt.clf()


		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.EDc3Hitx, expSampleOutb.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
		ticks = [1, 10, 100, 1000, 5000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title(r"$e'$"+" DC Outmost Layer Hits, Pre-fiducial, Outbending")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidDCOutb.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.PDc3Hitx, expSampleOutb.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticks = [1, 10, 100, 300]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_title(r"$p'$"+" DC Outmost Layer Hits, Pre-fiducial")
		# axs.set_xticks([-400, -200, 0, 200, 400])
		# axs.set_yticks([-400, -200, 0, 200, 400])
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidDCOutb.pdf")
		plt.clf()

		#after the fiducial cuts
		dvcsSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/outb/dvcs/4907.pkl")
		expSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/outb/exp/dvcs.pkl")

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.EDc3Hitx, expSampleOutb.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
		ticks = [1, 10, 100, 1000, 5000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title(r"$e'$"+" DC Outmost Layer Hits, Post-fiducial, Outbending")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_efidDCOutb.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.PDc3Hitx, expSampleOutb.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 100))
		ticks = [1, 10, 100]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_title(r"$p'$"+" DC Outmost Layer Hits, Post-fiducial")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_pfidDCOutb.pdf")
		plt.clf()

	#photon
	if args.figureofmerit == "photon":
		expSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/exp/dvcs.pkl")
		expSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/exp/dvcs.pkl")

		columns_needed = ["GcX", "GcY", "config", "Gsector"]
		expSample = expSample.loc[:, columns_needed]
		expSampleOutb = expSampleOutb.loc[:, columns_needed]
		df_gammaRec = pd.concat([expSample, expSampleOutb])

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(df_gammaRec.loc[df_gammaRec.config==3, "GcX"], df_gammaRec.loc[df_gammaRec.config==3, "GcY"], bins = [np.linspace(-20, 20, 101), np.linspace(-20, 20, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
		ticks = [1, 10, 100, 1000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_ylabel(r"$y_{\mathrm{FT}}$" + " ["+degree+"]")
		axs.set_xlabel(r"$x_{\mathrm{FT}}$" + " ["+degree+"]")
		theta = np.linspace(0, 2*np.pi, 101)
		circleCenterX1 = -8.419
		circleCenterY1 = 9.889
		circleRadius1 = 1.6

		circleCenterX2 = -9.89
		circleCenterY2 = -5.327
		circleRadius2 = 1.6

		circleCenterX3 = -6.15
		circleCenterY3 = -13
		circleRadius3 = 2.3
		        
		circleCenterX4 = 3.7
		circleCenterY4 = -6.5
		circleRadius4 = 2

		plt.plot(circleRadius1*np.cos(theta) + circleCenterX1, circleRadius1*np.sin(theta) + circleCenterY1, color = 'k', linewidth = 1, linestyle = '--')
		plt.plot(circleRadius2*np.cos(theta) + circleCenterX2, circleRadius2*np.sin(theta) + circleCenterY2, color = 'k', linewidth = 1, linestyle = '--')
		plt.plot(circleRadius3*np.cos(theta) + circleCenterX3, circleRadius3*np.sin(theta) + circleCenterY3, color = 'k', linewidth = 1, linestyle = '--')
		plt.plot(circleRadius4*np.cos(theta) + circleCenterX4, circleRadius4*np.sin(theta) + circleCenterY4, color = 'k', linewidth = 1, linestyle = '--')
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_gfidFT.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(df_gammaRec.loc[df_gammaRec.config<3, "GcX"], df_gammaRec.loc[df_gammaRec.config<3, "GcY"], bins = np.linspace(-450, 450, 100), cmin = 1, cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 3000))
		ticks = [1, 10, 100, 1000, 3000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xlabel("$x_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_ylabel("$y_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_xlim([-450, 450])
		axs.set_ylim([-450, 450])
		axs.set_xticks([-450, -300, -150, 0, 150, 300, 450])
		axs.set_yticks([-450, -300, -150, 0, 150, 300, 450])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_gfidPCAL.pdf", bbox_inches = 'tight')
		plt.clf()

		#after the fiducial cuts
		expSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/inb/exp/dvcs.pkl")
		expSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/outb/exp/dvcs.pkl")

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(df_gammaRec.loc[df_gammaRec.config==3, "GcX"], df_gammaRec.loc[df_gammaRec.config==3, "GcY"], bins = [np.linspace(-20, 20, 101), np.linspace(-20, 20, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticks = [1, 10, 100, 300]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_title(r"$\gamma$"+" FT-Cal Hits, Post-fiducial")
		axs.set_ylabel(r"$y_{\mathrm{FT}}$" + " ["+degree+"]")
		axs.set_xlabel(r"$x_{\mathrm{FT}}$" + " ["+degree+"]")
		theta = np.linspace(0, 2*np.pi, 101)
		circleCenterX1 = -8.419
		circleCenterY1 = 9.889
		circleRadius1 = 1.6

		circleCenterX2 = -9.89
		circleCenterY2 = -5.327
		circleRadius2 = 1.6

		circleCenterX3 = -6.15
		circleCenterY3 = -13
		circleRadius3 = 2.3
		        
		circleCenterX4 = 3.7
		circleCenterY4 = -6.5
		circleRadius4 = 2

		plt.plot(circleRadius1*np.cos(theta) + circleCenterX1, circleRadius1*np.sin(theta) + circleCenterY1, color = 'k', linewidth = 1, linestyle = '--')
		plt.plot(circleRadius2*np.cos(theta) + circleCenterX2, circleRadius2*np.sin(theta) + circleCenterY2, color = 'k', linewidth = 1, linestyle = '--')
		plt.plot(circleRadius3*np.cos(theta) + circleCenterX3, circleRadius3*np.sin(theta) + circleCenterY3, color = 'k', linewidth = 1, linestyle = '--')
		plt.plot(circleRadius4*np.cos(theta) + circleCenterX4, circleRadius4*np.sin(theta) + circleCenterY4, color = 'k', linewidth = 1, linestyle = '--')

		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_gfidFT.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		ang = -np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
		GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
		GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)
		h = axs.hist2d(GcX_rot, GcY_rot, bins = np.linspace(-450, 450, 100), cmin = 1, cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 3000))
		ticks = [1, 10, 100, 1000, 3000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticks)
		axs.set_xlabel("$x_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_ylabel("$y_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_xlim([-450, 450])
		axs.set_ylim([-450, 450])
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_gfidPCAL.pdf")
		plt.clf()


if chapter == 3:
	'''
	The chapter 3 is about the methodologies
	data set 2: with fiducial cut. no momenutm correction cut. allowing the same sectors
	data set 3: with fiducial cut. no momenutm correction cut. not allowing the same sectors.
	'''
	parent_dir = "/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr_samesectors"
	parent_exp = parent_dir + "/inb/exp/"
	epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
	parent_exp = parent_dir + "/outb/exp/"
	epgExpOutb = pd.read_pickle(parent_exp + "dvcs.pkl")

	epgExpInbCDFT = epgExpInb.loc[epgExpInb.config == 3]
	epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
	epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]

	epgExpOutbCDFT = epgExpOutb.loc[epgExpOutb.config == 3]
	epgExpOutbCD = epgExpOutb.loc[epgExpOutb.config == 2]
	epgExpOutbFD = epgExpOutb.loc[epgExpOutb.config == 1]

	epgExp = pd.concat([epgExpInb, epgExpOutb])

	fig, axs = plt.subplots(1, 1, figsize = (8, 5))
	h = axs.hist2d(epgExp.Gtheta, epgExp.Ptheta, bins = [np.linspace(0, 35, 101), np.linspace(0, 70, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
	ticks = [1, 10, 100, 1000, 5000]
	cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	axs.set_xlim([0, 35])
	axs.set_xticks([0, 5, 10, 15, 20, 25, 30, 35])
	axs.set_xticklabels([0, 5, 10, 15, 20, 25, 30, 35])
	axs.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
	axs.set_yticklabels([0, 10, 20, 30, 40, 50, 60, 70])
	axs.set_xlabel(r"$\theta_{\gamma}$"+ " ["+degree+"]")
	axs.set_ylabel(r"$\theta_{p'}$"+ " ["+degree+"]")
	plt.tight_layout()
	plt.savefig("plots/ch3/proton_correlation.pdf")
	plt.clf()


	fig, axs = plt.subplots(2, 3, figsize = (16, 10))
	h = axs[0, 0].hist2d(epgExpInbFD.xB, epgExpInbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	ticks = [1, 10, 100, 1000, 10000]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 1].hist2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 2].hist2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 0].hist2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 1].hist2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	titles = ["Inb, (FD, FD)", "Inb, (CD, FD)", "Inb, (CD, FT)", "Outb, (FD, FD)", "Outb, (CD, FD)", "Outb, (CD, FT)"]
	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$x_B$")
			axs[xind, yind].set_ylabel(r"$Q^2$" + " [" +GeVc+"]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/Q2xB_binning.pdf")
	plt.clf()

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))

	H1, xedges, yedges = np.histogram2d(epgExpInbFD.xB, epgExpInbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbFD.xB, epgExpInbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbFD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 0].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	ticks = [0.001, 0.01, 0.1, 1, 5]
	cbar = plt.colorbar(plot2d, ax = axs[0,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{min}$"+ " [ " + GeV2+"]\n", fontsize = 20)
	H1, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{min}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCDFT.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{min}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbFD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 0].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{min}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{min}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCDFT.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{min}$"+ " [ " + GeV2+"]\n", fontsize = 20)


	titles = ["Inb, (FD, FD)", "Inb, (CD, FD)", "Inb, (CD, FT)", "Outb, (FD, FD)", "Outb, (CD, FD)", "Outb, (CD, FT)"]
	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$x_B$")
			axs[xind, yind].set_ylabel(r"$Q^2$" + " [" +GeVc+"]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/Q2xB_binning_tmin.pdf")
	plt.clf()

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))

	H1, xedges, yedges = np.histogram2d(epgExpInbFD.xB, epgExpInbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbFD.xB, epgExpInbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbFD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 0].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	ticks = [0.1, 1, 10]
	cbar = plt.colorbar(plot2d, ax = axs[0,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{col}$"+ " [ " + GeV2+"]\n", fontsize = 20)
	H1, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{col}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCDFT.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{col}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbFD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 0].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{col}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{col}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCDFT.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	cbar.ax.set_title(r"$-t_{col}$"+ " [ " + GeV2+"]\n", fontsize = 20)

	titles = ["Inb, (FD, FD)", "Inb, (CD, FD)", "Inb, (CD, FT)", "Outb, (FD, FD)", "Outb, (CD, FD)", "Outb, (CD, FT)"]
	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$x_B$")
			axs[xind, yind].set_ylabel(r"$Q^2$" + " [" +GeVc+"]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/Q2xB_binning_tcol.pdf")
	plt.clf()

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))

	h = axs[0, 0].hist2d(epgExpInbFD.phi1, epgExpInbFD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	ticks = [1, 10, 100, 1000]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 1].hist2d(epgExpInbCD.phi1, epgExpInbCD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 2].hist2d(epgExpInbCDFT.phi1, epgExpInbCDFT.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 0].hist2d(epgExpOutbFD.phi1, epgExpOutbFD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 1].hist2d(epgExpOutbCD.phi1, epgExpOutbCD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.phi1, epgExpOutbCDFT.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)

	titles = ["Inb, (FD, FD)", "Inb, (CD, FD)", "Inb, (CD, FT)", "Outb, (FD, FD)", "Outb, (CD, FD)", "Outb, (CD, FT)"]
	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$\phi$" + " [" +degree+"]")
			axs[xind, yind].set_ylabel(r"$|t|$"+ " [" + GeV2 + "]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/tphi_binning.pdf")
	plt.clf()

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))
	h = axs[0, 0].hist2d(epgExpInbFD.Etheta, epgExpInbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	ticks = [1, 10, 100, 1000, 10000]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 1].hist2d(epgExpInbCD.Etheta, epgExpInbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 2].hist2d(epgExpInbCDFT.Etheta, epgExpInbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 0].hist2d(epgExpOutbFD.Etheta, epgExpOutbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 1].hist2d(epgExpOutbCD.Etheta, epgExpOutbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.Etheta, epgExpOutbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)

	titles = ["Inb, (FD, FD)", "Inb, (CD, FD)", "Inb, (CD, FT)", "Outb, (FD, FD)", "Outb, (CD, FD)", "Outb, (CD, FT)"]
	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$\theta_{e'}$" + " [" +degree+"]")
			axs[xind, yind].set_ylabel(r"$\theta_{e'\gamma}$"+ " [" + degree + "]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/coneAngle_etheta_samesector.pdf")
	plt.clf()

	fig, axs = plt.subplots(1, 1, figsize = (8, 5))
	(epgExpOutb.loc[epgExpOutb.Esector == epgExpOutb.Gsector, "phi1"]).hist( bins = np.linspace(0, 360, 361), ax = axs, color ='k')
	axs.set_xlabel(r"$\phi$" + " ["+degree+"]")
	plt.tight_layout()
	plt.savefig("plots/ch3/samesectors.pdf")
	plt.clf()


	# different sectors
	parent_dir = "/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr"
	parent_exp = parent_dir + "/inb/exp/"
	epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
	parent_exp = parent_dir + "/outb/exp/"
	epgExpOutb = pd.read_pickle(parent_exp + "dvcs.pkl")

	epgExpInbCDFT = epgExpInb.loc[epgExpInb.config == 3]
	epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
	epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]

	epgExpOutbCDFT = epgExpOutb.loc[epgExpOutb.config == 3]
	epgExpOutbCD = epgExpOutb.loc[epgExpOutb.config == 2]
	epgExpOutbFD = epgExpOutb.loc[epgExpOutb.config == 1]

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))
	h = axs[0, 0].hist2d(epgExpInbFD.Etheta, epgExpInbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	ticks = [1, 10, 100, 1000, 10000]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 1].hist2d(epgExpInbCD.Etheta, epgExpInbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[0, 2].hist2d(epgExpInbCDFT.Etheta, epgExpInbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 0].hist2d(epgExpOutbFD.Etheta, epgExpOutbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 1].hist2d(epgExpOutbCD.Etheta, epgExpOutbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.Etheta, epgExpOutbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticks)

	x = np.linspace(0, 50)

	axs[0, 0].plot(x, np.poly1d([0.0214, -0.379, 21.998])(x), color = 'k', linestyle = '--')
	axs[0, 0].plot(x, np.poly1d([0.0280, -1.001, 49.895])(x), color = 'k', linestyle = '--')
	axs[1, 0].plot(x, np.poly1d([0.0214, -0.379, 21.998])(x), color = 'k', linestyle = '--')
	axs[1, 0].plot(x, np.poly1d([0.0280, -1.001, 49.895])(x), color = 'k', linestyle = '--')

	axs[0, 1].plot(x, np.poly1d([0.0164, 0.408, 4.901])(x), color = 'k', linestyle = '--')
	axs[0, 1].plot(x, np.poly1d([0.0470, -1.677, 46.014])(x), color = 'k', linestyle = '--')
	axs[1, 1].plot(x, np.poly1d([0.0164, 0.408, 4.901])(x), color = 'k', linestyle = '--')
	axs[1, 1].plot(x, np.poly1d([0.0470, -1.677, 46.014])(x), color = 'k', linestyle = '--')

	axs[0, 2].plot(x, np.poly1d([0.0267, -0.0625, 7.730])(x), color = 'k', linestyle = '--')
	axs[0, 2].plot(x, np.poly1d([-0.00221, 0.863, 10.287])(x), color = 'k', linestyle = '--')
	axs[1, 2].plot(x, np.poly1d([0.0267, -0.0625, 7.730])(x), color = 'k', linestyle = '--')
	axs[1, 2].plot(x, np.poly1d([-0.00221, 0.863, 10.287])(x), color = 'k', linestyle = '--')

	x = np.linspace(0, 18)
	axs[0, 2].plot(x, np.poly1d([-0.000382, 0.777, 0.867])(x), color = 'b', linestyle = ':')
	axs[0, 2].plot(x, np.poly1d([0.0510, -0.0470, -0.492])(x), color = 'b', linestyle = ':')

	axs[1, 2].plot(x, np.poly1d([-0.000382, 0.777, 0.867])(x), color = 'b', linestyle = ':')
	axs[1, 2].plot(x, np.poly1d([0.0510, -0.0470, -0.492])(x), color = 'b', linestyle = ':')

	titles = ["Inb, (FD, FD)", "Inb, (CD, FD)", "Inb, (CD, FT)", "Outb, (FD, FD)", "Outb, (CD, FD)", "Outb, (CD, FT)"]
	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$\theta_{e'}$" + " [" +degree+"]")
			axs[xind, yind].set_ylabel(r"$\theta_{e'\gamma}$"+ " [" + degree + "]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/coneAngle_etheta_diffsector.pdf")
	plt.clf()

if chapter == 4:
	'''
	The chapter 4 is about the momentum correction.
	data set 3: with fiducial cut. no momenutm correction cut.
	data set 4: with fiducial cut. energy loss correction only (simulation only).
	'''

	inbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/inb/dvcs/4893.pkl")
	outbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_fid_noCorr/outb/dvcs/4907.pkl")
	columns_needed = ["Ep", "Etheta", "Ephi", "GenEp", "GenEtheta", "GenEphi", "Pp", "Ptheta", "Pphi", "GenPp", "GenPtheta", "GenPphi", "PDc1theta"]
	inbending = inbending.loc[:, columns_needed]
	outbending = outbending.loc[:, columns_needed]
	print("Drawing electron plots...")

	fig, axs = plt.subplots(2,5, figsize=(20,10))
	for row in range(2):
		for col in range(5):
			ind =col+5*row
			thetaCond = (inbending.Etheta >= 2*ind+7) & (inbending.Etheta < 2*ind+9)
			if ind == 9:
				thetaCond = inbending.Etheta>=25
			h = axs[row, col].hist2d(inbending.loc[thetaCond, "Ep"], inbending.loc[thetaCond, "GenEp"] - inbending.loc[thetaCond, "Ep"], bins = [np.linspace(2, 9, 101), np.linspace(-0.04, 0.04, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm())
			ticks = [1, 100]
			cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
			cbar.ax.set_yticklabels(ticks)
			axs[row, col].set_xlim([0, 9])
			axs[row, col].set_xticks([0, 3, 6, 9])
			axs[row, col].set_yticks([-0.04, -0.02, 0, 0.02, 0.04])
			axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
			axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
			axs[row, col].set_title(str(2*ind+7)+" "+degree + r" $\le\theta<$ " + str(2*ind+9)+" "+degree)
			if ind == 9:
				axs[row, col].set_title(str(2*ind+7)+" "+degree + r"$\le \theta<$")
	plt.tight_layout()
	# plt.show()
	plt.savefig("plots/ch4/electron_inb_mom.pdf")

	fig, axs = plt.subplots(2,5, figsize=(20,10))
	for row in range(2):
		for col in range(5):
			ind =col+5*row
			thetaCond = (outbending.Etheta >= 2*ind+5) & (outbending.Etheta < 2*ind+9)
			if ind == 9:
				thetaCond = outbending.Etheta>=23
			h = axs[row, col].hist2d(outbending.loc[thetaCond, "Ep"], outbending.loc[thetaCond, "GenEp"] - outbending.loc[thetaCond, "Ep"], bins = [np.linspace(2, 9, 101), np.linspace(-0.04, 0.04, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm())
			ticks = [1, 100]
			cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
			cbar.ax.set_yticklabels(ticks)
			axs[row, col].set_xlim([0, 9])
			axs[row, col].set_xticks([0, 3, 6, 9])
			axs[row, col].set_yticks([-0.04, -0.02, 0, 0.02, 0.04])
			axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
			axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
			axs[row, col].set_title(str(2*ind+5)+" "+degree + r" $\le\theta<$ " + str(2*ind+7)+" "+degree)
			if ind == 9:
				axs[row, col].set_title(str(2*ind+5)+" "+degree + r"$\le \theta<$")
	plt.tight_layout()
	# plt.show()
	plt.savefig("plots/ch4/electron_outb_mom.pdf")

	fig, axs = plt.subplots(2,5, figsize=(20,10))
	for row in range(2):
		for col in range(5):
			ind =col+5*row
			thetaCond = (inbending.Etheta >= 2*ind+7) & (inbending.Etheta < 2*ind+9)
			if ind == 9:
				thetaCond = inbending.Etheta>=25
			h = axs[row, col].hist2d(inbending.loc[thetaCond, "Ep"], inbending.loc[thetaCond, "GenEtheta"] - inbending.loc[thetaCond, "Etheta"], bins = [np.linspace(2, 9, 101), np.linspace(-0.2, 0.2, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm())
			ticks = [1, 100]
			cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
			cbar.ax.set_yticklabels(ticks)
			axs[row, col].set_xlim([0, 9])
			axs[row, col].set_xticks([0, 3, 6, 9])
			axs[row, col].set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
			axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
			axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
			axs[row, col].set_title(str(2*ind+7)+" "+degree + r" $\le\theta<$ " + str(2*ind+9)+" "+degree)
			if ind == 9:
				axs[row, col].set_title(str(2*ind+7)+" "+degree + r"$\le \theta<$")
	plt.tight_layout()
	# plt.show()
	plt.savefig("plots/ch4/electron_inb_theta.pdf")

	fig, axs = plt.subplots(2,5, figsize=(20,10))
	for row in range(2):
		for col in range(5):
			ind =col+5*row
			thetaCond = (outbending.Etheta >= 2*ind+5) & (outbending.Etheta < 2*ind+9)
			if ind == 9:
				thetaCond = outbending.Etheta>=23
			h = axs[row, col].hist2d(outbending.loc[thetaCond, "Ep"], outbending.loc[thetaCond, "GenEtheta"] - outbending.loc[thetaCond, "Etheta"], bins = [np.linspace(2, 9, 101), np.linspace(-0.2, 0.2, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm())
			ticks = [1, 100]
			axs[row, col].set_xlim([0, 9])
			axs[row, col].set_xticks([0, 3, 6, 9])
			axs[row, col].set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
			cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
			cbar.ax.set_yticklabels(ticks)
			axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
			axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
			axs[row, col].set_title(str(2*ind+5)+" "+degree + r" $\le\theta<$ " + str(2*ind+7)+" "+degree)
			if ind == 9:
				axs[row, col].set_title(str(2*ind+5)+" "+degree + r"$\le \theta<$")
	plt.tight_layout()
	# plt.show()
	plt.savefig("plots/ch4/electron_outb_theta.pdf")

	fig, axs = plt.subplots(2,5, figsize=(20,10))
	for row in range(2):
		for col in range(5):
			ind =col+5*row
			thetaCond = (inbending.Etheta >= 2*ind+7) & (inbending.Etheta < 2*ind+9)
			if ind == 9:
				thetaCond = inbending.Etheta>=25
			h = axs[row, col].hist2d(inbending.loc[thetaCond, "Ep"], inbending.loc[thetaCond, "GenEphi"] - inbending.loc[thetaCond, "Ephi"], bins = [np.linspace(2, 9, 101), np.linspace(-1, 1, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm())
			ticks = [1, 100]
			axs[row, col].set_xlim([0, 9])
			axs[row, col].set_xticks([0, 3, 6, 9])
			axs[row, col].set_yticks([-1, -0.5, 0, 0.5, 1])
			cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
			cbar.ax.set_yticklabels(ticks)
			axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
			axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
			axs[row, col].set_title(str(2*ind+7)+" "+degree + r" $\le\theta<$ " + str(2*ind+9)+" "+degree)
			if ind == 9:
				axs[row, col].set_title(str(2*ind+7)+" "+degree + r"$\le \theta<$")
	plt.tight_layout()
	# plt.show()
	plt.savefig("plots/ch4/electron_inb_phi.pdf")

	fig, axs = plt.subplots(2,5, figsize=(20,10))
	for row in range(2):
		for col in range(5):
			ind =col+5*row
			thetaCond = (outbending.Etheta >= 2*ind+5) & (outbending.Etheta < 2*ind+9)
			if ind == 9:
				thetaCond = outbending.Etheta>=23
			h = axs[row, col].hist2d(outbending.loc[thetaCond, "Ep"], outbending.loc[thetaCond, "GenEphi"] - outbending.loc[thetaCond, "Ephi"], bins = [np.linspace(2, 9, 101), np.linspace(-1, 1, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm())
			ticks = [1, 100]
			axs[row, col].set_xlim([0, 9])
			axs[row, col].set_xticks([0, 3, 6, 9])
			axs[row, col].set_yticks([-1, -0.5, 0, 0.5, 1])
			cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
			cbar.ax.set_yticklabels(ticks)
			axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
			axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
			axs[row, col].set_title(str(2*ind+5)+" "+degree + r" $\le\theta<$ " + str(2*ind+7)+" "+degree)
			if ind == 9:
				axs[row, col].set_title(str(2*ind+5)+" "+degree + r"$\le \theta<$")
	plt.tight_layout()
	# plt.show()
	plt.savefig("plots/ch4/electron_outb_phi.pdf")

	inbending_1 = inbending[inbending.GenPp - inbending.Pp - 0.4*0.022/inbending.Pp**1.5<0]
	inbending_2 = inbending[inbending.GenPp - inbending.Pp - 0.4*0.022/inbending.Pp**1.5>0]

	fig, axs = plt.subplots(2,3, figsize = (15, 10))
	for row in range(2):
		for col in range(3):
			axs[row, col].set_xlabel(r"$p$"+" ["+GeVc+"]")
			if col == 0:
				axs[row, col].set_ylabel(r"$\delta p$"+" ["+GeVc+"]")
			elif col == 1:
				axs[row, col].set_ylabel(r"$\theta_{\mathrm{DC, region 1}}$"+" ["+degree+"]")
			else:
				axs[row, col].set_ylabel(r"$\theta_{rec}$"+" ["+degree+"]")
	            
	axs[0, 0].hist2d(inbending_1.Pp, inbending_1.GenPp - inbending_1.Pp, bins = [np.linspace(0.3, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1)
	axs[1, 0].hist2d(inbending_2.Pp, inbending_2.GenPp - inbending_2.Pp, bins = [np.linspace(0.3, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1)
	axs[0, 1].hist2d(inbending_1.Pp, inbending_1.PDc1theta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	axs[1, 1].hist2d(inbending_2.Pp, inbending_2.PDc1theta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	axs[0, 2].hist2d(inbending_1.Pp, inbending_1.Ptheta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	axs[1, 2].hist2d(inbending_2.Pp, inbending_2.Ptheta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	t = np.linspace(0.03, 1.7, 101)
	axs[0, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)
	axs[1, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)
	t = np.linspace(0.3, 1.7, 101)
	axs[0, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)
	axs[1, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)

	plt.tight_layout()
	plt.savefig('separator1.pdf')

	params = [-53.14680163254601, 79.61307254040804, 0.3, 0.05739232362022314]#best_params#-52.99936209624629, 80.6709735338239, 0.3, 0.06899530845080828]#best_params#[-72.5, 100, 0.3, 0.055]#[-55.5, 80, 0.3, 0.04]
	inbending_check1 = inbending.loc[inbending.PDc1theta < corr(params, inbending.Pp), :]
	inbending_check2 = inbending.loc[inbending.PDc1theta > corr(params, inbending.Pp), :]

	fig, axs = plt.subplots(2,3, figsize = (15, 10))
	for row in range(2):
		for col in range(3):
			axs[row, col].set_xlabel(r"$p$"+" ["+GeVc+"]")
			if col == 0:
				axs[row, col].set_ylabel(r"$\delta p$"+" ["+GeVc+"]")
			elif col == 1:
				axs[row, col].set_ylabel(r"$\theta_{\mathrm{DC, region 1}}$"+" ["+degree+"]")
			else:
				axs[row, col].set_ylabel(r"$\theta_{rec}$"+" ["+degree+"]")
	            
	axs[0, 0].hist2d(inbending_check1.Pp, inbending_check1.GenPp - inbending_check1.Pp, bins = [np.linspace(0.3, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1)
	axs[1, 0].hist2d(inbending_check2.Pp, inbending_check2.GenPp - inbending_check2.Pp, bins = [np.linspace(0.3, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1)
	axs[0, 1].hist2d(inbending_check1.Pp, inbending_check1.PDc1theta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	axs[1, 1].hist2d(inbending_check2.Pp, inbending_check2.PDc1theta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	axs[0, 2].hist2d(inbending_check1.Pp, inbending_check1.Ptheta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)
	axs[1, 2].hist2d(inbending_check2.Pp, inbending_check2.Ptheta, bins = [np.linspace(0.3, 1.5, 101), np.linspace(0, 45, 101)], cmap = cmap, cmin = 1)

	t = np.linspace(0.03, 1.7, 101)
	axs[0, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)
	axs[1, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)

	t = np.linspace(0.3, 1.7, 101)
	axs[0, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)
	axs[1, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)

	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.4)
	plt.tight_layout()
	plt.savefig('separator2.pdf')


	# inbending_check1 = inbending.loc[inbending.PDc1theta < corr(params, inbending.Pp), :]
	# inbending_check2 = inbending.loc[inbending.PDc1theta > corr(params, inbending.Pp), :]

	# params_p = []
	# uncertainties_p = []
	# x0 = [-4.80389058e-05,  4.53222098e-03]
	# for i in range(12):
	# 	thetaCond = (inbending_check1.Ptheta >= 2.5*i+5) & (inbending_check1.Ptheta < 2.5*(i+1)+5)
	# 	dfi = copy(inbending_check1.loc[thetaCond, ["Pp", "GenPp"]])
	# 	dffit = copy(dfi[np.abs(dfi["GenPp"]-dfi["Pp"]-correction(x0, dfi["Pp"]))<0.01])
	# 	for i in range (0, 5):
	# 		res_lsq = least_squares(fun, x0, args=(dffit["Pp"], (dffit["GenPp"]-dffit["Pp"])))    
	# 		dffit = copy(dfi[np.abs(dfi["GenPp"]-dfi["Pp"]-correction(res_lsq.x, dfi["Pp"]))<0.01])
	# 		x0 = res_lsq.x

	# 	params_p.append(res_lsq.x)
	# 	# uncertainty
	# 	# https://github.com/scipy/scipy/blob/2526df72e5d4ca8bad6e2f4b3cbdfbc33e805865/scipy/optimize/minpack.py#L739
	# 	_, s, VT = np.linalg.svd(res_lsq.jac, full_matrices=False)
	# 	threshold = np.finfo(float).eps * max(res_lsq.jac.shape) * s[0]
	# 	s = s[s > threshold]
	# 	VT = VT[:s.size]
	# 	pcov = np.dot(VT.T / s**2, VT)
	# 	s_sq = np.sum((dfi["GenPp"]-dfi["Pp"]-correction(res_lsq.x, dfi["Pp"]))**2) / (len(dfi) - len(x0))
	# 	pcov = pcov * s_sq
	# 	uncertainties_p.append(np.sqrt(np.diag(pcov)))
	# 	params_p = np.array(params_p)
	# 	consts_p = params_p[:, 0]
	# 	coeffs_p = params_p[:, 1]

	# uncertainties_p = np.array(uncertainties_p)
	# consts_uncertainties_p = uncertainties_p[:, 0]
	# coeffs_uncertainties_p = uncertainties_p[:, 1]
	# fig, ax = plt.subplots(1,2, figsize=(15,5))
	# ax[0].errorbar(np.linspace(0, 11, 12)*2.5+ 5 + 1.25, consts_p, xerr= 1.25, yerr = consts_uncertainties_p, color='k', linestyle = '')
	# ax[1].errorbar(np.linspace(0, 11, 12)*2.5+ 5 + 1.25, coeffs_p, xerr= 1.25, yerr = coeffs_uncertainties_p, color='k', linestyle = '')
	# ax[0].plot(np.linspace(5, 35, 101), correction2(param1_p, np.linspace(5, 35, 101)), color = 'b')
	# ax[1].plot(np.linspace(5, 35, 101), correction3(param2_p, np.linspace(5, 35, 101)), color = 'b')
	# ax[0].set_xlabel(r"$\theta_{rec.}$"+" ["+degree+"]")
	# ax[0].set_ylabel(r"$A(\theta_{rec.})$"+" ["+GeVc+"]")
	# ax[0].set_xlim([5, 35])
	# ax[0].set_xticks(np.linspace(5, 35, 7))
	# ax[1].set_xlim([5, 35])
	# ax[1].set_xticks(np.linspace(5, 35, 7))
	# ax[1].set_xlabel(r"$\theta_{rec.}$"+" ["+degree+"]")
	# ax[1].set_ylabel(r"$B(\theta_{rec.})$"+" ["+GeVc2+"]")
	# plt.tight_layout()
	# # plt.show()
	# plt.savefig("plots/ch4/coeff_example.pdf")


if chapter == 5:
	'''
	data set 4: with fiducial cut. full momentum correction cuts. (nominal)
	data set 5: without fiducial cut. full momentum correction cuts. (eb)
	'''