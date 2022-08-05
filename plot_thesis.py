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

#chapter 2
if chapter == 2:
	'''
	The chapter 2 is about the fiducial cuts.
	data set 1: without fiducial cut. no momenutm correction cut.
	data set 2: with fiducial cut. no momenutm correction cut.
	'''

	#photon
	if args.figureofmerit == "photon":
		expSample = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/exp/dvcs.pkl")
		expSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/exp/dvcs.pkl")

		columns_needed = ["GcX", "GcY", "config", "Gsector"]
		# expSample = expSample.loc[:, columns_needed]
		# expSampleOutb = expSampleOutb.loc[:, columns_needed]
		df_gammaRec = pd.concat([expSample, expSampleOutb])

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(df_gammaRec.loc[df_gammaRec.config==3, "GcX"], df_gammaRec.loc[df_gammaRec.config==3, "GcY"], bins = [np.linspace(-20, 20, 101), np.linspace(-20, 20, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticks = [1, 10, 100, 1000, 10000]
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(a) " + r"$\gamma$"+" FT-Cal Hits, Pre-fiducial")
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
		h = axs.hist2d(df_gammaRec.loc[df_gammaRec.config<3, "GcX"], df_gammaRec.loc[df_gammaRec.config<3, "GcY"], bins = np.linspace(-450, 450, 100), cmin = 1, cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticks = [1, 10, 100, 1000, 10000]
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xlabel("$x_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_ylabel("$y_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_xlim([-450, 450])
		axs.set_ylim([-450, 450])
		axs.set_title("(a) " + r"$\gamma$"+" PCAL Hits, Pre-fiducial")
		axs.set_xticks([-450, -300, -150, 0, 150, 300, 450])
		axs.set_yticks([-450, -300, -150, 0, 150, 300, 450])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_gfidPCAL.pdf")
		plt.clf()

		#after the fiducial cuts
		expSample = electronFiducial(expSample, pol = "inbending", mc = False)
		expSample = protonFiducial(expSample, pol = "inbending")
		expSample = gammaFiducial(expSample)

		expSampleOutb = electronFiducial(expSampleOutb, pol = "outbending", mc = False)
		expSampleOutb = protonFiducial(expSampleOutb, pol = "outbending")
		expSampleOutb = gammaFiducial(expSampleOutb)

		df_gammaRec = pd.concat([expSample, expSampleOutb])

		#photon FD fiducial cuts by F.X. Girod
		df_gammaRec.loc[:, "GFid"] = 0

		#apply photon fiducial cuts
		sector_cond = [df_gammaRec.Gsector ==1, df_gammaRec.Gsector ==2, df_gammaRec.Gsector ==3, df_gammaRec.Gsector ==4, df_gammaRec.Gsector ==5, df_gammaRec.Gsector ==6]
		psplit = np.select(sector_cond, [87, 82, 85, 77, 78, 82])
		tleft = np.select(sector_cond, [58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822])
		tright = np.select(sector_cond, [58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914])
		sleft = np.select(sector_cond, [0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343])
		sright = np.select(sector_cond, [-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729])
		rleft = np.select(sector_cond, [64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997])
		rright = np.select(sector_cond, [65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641])
		qleft = np.select(sector_cond, [0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899])
		qright = np.select(sector_cond, [-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481])
		#first condition
		ang = np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
		GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
		GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)

		df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
		df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

		cond1_1 = df_gammaRec.GcX >= psplit
		cond1_2 = df_gammaRec.GcY < sleft * (df_gammaRec.GcX - tleft)
		cond1_3 = df_gammaRec.GcY > sright * (df_gammaRec.GcX - tright)
		cond1_4 = df_gammaRec.Gsector < 7
		cond1 = cond1_1 & cond1_2 & cond1_3 & cond1_4
		df_gammaRec.loc[cond1, "GFid"] = 1
		#second condition else if the first
		# cond2_0 = df_gammaRec.GFid == 0 # not necessary, because cond2_1 rules out the first (S. Lee)
		cond2_1 = df_gammaRec.GcX < psplit
		cond2_2 = df_gammaRec.GcY < qleft * (df_gammaRec.GcX - rleft)
		cond2_3 = df_gammaRec.GcY > qright * (df_gammaRec.GcX - rright)
		cond2_4 = df_gammaRec.Gsector < 7
		cond2 = cond2_1 & cond2_2 & cond2_3 & cond2_4
		df_gammaRec.loc[cond2, "GFid"] = 1

		df_gammaRec.loc[df_gammaRec.Gsector > 7, "GFid"] = 1

		#FT fiducial cuts
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

		circle1 = (df_gammaRec.GcX - circleCenterX1)**2 + (df_gammaRec.GcY - circleCenterY1)**2 < circleRadius1**2
		circle2 = (df_gammaRec.GcX - circleCenterX2)**2 + (df_gammaRec.GcY - circleCenterY2)**2 < circleRadius2**2
		circle3 = (df_gammaRec.GcX - circleCenterX3)**2 + (df_gammaRec.GcY - circleCenterY3)**2 < circleRadius3**2
		circle4 = (df_gammaRec.GcX - circleCenterX4)**2 + (df_gammaRec.GcY - circleCenterY4)**2 < circleRadius4**2

		df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle1, "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle2, "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle3, "GFid"] = 0
		df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle4, "GFid"] = 0

		df_gammaRec = df_gammaRec.loc[df_gammaRec.GFid==1, :]

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(df_gammaRec.loc[df_gammaRec.config==3, "GcX"], df_gammaRec.loc[df_gammaRec.config==3, "GcY"], bins = [np.linspace(-20, 20, 101), np.linspace(-20, 20, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticks = [1, 10, 100, 1000, 10000]
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(b) " + r"$\gamma$"+" FT-Cal Hits, Post-fiducial")
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
		h = axs.hist2d(GcX_rot, GcY_rot, bins = np.linspace(-450, 450, 100), cmin = 1, cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticks = [1, 10, 100, 1000, 10000]
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(b) " + r"$\gamma$"+" PCAL Hits, Post-fiducial")
		axs.set_xlabel("$x_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_ylabel("$y_{\mathrm{PCAL}}$"+" (cm)")
		axs.set_xlim([-450, 450])
		axs.set_ylim([-450, 450])
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_gfidPCAL.pdf")
		plt.clf()
		exit()
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

		axs.hist(expSample.Enphe, bins = np.linspace(0, 50, 51), histtype = 'step', color = 'k', density = True, label = 'Experimental Data')
		axs.hist(dvcsSample.Enphe, bins = np.linspace(0, 50, 51), histtype = 'step', color = 'r', density = True, label = 'Simulation')
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
		partp = np.linspace(2, 8, 101)
		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*(xind) + yind + 1
				h = axs[xind, yind].hist2d(expSample.loc[expSample.Esector==Esector, "Ep"],  expSample.loc[expSample.Esector==Esector, "ESamplFrac"], cmin =1, bins = [np.linspace(2, 8, 100), np.linspace(0, 0.35, 100)], cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
				ticklabels = [one, ten, hundred, thousand]
				ticks = [1, 10, 100, 1000]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				mean = ecal_e_sampl_mu[0][Esector-1] + ecal_e_sampl_mu[1][Esector-1]/1000*pow(partp-ecal_e_sampl_mu[2][Esector-1],2);
				sigma = ecal_e_sampl_sigm[0][Esector-1] + ecal_e_sampl_sigm[1][Esector-1]/(10*(partp-ecal_e_sampl_sigm[2][Esector-1]));
				axs[xind, yind].plot(partp, mean+3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].plot(partp, mean-3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].set_title(r"$e'$" + " Samp. Frac. Sector {}".format(Esector))
				axs[xind, yind].set_xlim([2, 8])
				axs[xind, yind].set_ylim([0, 0.37])
				#         axs[xind, yind].set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
				#         axs[xind, yind].set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_xticks([2, 4, 6, 8])
				axs[xind, yind].set_xticklabels([2, 4, 6, 8])
				axs[xind, yind].set_xlabel(r"$p_{e'}$" + " ["+GeVc+"]")
				axs[xind, yind].set_ylabel(r"$E_{dep.}/p_{e'}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutSampling.pdf")
		plt.clf()

		fig, axs = plt.subplots(2, 3, figsize = (15,10))
		partp = np.linspace(2, 8, 101)
		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*xind + yind  +1
				h = axs[xind, yind].hist2d(dvcsSample.loc[dvcsSample.Esector==Esector, "Ep"],  dvcsSample.loc[dvcsSample.Esector==Esector, "ESamplFrac"], cmin =1, bins = [np.linspace(2, 8, 100), np.linspace(0, 0.35, 100)], cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
				ticklabels = [one, ten, hundred, thousand, tenthousands]
				ticks = [1, 10, 100, 1000, 10000]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				mean = ecal_e_sampl_mu_mc[0][Esector-1] + ecal_e_sampl_mu_mc[1][Esector-1]/1000*pow(partp-ecal_e_sampl_mu_mc[2][Esector-1],2);
				sigma = ecal_e_sampl_sigm_mc[0][Esector-1] + ecal_e_sampl_sigm_mc[1][Esector-1]/(10*(partp-ecal_e_sampl_sigm_mc[2][Esector-1]));
				axs[xind, yind].plot(partp, mean+3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].plot(partp, mean-3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
				axs[xind, yind].set_title(r"$e'$" + " Samp. Frac. Sector {}".format(Esector))
				axs[xind, yind].set_xlim([2, 8])
				axs[xind, yind].set_ylim([0, 0.37])
				axs[xind, yind].set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
				axs[xind, yind].set_xticks([2, 4, 6, 8])
				axs[xind, yind].set_xticklabels([2, 4, 6, 8])
				axs[xind, yind].set_xlabel(r"$p_{e'}$" + " ["+GeVc+"]")
				axs[xind, yind].set_ylabel(r"$E_{dep.}/p_{e'}$")
		plt.tight_layout()
		plt.savefig("plots/ch2/precutSamplingMC.pdf")
		plt.clf()

		fig, axs = plt.subplots(2, 3, figsize = (15,10))
		for xind in range(0, 2):
			for yind in range(0,3):
				Esector = 3*(xind) + yind + 1
				h  = axs[xind, yind].hist2d(expSample.loc[expSample.Esector==Esector, "Eedep1"],  expSample.loc[expSample.Esector==Esector, "Eedep2"]+expSample.loc[expSample.Esector==Esector, "Eedep3"], cmin =1, bins = np.linspace(0, 0.8, 100), cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 200))
				ticklabels = [one, ten, hundred, two+times+hundred]
				ticks = [1, 10, 100, 200]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
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

		axs.hist(expSample.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'k', density = True, label = 'Experimental Data')
		axs.hist(dvcsSample.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'r', density = True, label = 'Simulation')
		axs.axvline(-13, color = 'k', linestyle = '--', linewidth = 5)
		axs.axvline(12, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlim([-20, 15])
		axs.set_xticks([-20, -15, -13, -10, -5, 0, 5, 10, 12, 15])
		axs.set_xticklabels([-20, "", -13, -10, -5, 0, 5, "",12, 15])
		axs.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
		axs.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2])
		# plt.hist(dvcsSample.Enphe, bins = np.linspace(0, 50, 51), density = True, histtype = 'step')
		plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 1)
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
				ticklabels = [one, ten, hundred, thousand]
				ticks = [1, 10, 100, 1000]
				cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
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
		h = axs.hist2d(expSample.EDc3Hitx, expSample.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		ticks = [1, 10, 100, 1000, 10000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title("(a) "+r"$e'$"+" DC Outmost Layer Hits, Pre-fiducial (Inb.)")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidDC.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EcalV1, expSample.ESamplFrac, bins = [np.linspace(0, 30, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 500))
		ticklabels = [one, ten, hundred, five+times+hundred]
		ticks = [1, 10, 100, 500]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.axvline(9, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlabel(r"$l_{V}$"+ " [cm]")
		axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
		axs.set_xticks([0, 5, 9, 10, 15, 20, 25, 30])
		axs.set_xticklabels([0, 5, "", 10, 15, 20, 25, 30])
		axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
		axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
		axs.set_title("(a) "+r"$e'$"+" PCAL "+r"$l_{V}$" +" Pre-fiducial (Inb.)")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidV.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EcalW1, expSample.ESamplFrac, bins = [np.linspace(0, 30, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 500))
		ticklabels = [one, ten, hundred, five+times+hundred]
		ticks = [1, 10, 100, 500]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.axvline(9, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlabel(r"$l_{W}$"+ " [cm]")
		axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
		axs.set_xticks([0, 5, 9, 10, 15, 20, 25, 30])
		axs.set_xticklabels([0, 5, "", 10, 15, 20, 25, 30])
		axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
		axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
		axs.set_title("(b) "+r"$e'$"+" PCAL "+r"$l_{W}$" +" Pre-fiducial (Inb.)")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidW.pdf")
		plt.clf()


		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config==1].Pp, expSample.loc[expSample.config==1].Pchi2pid, bins = [np.linspace(0, 1.6, 101), np.linspace(-6, 6, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
		ticklabels = [one, ten, hundred, thousand]
		ticks = [1, 10, 100, 1000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xticks([0, 0.4, 0.8, 1.2, 1.6])
		axs.set_xticklabels([0, 0.4, 0.8, 1.2, 1.6])
		axs.set_yticks([-6, -4, -2 ,0, 2, 4, 6])
		axs.set_yticklabels([-6, -4, -2 ,0, 2, 4, 6])
		axs.set_xlabel(r"$p_{p'}$"+ " ["+GeVc+"]")
		axs.set_ylabel(r"$\chi$")
		axs.set_title("(a) FD "+r"$p',~\chi-p_{p'}$"+" (Inb.)")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidchi_FD.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config>1].Pp, expSample.loc[expSample.config>1].Pchi2pid, bins = [np.linspace(0, 1.6, 101), np.linspace(-6, 6, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 1000))
		ticklabels = [one, ten, hundred, thousand]
		ticks = [1, 10, 100, 1000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xticks([0, 0.4, 0.8, 1.2, 1.6])
		axs.set_xticklabels([0, 0.4, 0.8, 1.2, 1.6])
		axs.set_yticks([-6, -4, -2 ,0, 2, 4, 6])
		axs.set_yticklabels([-6, -4, -2 ,0, 2, 4, 6])
		axs.set_xlabel(r"$p_{p'}$"+ " ["+GeVc+"]")
		axs.set_ylabel(r"$\chi$")
		axs.set_title("(b) CD "+r"$p',~\chi-p_{p'}$"+" (Inb.)")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidchi_CD.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.PDc3Hitx, expSample.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticklabels = [one, ten, hundred, three+times+hundred]
		ticks = [1, 10, 100, 300]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(a) " + r"$p'$"+" DC Outmost Layer Hits, Pre-fiducial (Inb.)")
		# axs.set_xticks([-400, -200, 0, 200, 400])
		# axs.set_yticks([-400, -200, 0, 200, 400])
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidDC.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		expSample.loc[(expSample.config>1), "Ptheta"].hist(bins = np.linspace(30, 80, 101), histtype = 'step', density = True, ax = axs, color = 'k', label = 'Experimental Data')
		dvcsSample.loc[(dvcsSample.config>1), "Ptheta"].hist(bins = np.linspace(30, 80, 101), histtype = 'step', density = True, ax = axs, color = 'r', label = 'Simulation')
		plt.axvline(64.23, color = 'k', linestyle = '--', linewidth = 2)
		plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 0.5)
		plt.xlabel(r"$\theta_{p'}$" + " ["+degree+"]")
		axs.set_title("(a)", loc = 'left')
		axs.set_xlim([30, 85])
		axs.set_xticks([30, 40, 50, 60, 64.23, 70, 80])
		axs.set_xticklabels([30, 40, 50, 60, "", 70, 80])
		axs.set_yticks([0, 0.03, 0.06])
		axs.set_ylim([0, 0.06])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		expSample.loc[(expSample.config>1), "PCvt12theta"].hist(bins = np.linspace(40, 90, 101), histtype = 'step', density = True, ax = axs, color = 'k', label = 'Experimental Data')
		dvcsSample.loc[(dvcsSample.config>1), "PCvt12theta"].hist(bins = np.linspace(40, 90, 101), histtype = 'step', density = True, ax = axs, color = 'r', label = 'Simulation')
		plt.axvline(44.5, color = 'k', linestyle = '--', linewidth = 2)
		plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 0.5)
		plt.xlabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
		axs.set_title("(a)", loc = 'left')
		axs.set_xlim([40, 90])
		axs.set_xticks([40, 44.5, 50, 60, 70, 80, 90])
		axs.set_xticklabels([40, "", 50, 60, 70, 80, 90])
		axs.set_yticks([0, 0.03, 0.06])
		axs.set_ylim([0, 0.06])
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd_cvt.pdf")
		plt.clf()

		def linearfit(args, x):
			x = np.array(x)
			a, b = args
			return a + b * x
		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config>1, "Ptheta"], expSample.loc[expSample.config>1, "PCvt12theta"], bins = [np.linspace(30,70, 101), np.linspace(40, 80, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 3000))
		ticks = [1, 10, 100, 1000, 3000]
		ticklabels = [one, ten, hundred, thousand, three + times + thousand]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
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
		axs.set_title("(a)", loc = 'left')
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd2.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.loc[expSample.config>1, "PCvt12phi"], expSample.loc[expSample.config>1, "PCvt12theta"], bins = [np.linspace(-180,180, 361), np.linspace(40, 80, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 100))
		ticks = [1, 10, 100]
		ticklabels = [one, ten, hundred]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
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
		axs.set_title("(b)", loc = 'left')
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
		axs.set_title("(b)", loc = 'left')
		axs.set_yticks([0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5])
		axs.set_ylim([0.75, 2.5])
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
		axs.set_yticks([0.5, 1, 1.5, 2, 2.5])
		axs.set_ylim([0.5, 2.7])
		axs.set_title("(b)", loc = 'left')
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd5.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		expSample.loc[(expSample.config>1), "PCvt12phi"].hist(bins = np.linspace(-180, 180, 181), histtype = 'step', density = True, ax = axs, color = 'k', label = 'Experimental Data')
		dvcsSample.loc[(dvcsSample.config>1), "PCvt12phi"].hist(bins = np.linspace(-180, 180, 181), histtype = 'step', density = True, ax = axs, color = 'r', label = 'Simulation')
		plt.xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
		axs.set_xticks([-180, -90, 0, 90, 180])
		axs.set_xticklabels([-180, -90, 0, 90, 180])
		plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 1)
		axs.set_xlim([-180, 180])
		axs.set_ylim([0, 0.006])
		axs.set_yticks([0, 0.002, 0.004, 0.006])
		axs.set_title("(a)", loc = 'left')
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
		axs.set_yticks([0, 0.5, 1, 1.5])
		axs.set_ylim([0, 1.5])
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
		axs.set_title("(b)", loc = 'left')
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfid_cd7.pdf")
		plt.clf()

		#after the fiducial cuts
		dvcsSample = electronFiducial(dvcsSample, pol = polarity, mc = True)
		dvcsSample = protonFiducial(dvcsSample, pol = polarity)
		dvcsSample = gammaFiducial(dvcsSample)

		expSample = electronFiducial(expSample, pol = polarity, mc = False)
		expSample = protonFiducial(expSample, pol = polarity)
		expSample = gammaFiducial(expSample)

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.EDc3Hitx, expSample.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		ticks = [1, 10, 100, 1000, 10000]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title("(b) " + r"$e'$"+" DC Outmost Layer Hits, Post-fiducial (Inb.)")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_efidDC.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSample.PDc3Hitx, expSample.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticklabels = [one, ten, hundred, three+times+hundred]
		ticks = [1, 10, 100, 300]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(b) " + r"$p'$"+" DC Outmost Layer Hits, Post-fiducial (Inb.)")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_pfidDC.pdf")
		plt.clf()

	if polarity == "outbending":
		dvcsSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/dvcs/4907.pkl")
		expSampleOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/exp/dvcs.pkl")


		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		axs.hist(expSampleOutb.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'k', density = True, label = 'Experimental Data')
		axs.hist(dvcsSampleOutb.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'r', density = True, label = 'Simulation')
		axs.axvline(-18, color = 'k', linestyle = '--', linewidth = 5)
		axs.axvline(10, color = 'k', linestyle = '--', linewidth = 5)
		axs.set_xlim([-20, 15])
		axs.set_xticks([-20, -18, -15, -10, -5, 0, 5, 10, 15])
		axs.set_xticklabels(["", -18,  -15, -10, -5, 0, 5, 10, 15])
		axs.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
		axs.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2])
		plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), title = 'Outbending', framealpha = 0.5)
		# axs.set_xlabel("Number of Photoelectrons in HTCC (" + r"$n_{phe.})$")
		axs.set_xlabel("$vz_{e'}$"+ " [cm]")
		plt.tight_layout()	
		plt.savefig("plots/ch2/precutVz_outb.pdf")
		plt.clf()


		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.EDc3Hitx, expSampleOutb.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticks = [1, 10, 100, 1000, 10000]
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title("(a) " + r"$e'$"+" DC Outmost Layer Hits, Pre-fiducial (Outb.)")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_efidDCOutb.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.PDc3Hitx, expSampleOutb.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticks = [1, 10, 100, 300]
		ticklabels = [one, ten, hundred, three + times + hundred]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(a) " + r"$p'$"+" DC Outmost Layer Hits, Pre-fiducial (Outb.)")
		# axs.set_xticks([-400, -200, 0, 200, 400])
		# axs.set_yticks([-400, -200, 0, 200, 400])
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/precut_pfidDCOutb.pdf")
		plt.clf()

		#after the fiducial cuts
		dvcsSampleOutb = electronFiducial(dvcsSampleOutb, pol = polarity, mc = True)
		dvcsSampleOutb = protonFiducial(dvcsSampleOutb, pol = polarity)
		dvcsSampleOutb = gammaFiducial(dvcsSampleOutb)

		expSampleOutb = electronFiducial(expSampleOutb, pol = polarity, mc = False)
		expSampleOutb = protonFiducial(expSampleOutb, pol = polarity)
		expSampleOutb = gammaFiducial(expSampleOutb)

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.EDc3Hitx, expSampleOutb.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 10000))
		ticks = [1, 10, 100, 1000, 10000]
		ticklabels = [one, ten, hundred, thousand, tenthousands]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_xlim([-300, 300])
		axs.set_ylim([-300, 300])
		axs.set_xticks([-300, -150, 0, 150, 300])
		axs.set_yticks([-300, -150, 0, 150, 300])
		axs.set_title("(b) " + r"$e'$"+" DC Outmost Layer Hits, Post-fiducial (Outb.)")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_efidDCOutb.pdf")
		plt.clf()

		fig, axs = plt.subplots(1, 1, figsize = (8, 5))
		h = axs.hist2d(expSampleOutb.PDc3Hitx, expSampleOutb.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm(vmin = 1, vmax = 300))
		ticks = [1, 10, 100, 300]
		ticklabels = [one, ten, hundred, three + times + hundred]
		cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		axs.set_title("(b) " + r"$p'$"+" DC Outmost Layer Hits, Post-fiducial (Outb.)")
		axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
		axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
		plt.tight_layout()
		plt.savefig("plots/ch2/postcut_pfidDCOutb.pdf")
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

	titles = ["Inb., (FD, FD)", "Inb., (CD, FD)", "Inb., (CD, FT)", "Outb., (FD, FD)", "Outb., (CD, FD)", "Outb., (CD, FT)"]

	fig, axs = plt.subplots(1, 1, figsize = (8, 5))
	h = axs.hist2d(epgExp.Gtheta, epgExp.Ptheta, bins = [np.linspace(0, 35, 101), np.linspace(0, 70, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 50000), rasterized = True)
	ticks = [1, 10, 100, 1000, 10000, 50000]
	ticklabels = [one, ten, hundred, thousand, tenthousands, five + times + tenthousands]
	cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
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
	ticklabels = [one, ten, hundred, thousand, tenthousands]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 1].hist2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 2].hist2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 0].hist2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 1].hist2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
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
	ticklabels = [thousandth, hundredth, tenth, one, five]
	cbar = plt.colorbar(plot2d, ax = axs[0,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{min}$"+ " [" + GeV2+"]", fontsize = 20)
	H1, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{min}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCDFT.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{min}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbFD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 0].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{min}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCD.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{min}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCDFT.tmin1)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.001, vmax =5), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{min}$"+ " [" + GeV2+"]", fontsize = 20)


	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$x_B$")
			axs[xind, yind].set_ylabel(r"$Q^2$" + " [" +GeVc+"]")
			axs[xind, yind].set_title(titles[xind*3+yind], loc = 'left')
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
	ticklabels = [tenth, one, ten]
	cbar = plt.colorbar(plot2d, ax = axs[0,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{col}$"+ " [" + GeV2+"]", fontsize = 20)
	H1, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCD.xB, epgExpInbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{col}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpInbCDFT.xB, epgExpInbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpInbCDFT.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[0, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[0,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{col}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbFD.xB, epgExpOutbFD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbFD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 0].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{col}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCD.xB, epgExpOutbCD.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCD.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 1].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{col}$"+ " [" + GeV2+"]", fontsize = 20)

	H1, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)])
	H1 = H1.T
	H2, xedges, yedges = np.histogram2d(epgExpOutbCDFT.xB, epgExpOutbCDFT.Q2, bins = [np.linspace(0, 1, 101), np.linspace(0, 10, 101)], weights = epgExpOutbCDFT.tcol)
	H2 = H2.T
	H = np.divide(H2, H1, out=np.zeros_like(H2), where=H1!=0)
	plot2d = axs[1, 2].imshow(H, interpolation='none', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto',cmap=cmap, norm = LogNorm(vmin = 0.1, vmax =10), rasterized = True)
	cbar = plt.colorbar(plot2d, ax = axs[1,2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	cbar.ax.set_title(r"$-t_{col}$"+ " [" + GeV2+"]", fontsize = 20)

	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$x_B$")
			axs[xind, yind].set_ylabel(r"$Q^2$" + " [" +GeVc+"]")
			axs[xind, yind].set_title(titles[xind*3+yind], loc = 'left')
	plt.tight_layout()
	plt.savefig("plots/ch3/Q2xB_binning_tcol.pdf")
	plt.clf()

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))

	h = axs[0, 0].hist2d(epgExpInbFD.phi1, epgExpInbFD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	ticks = [1, 10, 100, 1000]
	ticklabels = [one, ten, hundred, thousand]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 1].hist2d(epgExpInbCD.phi1, epgExpInbCD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 2].hist2d(epgExpInbCDFT.phi1, epgExpInbCDFT.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 0].hist2d(epgExpOutbFD.phi1, epgExpOutbFD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 1].hist2d(epgExpOutbCD.phi1, epgExpOutbCD.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.phi1, epgExpOutbCDFT.t1, bins = [np.linspace(0, 360, 181), np.linspace(0, 1.8, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 1000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)

	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$\phi$" + " [" +degree+"]")
			axs[xind, yind].set_ylabel(r"$|t|$"+ " [" + GeV2 + "]")
			axs[xind, yind].set_title(titles[xind*3+yind])
			axs[xind, yind].set_ylim([0, 2])
			axs[xind, yind].set_xlim([0, 360])
			axs[xind, yind].set_xticks([0, 90, 180, 270, 360])
	plt.tight_layout()
	plt.savefig("plots/ch3/tphi_binning.pdf")
	plt.clf()

	fig, axs = plt.subplots(2, 3, figsize = (16, 10))
	h = axs[0, 0].hist2d(epgExpInbFD.Etheta, epgExpInbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	ticks = [1, 10, 100, 1000, 10000]
	ticklabels = [one, ten, hundred, thousand, tenthousands]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 1].hist2d(epgExpInbCD.Etheta, epgExpInbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 2].hist2d(epgExpInbCDFT.Etheta, epgExpInbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 0].hist2d(epgExpOutbFD.Etheta, epgExpOutbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 1].hist2d(epgExpOutbCD.Etheta, epgExpOutbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.Etheta, epgExpOutbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)

	for xind in range(2):
		for yind in range(3):
			axs[xind, yind].set_xlabel(r"$\theta_{e'}$" + " [" +degree+"]")
			axs[xind, yind].set_ylabel(r"$\theta_{e'\gamma}$"+ " [" + degree + "]")
			axs[xind, yind].set_title(titles[xind*3+yind])
	plt.tight_layout()
	plt.savefig("plots/ch3/coneAngle_etheta_samesector.pdf")
	plt.clf()

	fig, axs = plt.subplots(1, 1, figsize = (8, 5))
	(epgExpOutb.loc[epgExpOutb.Esector == epgExpOutb.Gsector, "phi1"]).hist( bins = np.linspace(0, 360, 361), ax = axs, color ='k', histtype = 'stepfilled')
	axs.set_xlabel(r"$\phi$" + " ["+degree+"]")
	axs.set_xlim([0, 360])
	axs.set_xticks([0, 90, 180, 270, 360])
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
	ticklabels = [one, ten, hundred, thousand, tenthousands]
	cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 1].hist2d(epgExpInbCD.Etheta, epgExpInbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[0, 2].hist2d(epgExpInbCDFT.Etheta, epgExpInbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 0].hist2d(epgExpOutbFD.Etheta, epgExpOutbFD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 1].hist2d(epgExpOutbCD.Etheta, epgExpOutbCD.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)
	h = axs[1, 2].hist2d(epgExpOutbCDFT.Etheta, epgExpOutbCDFT.coneAngle, bins = [np.linspace(0, 40, 101), np.linspace(0, 50, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 10000), rasterized = True)
	cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
	cbar.ax.set_yticklabels(ticklabels)

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
	axs[0, 2].plot(x, np.poly1d([-0.000382, 0.777, 0.867])(x), color = 'k', linestyle = ':', linewidth = 3)
	axs[0, 2].plot(x, np.poly1d([0.0510, -0.0470, -0.492])(x), color = 'k', linestyle = ':', linewidth = 3)

	axs[1, 2].plot(x, np.poly1d([-0.000382, 0.777, 0.867])(x), color = 'k', linestyle = ':', linewidth = 3)
	axs[1, 2].plot(x, np.poly1d([0.0510, -0.0470, -0.492])(x), color = 'k', linestyle = ':', linewidth = 3)

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
	data set 4: with fiducial cut. energy loss correction only (Simulation only).
	'''
	def corr(x, t):
		x0, x1, x2, x3 = x
		return x0 + x1*np.power(t-np.ones(len(t))*0.3, x3)

	def fun(x, t, y):
		return x[0] + x[1]/t - y
	def correction(x, t):
		x0, x1 = x
		return x0 + x1/t

	def fun2(x, t, y):
		return x[0] + x[1]*t- y
	def correction2(x, t):
		x0, x1 = x
		return x0 + x1*t

	def fun3(x, t, y):
		return x[0] + x[1]*t + x[2]*t*t - y
	def correction3(x, t):
		x0, x1, x2 = x
		return x0 + x1*t + x2*t*t

	def fun4(x, t, y):
		return x[0] + x[1]*np.exp(x[2]*t) - y
	def correction4(x, t):
		x0, x1, x2 = x
		return x0 + x1*np.exp(x2*t)

	def fun5(x, t, y):
		return correction5(x, t) - y
	def correction5(x,t):
		x0, x1, x2, x3 = x
		return x0 + x1*t + x2*t*t +x3*t*t*t

	def fun6(x, t, y):
		return correction6(x,t) - y
	def correction6(x, t):
		x0, x1 = x
		return x0 + x1/t**2

	def fun7(x, t, y):
		return correction7(x,t) - y
	def correction7(x, t):
		x0, x1 = x
		return x0 + x1/t

	if args.figureofmerit == "electron":
		print("Drawing electron plots...")

		inbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/dvcs/4893.pkl")
		outbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/dvcs/4907.pkl")
		columns_needed = ["Ep", "Etheta", "Ephi", "GenEp", "GenEtheta", "GenEphi", "Pp", "Ptheta", "Pphi", "GenPp", "GenPtheta", "GenPphi", "PDc1theta", "Psector"]
		inbending = inbending.loc[:, columns_needed]
		outbending = outbending.loc[:, columns_needed]

		fig, axs = plt.subplots(2,5, figsize=(20,10))
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (inbending.Etheta >= 2*ind+7) & (inbending.Etheta < 2*ind+9)
				if ind == 9:
					thetaCond = inbending.Etheta>=25
				h = axs[row, col].hist2d(inbending.loc[thetaCond, "Ep"], inbending.loc[thetaCond, "GenEp"] - inbending.loc[thetaCond, "Ep"], bins = [np.linspace(2, 9, 101), np.linspace(-0.04, 0.04, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				axs[row, col].set_xlim([0, 9])
				axs[row, col].set_xticks([0, 3, 6, 9])
				axs[row, col].set_yticks([-0.04, -0.02, 0, 0.02, 0.04])
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
				axs[row, col].set_title(str(2*ind+7)+degree + r" $\le\theta<$ " + str(2*ind+9)+degree)
				if ind == 9:
					axs[row, col].set_title(r"$\theta\ge$" + str(2*ind+7)+degree)
		plt.tight_layout()
		plt.savefig("plots/ch4/electron_inb_mom.pdf")

		fig, axs = plt.subplots(2,5, figsize=(20,10))
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (outbending.Etheta >= 2*ind+5) & (outbending.Etheta < 2*ind+9)
				if ind == 9:
					thetaCond = outbending.Etheta>=23
				h = axs[row, col].hist2d(outbending.loc[thetaCond, "Ep"], outbending.loc[thetaCond, "GenEp"] - outbending.loc[thetaCond, "Ep"], bins = [np.linspace(2, 9, 101), np.linspace(-0.04, 0.04, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				axs[row, col].set_xlim([0, 9])
				axs[row, col].set_xticks([0, 3, 6, 9])
				axs[row, col].set_yticks([-0.04, -0.02, 0, 0.02, 0.04])
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
				axs[row, col].set_title(str(2*ind+5)+degree + r" $\le\theta<$ " + str(2*ind+7)+degree)
				if ind == 9:
					axs[row, col].set_title(r"$\theta\ge$" + str(2*ind+5)+degree)
		plt.tight_layout()
		plt.savefig("plots/ch4/electron_outb_mom.pdf")

		fig, axs = plt.subplots(2,5, figsize=(20,10))
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (inbending.Etheta >= 2*ind+7) & (inbending.Etheta < 2*ind+9)
				if ind == 9:
					thetaCond = inbending.Etheta>=25
				h = axs[row, col].hist2d(inbending.loc[thetaCond, "Ep"], inbending.loc[thetaCond, "GenEtheta"] - inbending.loc[thetaCond, "Etheta"], bins = [np.linspace(2, 9, 101), np.linspace(-0.2, 0.2, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				axs[row, col].set_xlim([0, 9])
				axs[row, col].set_xticks([0, 3, 6, 9])
				axs[row, col].set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
				axs[row, col].set_title(str(2*ind+7)+degree + r" $\le\theta<$ " + str(2*ind+9)+degree)
				if ind == 9:
					axs[row, col].set_title(r"$\theta\ge$" + str(2*ind+7)+degree)
		plt.tight_layout()
		plt.savefig("plots/ch4/electron_inb_theta.pdf")

		fig, axs = plt.subplots(2,5, figsize=(20,10))
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (outbending.Etheta >= 2*ind+5) & (outbending.Etheta < 2*ind+9)
				if ind == 9:
					thetaCond = outbending.Etheta>=23
				h = axs[row, col].hist2d(outbending.loc[thetaCond, "Ep"], outbending.loc[thetaCond, "GenEtheta"] - outbending.loc[thetaCond, "Etheta"], bins = [np.linspace(2, 9, 101), np.linspace(-0.2, 0.2, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]
				axs[row, col].set_xlim([0, 9])
				axs[row, col].set_xticks([0, 3, 6, 9])
				axs[row, col].set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
				axs[row, col].set_title(str(2*ind+5)+degree + r" $\le\theta<$ " + str(2*ind+7)+degree)
				if ind == 9:
					axs[row, col].set_title(r"$\theta\ge$" + str(2*ind+5)+degree)
		plt.tight_layout()
		plt.savefig("plots/ch4/electron_outb_theta.pdf")

		fig, axs = plt.subplots(2,5, figsize=(20,10))
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (inbending.Etheta >= 2*ind+7) & (inbending.Etheta < 2*ind+9)
				if ind == 9:
					thetaCond = inbending.Etheta>=25
				h = axs[row, col].hist2d(inbending.loc[thetaCond, "Ep"], inbending.loc[thetaCond, "GenEphi"] - inbending.loc[thetaCond, "Ephi"], bins = [np.linspace(2, 9, 101), np.linspace(-1, 1, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]
				axs[row, col].set_xlim([0, 9])
				axs[row, col].set_xticks([0, 3, 6, 9])
				axs[row, col].set_yticks([-1, -0.5, 0, 0.5, 1])
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
				axs[row, col].set_title(str(2*ind+7)+degree + r" $\le\theta<$ " + str(2*ind+9)+degree)
				if ind == 9:
					axs[row, col].set_title(r"$\theta\ge$" + str(2*ind+7)+degree)
		plt.tight_layout()
		plt.savefig("plots/ch4/electron_inb_phi.pdf")

		fig, axs = plt.subplots(2,5, figsize=(20,10))
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (outbending.Etheta >= 2*ind+5) & (outbending.Etheta < 2*ind+9)
				if ind == 9:
					thetaCond = outbending.Etheta>=23
				h = axs[row, col].hist2d(outbending.loc[thetaCond, "Ep"], outbending.loc[thetaCond, "GenEphi"] - outbending.loc[thetaCond, "Ephi"], bins = [np.linspace(2, 9, 101), np.linspace(-1, 1, 101)], cmap = cmap, cmin =1, rasterized = True, norm = LogNorm(vmin = 1, vmax = 5000))
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]
				axs[row, col].set_xlim([0, 9])
				axs[row, col].set_xticks([0, 3, 6, 9])
				axs[row, col].set_yticks([-1, -0.5, 0, 0.5, 1])
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
				axs[row, col].set_title(str(2*ind+5)+degree + r" $\le\theta<$ " + str(2*ind+7)+degree)
				if ind == 9:
					axs[row, col].set_title(r"$\theta\ge$" + str(2*ind+5)+degree)
		plt.tight_layout()
		plt.savefig("plots/ch4/electron_outb_phi.pdf")
		exit()

	elif args.figureofmerit == "benchmark":
		columns_needed = ["Pp", "Ptheta", "Pphi", "GenPp", "GenPtheta", "GenPphi", "Psector", "MM2_ep"]

		print("reading inbending no corr..")
		inbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/dvcs/4893.pkl")
		print("reading outbending no corr..")
		outbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/dvcs/4907.pkl")

		inb_FD = inbending.loc[inbending.Psector<7, columns_needed]
		outb_FD = outbending.loc[outbending.Psector<7, columns_needed]
		inb_CD = inbending.loc[inbending.Psector>7, columns_needed]
		outb_CD = outbending.loc[outbending.Psector>7, columns_needed]

		print("reading inbending with corr..")
		inbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_Corr/inb/dvcs/4893.pkl")
		print("reading outbending with corr..")
		outbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_Corr/outb/dvcs/4907.pkl")

		inb_FD_corr = inbending.loc[inbending.Psector<7, columns_needed]
		outb_FD_corr = outbending.loc[outbending.Psector<7, columns_needed]
		inb_CD_corr = inbending.loc[inbending.Psector>7, columns_needed]
		outb_CD_corr = outbending.loc[outbending.Psector>7, columns_needed]

		fig, axs = plt.subplots(1,4, figsize = (20, 8))
		for row in range(4):
			axs[row].set_xlabel(r"$MM^2_{ep}$"+" ["+GeV2+"]")
			axs[row].set_xlim([-0.4, 0.4])
			axs[row].set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
		bins = np.linspace(-0.4, 0.4, 101)

		inb_FD.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'r', ax = axs[0])
		inb_FD_corr.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'b', ax = axs[0])
		axs[0].set_title("Inb. FD  " + r"$MM^2_{ep}$")
		axs[0].set_ylim([0, 60000])
		axs[0].set_yticks([0, 20000, 40000, 60000])
		axs[0].set_yticklabels([zero,two+times+tenthousands, four+times+tenthousands, six+times+tenthousands])

		outb_FD.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'r', ax = axs[1])
		outb_FD_corr.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'b', ax = axs[1])
		axs[1].set_title("Outb. FD  " + r"$MM^2_{ep}$")
		axs[1].set_ylim([0, 20000])
		axs[1].set_yticks([0, 5000, 10000, 15000, 20000])
		axs[1].set_yticklabels([zero,five+times+thousand, tenthousands, one+point+five+times+tenthousands, two +times+tenthousands])

		inb_CD.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'r', ax = axs[2])
		inb_CD_corr.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'b', ax = axs[2])
		axs[2].set_title("Inb. CD  " + r"$MM^2_{ep}$")
		axs[2].set_ylim([0, 400000])
		axs[2].set_yticks([0, 100000, 200000, 300000, 400000])
		axs[2].set_yticklabels([zero,hundredthousands, two+times+hundredthousands, three+times+hundredthousands, four+times+hundredthousands])

		outb_CD.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'r', ax = axs[3])
		outb_CD_corr.MM2_ep.hist(bins = bins, histtype = 'step', edgecolor = 'b', ax = axs[3])
		axs[3].set_title("Outb. CD  " + r"$MM^2_{ep}$")
		axs[3].set_ylim([0, 300000])
		axs[3].set_yticks([0, 100000, 200000, 300000])
		axs[3].set_yticklabels([zero,hundredthousands, two+times+hundredthousands, three+times+hundredthousands])
		plt.suptitle("(a)", x = 0.05, size = 40)
		plt.tight_layout()
		plt.savefig("plots/ch4/benchmarks_excl.pdf")	

		fig, axs = plt.subplots(1,4, figsize = (20, 8))
		for row in range(4):
			axs[row].set_xlabel(r"$\delta p$"+" ["+GeVc+"]")
			axs[row].set_xlim([-0.05, 0.05])
		bins_p = np.linspace(-0.05, 0.05, 101)

		(inb_FD.GenPp - inb_FD.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'r', ax = axs[0])
		(inb_FD_corr.GenPp - inb_FD_corr.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'b', ax = axs[0])
		axs[0].set_title("Inb. FD  " + r"$\delta p$")
		axs[0].set_ylim([0, 100000])
		axs[0].set_yticks([0, 20000, 40000, 60000, 80000, 100000])
		axs[0].set_yticklabels([zero,two+times+tenthousands, four+times+tenthousands, six+times+tenthousands, eight+times+tenthousands, hundredthousands])

		(outb_FD.GenPp - outb_FD.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'r', ax = axs[1])
		(outb_FD_corr.GenPp - outb_FD_corr.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'b', ax = axs[1])
		axs[1].set_title("Outb. FD  " + r"$\delta p$")
		axs[1].set_ylim([0, 25000])
		axs[1].set_yticks([0, 5000, 10000, 15000, 20000, 25000])
		axs[1].set_yticklabels([zero,five+times+tenthousands, tenthousands, one+point+five+times+tenthousands, two+times+tenthousands, two+point+five+times+tenthousands])

		(inb_CD.GenPp - inb_CD.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'r', ax = axs[2])
		(inb_CD_corr.GenPp - inb_CD_corr.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'b', ax = axs[2])
		axs[2].set_title("Inb. CD  " + r"$\delta p$")
		axs[2].set_ylim([0, 200000])
		axs[2].set_yticks([0, 50000, 100000, 150000, 200000])
		axs[2].set_yticklabels([zero,five+times+tenthousands, tenthousands, one+point+five+times+tenthousands, two+times+tenthousands])

		(outb_CD.GenPp - outb_CD.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'r', ax = axs[3])
		(outb_CD_corr.GenPp - outb_CD_corr.Pp).hist(bins = bins_p, histtype = 'step', edgecolor = 'b', ax = axs[3])
		axs[3].set_title("Outb. CD  " + r"$\delta p$")
		axs[3].set_ylim([0, 150000])
		axs[3].set_yticks([0, 50000, 100000, 150000])
		axs[3].set_yticklabels([zero,five+times+tenthousands, tenthousands, one+point+five+times+tenthousands])
		plt.suptitle("(b)", x = 0.05, size = 40)
		plt.tight_layout()
		plt.savefig("plots/ch4/benchmarks_p.pdf")

		fig, axs = plt.subplots(1,4, figsize = (20, 8))
		for row in range(4):
			axs[row].set_xlabel(r"$\delta \theta$"+" ["+degree+"]")
			axs[row].set_xlim([-2, 2])
		bins_theta = np.linspace(-2, 2, 101)

		(inb_FD.GenPtheta - inb_FD.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'r', ax = axs[0])
		(inb_FD_corr.GenPtheta - inb_FD_corr.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'b', ax = axs[0])
		axs[0].set_title("Inb. FD  " + r"$\delta \theta$")
		axs[0].set_ylim([0, 50000])
		axs[0].set_yticks([0, 10000, 20000, 30000, 40000, 50000])
		axs[0].set_yticklabels([zero, tenthousands, two+times+tenthousands, three+times+tenthousands, four+times+tenthousands, five+times+tenthousands])

		(outb_FD.GenPtheta - outb_FD.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'r', ax = axs[1])
		(outb_FD_corr.GenPtheta - outb_FD_corr.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'b', ax = axs[1])
		axs[1].set_title("Outb. FD  " + r"$\delta \theta$")
		axs[1].set_ylim([0, 20000])
		axs[1].set_yticks([0, 5000, 10000, 15000, 20000])
		axs[1].set_yticklabels([zero,five+times+thousand, tenthousands, one+point+five+times+tenthousands, two +times+tenthousands])

		(inb_CD.GenPtheta - inb_CD.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'r', ax = axs[2])
		(inb_CD_corr.GenPtheta - inb_CD_corr.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'b', ax = axs[2])
		axs[2].set_title("Inb. CD  " + r"$\delta \theta$")
		axs[2].set_ylim([0, 300000])
		axs[2].set_yticks([0, 100000, 200000, 300000])
		axs[2].set_yticklabels([zero,hundredthousands, two+times+hundredthousands, three+times+hundredthousands])

		(outb_CD.GenPtheta - outb_CD.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'r', ax = axs[3])
		(outb_CD_corr.GenPtheta - outb_CD_corr.Ptheta).hist(bins = bins_theta, histtype = 'step', edgecolor = 'b', ax = axs[3])
		axs[3].set_title("Outb. CD  " + r"$\delta \theta$")
		axs[3].set_ylim([0, 250000])
		axs[3].set_yticks([0, 50000, 100000, 150000, 200000, 250000])
		axs[3].set_yticklabels([zero,five+times+tenthousands, tenthousands, one+point+five+times+tenthousands, two+times+tenthousands, two+point+five+times+tenthousands])
		plt.suptitle("(c)", x = 0.05, size = 40)
		plt.tight_layout()
		plt.savefig("plots/ch4/benchmarks_theta.pdf")

		fig, axs = plt.subplots(1,4, figsize = (20, 8))
		for row in range(4):
			axs[row].set_xlabel(r"$\delta \phi$"+" ["+degree+"]")
			axs[row].set_xlim([-2, 2])
		bins_phi = np.linspace(-2, 2, 101)

		(inb_FD.GenPphi - inb_FD.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'r', ax = axs[0])
		(inb_FD_corr.GenPphi - inb_FD_corr.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'b', ax = axs[0])
		axs[0].set_title("Inb. FD  " + r"$\delta \phi$")
		axs[0].set_ylim([0, 25000])
		axs[0].set_yticks([0, 5000, 10000, 15000, 20000, 25000])
		axs[0].set_yticklabels([zero,five+times+thousand, thousand, one+point+five+times+thousand, two+times+thousand, two+point+five+times+thousand])

		(outb_FD.GenPphi - outb_FD.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'r', ax = axs[1])
		(outb_FD_corr.GenPphi - outb_FD_corr.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'b', ax = axs[1])
		axs[1].set_title("Outb. FD  " + r"$\delta \phi$")
		axs[1].set_ylim([0, 12000])
		axs[1].set_yticks([0, 2000, 4000, 6000, 8000, 10000, 12000])
		axs[1].set_yticklabels([zero,two+times+thousand, four+times+thousand, six+times+thousand, eight+times+thousand, tenthousands, one+point+two+times+tenthousands])

		(inb_CD.GenPphi - inb_CD.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'r', ax = axs[2])
		(inb_CD_corr.GenPphi - inb_CD_corr.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'b', ax = axs[2])
		axs[2].set_title("Inb. CD  " + r"$\delta \phi$")
		axs[2].set_ylim([0, 600000])
		axs[2].set_yticks([0, 200000, 400000, 600000])
		axs[2].set_yticklabels([zero,two+times+hundredthousands, four+times+hundredthousands, six+times+hundredthousands])

		(outb_CD.GenPphi - outb_CD.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'r', ax = axs[3])
		(outb_CD_corr.GenPphi - outb_CD_corr.Pphi).hist(bins = bins_phi, histtype = 'step', edgecolor = 'b', ax = axs[3])
		axs[3].set_title("Outb. CD  " + r"$\delta \phi$")
		axs[3].set_ylim([0, 300000])
		axs[3].set_yticks([0, 100000, 200000, 300000, 400000])
		axs[3].set_yticklabels([zero,hundredthousands, two+times+hundredthousands, three+times+hundredthousands, four+times+hundredthousands])
		plt.suptitle("(d)", x = 0.05, size = 40)
		plt.tight_layout()
		plt.savefig("plots/ch4/benchmarks_phi.pdf")

	elif args.figureofmerit == "fitting":
		inbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/dvcs/4893.pkl")
		outbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/dvcs/4907.pkl")
		columns_needed = ["Ep", "Etheta", "Ephi", "GenEp", "GenEtheta", "GenEphi", "Pp", "Ptheta", "Pphi", "GenPp", "GenPtheta", "GenPphi", "PDc1theta", "Psector"]
		inbending = inbending.loc[:, columns_needed]
		outbending = outbending.loc[:, columns_needed]

		inb_FD = inbending.loc[inbending.Psector<7, :]
		inb_CD = inbending.loc[inbending.Psector>7, :]
		outb_FD = outbending.loc[outbending.Psector<7, :]
		outb_CD = outbending.loc[outbending.Psector>7, :]
		fig, axs = plt.subplots(1, 4, figsize = (20, 5))
		for row in range(4):
			axs[row].set_xlim([0, 1.5])
			axs[row].set_xticks([0, 0.5, 1, 1.5])
			axs[row].set_yticks([-0.05, 0, 0.05, 0.1])
			axs[row].set_xlabel(r"$p$"+" ["+GeVc+"]")
			axs[row].set_ylabel(r"$\delta p$"+" ["+GeVc+"]")
		ticks = [1, 10, 100, 1000, 2000]
		ticklabels = [one, ten, hundred, thousand, two+times+thousand]
		h = axs[0].hist2d(inb_FD.Pp, inb_FD.GenPp - inb_FD.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin =1, vmax = 2000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1].hist2d(inb_CD.Pp, inb_CD.GenPp - inb_CD.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin =1, vmax = 2000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[2].hist2d(outb_FD.Pp, outb_FD.GenPp - outb_FD.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin =1, vmax = 2000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[2], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[3].hist2d(outb_CD.Pp, outb_CD.GenPp - outb_CD.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin =1, vmax = 2000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[3], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)

		axs[0].set_title("(a) Inb., FD " + r"$p'$")
		axs[1].set_title("(b) Inb., CD " + r"$p'$")
		axs[2].set_title("(c) Outb., FD " + r"$p'$")
		axs[3].set_title("(d) Outb., CD " + r"$p'$")

		plt.tight_layout()
		plt.savefig("plots/ch4/protons_all_conf.pdf")
		plt.clf()

		inbending_1 = inb_FD.loc[inb_FD.GenPp - inb_FD.Pp - 0.4*0.022/inb_FD.Pp**1.5<0, :]
		inbending_2 = inb_FD.loc[inb_FD.GenPp - inb_FD.Pp - 0.4*0.022/inb_FD.Pp**1.5>0, :]

		params = [-53.14680163254601, 79.61307254040804, 0.3, 0.05739232362022314]#best_params#-52.99936209624629, 80.6709735338239, 0.3, 0.06899530845080828]#best_params#[-72.5, 100, 0.3, 0.055]#[-55.5, 80, 0.3, 0.04]

		fig, axs = plt.subplots(2,3, figsize = (15, 10))
		for row in range(2):
			for col in range(3):
				axs[row, col].set_xticks([0, 0.5, 1, 1.5])
				axs[row, col].set_xlabel(r"$p$"+" ["+GeVc+"]")
				if col == 0:
					axs[row, col].set_ylabel(r"$\delta p$"+" ["+GeVc+"]")
				elif col == 1:
					axs[row, col].set_ylabel(r"$\theta_{\mathrm{DC, region 1}}$"+" ["+degree+"]")
				else:
					axs[row, col].set_ylabel(r"$\theta_{rec}$"+" ["+degree+"]")
				
		ticks = [1, 10, 100, 1000, 3000]
		ticklabels = [one, ten, hundred, thousand, three+times+thousand]
		h = axs[0, 0].hist2d(inbending_1.Pp, inbending_1.GenPp - inbending_1.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1, 0].hist2d(inbending_2.Pp, inbending_2.GenPp - inbending_2.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[0, 1].hist2d(inbending_1.Pp, inbending_1.PDc1theta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1, 1].hist2d(inbending_2.Pp, inbending_2.PDc1theta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[0, 2].hist2d(inbending_1.Pp, inbending_1.Ptheta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1, 2].hist2d(inbending_2.Pp, inbending_2.Ptheta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		t = np.linspace(0.03, 1.7, 101)
		axs[0, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)
		axs[1, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)
		t = np.linspace(0.3, 1.7, 101)
		axs[0, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)
		axs[1, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)
		axs[0, 0].set_title("(a)", loc = 'left')
		axs[0, 1].set_title("(b)", loc = 'left')
		axs[0, 2].set_title("(c)", loc = 'left')
		#plt.suptitle("FD "+rx = 0.05ided in "+r"$\delta p-p$"+".", size = 40)
		plt.tight_layout()
		plt.savefig('plots/ch4/separator1.pdf')
		plt.clf()

		inbending_check1 = inb_FD.loc[inb_FD.PDc1theta < corr(params, inb_FD.Pp), :]
		inbending_check2 = inb_FD.loc[inb_FD.PDc1theta > corr(params, inb_FD.Pp), :]

		fig, axs = plt.subplots(2,3, figsize = (15, 10))
		for row in range(2):
			for col in range(3):
				axs[row, col].set_xticks([0, 0.5, 1, 1.5])
				axs[row, col].set_xlabel(r"$p$"+" ["+GeVc+"]")
				if col == 0:
					axs[row, col].set_ylabel(r"$\delta p$"+" ["+GeVc+"]")
				elif col == 1:
					axs[row, col].set_ylabel(r"$\theta_{\mathrm{DC, region 1}}$"+" ["+degree+"]")
				else:
					axs[row, col].set_ylabel(r"$\theta_{rec}$"+" ["+degree+"]")
				
		ticks = [1, 10, 100, 1000, 3000]
		ticklabels = [one, ten, hundred, thousand, three+times+thousand]
		h = axs[0, 0].hist2d(inbending_check1.Pp, inbending_check1.GenPp - inbending_check1.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0, 0], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1, 0].hist2d(inbending_check2.Pp, inbending_check2.GenPp - inbending_check2.Pp, bins = [np.linspace(0, 1.5, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1, 0], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[0, 1].hist2d(inbending_check1.Pp, inbending_check1.PDc1theta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0, 1], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1, 1].hist2d(inbending_check2.Pp, inbending_check2.PDc1theta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1, 1], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[0, 2].hist2d(inbending_check1.Pp, inbending_check1.Ptheta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[0, 2], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)
		h = axs[1, 2].hist2d(inbending_check2.Pp, inbending_check2.Ptheta, bins = [np.linspace(0, 1.5, 101), np.linspace(0, 40, 101)], cmap = cmap, cmin = 1, norm = LogNorm(vmin = 1, vmax = 3000), rasterized = True)
		cbar = plt.colorbar(h[3], ax = axs[1, 2], ticks = ticks)
		cbar.ax.set_yticklabels(ticklabels)

		t = np.linspace(0.03, 1.7, 101)
		axs[0, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)
		axs[1, 0].plot(t, 0.4*0.022/t**1.5, color = 'k', linewidth = 3)

		t = np.linspace(0.3, 1.7, 101)
		axs[0, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)
		axs[1, 1].plot(t, corr(params, t), color = 'k', linewidth = 3)
		axs[0, 0].set_title("(a)", loc = 'left')
		axs[0, 1].set_title("(b)", loc = 'left')
		axs[0, 2].set_title("(c)", loc = 'left')
		#plt.suptitle("FD "+rx = 0.05ided in "+r"$\theta_{\mathrm{DC, region 1}}-p$"+".", fontsize = 40)
		plt.tight_layout()
		plt.savefig('plots/ch4/separator2.pdf')
		plt.clf()

		params_p = []
		uncertainties_p = []
		for i in range(12):
			x0 = [correction2([-0.00051894, -0.00018104],2.5*i+5+1.25), correction3([ 3.29466917e-03,  5.73663160e-04, -1.40807209e-05], 2.5*i+5+1.25)]
			thetaCond = (inbending_check1.Ptheta >= 2.5*i+5) & (inbending_check1.Ptheta < 2.5*(i+1)+5)
			dfi = copy(inbending_check1.loc[thetaCond, ["Pp", "GenPp"]])
			dffit = copy(dfi[np.abs(dfi["GenPp"]-dfi["Pp"]-correction(x0, dfi["Pp"]))<0.01])
			for i in range (0, 5):
				res_lsq = least_squares(fun, x0, args=(dffit["Pp"], (dffit["GenPp"]-dffit["Pp"])))    
				dffit = copy(dfi[np.abs(dfi["GenPp"]-dfi["Pp"]-correction(res_lsq.x, dfi["Pp"]))<0.01])
				x0 = res_lsq.x

			params_p.append(res_lsq.x)
			# uncertainty
			# https://github.com/scipy/scipy/blob/2526df72e5d4ca8bad6e2f4b3cbdfbc33e805865/scipy/optimize/minpack.py#L739
			_, s, VT = np.linalg.svd(res_lsq.jac, full_matrices=False)
			threshold = np.finfo(float).eps * max(res_lsq.jac.shape) * s[0]
			s = s[s > threshold]
			VT = VT[:s.size]
			pcov = np.dot(VT.T / s**2, VT)
			s_sq = np.sum((dfi["GenPp"]-dfi["Pp"]-correction(res_lsq.x, dfi["Pp"]))**2) / (len(dfi) - len(x0))
			pcov = pcov * s_sq
			uncertainties_p.append(np.sqrt(np.diag(pcov)))
		params_p = np.array(params_p)
		consts_p = params_p[:, 0]
		coeffs_p = params_p[:, 1]
		uncertainties_p = np.array(uncertainties_p)
		consts_uncertainties_p = uncertainties_p[:, 0]
		coeffs_uncertainties_p = uncertainties_p[:, 1]
		x = np.linspace(0, 11, 12)*2.5+ 5 + 1.25
		x = x[:-4]
		res_lsq = least_squares(fun2, [0.005, -0.002], args=(x, consts_p[:-4]))    
		param1_p = res_lsq.x
		res_lsq = least_squares(fun3, [0, 0, 0], args=(x, coeffs_p[:-4]))    
		param2_p = res_lsq.x
		fig, ax = plt.subplots(1,2, figsize=(15,5))
		ax[0].errorbar(np.linspace(0, 11, 12)*2.5+ 5 + 1.25, consts_p, xerr= 1.25, yerr = consts_uncertainties_p, color='k', linestyle = '')
		ax[1].errorbar(np.linspace(0, 11, 12)*2.5+ 5 + 1.25, coeffs_p, xerr= 1.25, yerr = coeffs_uncertainties_p, color='k', linestyle = '')
		ax[0].plot(np.linspace(5, 35, 101), correction2(param1_p, np.linspace(5, 35, 101)), color = 'r', linestyle = '--')
		ax[1].plot(np.linspace(5, 35, 101), correction3(param2_p, np.linspace(5, 35, 101)), color = 'r', linestyle = '--')
		# ax[0].plot(np.linspace(5, 35, 101), correction2([-0.00051894, -0.00018104], np.linspace(5, 35, 101)), color = 'b')
		# ax[1].plot(np.linspace(5, 35, 101), correction3([ 3.29466917e-03,  5.73663160e-04, -1.40807209e-05], np.linspace(5, 35, 101)), color = 'b')
		ax[0].set_xlabel(r"$\theta_{rec.}$"+" ["+degree+"]")
		ax[0].set_ylabel(r"$A(\theta_{rec.})$"+" ["+GeVc+"]")
		ax[0].set_xlim([5, 35])
		ax[0].set_ylim([-0.05, 0.01])
		ax[0].set_yticks([-0.05, -0.04, -0.03,-0.02, -0.01, 0, 0.01])
		ax[0].set_xticks(np.linspace(5, 35, 7))
		ax[1].set_xlim([5, 35])
		ax[1].set_ylim([0, 0.03])
		ax[1].set_xticks(np.linspace(5, 35, 7))
		ax[1].set_xlabel(r"$\theta_{rec.}$"+" ["+degree+"]")
		ax[1].set_ylabel(r"$B(\theta_{rec.})$"+" ["+GeVc2+"]")
		plt.suptitle("(b)", x = 0.05, size = 40)
		plt.tight_layout()
		plt.savefig("plots/ch4/coeff_example.pdf")
		plt.clf()

		x = np.linspace(0, 11, 12)*2.5+ 5 + 1.25
		fig, axs = plt.subplots(2, 5, figsize=(20,10))
		ticks = [1, 10, 100, 200]
		ticklabels = [one, ten, hundred, two+times+hundred]
		for row in range(2):
			for col in range(5):
				ind =col+5*row
				thetaCond = (inbending_check1.Ptheta >= 2.5*ind+5) & (inbending_check1.Ptheta < 2.5*(ind+1)+5)
				h = axs[row, col].hist2d(inbending_check1.loc[thetaCond, "Pp"], inbending_check1.loc[thetaCond, "GenPp"] - inbending_check1.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.05, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax= 200), rasterized = True)
				cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
				cbar.ax.set_yticklabels(ticklabels)
				param = params_p[ind]
				axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
				axs[row, col].set_xlim([0, 2])
				axs[row, col].set_ylim([-0.05, 0.05])
				axs[row, col].set_yticks([-0.05, 0, 0.05])
				axs[row, col].set_xticks([0, 1, 2])
				axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
				axs[row, col].set_xticklabels(['', '', '', '', '', '', '', '', ''], minor = True)
				axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
				axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
				axs[row, col].set_title(str(2.5*ind+5)+degree + r" $\le\theta_{rec.}<$ " + str(2.5*(ind+1)+5)+degree)
		plt.suptitle("(a)", x = 0.05, size = 40)
		plt.tight_layout()
		plt.savefig("plots/ch4/protonFD_inb_mom_1_procedure.pdf")
		plt.clf()
		exit()
	elif args.figureofmerit == "proton":
		columns_needed = ["Pp", "Ptheta", "Pphi", "GenPp", "GenPtheta", "GenPphi", "PDc1theta", "Psector"]
		params = [-53.14680163254601, 79.61307254040804, 0.3, 0.05739232362022314]#best_params#-52.99936209624629, 80.6709735338239, 0.3, 0.06899530845080828]#best_params#[-72.5, 100, 0.3, 0.055]#[-55.5, 80, 0.3, 0.04]

		if polarity == "inbending":
			inbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/inb/dvcs/4893.pkl")
			inbending = inbending.loc[:, columns_needed]
			inb_FD = inbending.loc[inbending.Psector<7, :]
			inbendingCD = inbending.loc[inbending.Psector>7, :]
			if args.detector == "FD1":
				inbending_check1 = inb_FD.loc[inb_FD.PDc1theta < corr(params, inb_FD.Pp), :]
				param1_p, param2_p = [-0.00051894, -0.00018104], [ 3.29466917e-03,  5.73663160e-04, -1.40807209e-05]
				x = np.linspace(0, 11, 12)*2.5+ 5 + 1.25
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				ticks = [1, 10, 100, 200]
				ticklabels = [one, ten, hundred, two+times+hundred]
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbending_check1.Ptheta >= 2.5*ind+5) & (inbending_check1.Ptheta < 2.5*(ind+1)+5)
						h = axs[row, col].hist2d(inbending_check1.loc[thetaCond, "Pp"], inbending_check1.loc[thetaCond, "GenPp"] - inbending_check1.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.05, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 200), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = [correction2(param1_p, x[ind]), correction3(param2_p, x[ind])]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_ylim([-0.05, 0.05])
						axs[row, col].set_yticks([-0.05, 0, 0.05])
						axs[row, col].set_xticklabels(['', '', '', '', '', '', '', '', ''], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
						axs[row, col].set_title(str(2.5*ind+5)+degree + r" $\le\theta_{rec.}<$ " + str(2.5*(ind+1)+5)+degree)
				plt.suptitle("(c)", x = 0.05, size = 40)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_inb_mom_1.pdf")
				plt.clf()

				fig, axs = plt.subplots(2,5, figsize=(20,10))
				param1_theta, param2_theta = [-0.16742969,  0.00697925], [ 0.23352115, -0.01338697]
				ticks = [1, 10, 50]
				ticklabels = [one, ten, r'$50$']
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbending_check1.Ptheta >= 2.5*ind+5) & (inbending_check1.Ptheta < 2.5*(ind+1)+5)
						h = axs[row, col].hist2d(inbending_check1.loc[thetaCond, "Pp"], inbending_check1.loc[thetaCond, "GenPtheta"] - inbending_check1.loc[thetaCond, "Ptheta"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-1, 1, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin =1, vmax = 50), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = [correction2(param1_theta, x[ind]), correction2(param2_theta, x[ind])]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_ylim(-2, 2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(2.5*ind+5)+degree + r" $\le\theta_{rec.}<$ " + str(2.5*(ind+1)+5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_inb_theta_1.pdf")
				plt.clf()

				fig, axs = plt.subplots(2,5, figsize=(20,10))
				ticks = [1, 10, 100, 200]
				ticklabels = [one, ten, hundred, two+times+hundred]
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbending_check1.Ptheta >= 2.5*ind+5) & (inbending_check1.Ptheta < 2.5*(ind+1)+5)
						h = axs[row, col].hist2d(inbending_check1.loc[thetaCond, "Pp"], inbending_check1.loc[thetaCond, "GenPphi"] - inbending_check1.loc[thetaCond, "Pphi"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-5, 5, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 200), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						const_FD = 0.21192125 -0.0115175 * x[ind]
						coeff_FD = -8.94307411*0.1 + 1.66349766*0.1 * x[ind] -8.90617559*0.001 * x[ind] * x[ind] + 1.64803754*0.0001 * x[ind] * x[ind] * x[ind]
						param = [const_FD, coeff_FD]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_ylim(-2, 2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
						axs[row, col].set_title(str(2.5*ind+5)+degree + r" $\le\theta_{rec.}<$ " + str(2.5*(ind+1)+5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_inb_phi_1.pdf")
				plt.clf()

			if args.detector == "FD2":
				inbending_check2 = inb_FD.loc[inb_FD.PDc1theta > corr(params, inb_FD.Pp), :]
				ticks = [1, 10, 100, 300]
				ticklabels = [one, ten, hundred, three+times+hundred]

				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbending_check2.Ptheta >= 2*ind+23) & (inbending_check2.Ptheta < 2*(ind+1)+23)
						h = axs[row, col].hist2d(inbending_check2.loc[thetaCond, "Pp"], inbending_check2.loc[thetaCond, "GenPp"] - inbending_check2.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.1, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin =1, vmax = 300), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						theta = 2*ind + 23 + 1
						const_FD = -3.03346359*10**(-1) + 1.83368163*10**(-2)*theta - 2.86486404*10**(-4)*theta*theta
						coeff_FD =  2.01023276*10**(-1) - 1.13312215*10**(-2)*theta + 1.82487916*10**(-4)*theta*theta
						param = [const_FD, coeff_FD]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
						axs[row, col].set_title(str(2*ind+23)+degree + r" $\le\theta_{rec.}<$ " + str(2*(ind+1)+23)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_inb_mom_2.pdf")
				plt.clf()

				fig, axs = plt.subplots(2,5, figsize=(20,10))
				param1_theta, param2_theta = [ 2.04334532e+01, -1.81052405e+00,  5.32556360e-02, -5.23157558e-04], [ 8.74233279e+00, -7.63869344e-01,  2.22376362e-02, -2.16457260e-04]
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbending_check2.Ptheta >= 2*ind+23) & (inbending_check2.Ptheta < 2*(ind+1)+23)
						h = axs[row, col].hist2d(inbending_check2.loc[thetaCond, "Pp"], inbending_check2.loc[thetaCond, "GenPtheta"] - inbending_check2.loc[thetaCond, "Ptheta"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-2, 2, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin =1, vmax = 300), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						theta = 2*ind + 23 + 1
						param = [correction5(param1_theta, theta), correction5(param2_theta, theta)]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(2*ind+23)+degree + r" $\le\theta_{rec.}<$ " + str(2*(ind+1)+23)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_inb_theta_2.pdf")
				plt.clf()

				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbending_check2.Ptheta >= 2*ind+23) & (inbending_check2.Ptheta < 2*(ind+1)+23)
						h = axs[row, col].hist2d(inbending_check2.loc[thetaCond, "Pp"], inbending_check2.loc[thetaCond, "GenPphi"] - inbending_check2.loc[thetaCond, "Pphi"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-10, 10, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin =1, vmax = 300), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						theta = 2*ind + 23 + 1
						const_FD = 0.54697831 -0.04896981*theta +  0.00111376*theta*theta
						coeff_FD = -4.06733541*10**2 + 2.43696202*10*theta -3.36144736*10**(-1)*theta*theta
						coeff2_FD = 2.06378660*10 - 1.42866062*theta + 2.01085440*10**(-2)*theta*theta
						param = [const_FD, coeff_FD, coeff2_FD]#params_phi[ind]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction4(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
						axs[row, col].set_title(str(2*ind+23)+degree + r" $\le\theta_{rec.}<$ " + str(2*(ind+1)+23)+degree)
						plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_inb_phi_2.pdf")
				plt.clf()

			if args.detector == "CD":
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five + times + thousand]

				theta = np.linspace(0, 9, 10)*2.5+40 + 1.25
				const_CD = 1.93686914 - 0.116288824*theta + 0.00223685833*theta**2 - 1.40771969 * 10**(-5)*theta**3
				coeff_CD = -0.738047800 + 0.0443343685*theta - 8.50985972*10**(-4)*theta*theta + 5.36810280 * 10**(-6) * theta**3
				params_p = np.array([const_CD, coeff_CD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						if ind >12:
							continue
						thetaCond = (inbendingCD.Ptheta >= 2.5*(ind)+40) & (inbendingCD.Ptheta < 2.5*(ind+1)+40)
						h = axs[row, col].hist2d(inbendingCD.loc[thetaCond, "Pp"], inbendingCD.loc[thetaCond, "GenPp"] - inbendingCD.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.2, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_p[ind]
						axs[row, col].plot(np.linspace(0.25, 1.7, 101), correction(param, np.linspace(0.25, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
						axs[row, col].set_title(str(2.5*ind+40)+degree + r" $\le\theta<$ " + str(2.5*(ind+1)+40)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonCD_inb_mom.pdf")
				plt.clf()

				const_CD = -1.09849291*100 + 8.86664014 * theta - 0.26643881 * theta**2 + 3.53814210 * 10**(-3) * theta**3 - 1.75297107 * 10**(-5) * theta**4
				coeff_CD = 9.52034523*100 -5.74808292 * 10 * theta + 1.15386949 * theta**2 - 7.57970373 * 0.001 * theta**3
				coeff2_CD = -2.00387313*100 + 1.18979079 * 10 * theta - 2.37730217*0.1 * theta**2 + 1.55153003*0.001*theta**3
				params_theta = np.array([const_CD, coeff_CD, coeff2_CD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						if ind >12:
							continue
						thetaCond = (inbendingCD.Ptheta >= 2.5*(ind)+40) & (inbendingCD.Ptheta < 2.5*(ind+1)+40)
						h = axs[row, col].hist2d(inbendingCD.loc[thetaCond, "Pp"], inbendingCD.loc[thetaCond, "GenPtheta"] - inbendingCD.loc[thetaCond, "Ptheta"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-5, 5, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_theta[ind]
						axs[row, col].plot(np.linspace(0.25, 1.7, 101), correction4(param, np.linspace(0.25, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(2.5*ind+40)+degree + r" $\le\theta<$ " + str(2.5*(ind+1)+40)+degree)
						axs[row, col].set_yticks(np.linspace(-5, 5, 5))
				plt.tight_layout()
				plt.savefig("plots/ch4/protonCD_inb_theta.pdf")
				plt.clf()

				const_CD = 4.94546178 -3.26662886*0.1 * theta +  7.39069603 * 0.001 * theta**2 -6.83599356*10**(-5) * theta**3 + 2.12303103*10**(-7) * theta**4
				coeff_CD = 1.72181613*10**(5) -1.36827111*10**(4) * theta + 4.00923146*10**(2) * theta**2 - 5.12792347 * theta**3 + 2.41793167*10**(-2) * theta**4
				coeff2_CD =  1.20477219*10**(2) -5.86630228 * theta + 7.44007875*10**(-2) * theta**2 -2.42652473*10**(-4) * theta**3
				params_theta = np.array([const_CD, coeff_CD, coeff2_CD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						thetaCond = (inbendingCD.Ptheta >= 2.5*(ind)+40) & (inbendingCD.Ptheta < 2.5*(ind+1)+40)
						h = axs[row, col].hist2d(inbendingCD.loc[thetaCond, "Pp"], inbendingCD.loc[thetaCond, "GenPphi"] - inbendingCD.loc[thetaCond, "Pphi"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-5, 5, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_theta[ind]
						axs[row, col].plot(np.linspace(0.25, 1.7, 101), correction4(param, np.linspace(0.25, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						# axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
						axs[row, col].set_title(str(2.5*ind+40)+degree + r" $\le\theta<$ " + str(2.5*(ind+1)+40)+degree)
						axs[row, col].set_yticks(np.linspace(-5, 5, 5))
				plt.tight_layout()
				plt.savefig("plots/ch4/protonCD_inb_phi.pdf")
				plt.clf()

		elif polarity == "outbending":
			outbending = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_nofid_noCorr/outb/dvcs/4907.pkl")
			outbending = outbending.loc[:, columns_needed]
			outb_FD = outbending.loc[outbending.Psector<7, :]
			outbendingCD = outbending.loc[outbending.Psector>7, :]

			if args.detector == "FD1":
				outb_FD_1 = outb_FD.loc[outb_FD.PDc1theta < corr(params, outb_FD.Pp), :]
				ticks = [1, 10, 50]
				ticklabels = [one, ten, r'$50$']

				x = np.linspace(0, 9, 10)*1.5 + 16.5 + 1.5/2
				const_FD = 0.05083242 -0.00469777*x + 0.0001082*x*x
				coeff_FD = -1.47443264*0.01 + 1.58220893*0.001*x -3.19490013*0.00001*x*x
				params_p = np.array([const_FD, coeff_FD]).T
				fig, axs = plt.subplots(2, 5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+2*row
						thetaCond = (outb_FD_1.Ptheta >= 1.5*ind+16.5) & (outb_FD_1.Ptheta < 1.5*(ind+1)+16.5)
						h = axs[row, col].hist2d(outb_FD_1.loc[thetaCond, "Pp"], outb_FD_1.loc[thetaCond, "GenPp"] - outb_FD_1.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.05, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 50), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_p[ind]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
						axs[row, col].set_title(str(1.5*ind+16.5)+degree + r" $\le\theta_{rec.}<$ " + str(1.5*(ind+1)+16.5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_outb_mom_1.pdf")
				plt.clf()

				const_FD = -2.56460305*10 + 3.29877542*x -1.43106886*0.1*x*x + 2.08341898*0.001*x*x*x
				coeff_FD =  9.12532740*10 -1.20100762*10*x + 5.27654711*0.1*x*x -7.72656759*0.001*x*x*x
				params_theta = np.array([const_FD, coeff_FD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+2*row
						thetaCond = (outb_FD_1.Ptheta >= 1.5*ind+16.5) & (outb_FD_1.Ptheta < 1.5*(ind+1)+16.5)
						param = params_theta[ind]
						h = axs[row, col].hist2d(outb_FD_1.loc[thetaCond, "Pp"], outb_FD_1.loc[thetaCond, "GenPtheta"] - outb_FD_1.loc[thetaCond, "Ptheta"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-1, 1, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 50), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color = 'k', linewidth=2, linestyle = '--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+degree+"]")
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(1.5*ind+16.5)+degree + r" $\le\theta_{rec.}<$ " + str(1.5*(ind+1)+16.5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_outb_theta_1.pdf")
				plt.clf()

				const_FD = -20.4780893 + 1.67020488*x - 0.03419348*x*x
				coeff_FD = 35.02807194 - 2.9098043*x +  0.06037906*x*x
				params_phi = np.array([const_FD, coeff_FD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+2*row
						thetaCond = (outb_FD_1.Ptheta >= 1.5*ind+16.5) & (outb_FD_1.Ptheta < 1.5*(ind+1)+16.5)
						param = params_phi[ind]
						h = axs[row, col].hist2d(outb_FD_1.loc[thetaCond, "Pp"], outb_FD_1.loc[thetaCond, "GenPphi"] - outb_FD_1.loc[thetaCond, "Pphi"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-2, 2, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 50), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
						axs[row, col].set_title(str(1.5*ind+16.5)+degree + r" $\le\theta_{rec.}<$ " + str(1.5*(ind+1)+16.5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_outb_phi_1.pdf")
				plt.clf()

			if args.detector == "FD2":
				outb_FD_2 = outb_FD.loc[outb_FD.PDc1theta > corr(params, outb_FD.Pp), :]
				ticks = [1, 10, 100]
				ticklabels = [one, ten, hundred]

				x = np.linspace(0, 9, 10)*1.5 + 27.5 + 1.5/2
				const_FD = 0.09832589 -0.0066463*x + 0.00010312*x*x
				coeff_FD = -9.61421691*0.01 + 6.85638807*0.001*x -9.75766427*0.00001*x*x
				params_p = np.array([const_FD, coeff_FD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+2*row
						thetaCond = (outb_FD_2.Ptheta >= 1.5*ind+27.5) & (outb_FD_2.Ptheta < 1.5*(ind+1)+27.5)
						h = axs[row, col].hist2d(outb_FD_2.loc[thetaCond, "Pp"], outb_FD_2.loc[thetaCond, "GenPp"] - outb_FD_2.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.05, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_p[ind]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
						axs[row, col].set_title(str(1.25*ind+27.5)+degree + r" $\le\theta_{rec.}<$ " + str(1.25*(ind+1)+27.5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_outb_mom_2.pdf")
				plt.clf()

				const_FD = -1.68873940 + 9.56867163*0.01*x -1.43741464*0.001*x*x
				coeff_FD = 1.49978357*10 -1.40137094*x + 4.38501543*0.01*x*x -4.57982872*0.0001*x*x*x
				params_theta = np.array([const_FD, coeff_FD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+2*row
						thetaCond = (outb_FD_2.Ptheta >= 1.5*ind+27.5) & (outb_FD_2.Ptheta < 1.5*(ind+1)+27.5)
						h = axs[row, col].hist2d(outb_FD_2.loc[thetaCond, "Pp"], outb_FD_2.loc[thetaCond, "GenPtheta"] - outb_FD_2.loc[thetaCond, "Ptheta"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-2, 2, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_theta[ind]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(1.25*ind+27.5)+degree + r" $\le\theta_{rec.}<$ " + str(1.25*(ind+1)+27.5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_outb_theta_2.pdf")
				plt.clf()

				const_FD = 6.75359137 - 0.43199851*x + 0.0068995*x*x
				coeff_FD = -1.68588219 + 1.05609627*0.1*x -1.50452832*0.001*x*x
				params_phi = np.array([const_FD, coeff_FD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+2*row
						thetaCond = (outb_FD_2.Ptheta >= 1.5*ind+27.5) & (outb_FD_2.Ptheta < 1.5*(ind+1)+27.5)
						h = axs[row, col].hist2d(outb_FD_2.loc[thetaCond, "Pp"], outb_FD_2.loc[thetaCond, "GenPphi"] - outb_FD_2.loc[thetaCond, "Pphi"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-2, 2, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 100), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_phi[ind]
						axs[row, col].plot(np.linspace(0.3, 1.7, 101), correction6(param, np.linspace(0.3, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(1.25*ind+27.5)+degree + r" $\le\theta_{rec.}<$ " + str(1.25*(ind+1)+27.5)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonFD_outb_phi_2.pdf")
				plt.clf()

			if args.detector == "CD":
				ticks = [1, 10, 100, 1000, 5000]
				ticklabels = [one, ten, hundred, thousand, five+times+thousand]

				theta = np.linspace(0, 9, 10)*2.5+40 + 1.25
				const_CD = 1.92657376 - 0.113836734*theta + 0.00215038526*theta**2 - 1.32525053 * 10**(-5)*theta**3
				coeff_CD = -0.755650043 + 0.0445538936*theta - 8.38241864*10**(-4)*theta*theta + 5.16887255 * 10**(-6) * theta**3
				params_p = np.array([const_CD, coeff_CD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						if ind >12:
							continue
						thetaCond = (outbendingCD.Ptheta >= 2.5*(ind)+40) & (outbendingCD.Ptheta < 2.5*(ind+1)+40)
						h = axs[row, col].hist2d(outbendingCD.loc[thetaCond, "Pp"], outbendingCD.loc[thetaCond, "GenPp"] - outbendingCD.loc[thetaCond, "Pp"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-0.05, 0.2, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_p[ind]
						axs[row, col].plot(np.linspace(0.25, 1.7, 101), correction(param, np.linspace(0.25, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta p$" + " ["+GeVc+"]")
						axs[row, col].set_title(str(2.5*ind+40)+degree + r" $\le\theta<$ " + str(2.5*(ind+1)+40)+degree)
				plt.tight_layout()
				plt.savefig("plots/ch4/protonCD_outb_mom.pdf")
				plt.clf()

				const_CD = -5.79024055*10 + 4.67197531 * theta - 0.140156897 * theta**2 + 1.85853057 * 10**(-3) * theta**3 - 9.19989908 * 10**(-6) * theta**4
				coeff_CD = 2.99700765*1000 - 2.18027982 * 10**2 * theta + 5.84757503 * theta**2 - 6.80409195 * 0.01 * theta**3 + 2.89244618 * 0.0001 * theta**4
				coeff2_CD = -1.82237904*100 + 1.10153549 * 10 * theta - 2.24699931*0.1 * theta**2 + 1.49390960*0.001*theta**3
				params_theta = np.array([const_CD, coeff_CD, coeff2_CD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						if ind >12:
							continue
						thetaCond = (outbendingCD.Ptheta >= 2.5*(ind)+40) & (outbendingCD.Ptheta < 2.5*(ind+1)+40)
						h = axs[row, col].hist2d(outbendingCD.loc[thetaCond, "Pp"], outbendingCD.loc[thetaCond, "GenPtheta"] - outbendingCD.loc[thetaCond, "Ptheta"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-5, 5, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_theta[ind]
						axs[row, col].plot(np.linspace(0.25, 1.7, 101), correction4(param, np.linspace(0.25, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \theta$" + " ["+degree+"]")
						axs[row, col].set_title(str(2.5*ind+40)+degree + r" $\le\theta<$ " + str(2.5*(ind+1)+40)+degree)
						axs[row, col].set_yticks(np.linspace(-5, 5, 5))
				plt.tight_layout()
				plt.savefig("plots/ch4/protonCD_outb_theta.pdf")
				plt.clf()

				const_CD = 7.58761670 - 5.28224578*0.1 * theta +  1.31580117 * 0.01 * theta**2 -1.41738951*10**(-4) * theta**3 + 5.62884363*10**(-7) * theta**4
				coeff_CD = 1.07644097*10**(5) - 8.67994639*10**(3) * theta + 2.57187193*10**(2) * theta**2 - 3.31379317 * theta**3 + 1.56896621*10**(-2) * theta**4
				coeff2_CD =  1.92263184*10**(2) -1.00870704 * 10 * theta + 1.56575252*10**(-1) * theta**2 -7.71489734*10**(-4) * theta**3
				params_phi = np.array([const_CD, coeff_CD, coeff2_CD]).T
				fig, axs = plt.subplots(2,5, figsize=(20,10))
				for row in range(2):
					for col in range(5):
						ind =col+5*row
						if ind >12:
							continue
						thetaCond = (outbendingCD.Ptheta >= 2.5*(ind)+40) & (outbendingCD.Ptheta < 2.5*(ind+1)+40)
						h = axs[row, col].hist2d(outbendingCD.loc[thetaCond, "Pp"], outbendingCD.loc[thetaCond, "GenPphi"] - outbendingCD.loc[thetaCond, "Pphi"], bins = [np.linspace(0.3, 1.7, 101), np.linspace(-5, 5, 101)], cmap = cmap, cmin =1, norm = LogNorm(vmin = 1, vmax = 5000), rasterized = True)
						cbar = plt.colorbar(h[3], ax = axs[row, col], ticks = ticks)
						cbar.ax.set_yticklabels(ticklabels)
						param = params_phi[ind]
						axs[row, col].plot(np.linspace(0.25, 1.7, 101), correction4(param, np.linspace(0.25, 1.7, 101)), color='k', linewidth=5, linestyle='--')
						#axs[row, col].set_xlim(0.3, 1.2)
						axs[row, col].set_xlim([0, 2])
						axs[row, col].set_xticks([0, 1, 2])
						axs[row, col].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], minor = True)
						axs[row, col].set_xlabel(r"$p$" + " ["+GeVc+"]")
						axs[row, col].set_ylabel(r"$\delta \phi$" + " ["+degree+"]")
						axs[row, col].set_title(str(2.5*ind+40)+degree + r" $\le\theta<$ " + str(2.5*(ind+1)+40)+degree)
						axs[row, col].set_yticks(np.linspace(-5, 5, 5))
				plt.tight_layout()
				plt.savefig("plots/ch4/protonCD_outb_phi.pdf")
				plt.clf()

if chapter == 5:
	'''
	data set 4: without fiducial cut. full momentum correction cuts. (nominal)
	'''
	if polarity == "inbending":
		parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/"
		parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_1g/"
		parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_2g/"
		parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/exp/"

		#epg Exp
		epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
		dvcsSimInb = pd.read_pickle(parent_MC + "4893.pkl")
		bkgSimInb = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
		pi0ExpInb = pd.read_pickle(parent_exp + "pi0.pkl")
		pi0SimInb = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")

		epgExpInbCDFT = epgExpInb.loc[epgExpInb.config == 3]
		dvcsSimInbCDFT = dvcsSimInb.loc[dvcsSimInb.config == 3]
		bkgSimInbCDFT = bkgSimInb.loc[bkgSimInb.config == 3]
		pi0ExpInbCDFT = pi0ExpInb.loc[(pi0ExpInb.config == 3)]
		pi0SimInbCDFT = pi0SimInb.loc[(pi0SimInb.config == 3)]

		epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
		dvcsSimInbCD = dvcsSimInb.loc[dvcsSimInb.config == 2]
		bkgSimInbCD = bkgSimInb.loc[bkgSimInb.config == 2]
		pi0ExpInbCD = pi0ExpInb.loc[(pi0ExpInb.config == 2)]
		pi0SimInbCD = pi0SimInb.loc[(pi0SimInb.config == 2)]

		epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]
		dvcsSimInbFD = dvcsSimInb.loc[dvcsSimInb.config == 1]
		bkgSimInbFD = bkgSimInb.loc[bkgSimInb.config == 1]
		pi0ExpInbFD = pi0ExpInb.loc[(pi0ExpInb.config == 1)]
		pi0SimInbFD = pi0SimInb.loc[(pi0SimInb.config == 1)]

		contInbCD = len(pi0ExpInbCD)*len(bkgSimInbCD)/len(pi0SimInbCD)/len(epgExpInbCD)
		contInbFD = len(pi0ExpInbFD)*len(bkgSimInbFD)/len(pi0SimInbFD)/len(epgExpInbFD)
		contInbCDFT = len(pi0ExpInbCDFT)*len(bkgSimInbCDFT)/len(pi0SimInbCDFT)/len(epgExpInbCDFT)

		varstoplot = ["coneAngle", "MM2_eg", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
		title = [r"$\theta_{e'\gamma}$", r"$MM^2_{e'\gamma}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi_{\vec{L}\vec{\Gamma}}$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
		unit = [degree, GeV2, degree, degree, GeV, GeV2, GeV2, GeVc]
		df3 = epgExpInbCDFT
		df1 = dvcsSimInbCDFT
		df2 = bkgSimInbCDFT
		df4 = pi0ExpInbCDFT
		df5 = pi0SimInbCDFT
		xlb = [10, 0.3, 0, 0, -.3, -0.01, -0.3, 0]
		xub = [35, 1.5, 0.75, 6, .3, 0.01, 0.3, 0.1]
		yub = [0.15, 2, 3, 0.4, 3, 250, 4, 25]
		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
				simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
				simDist = (1-contInbCDFT)*simDist_dvcs + contInbCDFT*simDist_dvpi0
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
		axs[0, 0].set_xlim([10, 30])
		axs[0, 0].set_xticks([10, 20, 30])
		axs[0, 1].set_xlim([0.5, 1.5])
		axs[0, 1].set_xticks([0.5, 1, 1.5])
		axs[0, 2].set_xticks([0, 0.2, 0.4, 0.6])
		axs[0, 3].set_xticks([0, 2, 4, 6])
		axs[1, 0].set_xlim([-0.3, 0.3])
		axs[1, 0].set_xticks([-0.3, 0, 0.3])
		axs[1, 1].set_xlim([-0.01, 0.01])
		axs[1, 1].set_xticks([-0.01, 0, 0.01])
		axs[1, 2].set_xlim([-0.3, 0.3])
		axs[1, 2].set_xticks([-0.3, 0, 0.3])
		axs[1, 3].set_xlim([0, 0.08])
		axs[1, 3].set_xticks([0, 0.04, 0.08])
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/dvcsInbCDFTexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		df3 = epgExpInbCD
		df1 = dvcsSimInbCD
		df2 = bkgSimInbCD
		df4 = pi0ExpInbCD
		df5 = pi0SimInbCD
		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
				simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
				simDist = (1-contInbCD)*simDist_dvcs + contInbCD*simDist_dvpi0
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/dvcsInbCDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		df3 = epgExpInbFD
		df1 = dvcsSimInbFD
		df2 = bkgSimInbFD
		df4 = pi0ExpInbFD
		df5 = pi0SimInbFD
		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
				simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
				simDist = (1-contInbFD)*simDist_dvcs + contInbFD*simDist_dvpi0
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
		axs[0, 0].set_xticks([25, 30, 35, 40, 45])
		axs[0, 0].set_xlim([25, 40])
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/dvcsInbFDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		varstoplot = ["Mpi0", "MM2_egg", "reconPi", "coplanarity", "ME_epgg", "MM2_epgg", "MM2_ep", "MPt"]
		title = [r"$IM_{\pi^0}$", r"$MM^2_{e'\pi^0}$", r"$\theta_{\pi^0_{det.}\pi^0_{rec.}}$", r"$\Delta\phi$" , "ME"+r"${}_{ep\pi^0}$", "MM"+r"${}^{2}_{ep\pi^0}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{ep\pi^0}$"]
		unit = [GeV, GeV2, degree, degree, GeV, GeV2, GeV2, GeVc]

		df4 = pi0ExpInbCDFT
		df5 = pi0SimInbCDFT

		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/pi0InbCDFTexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()


		df4 = pi0ExpInbCD
		df5 = pi0SimInbCD

		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/pi0InbCDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()


		df4 = pi0ExpInbFD
		df5 = pi0SimInbFD

		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/pi0InbFDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

	else:
		parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/dvcs/"
		parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/bkg_1g/"
		parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/bkg_2g/"
		parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/outb/exp/"

		#epg Exp
		epgExpOutb = pd.read_pickle(parent_exp + "dvcs.pkl")
		dvcsSimOutb= pd.read_pickle(parent_MC + "4907.pkl")
		bkgSimOutb = pd.read_pickle(parent_MC_bkg1g + "4243.pkl")
		pi0ExpOutb = pd.read_pickle(parent_exp + "pi0.pkl")
		pi0SimOutb = pd.read_pickle(parent_MC_bkg2g + "4243.pkl")

		epgExpOutbCDFT = epgExpOutb.loc[epgExpOutb.config == 3]
		dvcsSimOutbCDFT = dvcsSimOutb.loc[dvcsSimOutb.config == 3]
		bkgSimOutbCDFT = bkgSimOutb.loc[bkgSimOutb.config == 3]
		pi0ExpOutbCDFT = pi0ExpOutb.loc[(pi0ExpOutb.config == 3)]
		pi0SimOutbCDFT = pi0SimOutb.loc[(pi0SimOutb.config == 3)]

		epgExpOutbCD = epgExpOutb.loc[epgExpOutb.config == 2]
		dvcsSimOutbCD = dvcsSimOutb.loc[dvcsSimOutb.config == 2]
		bkgSimOutbCD = bkgSimOutb.loc[bkgSimOutb.config == 2]
		pi0ExpOutbCD = pi0ExpOutb.loc[(pi0ExpOutb.config == 2)]
		pi0SimOutbCD = pi0SimOutb.loc[(pi0SimOutb.config == 2)]

		epgExpOutbFD = epgExpOutb.loc[epgExpOutb.config == 1]
		dvcsSimOutbFD = dvcsSimOutb.loc[dvcsSimOutb.config == 1]
		bkgSimOutbFD = bkgSimOutb.loc[bkgSimOutb.config == 1]
		pi0ExpOutbFD = pi0ExpOutb.loc[(pi0ExpOutb.config == 1)]
		pi0SimOutbFD = pi0SimOutb.loc[(pi0SimOutb.config == 1)]

		contOutbCD = len(pi0ExpOutbCD)*len(bkgSimOutbCD)/len(pi0SimOutbCD)/len(epgExpOutbCD)
		contOutbFD = len(pi0ExpOutbFD)*len(bkgSimOutbFD)/len(pi0SimOutbFD)/len(epgExpOutbFD)
		contOutbCDFT = len(pi0ExpOutbCDFT)*len(bkgSimOutbCDFT)/len(pi0SimOutbCDFT)/len(epgExpOutbCDFT)

		varstoplot = ["coneAngle", "MM2_eg", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
		title = [r"$\theta_{e'\gamma}$", r"$MM^2_{e'\gamma}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi_{\vec{L}\vec{\Gamma}}$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
		unit = [degree, GeV2, degree, degree, GeV, GeV2, GeV2, GeVc]

		df3 = epgExpOutbCDFT
		df1 = dvcsSimOutbCDFT
		df2 = bkgSimOutbCDFT
		df4 = pi0ExpOutbCDFT
		df5 = pi0SimOutbCDFT
		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
				simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
				simDist = (1-contOutbCDFT)*simDist_dvcs + contOutbCDFT*simDist_dvpi0
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
		axs[0, 0].set_xlim([5, 35])
		axs[0, 0].set_xticks([5, 10, 15, 20, 25, 30, 35])
		axs[0, 1].set_xlim([0.5, 1.5])
		axs[0, 1].set_xticks([0.5, 1, 1.5])
		axs[0, 2].set_xticks([0, 0.2, 0.4, 0.6])
		axs[0, 3].set_xticks([0, 2.5, 5])
		axs[1, 0].set_xlim([-0.3, 0.3])
		axs[1, 0].set_xticks([-0.3, 0, 0.3])
		axs[1, 1].set_xlim([-0.01, 0.01])
		axs[1, 1].set_xticks([-0.01, 0, 0.01])
		axs[1, 2].set_xlim([-0.3, 0.3])
		axs[1, 2].set_xticks([-0.3, 0, 0.3])
		axs[1, 3].set_xlim([0, 0.08])
		axs[1, 3].set_xticks([0, 0.04, 0.08])
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/dvcsOutbCDFTexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		df3 = epgExpOutbCD
		df1 = dvcsSimOutbCD
		df2 = bkgSimOutbCD
		df4 = pi0ExpOutbCD
		df5 = pi0SimOutbCD

		contOutbCD = len(df2)*len(df4)/len(df3)/len(df5)
		fig, axs = plt.subplots(2, 4, figsize = (16, 8))

		# contOutbCD = 0.0
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
				simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
				simDist = (1-contOutbCD)*simDist_dvcs + contOutbCD*simDist_dvpi0
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
		axs[0, 0].set_xlim([10, 30])
		axs[0, 0].set_xticks([10, 20, 30])
		axs[0, 1].set_xlim([0, 2])
		axs[0, 1].set_xticks([0, 0.5, 1, 1.5, 2])
		axs[0, 2].set_xlim([0, 0.8])
		axs[0, 2].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
		axs[0, 3].set_xticks([0, 2.5, 5])
		axs[1, 0].set_xlim([-0.6, 0.6])
		axs[1, 0].set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
		axs[1, 1].set_xlim([-0.015, 0.015])
		axs[1, 1].set_xticks([-0.015, 0, 0.015])
		axs[1, 2].set_xlim([-0.25, 0.25])
		axs[1, 2].set_xticks([-0.25, 0, 0.25])
		axs[1, 3].set_xlim([0, 0.15])
		axs[1, 3].set_xticks([0, 0.05, 0.1, 0.15])
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/dvcsOutbCDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		df3 = epgExpOutbFD
		df1 = dvcsSimOutbFD
		df2 = bkgSimOutbFD
		df4 = pi0ExpOutbFD
		df5 = pi0SimOutbFD
		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
				simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
				simDist = (1-contOutbFD)*simDist_dvcs + contOutbFD*simDist_dvpi0
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
		axs[0, 0].set_xticks([25, 30, 35, 40, 45])
		axs[0, 0].set_xlim([25, 40])
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/dvcsOutbFDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		varstoplot = ["Mpi0", "MM2_egg", "reconPi", "coplanarity", "ME_epgg", "MM2_epgg", "MM2_ep", "MPt"]
		title = [r"$IM_{\pi^0}$", r"$MM^2_{e'\pi^0}$", r"$\theta_{\pi^0_{det.}\pi^0_{rec.}}$", r"$\Delta\phi$" , "ME"+r"${}_{ep\pi^0}$", "MM"+r"${}^{2}_{ep\pi^0}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{ep\pi^0}$"]
		unit = [GeV, GeV2, degree, degree, GeV, GeV2, GeV2, GeVc]

		df4 = pi0ExpOutbCDFT
		df5 = pi0SimOutbCDFT

		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
		axs[0, 0].set_xticks([0.1, 0.13, 0.16])
		axs[0, 0].set_xlim([0.1, 0.16])
		axs[0, 1].set_xticks([0, 0.4, 0.8, 1.2, 1.6])
		axs[0, 2].set_xticks([0, 0.4, 0.8, 1.2])
		axs[0, 3].set_xlim([0, 10])
		axs[0, 3].set_xticks([0, 5, 10])
		axs[1, 0].set_xlim([-0.5, 0.5])
		axs[1, 0].set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
		axs[1, 1].set_xlim([-0.03, 0.03])
		axs[1, 1].set_xticks([-0.03, 0, 0.03])
		axs[1, 2].set_xlim([-0.5, 0.5])
		axs[1, 2].set_xticks([-0.5, 0, 0.5])
		axs[1, 3].set_xlim([0, 0.15])
		axs[1, 3].set_xticks([0, 0.05, 0.1, 0.15])
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/pi0OutbCDFTexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		df4 = pi0ExpOutbCD
		df5 = pi0SimOutbCD

		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
axes.yaxis.get_major_ticks()[1].label1.set_visible(False)
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
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/pi0OutbCDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()

		df4 = pi0ExpOutbFD
		df5 = pi0SimOutbFD

		fig, axs = plt.subplots(2, 4, figsize = (16, 8))
		for yind in range(0, 2):
			for xind in range(0, 4):
				ind = 4*yind + xind
				simDist, bins = np.histogram(df5[varstoplot[ind]], 100, density = True)
				bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='b', linewidth=1, label = 'Experimental Data')
				axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='r', density=True, linewidth=1, label = "Simulation")
				axs[yind, xind].set_title(title[ind])
				if (unit[ind]):
					axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				else:
					axs[yind, xind].set_xlabel(title[ind])
				axs[yind, xind].set_xlim([bins[0], bins[-1]])
				axs[yind, xind].yaxis.get_major_ticks()[0].label1.set_visible(False)
		axs[0, 0].set_xticks([0.1, 0.13, 0.16])
		axs[0, 0].set_xlim([0.1, 0.16])
		axs[0, 1].set_xticks([0, 0.5, 1, 1.5, 2])
		axs[0, 2].set_xticks([0, 0.5, 1, 1.5, 2])
		axs[0, 3].set_xticks([0, 4, 8])
		axs[1, 0].set_xlim([-0.6, 0.6])
		axs[1, 0].set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
		axs[1, 1].set_xlim([-0.02, 0.02])
		axs[1, 1].set_xticks([-0.02, 0, 0.02])
		axs[1, 2].set_xlim([-0.3, 0.3])
		axs[1, 2].set_xticks([-0.3, 0, 0.3])
		axs[1, 3].set_xlim([0, 0.2])
		axs[1, 3].set_xticks([0, 0.05, 0.1, 0.15, 0.2])
		handles, labels = axs[0, 0].get_legend_handles_labels()
		lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
		plt.tight_layout()
		plt.savefig("plots/ch5/pi0OutbFDexcl.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		plt.clf()