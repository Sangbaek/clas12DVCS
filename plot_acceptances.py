#!/usr/bin/env python3
"""
Script to study the acceptances.
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

# # Electron reconstruction
# parent_MC = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/"
# parent_MC_BH = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bh/"
# parent_MC_bkg1g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_1g/"
# parent_MC_bkg2g = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/bkg_2g/"
# parent_exp = "/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/exp/"

# #inbending
# epgExpInb = pd.read_pickle(parent_exp + "dvcs.pkl")
# pi0ExpInb = pd.read_pickle(parent_exp + "pi0.pkl")
# dvcsSimInb = pd.read_pickle(parent_MC + "4893.pkl")
# bhSimInb = pd.read_pickle(parent_MC + "4892.pkl")
# bkgSimInb = pd.read_pickle(parent_MC_bkg1g + "4076.pkl")
# pi0SimInb = pd.read_pickle(parent_MC_bkg2g + "4076.pkl")

print("read exp...")
basedir = "/volatile/clas12/sangbaek/clas12DVCS/"
epgExp = pd.read_pickle(basedir + "/nphistograms/epgExp.pkl")

k = 3
i = 3
xBbins  = collection_xBbins[k]
Q2bins  = collection_Q2bins[k]
tbins   = collection_tbins [k]
phibins = collection_phibins[k]

histExpInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histExpInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histExpInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histExpInbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHDVCSInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)])
histBHDVCSInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)])
histBHDVCSInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)])
histBHDVCSInbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "cont{}".format(i)])


histExpInb = histExpInbFD + histExpInbCD + histExpInbCDFT + histExpInbCR
histBHDVCSInb = histBHDVCSInbFD + histBHDVCSInbCD + histBHDVCSInbCDFT + histBHDVCSInbCR

histExpOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histExpOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histExpOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histExpOutbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHDVCSOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)])
histBHDVCSOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)])
histBHDVCSOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)])
histBHDVCSOutbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "cont{}".format(i)])

histExpOutb = histExpOutbFD + histExpOutbCD + histExpOutbCDFT
histBHDVCSOutb = histBHDVCSOutbFD + histBHDVCSOutbCD + histBHDVCSOutbCDFT


print("reading bhs - inbending ")

# Count Rec BH
histBHInb45nA, histBHInbFD45nA, histBHInbCD45nA, histBHInbCDFT45nA, histBHInbCR45nA = 0, 0, 0, 0, 0
for jobNum in [*runs_inb_bh45nA, *runs_inb_bh50nA, *runs_inb_bh55nA]:
	histBHInb45nA = histBHInb45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec.npz".format('', k, jobNum))["hist"]
	histBHInbFD45nA = histBHInbFD45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec1.npz".format('', k, jobNum))["hist"]
	histBHInbCD45nA = histBHInbCD45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec2.npz".format('', k, jobNum))["hist"]
	histBHInbCDFT45nA = histBHInbCDFT45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec3.npz".format('', k, jobNum))["hist"]
	histBHInbCR45nA = histBHInbCR45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec4.npz".format('', k, jobNum))["hist"]

# Count Gen BH
histBHGenInb45nA, histBHGenInbFD45nA, histBHGenInbCD45nA, histBHGenInbCDFT45nA, histBHGenInbCR45nA = 0, 0, 0, 0, 0
for jobNum in [*runs_inb_bh45nA, *runs_inb_bh50nA, *runs_inb_bh55nA]:
	histBHGenInb45nA = histBHGenInb45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
	histBHGenInbFD45nA = histBHGenInbFD45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
	histBHGenInbCD45nA = histBHGenInbCD45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
	histBHGenInbCDFT45nA = histBHGenInbCDFT45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
	histBHGenInbCR45nA = histBHGenInbCR45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]

print("reading vggs  - inbending")
# Count Rec VGG
histVGGInb45nA, histVGGInbFD45nA, histVGGInbCD45nA, histVGGInbCDFT45nA, histVGGInbCR45nA = 0, 0, 0, 0, 0
for jobNum in [*runs_inb_vgg45nA, *runs_inb_vgg50nA, *runs_inb_vgg55nA]:
	histVGGInb45nA = histVGGInb45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec.npz".format('', k, jobNum))["hist"]
	histVGGInbFD45nA = histVGGInbFD45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec1.npz".format('', k, jobNum))["hist"]
	histVGGInbCD45nA = histVGGInbCD45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec2.npz".format('', k, jobNum))["hist"]
	histVGGInbCDFT45nA = histVGGInbCDFT45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec3.npz".format('', k, jobNum))["hist"]
	histVGGInbCR45nA = histVGGInbCR45nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec4.npz".format('', k, jobNum))["hist"]
#Count Gen VGG
histVGGGenInb45nA, histVGGGenInbFD45nA, histVGGGenInbCD45nA, histVGGGenInbCDFT45nA, histVGGGenInbCR45nA = 0, 0, 0, 0, 0
histVGGGenInb45nA_plus, histVGGGenInbFD45nA_plus, histVGGGenInbCD45nA_plus, histVGGGenInbCDFT45nA_plus, histVGGGenInbCR45nA_plus = 0, 0, 0, 0, 0
histVGGGenInb45nA_minus, histVGGGenInbFD45nA_minus, histVGGGenInbCD45nA_minus, histVGGGenInbCDFT45nA_minus, histVGGGenInbCR45nA_minus = 0, 0, 0, 0, 0
for jobNum in [*runs_inb_vgg45nA, *runs_inb_vgg50nA, *runs_inb_vgg55nA]:
	histVGGGenInb45nA = histVGGGenInb45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
	histVGGGenInbFD45nA = histVGGGenInbFD45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
	histVGGGenInbCD45nA = histVGGGenInbCD45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
	histVGGGenInbCDFT45nA = histVGGGenInbCDFT45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
	histVGGGenInbCR45nA = histVGGGenInbCR45nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]

print("reading bhs  - outbending")
# bh rec
histBHOutb50nA, histBHOutbFD50nA, histBHOutbCD50nA, histBHOutbCDFT50nA, histBHOutbCR50nA = 0, 0, 0, 0, 0
histBHOutb50nA_plus, histBHOutbFD50nA_plus, histBHOutbCD50nA_plus, histBHOutbCDFT50nA_plus, histBHOutbCR50nA_plus = 0, 0, 0, 0, 0
histBHOutb50nA_minus, histBHOutbFD50nA_minus, histBHOutbCD50nA_minus, histBHOutbCDFT50nA_minus, histBHOutbCR50nA_minus = 0, 0, 0, 0, 0
for jobNum in [*runs_outb_bh50nA, *runs_outb_bh40nA, *runs_outb_bh40nAT]:
	histBHOutb50nA = histBHOutb50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec.npz".format('', k, jobNum))["hist"]
	histBHOutbFD50nA = histBHOutbFD50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec1.npz".format('', k, jobNum))["hist"]
	histBHOutbCD50nA = histBHOutbCD50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec2.npz".format('', k, jobNum))["hist"]
	histBHOutbCDFT50nA = histBHOutbCDFT50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec3.npz".format('', k, jobNum))["hist"]
	histBHOutbCR50nA = histBHOutbCR50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec4.npz".format('', k, jobNum))["hist"]

# bh gen
histBHGenOutb50nA, histBHGenOutbFD50nA, histBHGenOutbCD50nA, histBHGenOutbCDFT50nA, histBHGenOutbCR50nA = 0, 0, 0, 0, 0
histBHGenOutb50nA_plus, histBHGenOutbFD50nA_plus, histBHGenOutbCD50nA_plus, histBHGenOutbCDFT50nA_plus, histBHGenOutbCR50nA_plus = 0, 0, 0, 0, 0
histBHGenOutb50nA_minus, histBHGenOutbFD50nA_minus, histBHGenOutbCD50nA_minus, histBHGenOutbCDFT50nA_minus, histBHGenOutbCR50nA_minus = 0, 0, 0, 0, 0
for jobNum in [*runs_outb_bh50nA, *runs_outb_bh40nA, *runs_outb_bh40nAT]:
	histBHGenOutb50nA = histBHGenOutb50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
	histBHGenOutbFD50nA = histBHGenOutbFD50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
	histBHGenOutbCD50nA = histBHGenOutbCD50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
	histBHGenOutbCDFT50nA = histBHGenOutbCDFT50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
	histBHGenOutbCR50nA = histBHGenOutbCR50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]

print("reading vggs  - outbending")
# vgg rec
histVGGOutb50nA, histVGGOutbFD50nA, histVGGOutbCD50nA, histVGGOutbCDFT50nA, histVGGOutbCR50nA = 0, 0, 0, 0, 0
histVGGOutb50nA_plus, histVGGOutbFD50nA_plus, histVGGOutbCD50nA_plus, histVGGOutbCDFT50nA_plus, histVGGOutbCR50nA_plus = 0, 0, 0, 0, 0
histVGGOutb50nA_minus, histVGGOutbFD50nA_minus, histVGGOutbCD50nA_minus, histVGGOutbCDFT50nA_minus, histVGGOutbCR50nA_minus = 0, 0, 0, 0, 0
for jobNum in [*runs_outb_vgg50nA, *runs_outb_vgg40nA, *runs_outb_vgg40nAT]:
	histVGGOutb50nA = histVGGOutb50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec.npz".format('', k, jobNum))["hist"]
	histVGGOutbFD50nA = histVGGOutbFD50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec1.npz".format('', k, jobNum))["hist"]
	histVGGOutbCD50nA = histVGGOutbCD50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec2.npz".format('', k, jobNum))["hist"]
	histVGGOutbCDFT50nA = histVGGOutbCDFT50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec3.npz".format('', k, jobNum))["hist"]
	histVGGOutbCR50nA = histVGGOutbCR50nA + np.load(basedir + "/nphistograms{}/binscheme{}/{}Rec4.npz".format('', k, jobNum))["hist"]

#vgg gen
histVGGGenOutb50nA, histVGGGenOutbFD50nA, histVGGGenOutbCD50nA, histVGGGenOutbCDFT50nA, histVGGGenOutbCR50nA = 0, 0, 0, 0, 0
histVGGGenOutb50nA_plus, histVGGGenOutbFD50nA_plus, histVGGGenOutbCD50nA_plus, histVGGGenOutbCDFT50nA_plus, histVGGGenOutbCR50nA_plus = 0, 0, 0, 0, 0
histVGGGenOutb50nA_minus, histVGGGenOutbFD50nA_minus, histVGGGenOutbCD50nA_minus, histVGGGenOutbCDFT50nA_minus, histVGGGenOutbCR50nA_minus = 0, 0, 0, 0, 0
for jobNum in [*runs_outb_vgg50nA, *runs_outb_vgg40nA, *runs_outb_vgg40nAT]:
	histVGGGenOutb50nA = histVGGGenOutb50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
	histVGGGenOutbFD50nA = histVGGGenOutbFD50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
	histVGGGenOutbCD50nA = histVGGGenOutbCD50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
	histVGGGenOutbCDFT50nA = histVGGGenOutbCDFT50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
	histVGGGenOutbCR50nA = histVGGGenOutbCR50nA + np.load(basedir + "/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]


#acc correction with stat error - inbending
#with VGG
accCorrectedInbFD_VGG = divideHist(histBHDVCSInbFD*histVGGGenInbFD45nA , histVGGInbFD45nA)
accCorrectedInbCD_VGG = divideHist(histBHDVCSInbCD*histVGGGenInbCD45nA , histVGGInbCD45nA)
accCorrectedInbCDFT_VGG = divideHist(histBHDVCSInbCDFT*histVGGGenInbCDFT45nA , histVGGInbCDFT45nA)
accCorrectedInbCR_VGG = divideHist(histBHDVCSInbCR*histVGGGenInbCR45nA , histVGGInbCR45nA)
accCorrectedInb_VGG = accCorrectedInbFD_VGG + accCorrectedInbCD_VGG + accCorrectedInbCDFT_VGG + accCorrectedInbCR_VGG
accCorrectedInb_VGG2 = divideHist((histBHDVCSInbFD+histBHDVCSInbCD+histBHDVCSInbCDFT+histBHDVCSInbCR)*histVGGGenInb45nA , histVGGInbFD45nA+histVGGInbCD45nA+histVGGInbCDFT45nA+histVGGInbCR45nA)
#FD inbending threshold [:, :, 0, :]
accCorrectedInb_VGG[:, :, 0, :] = accCorrectedInbCD_VGG[:, :, 0, :] + accCorrectedInbCDFT_VGG[:, :, 0, :] + accCorrectedInbCR_VGG[:, :, 0, :]
accCorrectedInb_VGG2[:, :, 0, :] = divideHist((histBHDVCSInbCD+histBHDVCSInbCDFT+histBHDVCSInbCR)*histVGGGenInb45nA , histVGGInbCD45nA+histVGGInbCDFT45nA+histVGGInbCR45nA)[:, :, 0, :]

#Correct by the ratio of covered/uncovered
accCorrectedInb_VGG = divideHist(accCorrectedInb_VGG*histVGGGenInb45nA, histVGGGenInbFD45nA + histVGGGenInbCD45nA + histVGGGenInbCDFT45nA + histVGGGenInbCR45nA)
accCorrectedInb_VGG[:, :, 0, :] = divideHist(accCorrectedInb_VGG[:, :, 0, :]*histVGGGenInb45nA[:, :, 0, :], histVGGGenInbCD45nA[:, :, 0, :] + histVGGGenInbCDFT45nA[:, :, 0, :] + histVGGGenInbCR45nA[:, :, 0, :])


# acc correction inbending, with BH
accCorrectedInbFD_BH = divideHist(histBHDVCSInbFD*histBHGenInbFD45nA , histBHInbFD45nA)
accCorrectedInbCD_BH = divideHist(histBHDVCSInbCD*histBHGenInbCD45nA , histBHInbCD45nA)
accCorrectedInbCDFT_BH = divideHist(histBHDVCSInbCDFT*histBHGenInbCDFT45nA , histBHInbCDFT45nA)
accCorrectedInbCR_BH = divideHist(histBHDVCSInbCR*histBHGenInbCR45nA , histBHInbCR45nA)
accCorrectedInb_BH = accCorrectedInbFD_BH + accCorrectedInbCD_BH + accCorrectedInbCDFT_BH + accCorrectedInbCR_BH
accCorrectedInb_BH2 = divideHist((histBHDVCSInbFD+histBHDVCSInbCD+histBHDVCSInbCDFT+histBHDVCSInbCR)*histBHGenInb45nA , histBHInbFD45nA+histBHInbCD45nA+histBHInbCDFT45nA+histBHInbCR45nA)
#FD inbending threshold [:, :, 0, :]
accCorrectedInb_BH[:, :, 0, :] = accCorrectedInbCD_BH[:, :, 0, :] + accCorrectedInbCDFT_BH[:, :, 0, :] + accCorrectedInbCR_BH[:, :, 0, :]
accCorrectedInb_BH2[:, :, 0, :] = divideHist((histBHDVCSInbCD+histBHDVCSInbCDFT+histBHDVCSInbCR)*histBHGenInb45nA , histBHInbCD45nA+histBHInbCDFT45nA+histBHInbCR45nA)[:, :, 0, :]

#Correct by the ratio of covered/uncovered
accCorrectedInb_BH = divideHist(accCorrectedInb_BH*histBHGenInb45nA, histBHGenInbFD45nA + histBHGenInbCD45nA + histBHGenInbCDFT45nA + histBHGenInbCR45nA)
accCorrectedInb_BH[:, :, 0, :] = divideHist(accCorrectedInb_BH[:, :, 0, :]*histBHGenInb45nA[:, :, 0, :], histBHGenInbCD45nA[:, :, 0, :] + histBHGenInbCDFT45nA[:, :, 0, :] + histBHGenInbCR45nA[:, :, 0, :])

# acceptance correction with stat error - outbending
# with VGG
accCorrectedOutbFD_VGG = divideHist(histBHDVCSOutbFD*histVGGGenOutbFD50nA , histVGGOutbFD50nA)
accCorrectedOutbCD_VGG = divideHist(histBHDVCSOutbCD*histVGGGenOutbCD50nA , histVGGOutbCD50nA)
accCorrectedOutbCDFT_VGG = divideHist(histBHDVCSOutbCDFT*histVGGGenOutbCDFT50nA , histVGGOutbCDFT50nA)
accCorrectedOutbCR_VGG = divideHist(histBHDVCSOutbCR*histVGGGenOutbCR50nA , histVGGOutbCR50nA)
accCorrectedOutb_VGG = accCorrectedOutbFD_VGG + accCorrectedOutbCD_VGG + accCorrectedOutbCDFT_VGG + accCorrectedOutbCR_VGG
accCorrectedOutb_VGG2 = divideHist((histBHDVCSOutbFD+histBHDVCSOutbCD+histBHDVCSOutbCDFT+histBHDVCSOutbCR)*histVGGGenOutb50nA , histVGGOutbFD50nA+histVGGOutbCD50nA+histVGGOutbCDFT50nA+histVGGOutbCR50nA)
#FD inbending threshold [:, :, [0, 1], :]
accCorrectedOutb_VGG[:, :, [0, 1], :] = accCorrectedOutbCD_VGG[:, :, [0, 1], :] + accCorrectedOutbCDFT_VGG[:, :, [0, 1], :] + accCorrectedOutbCR_VGG[:, :, [0, 1], :]
accCorrectedOutb_VGG2[:, :, [0, 1], :] = divideHist((histBHDVCSOutbCD+histBHDVCSOutbCDFT+histBHDVCSOutbCR)*histVGGGenOutb50nA , histVGGOutbCD50nA+histVGGOutbCDFT50nA+histVGGOutbCR50nA)[:, :, [0, 1], :]

#Correct by the ratio of covered/uncovered
accCorrectedOutb_VGG = divideHist(accCorrectedOutb_VGG*histVGGGenOutb50nA, histVGGGenOutbFD50nA + histVGGGenOutbCD50nA + histVGGGenOutbCDFT50nA + histVGGGenOutbCR50nA)
accCorrectedOutb_VGG[:, :, [0, 1], :] = divideHist(accCorrectedOutb_VGG[:, :, [0, 1], :]*histVGGGenOutb50nA[:, :, [0, 1], :], histVGGGenOutbCD50nA[:, :, [0, 1], :] + histVGGGenOutbCDFT50nA[:, :, [0, 1], :] + histVGGGenOutbCR50nA[:, :, [0, 1], :])

# acc correction with BH, outbending
accCorrectedOutbFD_BH = divideHist(histBHDVCSOutbFD*histBHGenOutbFD50nA , histBHOutbFD50nA)
accCorrectedOutbCD_BH = divideHist(histBHDVCSOutbCD*histBHGenOutbCD50nA , histBHOutbCD50nA)
accCorrectedOutbCDFT_BH = divideHist(histBHDVCSOutbCDFT*histBHGenOutbCDFT50nA , histBHOutbCDFT50nA)
accCorrectedOutbCR_BH = divideHist(histBHDVCSOutbCR*histBHGenOutbCR50nA , histBHOutbCR50nA)
accCorrectedOutb_BH = accCorrectedOutbFD_BH + accCorrectedOutbCD_BH + accCorrectedOutbCDFT_BH + accCorrectedOutbCR_BH
accCorrectedOutb_BH2 = divideHist((histBHDVCSOutbFD+histBHDVCSOutbCD+histBHDVCSOutbCDFT+histBHDVCSOutbCR)*histBHGenOutb50nA , histBHOutbFD50nA+histBHOutbCD50nA+histBHOutbCDFT50nA+histBHOutbCR50nA)
#FD inbending threshold [:, :, [0, 1], :]
accCorrectedOutb_BH[:, :, [0, 1], :] = accCorrectedOutbCD_BH[:, :, [0, 1], :] + accCorrectedOutbCDFT_BH[:, :, [0, 1], :] + accCorrectedOutbCR_BH[:, :, [0, 1], :]
accCorrectedOutb_BH2[:, :, [0, 1], :] = divideHist((histBHDVCSOutbCD+histBHDVCSOutbCDFT+histBHDVCSOutbCR)*histBHGenOutb50nA , histBHOutbCD50nA+histBHOutbCDFT50nA+histBHOutbCR50nA)[:, :, [0, 1], :]

#Correct by the ratio of covered/uncovered
accCorrectedOutb_BH = divideHist(accCorrectedOutb_BH*histBHGenOutb50nA, histBHGenOutbFD50nA + histBHGenOutbCD50nA + histBHGenOutbCDFT50nA + histBHGenOutbCR50nA)
accCorrectedOutb_BH[:, :, [0, 1], :] = divideHist(accCorrectedOutb_BH[:, :, [0, 1], :]*histBHGenOutb50nA[:, :, [0, 1], :], histBHGenOutbCD50nA[:, :, [0, 1], :] + histBHGenOutbCDFT50nA[:, :, [0, 1], :] + histBHGenOutbCR50nA[:, :, [0, 1], :])

ActiveAll       = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll.npz".format('', k, i))["hist"]
ActiveAny       = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format('', k, i))["hist"]
ActiveInb       = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb.npz".format('', k, i))["hist"]
ActiveOutb      = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb.npz".format('', k, i))["hist"]

ActiveAll_int       = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll_int.npz".format('', k, i))["hist"]
ActiveAny_int       = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny_int.npz".format('', k, i))["hist"]
ActiveInb_int       = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb_int.npz".format('', k, i))["hist"]
ActiveOutb_int      = np.load(basedir + "/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format('', k, i))["hist"]

phi1avg_VGG = np.load(basedir + "/nphistograms/binscheme{}/phi1avg_VGG.npz".format(k))["hist"]
xBavg_VGG   = np.load(basedir + "/nphistograms/binscheme{}/xBavg_VGG.npz".format(k))["hist"]
Q2avg_VGG   = np.load(basedir + "/nphistograms/binscheme{}/Q2avg_VGG.npz".format(k))["hist"]
t1avg_VGG   = np.load(basedir + "/nphistograms/binscheme{}/t1avg_VGG.npz".format(k))["hist"]

phi1avg_BH  = np.load(basedir + "/nphistograms/binscheme{}/phi1avg_BH.npz".format(k))["hist"]
xBavg_BH    = np.load(basedir + "/nphistograms/binscheme{}/xBavg_BH.npz".format(k))["hist"]
Q2avg_BH    = np.load(basedir + "/nphistograms/binscheme{}/Q2avg_BH.npz".format(k))["hist"]
t1avg_BH    = np.load(basedir + "/nphistograms/binscheme{}/t1avg_BH.npz".format(k))["hist"]


def badBinCondxBQ2t(xBbin, Q2bin, tbin, k = 0):
	# if k ==0:
	# 	return (xBbin==1 and Q2bin == 0) or (xBbin==0 and Q2bin==4) or (tbin==0 and xBbin==1)
	# else:
	return ~ActiveAny_int[xBbin, Q2bin, tbin, :].any()

def badBinCondxBQ2(xBbin, Q2bin, k = 0):
	# if k ==0:
	# 	return (xBbin==1 and Q2bin == 0) or (xBbin==0 and Q2bin==4)
	# else:
	return ~ActiveAny_int[xBbin, Q2bin, :, :].any()

def badBinCondxBt(xBbin, tbin, k = 0):
	# if k ==0:
	# 	return (xBbin==1 and tbin == 0)
	# else:
	return ~ActiveAny_int[xBbin, :, tbin, :].any()

num_plotQ2 = len(Q2bins) - 1
num_plotxB = len(xBbins) - 1
num_plott = len(tbins) - 1

if k == 0:
	num_plotQ2 = 6 - 1
	num_plotxB = 3 - 1
if k == 3:
	num_plotQ2 = 6
	num_plotxB = 5
	num_plott = 5


# '''start of unpol landscape'''
for tbin in range(6):
	active = 0
	ttitle = "{:.3f} GeV".format(tbins[tbin])+r"$^2<|t|<$"+"{:.3f} GeV".format(tbins[tbin+1])+r"$^2$"
	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (32.5, 45))
	for xBbin in range(num_plotxB):
		for Q2bin in range(num_plotQ2):
			#skip inactive bins
			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
				axs[num_plotQ2-Q2bin-1 , xBbin].axis('off')
				continue
			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			if len(phibin):
				pass
			else:
				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
				axs[num_plotQ2-Q2bin-1 , xBbin].axis('off')
				continue
			Nplot = 40
			xBs = np.ones(Nplot)*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
			Q2s = np.ones(Nplot)*Q2avg_BH[xBbin, Q2bin, tbin, phibin][0]
			t1s = np.ones(Nplot)*t1avg_BH[xBbin, Q2bin, tbin, phibin][0]
			phi1s = np.linspace(0, 360, Nplot)

			axs[num_plotQ2-Q2bin-1 , xBbin].hist(phibins[:-1], phibins, weights = histBHDVCSInbFD[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'k', label = "(FD, FD)")
			axs[num_plotQ2-Q2bin-1 , xBbin].hist(phibins[:-1], phibins, weights = histBHDVCSInbCD[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'r', label = "(CD, FD)")
			axs[num_plotQ2-Q2bin-1 , xBbin].hist(phibins[:-1], phibins, weights = (histBHDVCSInbCDFT+histBHDVCSInbCR)[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'g', label = "(CD, FT)")

			xBheader = "{}. ".format(xBbin)+ r"$<x_B>=$"+"{:.3f}\n".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			Q2header = "{}. ".format(Q2bin) + r"$<Q^2>=$"+"{:.3f} (GeV/c)".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"$^2$" +"\n"
			theader = "{}. ".format(tbin) + r"$<|t|>=$"+"{:.3f} GeV".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"$^2$"
			header = xBheader +Q2header + theader
			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, loc = 'left', fontsize = 20)
			axs[num_plotQ2-Q2bin-1, xBbin].set_xlim([0, 360])
			# axs[num_plotQ2-Q2bin-1, xBbin].set_ylim([0, 300])
			# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			axs[num_plotQ2-Q2bin-1, xBbin].set_xticklabels([0, 90, 180, 270, 360], fontsize = 25)
			axs[num_plotQ2-Q2bin-1, xBbin].set_yticks([0, 50, 100, 150, 200, 250, 300])
			axs[num_plotQ2-Q2bin-1, xBbin].set_yticklabels([0, 50, 100, 150, 200, 250, 300])
			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " ["+degree+"]", fontsize = 30)
			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel("Counts (bkg. subtracted)", fontsize = 30)

			# axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, loc = 'left', fontsize = 20)
			# axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")
			# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# axs[num_plotQ2-Q2bin-1, xBbin].set_xlim([0, 360])
			# axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# axs[num_plotQ2-Q2bin-1, xBbin].set_xticklabels([0, 90, 180, 270, 360], fontsize = 24)
			# axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]", fontsize = 24)
			if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
				active = 1
	# lgd = plt.figlegend([handles[idx] for idx in order_unpol],[labels[idx] for idx in order_unpol], loc='upper left', fontsize= 30, title_fontsize = 30, title = ttitle, bbox_to_anchor = (0.2, 0.8))
	lgd = plt.figlegend([handles, labels], loc='upper left', fontsize= 30, title_fontsize = 30, title = ttitle, bbox_to_anchor = (0.2, 0.8))
	fig.subplots_adjust(wspace = 0.7, hspace = 0.5)
	plt.savefig(basedir + "/plots/rawyields_inb_tbin{}.pdf".format(tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
	plt.clf()
# '''end of unpol landscape'''

# fig, axs = plt.subplots(1, 1, figsize = (8, 6))

# axs.hist(phibins[:-1], phibins, weights = accCorrectedInb_BH[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'k', label = "Acc. Separately, BH")
# axs.hist(phibins[:-1], phibins, weights = accCorrectedInb_BH2[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'r', label = "Acc. Entirely, BH")
# axs.hist(phibins[:-1], phibins, weights = accCorrectedInb_VGG[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'g', label = "Acc. Separately, BH-DVCS")

# axs.set_xlim([0, 360])
# axs.set_ylim([0, 6000])
# # axs.set_yscale('log')
# axs.set_xticks([0, 90, 180, 270, 360])
# axs.set_xticklabels([0, 90, 180, 270, 360], fontsize = 25)
# axs.set_yticks([0, 1000, 2000, 3000, 4000, 5000, 6000])
# axs.set_yticklabels([0, 1000, 2000, 3000, 4000, 5000, 6000])
# axs.set_xlabel(r"$\phi$" + " ["+degree+"]", fontsize = 30)
# axs.set_ylabel("Acc. Corrected Counts (bkg. subtracted)", fontsize = 30)
# axs.set_title("(b)  Acc. Corrected Yields", fontsize = 30)

# xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}".format(xBbins[xBbin+1]) + "\n"
# Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f} ".format(Q2bins[Q2bin+1])+ "\n"
# theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f} ".format(tbins[tbin+1])
# header = xBheader + Q2header + theader
# handles, labels = axs.get_legend_handles_labels()
# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, bbox_to_anchor = (0.8, 1.1), title = header)

# plt.savefig(basedir+"/plots/accCorrectedYields{}{}{}.pdf".format(xBbin, Q2bin, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
# plt.clf()

# fig, axs = plt.subplots(1, 1, figsize = (8, 6))

# axs.hist(phibins[:-1], phibins, weights = divideHist(histBHDVCSInb, accCorrectedInb_BH)[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'k', label = "Acc. Separately, BH")
# axs.hist(phibins[:-1], phibins, weights = divideHist(histBHDVCSInb, accCorrectedInb_BH2)[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'r', label = "Acc. Entirely, BH")
# axs.hist(phibins[:-1], phibins, weights = divideHist(histBHDVCSInb, accCorrectedInb_VGG)[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'g', label = "Acc. Separately, BH-DVCS")

# axs.set_xlim([0, 360])
# # axs.set_ylim([0, 200])
# # axs.set_yscale('log')
# axs.set_xticks([0, 90, 180, 270, 360])
# axs.set_xticklabels([0, 90, 180, 270, 360], fontsize = 25)
# axs.set_yticks([0, 0.02, 0.04, 0.06, 0.08])
# axs.set_yticklabels([0, 0.02, 0.04, 0.06, 0.08])
# axs.set_xlabel(r"$\phi$" + " ["+degree+"]", fontsize = 30)
# axs.set_ylabel("Effective Acceptances", fontsize = 30)
# axs.set_title("(c)  Effective Acceptances", fontsize = 30)

# xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}".format(xBbins[xBbin+1]) + "\n"
# Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f} ".format(Q2bins[Q2bin+1])+ "\n"
# theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f} ".format(tbins[tbin+1])
# header = xBheader + Q2header + theader
# handles, labels = axs.get_legend_handles_labels()
# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, bbox_to_anchor = (0.8, 1.1), title = header)

# plt.savefig(basedir+"/plots/EffectiveAcc{}{}{}.pdf".format(xBbin, Q2bin, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
# plt.clf()
