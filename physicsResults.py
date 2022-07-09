#!/usr/bin/env python3
"""
Main Script to save the cross sections
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

parser.add_argument("-sc","--savecont", help="save cont", action = "store_true")
parser.add_argument("-sx","--savexsec", help="save xsec", action = "store_true")
parser.add_argument("-ct","--contplot", help="save cont plots", action = "store_true")
parser.add_argument("-sp","--saveplot", help="save plots", action = "store_true")
parser.add_argument("-pr","--parent", help="parent", default = "/volatile/clas12/sangbaek/nov2021/")
parser.add_argument("-k","--binscheme", help="binning scheme number", default=None)
parser.add_argument("-kstart","--kstart", help="binning scheme start", default=None)
parser.add_argument("-kend","--kend", help="binning scheme end", default=None)
parser.add_argument("-opt","--optionaltag", help="optional tag: none, eb, 2sigma, 3sigma, sm09, sm11", default=None)
parser.add_argument("-sb", "--savebinVolume", help = "save binVolume", action = "store_true")
parser.add_argument("-sy", "--saveyields", help = "save yields", action = "store_true")
parser.add_argument("-sk", "--savekine", help = "save kinematic variables", action = "store_true")
parser.add_argument("-ss", "--savesyst", help = "save systematic uncertainties", action = "store_true")
parser.add_argument("-ss2", "--savesyst2", help = "save systematic uncertainties", action = "store_true")
parser.add_argument("-rp", "--radplot", help = "save rad plots", action = "store_true")

args = parser.parse_args()

if args.optionaltag:
	optionaltag = args.optionaltag
	istart = len(collection_cont_xBbins) - 1
	kstart = len(collection_xBbins) - 1
	kend = len(collection_xBbins)
else:
	optionaltag = ''
	istart = 0
	kstart = 0

	if args.kstart:
		kstart = int(args.kstart)
	else:
		kstart = len(collection_xBbins) - 1 
	if args.kend:
		kend = int(args.kend)
	else:
		kend = len(collection_xBbins)

inbcharge_epg, outbcharge_epg = 30398709.057119943, 32085430.046131916
charge_epg = inbcharge_epg + outbcharge_epg

# def FourierSeries(args, x):
def FourierSeries(x, a, b, c):
#     df = args
    # a, b, c = args
    return a + b * np.cos(np.radians(x)) + c* np.cos(2*np.radians(x))

def lstsq_FourierSeries(args, x, y):
#     print(args, x, y)
    return FourierSeries(args, x) - y

if args.savecont:
	# read exp
	parent = args.parent

	parent_MC_bkg1g_inb = parent + "convPkl_full{}/inb/bkg_1g/".format(optionaltag)
	parent_MC_bkg2g_inb = parent + "convPkl_full{}/inb/bkg_2g/".format(optionaltag)
	parent_exp_inb = parent + "convPkl_full{}/inb/exp/".format(optionaltag)

	parent_MC_bkg1g_outb = parent + "convPkl_full{}/outb/bkg_1g/".format(optionaltag)
	parent_MC_bkg2g_outb = parent + "convPkl_full{}/outb/bkg_2g/".format(optionaltag)
	parent_exp_outb = parent + "convPkl_full{}/outb/exp/".format(optionaltag)

	print("read exp...")
	epgExpInb = pd.read_pickle(parent_exp_inb + "dvcs.pkl")
	pi0ExpInb = pd.read_pickle(parent_exp_inb + "pi0.pkl")

	epgExpOutb = pd.read_pickle(parent_exp_outb + "dvcs.pkl")
	pi0ExpOutb = pd.read_pickle(parent_exp_outb + "pi0.pkl")

	runConfig =  pd.read_csv(parent + "rga_runconfigs.txt", sep = '\t')
	beamCurrent_dict = {runConfig.RunNum[i]:runConfig.BeamCurrent_request[i] for i in runConfig.index}


	epgExpInb.loc[:, "polarity"] = -1
	epgExpOutb.loc[:, "polarity"] = +1
	pi0ExpInb.loc[:, "polarity"] = -1
	pi0ExpOutb.loc[:, "polarity"] = +1

	epgExpInb.loc[:, "beamCurrent"] = epgExpInb.RunNum.map(beamCurrent_dict)
	epgExpOutb.loc[:, "beamCurrent"] = epgExpOutb.RunNum.map(beamCurrent_dict)
	pi0ExpInb.loc[:, "beamCurrent"] = pi0ExpInb.RunNum.map(beamCurrent_dict)
	pi0ExpOutb.loc[:, "beamCurrent"] = pi0ExpOutb.RunNum.map(beamCurrent_dict)

	# runNum = np.sort(epgExpInb.RunNum.unique())
	# charges = []
	# charges_QA = []
	# for run in runNum:
	#     charges.append(epgExpInb.loc[epgExpInb.RunNum==run, "beamQ"].max() - epgExpInb.loc[epgExpInb.RunNum==run, "beamQ"].min())
	#     charges_QA.append(epgExpInb.loc[epgExpInb.RunNum==run, "beamQ_QA"].max())
	# charges = np.array(charges)
	# charges_QA = np.array(charges_QA)
	# inbcharge_epg = np.sum(charges_QA)

	# runNum = np.sort(epgExpOutb.RunNum.unique())
	# charges = []
	# charges_QA = []
	# for run in runNum:
	#     charges.append(epgExpOutb.loc[epgExpOutb.RunNum==run, "beamQ"].max() - epgExpOutb.loc[epgExpOutb.RunNum==run, "beamQ"].min())
	#     charges_QA.append(epgExpOutb.loc[epgExpOutb.RunNum==run, "beamQ_QA"].max())
	# charges = np.array(charges)
	# charges_QA = np.array(charges_QA)
	# outbcharge_epg = np.sum(charges_QA)

	epgExpInb = makeReduced(epgExpInb)
	pi0ExpInb = makeReduced(pi0ExpInb)
	epgExpOutb = makeReduced(epgExpOutb)
	pi0ExpOutb = makeReduced(pi0ExpOutb)

	epgExp = pd.concat([epgExpInb, epgExpOutb])
	pi0Exp = pd.concat([pi0ExpInb, pi0ExpOutb])

	# read pi0 simulations
	print("read bkg...")
	df_bkg1gs_inb = []
	for jobNum in runs_inb_bkg50nA:
		df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g_inb, jobNum, -1, 50))
	for jobNum in runs_inb_bkg55nA:
		df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g_inb, jobNum, -1, 55))
	for jobNum in runs_inb_bkg45nA:
		df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g_inb, jobNum, -1, 45))
	for jobNum in runs_inb_bkg0nA:
		df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g_inb, jobNum, -1, 0))

	df_bkg2gs_inb = []
	for jobNum in runs_inb_bkg50nA:
		df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g_inb, jobNum, -1, 50))
	for jobNum in runs_inb_bkg55nA:
		df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g_inb, jobNum, -1, 55))
	for jobNum in runs_inb_bkg45nA:
		df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g_inb, jobNum, -1, 45))
	for jobNum in runs_inb_bkg0nA:
		df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g_inb, jobNum, -1, 0))
		
	df_bkg1gs_inb = pd.concat(df_bkg1gs_inb)
	df_bkg2gs_inb = pd.concat(df_bkg2gs_inb)

	df_bkg1gs_outb = []
	for jobNum in runs_outb_bkg50nA:
		df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g_outb, jobNum, -1, 50))
	for jobNum in runs_outb_bkg40nA:
		df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g_outb, jobNum, -1, 40))
	for jobNum in runs_outb_bkg0nA:
		df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g_outb, jobNum, -1, 0))
	for jobNum in runs_outb_bkg40nAT:
		df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g_outb, jobNum, -1, 40))

	df_bkg2gs_outb = []
	for jobNum in runs_outb_bkg50nA:
		df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g_outb, jobNum, -1, 50))
	for jobNum in runs_outb_bkg40nA:
		df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g_outb, jobNum, -1, 40))
	for jobNum in runs_outb_bkg0nA:
		df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g_outb, jobNum, -1, 0))
	for jobNum in runs_outb_bkg40nAT:
		df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g_outb, jobNum, -1, 40))
		
	df_bkg1gs_outb = pd.concat(df_bkg1gs_outb)
	df_bkg2gs_outb = pd.concat(df_bkg2gs_outb)

	i = 0	# trial
	for i in range(istart, len(collection_cont_xBbins)):

		xBbins = collection_cont_xBbins  [i]
		Q2bins = collection_cont_Q2bins  [i]
		tbins = collection_cont_tbins    [i]
		phibins = collection_cont_phibins[i]

		histExpInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbFD, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.polarity == -1) & (pi0Exp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgInbFD, bins = np.histogramdd(df_bkg1gs_inb.loc[df_bkg1gs_inb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histRefInbFD, bins = np.histogramdd(df_bkg2gs_inb.loc[df_bkg2gs_inb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbCD, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.polarity == -1) & (pi0Exp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgInbCD, bins = np.histogramdd(df_bkg1gs_inb.loc[df_bkg1gs_inb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histRefInbCD, bins = np.histogramdd(df_bkg2gs_inb.loc[df_bkg2gs_inb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbCDFT, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.polarity == -1) & (pi0Exp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgInbCDFT, bins = np.histogramdd(df_bkg1gs_inb.loc[df_bkg1gs_inb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histRefInbCDFT, bins = np.histogramdd(df_bkg2gs_inb.loc[df_bkg2gs_inb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgInbCR, bins = np.histogramdd(df_bkg1gs_inb.loc[df_bkg1gs_inb.config == 4 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

		histExpOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbFD, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.polarity == 1) & (pi0Exp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgOutbFD, bins = np.histogramdd(df_bkg1gs_outb.loc[df_bkg1gs_outb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histRefOutbFD, bins = np.histogramdd(df_bkg2gs_outb.loc[df_bkg2gs_outb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbCD, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.polarity == 1) & (pi0Exp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgOutbCD, bins = np.histogramdd(df_bkg1gs_outb.loc[df_bkg1gs_outb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histRefOutbCD, bins = np.histogramdd(df_bkg2gs_outb.loc[df_bkg2gs_outb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbCDFT, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.polarity == 1) & (pi0Exp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgOutbCDFT, bins = np.histogramdd(df_bkg1gs_outb.loc[df_bkg1gs_outb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histRefOutbCDFT, bins = np.histogramdd(df_bkg2gs_outb.loc[df_bkg2gs_outb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histBkgOutbCR, bins = np.histogramdd(df_bkg1gs_outb.loc[df_bkg1gs_outb.config == 4 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

		histExpInbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbFD_plus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == -1) & (pi0Exp.polarity == -1) & (pi0Exp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbCD_plus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == -1) & (pi0Exp.polarity == -1) & (pi0Exp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbCDFT_plus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == -1) & (pi0Exp.polarity == -1) & (pi0Exp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

		histExpOutbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbFD_plus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == -1) & (pi0Exp.polarity == 1) & (pi0Exp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbCD_plus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == -1) & (pi0Exp.polarity == 1) & (pi0Exp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbCDFT_plus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == -1) & (pi0Exp.polarity == 1) & (pi0Exp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

		histExpInbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbFD_minus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == 1) & (pi0Exp.polarity == -1) & (pi0Exp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbCD_minus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == 1) & (pi0Exp.polarity == -1) & (pi0Exp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0InbCDFT_minus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == 1) & (pi0Exp.polarity == -1) & (pi0Exp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpInbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

		histExpOutbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbFD_minus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == 1) & (pi0Exp.polarity == 1) & (pi0Exp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbCD_minus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == 1) & (pi0Exp.polarity == 1) & (pi0Exp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histPi0OutbCDFT_minus, bins = np.histogramdd(pi0Exp.loc[(pi0Exp.helicity == 1) & (pi0Exp.polarity == 1) & (pi0Exp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		histExpOutbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

		epgExp.loc[:, "newxBbin{}".format(i)] = (len(phibins)-1)*(len(tbins)-1)*(len(Q2bins)-1)*(np.digitize(epgExp.xB, xBbins)-1)
		epgExp.loc[:, "newQ2bin{}".format(i)] = (len(phibins)-1)*(len(tbins)-1)*(np.digitize(epgExp.Q2, Q2bins)-1)
		epgExp.loc[:, "newtbin{}".format(i)] = (len(phibins)-1)*(np.digitize(epgExp.t1, tbins)-1)
		epgExp.loc[:, "newphibin{}".format(i)] = np.digitize(epgExp.phi1, phibins)-1
		epgExp.loc[:, "newbin{}".format(i)] = np.sum(epgExp.loc[:, ["newxBbin{}".format(i), "newQ2bin{}".format(i), "newtbin{}".format(i), "newphibin{}".format(i)]], axis = 1)

		contInbFD = divideHist(histBkgInbFD*histPi0InbFD, histRefInbFD*histExpInbFD)
		contInbCD = divideHist(histBkgInbCD*histPi0InbCD, histRefInbCD*histExpInbCD)
		contInbCDFT = divideHist(histBkgInbCDFT*histPi0InbCDFT, histRefInbCDFT*histExpInbCDFT)
		contInbCR = divideHist(histBkgInbCR*histPi0InbCDFT, histRefInbCDFT*histExpInbCR)
		contOutbFD = divideHist(histBkgOutbFD*histPi0OutbFD, histRefOutbFD*histExpOutbFD)
		contOutbCD = divideHist(histBkgOutbCD*histPi0OutbCD, histRefOutbCD*histExpOutbCD)
		contOutbCDFT = divideHist(histBkgOutbCDFT*histPi0OutbCDFT, histRefOutbCDFT*histExpOutbCDFT)
		contOutbCR = divideHist(histBkgOutbCR*histPi0OutbCDFT, histRefOutbCDFT*histExpOutbCR)

		unccontInbFD = contInbFD*np.sqrt(inverseHist(histBkgInbFD)+inverseHist(histPi0InbFD)+inverseHist(histRefInbFD)+inverseHist(histExpInbFD))
		unccontInbCD = contInbCD*np.sqrt(inverseHist(histBkgInbCD)+inverseHist(histPi0InbCD)+inverseHist(histRefInbCD)+inverseHist(histExpInbCD))
		unccontInbCDFT = contInbCDFT*np.sqrt(inverseHist(histBkgInbCDFT)+inverseHist(histPi0InbCDFT)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCDFT))
		unccontInbCR = contInbCR*np.sqrt(inverseHist(histBkgInbCR)+inverseHist(histPi0InbCDFT)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCR))
		unccontOutbFD = contOutbFD*np.sqrt(inverseHist(histBkgOutbFD)+inverseHist(histPi0OutbFD)+inverseHist(histRefOutbFD)+inverseHist(histExpOutbFD))
		unccontOutbCD = contOutbCD*np.sqrt(inverseHist(histBkgOutbCD)+inverseHist(histPi0OutbCD)+inverseHist(histRefOutbCD)+inverseHist(histExpOutbCD))
		unccontOutbCDFT = contOutbCDFT*np.sqrt(inverseHist(histBkgOutbCDFT)+inverseHist(histPi0OutbCDFT)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCDFT))
		unccontOutbCR = contOutbCR*np.sqrt(inverseHist(histBkgOutbCR)+inverseHist(histPi0OutbCDFT)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCR))

		contInbFD_plus = divideHist(histBkgInbFD*histPi0InbFD_plus, histRefInbFD*histExpInbFD_plus)
		contInbCD_plus = divideHist(histBkgInbCD*histPi0InbCD_plus, histRefInbCD*histExpInbCD_plus)
		contInbCDFT_plus = divideHist(histBkgInbCDFT*histPi0InbCDFT_plus, histRefInbCDFT*histExpInbCDFT_plus)
		contInbCR_plus = divideHist(histBkgInbCR*histPi0InbCDFT_plus, histRefInbCDFT*histExpInbCR_plus)
		contOutbFD_plus = divideHist(histBkgOutbFD*histPi0OutbFD_plus, histRefOutbFD*histExpOutbFD_plus)
		contOutbCD_plus = divideHist(histBkgOutbCD*histPi0OutbCD_plus, histRefOutbCD*histExpOutbCD_plus)
		contOutbCDFT_plus = divideHist(histBkgOutbCDFT*histPi0OutbCDFT_plus, histRefOutbCDFT*histExpOutbCDFT_plus)
		contOutbCR_plus = divideHist(histBkgOutbCR*histPi0OutbCDFT_plus, histRefOutbCDFT*histExpOutbCR_plus)

		unccontInbFD_plus = contInbFD_plus*np.sqrt(inverseHist(histBkgInbFD)+inverseHist(histPi0InbFD_plus)+inverseHist(histRefInbFD)+inverseHist(histExpInbFD_plus))
		unccontInbCD_plus = contInbCD_plus*np.sqrt(inverseHist(histBkgInbCD)+inverseHist(histPi0InbCD_plus)+inverseHist(histRefInbCD)+inverseHist(histExpInbCD_plus))
		unccontInbCDFT_plus = contInbCDFT_plus*np.sqrt(inverseHist(histBkgInbCDFT)+inverseHist(histPi0InbCDFT_plus)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCDFT_plus))
		unccontInbCR_plus = contInbCR_plus*np.sqrt(inverseHist(histBkgInbCR)+inverseHist(histPi0InbCDFT_plus)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCR_plus))
		unccontOutbFD_plus = contOutbFD_plus*np.sqrt(inverseHist(histBkgOutbFD)+inverseHist(histPi0OutbFD_plus)+inverseHist(histRefOutbFD)+inverseHist(histExpOutbFD_plus))
		unccontOutbCD_plus = contOutbCD_plus*np.sqrt(inverseHist(histBkgOutbCD)+inverseHist(histPi0OutbCD_plus)+inverseHist(histRefOutbCD)+inverseHist(histExpOutbCD_plus))
		unccontOutbCDFT_plus = contOutbCDFT_plus*np.sqrt(inverseHist(histBkgOutbCDFT)+inverseHist(histPi0OutbCDFT_plus)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCDFT_plus))
		unccontOutbCR_plus = contOutbCR_plus*np.sqrt(inverseHist(histBkgOutbCR)+inverseHist(histPi0OutbCDFT_plus)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCR_plus))

		contInbFD_minus = divideHist(histBkgInbFD*histPi0InbFD_minus, histRefInbFD*histExpInbFD_minus)
		contInbCD_minus = divideHist(histBkgInbCD*histPi0InbCD_minus, histRefInbCD*histExpInbCD_minus)
		contInbCDFT_minus = divideHist(histBkgInbCDFT*histPi0InbCDFT_minus, histRefInbCDFT*histExpInbCDFT_minus)
		contInbCR_minus = divideHist(histBkgInbCR*histPi0InbCDFT_minus, histRefInbCDFT*histExpInbCR_minus)
		contOutbFD_minus = divideHist(histBkgOutbFD*histPi0OutbFD_minus, histRefOutbFD*histExpOutbFD_minus)
		contOutbCD_minus = divideHist(histBkgOutbCD*histPi0OutbCD_minus, histRefOutbCD*histExpOutbCD_minus)
		contOutbCDFT_minus = divideHist(histBkgOutbCDFT*histPi0OutbCDFT_minus, histRefOutbCDFT*histExpOutbCDFT_minus)
		contOutbCR_minus = divideHist(histBkgOutbCR*histPi0OutbCDFT_minus, histRefOutbCDFT*histExpOutbCR_minus)

		unccontInbFD_minus = contInbFD_minus*np.sqrt(inverseHist(histBkgInbFD)+inverseHist(histPi0InbFD_minus)+inverseHist(histRefInbFD)+inverseHist(histExpInbFD_minus))
		unccontInbCD_minus = contInbCD_minus*np.sqrt(inverseHist(histBkgInbCD)+inverseHist(histPi0InbCD_minus)+inverseHist(histRefInbCD)+inverseHist(histExpInbCD_minus))
		unccontInbCDFT_minus = contInbCDFT_minus*np.sqrt(inverseHist(histBkgInbCDFT)+inverseHist(histPi0InbCDFT_minus)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCDFT_minus))
		unccontInbCR_minus = contInbCR_minus*np.sqrt(inverseHist(histBkgInbCR)+inverseHist(histPi0InbCDFT_minus)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCR_minus))
		unccontOutbFD_minus = contOutbFD_minus*np.sqrt(inverseHist(histBkgOutbFD)+inverseHist(histPi0OutbFD_minus)+inverseHist(histRefOutbFD)+inverseHist(histExpOutbFD_minus))
		unccontOutbCD_minus = contOutbCD_minus*np.sqrt(inverseHist(histBkgOutbCD)+inverseHist(histPi0OutbCD_minus)+inverseHist(histRefOutbCD)+inverseHist(histExpOutbCD_minus))
		unccontOutbCDFT_minus = contOutbCDFT_minus*np.sqrt(inverseHist(histBkgOutbCDFT)+inverseHist(histPi0OutbCDFT_minus)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCDFT_minus))
		unccontOutbCR_minus = contOutbCR_minus*np.sqrt(inverseHist(histBkgOutbCR)+inverseHist(histPi0OutbCDFT_minus)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCR_minus))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(contInbCR.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(contOutbCR.flatten())))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontInbFD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontInbCD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontInbCDFT.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(unccontInbCR.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontOutbFD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCDFT.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCR.flatten())))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT_plus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(contInbCR_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "cont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(contOutbCR_plus.flatten())))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontInbFD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontInbCD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontInbCDFT_plus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(unccontInbCR_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontOutbFD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCD_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCDFT_plus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "unccont_plus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCR_plus.flatten())))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT_minus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(contInbCR_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "cont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(contOutbCR_minus.flatten())))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontInbFD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontInbCD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontInbCDFT_minus.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(unccontInbCR_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontOutbFD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCD_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCDFT_minus.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "unccont_minus{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCR_minus.flatten())))

		epgExp.loc[epgExp.loc[:, "cont{}".format(i)] > 1, "cont{}".format(i)] = 1
		epgExp.loc[epgExp.loc[:, "cont_plus{}".format(i)] > 1, "cont_plus{}".format(i)] = 1
		epgExp.loc[epgExp.loc[:, "cont_minus{}".format(i)] > 1, "cont_minus{}".format(i)] = 1

		print("saved contaminations {}".format(i))

	print("clear memory...")

	del df_bkg1gs_inb
	del df_bkg1gs_outb
	del df_bkg2gs_inb
	del df_bkg2gs_outb
	gc.collect()

	epgExp.to_pickle("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/epgExp.pkl".format(optionaltag))

if args.savebinVolume:
	for k in range(kstart, kend):
		print("saving bin scheme {}".format(k))
		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		binVolume = np.zeros((len(xBbins)-1, len(Q2bins)-1, len(tbins)-1, len(phibins)-1))
		N = 1000
		n = 1000

		for xBind in range(len(xBbins)-1):
			print(xBind)
			for Q2ind in range(len(Q2bins)-1):
				for tind in range(len(tbins)-1):
					for phiind in range(len(phibins)-1):
						xBi = xBbins[xBind]
						xBf = xBbins[xBind+1]
						Q2i = Q2bins[Q2ind]
						Q2f = Q2bins[Q2ind+1]
						ti = tbins[tind]
						tf = tbins[tind+1]
						phii = np.radians(phibins[phiind])
						phif = np.radians(phibins[phiind+1])

						ratio = []
						for trial in range(n):
							xB = np.random.uniform(xBi, xBf, N)
							Q2 = np.random.uniform(Q2i, Q2f, N)
							t1 = np.random.uniform(ti, tf, N)
							phi1 = np.random.uniform(phii, phif, N)
							nub  = nu(xB, Q2, t1, phi1)
							yb= y(xB, Q2, t1, phi1)
							tmin1 = np.abs(tmin(xB, Q2, t1, phi1))
							tmax1 = np.abs(tmax(xB, Q2, t1, phi1))
							cond = (nub>2 )&(nub<10.604 - 2 ) & (W2(xB, Q2, t1, phi1)>4) & (t1>tmin1) & (Kfac2(xB, Q2, t1, phi1)>0)
							ratio.append(np.sum(cond)/N)
						binVolume[xBind, Q2ind, tind, phiind] = (xBf-xBi)*(Q2f-Q2i)*(tf-ti)*(phif-phii)*np.mean(ratio)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/binVolume.npz".format(k), hist = binVolume)

if args.saveyields:
	print("read exp...")
	epgExp = pd.read_pickle("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/epgExp.pkl".format(optionaltag))
	for k in range(kstart, kend):
		print("bin scheme {}".format(k))
		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		for i in range(istart, len(collection_cont_xBbins)):
			#inbending
			#exp - all
			histExpInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)])
			histBHDVCSInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)])
			histBHDVCSInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)])
			histBHDVCSInbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "cont{}".format(i)])

			histBkgUncInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "unccont{}".format(i)])
			histBkgUncInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "unccont{}".format(i)])
			histBkgUncInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "unccont{}".format(i)])
			histBkgUncInbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 4), "unccont{}".format(i)])
			histExpUncInbFD = np.sqrt(histExpInbFD)
			histExpUncInbCD = np.sqrt(histExpInbCD)
			histExpUncInbCDFT = np.sqrt(histExpInbCDFT)
			histExpUncInbCR = np.sqrt(histExpInbCR)

			#exp - plus helicity (helcity flip issue)
			histExpInbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSInbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 1), "cont_plus{}".format(i)])
			histBHDVCSInbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 2), "cont_plus{}".format(i)])
			histBHDVCSInbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 3), "cont_plus{}".format(i)])
			histBHDVCSInbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 4), "cont_plus{}".format(i)])

			histBkgUncInbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 1), "unccont_plus{}".format(i)])
			histBkgUncInbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 2), "unccont_plus{}".format(i)])
			histBkgUncInbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 3), "unccont_plus{}".format(i)])
			histBkgUncInbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == -1) & (epgExp.config == 4), "unccont_plus{}".format(i)])
			histExpUncInbFD_plus = np.sqrt(histExpInbFD_plus)
			histExpUncInbCD_plus = np.sqrt(histExpInbCD_plus)
			histExpUncInbCDFT_plus = np.sqrt(histExpInbCDFT_plus)
			histExpUncInbCR_plus = np.sqrt(histExpInbCR_plus)

			#exp -minus helicity (helcity flip issue)
			histExpInbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSInbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 1), "cont_minus{}".format(i)])
			histBHDVCSInbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 2), "cont_minus{}".format(i)])
			histBHDVCSInbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 3), "cont_minus{}".format(i)])
			histBHDVCSInbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 4), "cont_minus{}".format(i)])

			histBkgUncInbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 1), "unccont_minus{}".format(i)])
			histBkgUncInbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 2), "unccont_minus{}".format(i)])
			histBkgUncInbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 3), "unccont_minus{}".format(i)])
			histBkgUncInbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == -1) & (epgExp.config == 4), "unccont_minus{}".format(i)])
			histExpUncInbFD_minus = np.sqrt(histExpInbFD_minus)
			histExpUncInbCD_minus = np.sqrt(histExpInbCD_minus)
			histExpUncInbCDFT_minus = np.sqrt(histExpInbCDFT_minus)
			histExpUncInbCR_minus = np.sqrt(histExpInbCR_plus)

			histExpInb = histExpInbFD + histExpInbCD + histExpInbCDFT + histExpInbCR
			histBHDVCSInb = histBHDVCSInbFD + histBHDVCSInbCD + histBHDVCSInbCDFT + histBHDVCSInbCR

			ActiveInb = histExpInb>20
			if k <  2:
				ActiveInb[:, 0, :, :] = False
			elif k == 2:
				ActiveInb[:, 0, :, :] = False
				ActiveInb[2, 1, 2, 4] = False
				ActiveInb[2, 1, 2, 5] = False
				ActiveInb[3, 1, 2, 22] = False
			elif k == 3:
				ActiveInb[:, 0, :, :] = False
				ActiveInb[:, 1, :, :] = False
			ActiveInb_int = np.stack([np.sum(ActiveInb, axis=-1)>8]*(len(phibins)-1), axis = -1)

			print("reading bhs - inbending ")

			# Count Rec BH
			histBHInb45nA, histBHInbFD45nA, histBHInbCD45nA, histBHInbCDFT45nA, histBHInbCR45nA = 0, 0, 0, 0, 0
			histBHInb45nA_plus, histBHInbFD45nA_plus, histBHInbCD45nA_plus, histBHInbCDFT45nA_plus, histBHInbCR45nA_plus = 0, 0, 0, 0, 0
			histBHInb45nA_minus, histBHInbFD45nA_minus, histBHInbCD45nA_minus, histBHInbCDFT45nA_minus, histBHInbCR45nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_inb_bh45nA:
				histBHInb45nA = histBHInb45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbFD45nA = histBHInbFD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCD45nA = histBHInbCD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCDFT45nA = histBHInbCDFT45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCR45nA = histBHInbCR45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInb45nA_plus = histBHInb45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recplus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbFD45nA_plus = histBHInbFD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCD45nA_plus = histBHInbCD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCDFT45nA_plus = histBHInbCDFT45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCR45nA_plus = histBHInbCR45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInb45nA_minus = histBHInb45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recminus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbFD45nA_minus = histBHInbFD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1minus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCD45nA_minus = histBHInbCD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2minus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCDFT45nA_minus = histBHInbCDFT45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3minus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHInbCR45nA_minus = histBHInbCR45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4minus.npz".format(optionaltag, k, jobNum))["hist"]

			# Count Gen BH
			histBHGenInb45nA, histBHGenInbFD45nA, histBHGenInbCD45nA, histBHGenInbCDFT45nA, histBHGenInbCR45nA = 0, 0, 0, 0, 0
			histBHGenInb45nA_plus, histBHGenInbFD45nA_plus, histBHGenInbCD45nA_plus, histBHGenInbCDFT45nA_plus, histBHGenInbCR45nA_plus = 0, 0, 0, 0, 0
			histBHGenInb45nA_minus, histBHGenInbFD45nA_minus, histBHGenInbCD45nA_minus, histBHGenInbCDFT45nA_minus, histBHGenInbCR45nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_inb_bh45nA:
				histBHGenInb45nA = histBHGenInb45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histBHGenInbFD45nA = histBHGenInbFD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histBHGenInbCD45nA = histBHGenInbCD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histBHGenInbCDFT45nA = histBHGenInbCDFT45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
				histBHGenInbCR45nA = histBHGenInbCR45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]
				histBHGenInb45nA_plus = histBHGenInb45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genplus.npz".format(k, jobNum))["hist"]
				histBHGenInbFD45nA_plus = histBHGenInbFD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1plus.npz".format(k, jobNum))["hist"]
				histBHGenInbCD45nA_plus = histBHGenInbCD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2plus.npz".format(k, jobNum))["hist"]
				histBHGenInbCDFT45nA_plus = histBHGenInbCDFT45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3plus.npz".format(k, jobNum))["hist"]
				histBHGenInbCR45nA_plus = histBHGenInbCR45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4plus.npz".format(k, jobNum))["hist"]
				histBHGenInb45nA_minus = histBHGenInb45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genminus.npz".format(k, jobNum))["hist"]
				histBHGenInbFD45nA_minus = histBHGenInbFD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1minus.npz".format(k, jobNum))["hist"]
				histBHGenInbCD45nA_minus = histBHGenInbCD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2minus.npz".format(k, jobNum))["hist"]
				histBHGenInbCDFT45nA_minus = histBHGenInbCDFT45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3minus.npz".format(k, jobNum))["hist"]
				histBHGenInbCR45nA_minus = histBHGenInbCR45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4minus.npz".format(k, jobNum))["hist"]

			print("reading vggs  - inbending")
			# Count Rec VGG
			histVGGInb45nA, histVGGInbFD45nA, histVGGInbCD45nA, histVGGInbCDFT45nA, histVGGInbCR45nA = 0, 0, 0, 0, 0
			histVGGInb45nA_plus, histVGGInbFD45nA_plus, histVGGInbCD45nA_plus, histVGGInbCDFT45nA_plus, histVGGInbCR45nA_plus = 0, 0, 0, 0, 0
			histVGGInb45nA_minus, histVGGInbFD45nA_minus, histVGGInbCD45nA_minus, histVGGInbCDFT45nA_minus, histVGGInbCR45nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_inb_vgg45nA:
				histVGGInb45nA = histVGGInb45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbFD45nA = histVGGInbFD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCD45nA = histVGGInbCD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCDFT45nA = histVGGInbCDFT45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCR45nA = histVGGInbCR45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInb45nA_plus = histVGGInb45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recplus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbFD45nA_plus = histVGGInbFD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCD45nA_plus = histVGGInbCD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCDFT45nA_plus = histVGGInbCDFT45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCR45nA_plus = histVGGInbCR45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInb45nA_minus = histVGGInb45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recminus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbFD45nA_minus = histVGGInbFD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1minus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCD45nA_minus = histVGGInbCD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2minus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCDFT45nA_minus = histVGGInbCDFT45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3minus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGInbCR45nA_minus = histVGGInbCR45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4minus.npz".format(optionaltag, k, jobNum))["hist"]

			#Count Gen VGG
			histVGGGenInb45nA, histVGGGenInbFD45nA, histVGGGenInbCD45nA, histVGGGenInbCDFT45nA, histVGGGenInbCR45nA = 0, 0, 0, 0, 0
			histVGGGenInb45nA_plus, histVGGGenInbFD45nA_plus, histVGGGenInbCD45nA_plus, histVGGGenInbCDFT45nA_plus, histVGGGenInbCR45nA_plus = 0, 0, 0, 0, 0
			histVGGGenInb45nA_minus, histVGGGenInbFD45nA_minus, histVGGGenInbCD45nA_minus, histVGGGenInbCDFT45nA_minus, histVGGGenInbCR45nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_inb_vgg45nA:
				histVGGGenInb45nA = histVGGGenInb45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD45nA = histVGGGenInbFD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD45nA = histVGGGenInbCD45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT45nA = histVGGGenInbCDFT45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
				histVGGGenInbCR45nA = histVGGGenInbCR45nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]
				histVGGGenInb45nA_plus = histVGGGenInb45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genplus.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD45nA_plus = histVGGGenInbFD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1plus.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD45nA_plus = histVGGGenInbCD45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2plus.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT45nA_plus = histVGGGenInbCDFT45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3plus.npz".format(k, jobNum))["hist"]
				histVGGGenInbCR45nA_plus = histVGGGenInbCR45nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4plus.npz".format(k, jobNum))["hist"]
				histVGGGenInb45nA_minus = histVGGGenInb45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genminus.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD45nA_minus = histVGGGenInbFD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1minus.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD45nA_minus = histVGGGenInbCD45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2minus.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT45nA_minus = histVGGGenInbCDFT45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3minus.npz".format(k, jobNum))["hist"]
				histVGGGenInbCR45nA_minus = histVGGGenInbCR45nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4minus.npz".format(k, jobNum))["hist"]

			#outbending
			#exp -all
			histExpOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)])
			histBHDVCSOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)])
			histBHDVCSOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)])
			histBHDVCSOutbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "cont{}".format(i)])

			histBkgUncOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "unccont{}".format(i)])
			histBkgUncOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "unccont{}".format(i)])
			histBkgUncOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "unccont{}".format(i)])
			histBkgUncOutbCR, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 4), "unccont{}".format(i)])
			histExpUncOutbFD = np.sqrt(histExpOutbFD)
			histExpUncOutbCD = np.sqrt(histExpOutbCD)
			histExpUncOutbCDFT = np.sqrt(histExpOutbCDFT)
			histExpUncOutbCR = np.sqrt(histExpOutbCR)

			#exp-positive helicity
			histExpOutbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSOutbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)])
			histBHDVCSOutbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)])
			histBHDVCSOutbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)])
			histBHDVCSOutbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 4), "cont{}".format(i)])

			histBkgUncOutbFD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 1), "unccont{}".format(i)])
			histBkgUncOutbCD_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 2), "unccont{}".format(i)])
			histBkgUncOutbCDFT_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 3), "unccont{}".format(i)])
			histBkgUncOutbCR_plus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == -1) & (epgExp.polarity == 1) & (epgExp.config == 4), "unccont{}".format(i)])
			histExpUncOutbFD_plus = np.sqrt(histExpOutbFD_plus)
			histExpUncOutbCD_plus = np.sqrt(histExpOutbCD_plus)
			histExpUncOutbCDFT_plus = np.sqrt(histExpOutbCDFT_plus)
			histExpUncOutbCR_plus = np.sqrt(histExpOutbCR_plus)

			#exp-negative helicity
			histExpOutbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSOutbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)])
			histBHDVCSOutbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)])
			histBHDVCSOutbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)])
			histBHDVCSOutbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 4), "cont{}".format(i)])

			histBkgUncOutbFD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 1), "unccont{}".format(i)])
			histBkgUncOutbCD_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 2), "unccont{}".format(i)])
			histBkgUncOutbCDFT_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 3), "unccont{}".format(i)])
			histBkgUncOutbCR_minus, bins = np.histogramdd(epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 4) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.helicity == 1) & (epgExp.polarity == 1) & (epgExp.config == 4), "unccont{}".format(i)])
			histExpUncOutbFD_minus = np.sqrt(histExpOutbFD_minus)
			histExpUncOutbCD_minus = np.sqrt(histExpOutbCD_minus)
			histExpUncOutbCDFT_minus = np.sqrt(histExpOutbCDFT_minus)
			histExpUncOutbCR_minus = np.sqrt(histExpOutbCR_minus)

			histExpOutb = histExpOutbFD + histExpOutbCD + histExpOutbCDFT
			histBHDVCSOutb = histBHDVCSOutbFD + histBHDVCSOutbCD + histBHDVCSOutbCDFT

			ActiveOutb = histExpOutb>20
			ActiveOutb_int = np.stack([np.sum(ActiveOutb, axis=-1)>8]*(len(phibins)-1), axis = -1)
			ActiveAll = ActiveInb & ActiveOutb
			ActiveAny = ActiveInb | ActiveOutb
			ActiveAll_int = ActiveInb_int & ActiveOutb_int
			ActiveAny_int = ActiveInb_int | ActiveOutb_int

			print("reading bhs  - outbending")
			# bh rec
			histBHOutb50nA, histBHOutbFD50nA, histBHOutbCD50nA, histBHOutbCDFT50nA, histBHOutbCR50nA = 0, 0, 0, 0, 0
			histBHOutb50nA_plus, histBHOutbFD50nA_plus, histBHOutbCD50nA_plus, histBHOutbCDFT50nA_plus, histBHOutbCR50nA_plus = 0, 0, 0, 0, 0
			histBHOutb50nA_minus, histBHOutbFD50nA_minus, histBHOutbCD50nA_minus, histBHOutbCDFT50nA_minus, histBHOutbCR50nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_outb_bh50nA:
				histBHOutb50nA = histBHOutb50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbFD50nA = histBHOutbFD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCD50nA = histBHOutbCD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCDFT50nA = histBHOutbCDFT50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCR50nA = histBHOutbCR50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutb50nA_plus = histBHOutb50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recplus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbFD50nA_plus = histBHOutbFD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCD50nA_plus = histBHOutbCD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCDFT50nA_plus = histBHOutbCDFT50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCR50nA_plus = histBHOutbCR50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4plus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutb50nA_minus = histBHOutb50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recminus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbFD50nA_minus = histBHOutbFD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1minus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCD50nA_minus = histBHOutbCD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2minus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCDFT50nA_minus = histBHOutbCDFT50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3minus.npz".format(optionaltag, k, jobNum))["hist"]
				histBHOutbCR50nA_minus = histBHOutbCR50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4minus.npz".format(optionaltag, k, jobNum))["hist"]

			# bh gen
			histBHGenOutb50nA, histBHGenOutbFD50nA, histBHGenOutbCD50nA, histBHGenOutbCDFT50nA, histBHGenOutbCR50nA = 0, 0, 0, 0, 0
			histBHGenOutb50nA_plus, histBHGenOutbFD50nA_plus, histBHGenOutbCD50nA_plus, histBHGenOutbCDFT50nA_plus, histBHGenOutbCR50nA_plus = 0, 0, 0, 0, 0
			histBHGenOutb50nA_minus, histBHGenOutbFD50nA_minus, histBHGenOutbCD50nA_minus, histBHGenOutbCDFT50nA_minus, histBHGenOutbCR50nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_outb_bh50nA:
				histBHGenOutb50nA = histBHGenOutb50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histBHGenOutbFD50nA = histBHGenOutbFD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histBHGenOutbCD50nA = histBHGenOutbCD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histBHGenOutbCDFT50nA = histBHGenOutbCDFT50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
				histBHGenOutbCR50nA = histBHGenOutbCR50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]
				histBHGenOutb50nA_plus = histBHGenOutb50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genplus.npz".format(k, jobNum))["hist"]
				histBHGenOutbFD50nA_plus = histBHGenOutbFD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1plus.npz".format(k, jobNum))["hist"]
				histBHGenOutbCD50nA_plus = histBHGenOutbCD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2plus.npz".format(k, jobNum))["hist"]
				histBHGenOutbCDFT50nA_plus = histBHGenOutbCDFT50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3plus.npz".format(k, jobNum))["hist"]
				histBHGenOutbCR50nA_plus = histBHGenOutbCR50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4plus.npz".format(k, jobNum))["hist"]
				histBHGenOutb50nA_minus = histBHGenOutb50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genminus.npz".format(k, jobNum))["hist"]
				histBHGenOutbFD50nA_minus = histBHGenOutbFD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1minus.npz".format(k, jobNum))["hist"]
				histBHGenOutbCD50nA_minus = histBHGenOutbCD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2minus.npz".format(k, jobNum))["hist"]
				histBHGenOutbCDFT50nA_minus = histBHGenOutbCDFT50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3minus.npz".format(k, jobNum))["hist"]
				histBHGenOutbCR50nA_minus = histBHGenOutbCR50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4minus.npz".format(k, jobNum))["hist"]

			print("reading vggs  - outbending")
			# vgg rec
			histVGGOutb50nA, histVGGOutbFD50nA, histVGGOutbCD50nA, histVGGOutbCDFT50nA, histVGGOutbCR50nA = 0, 0, 0, 0, 0
			histVGGOutb50nA_plus, histVGGOutbFD50nA_plus, histVGGOutbCD50nA_plus, histVGGOutbCDFT50nA_plus, histVGGOutbCR50nA_plus = 0, 0, 0, 0, 0
			histVGGOutb50nA_minus, histVGGOutbFD50nA_minus, histVGGOutbCD50nA_minus, histVGGOutbCDFT50nA_minus, histVGGOutbCR50nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_outb_vgg50nA:
				histVGGOutb50nA = histVGGOutb50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbFD50nA = histVGGOutbFD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCD50nA = histVGGOutbCD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCDFT50nA = histVGGOutbCDFT50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCR50nA = histVGGOutbCR50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutb50nA_plus = histVGGOutb50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recplus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbFD50nA_plus = histVGGOutbFD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCD50nA_plus = histVGGOutbCD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCDFT50nA_plus = histVGGOutbCDFT50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCR50nA_plus = histVGGOutbCR50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4plus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutb50nA_minus = histVGGOutb50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Recminus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbFD50nA_minus = histVGGOutbFD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec1minus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCD50nA_minus = histVGGOutbCD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec2minus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCDFT50nA_minus = histVGGOutbCDFT50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec3minus.npz".format(optionaltag, k, jobNum))["hist"]
				histVGGOutbCR50nA_minus = histVGGOutbCR50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/{}Rec4minus.npz".format(optionaltag, k, jobNum))["hist"]

			#vgg gen
			histVGGGenOutb50nA, histVGGGenOutbFD50nA, histVGGGenOutbCD50nA, histVGGGenOutbCDFT50nA, histVGGGenOutbCR50nA = 0, 0, 0, 0, 0
			histVGGGenOutb50nA_plus, histVGGGenOutbFD50nA_plus, histVGGGenOutbCD50nA_plus, histVGGGenOutbCDFT50nA_plus, histVGGGenOutbCR50nA_plus = 0, 0, 0, 0, 0
			histVGGGenOutb50nA_minus, histVGGGenOutbFD50nA_minus, histVGGGenOutbCD50nA_minus, histVGGGenOutbCDFT50nA_minus, histVGGGenOutbCR50nA_minus = 0, 0, 0, 0, 0
			for jobNum in runs_outb_vgg50nA:
				histVGGGenOutb50nA = histVGGGenOutb50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD50nA = histVGGGenOutbFD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD50nA = histVGGGenOutbCD50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT50nA = histVGGGenOutbCDFT50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCR50nA = histVGGGenOutbCR50nA + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4.npz".format(k, jobNum))["hist"]
				histVGGGenOutb50nA_plus = histVGGGenOutb50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genplus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD50nA_plus = histVGGGenOutbFD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1plus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD50nA_plus = histVGGGenOutbCD50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2plus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT50nA_plus = histVGGGenOutbCDFT50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3plus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCR50nA_plus = histVGGGenOutbCR50nA_plus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4plus.npz".format(k, jobNum))["hist"]
				histVGGGenOutb50nA_minus = histVGGGenOutb50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genminus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD50nA_minus = histVGGGenOutbFD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen1minus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD50nA_minus = histVGGGenOutbCD50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen2minus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT50nA_minus = histVGGGenOutbCDFT50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen3minus.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCR50nA_minus = histVGGGenOutbCR50nA_minus + np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gen4minus.npz".format(k, jobNum))["hist"]

			#acc correction with stat error - inbending
			#with VGG
			accCorrectedInbFD_VGG = divideHist(histBHDVCSInbFD*histVGGGenInbFD45nA , histVGGInbFD45nA)
			accCorrectedInbCD_VGG = divideHist(histBHDVCSInbCD*histVGGGenInbCD45nA , histVGGInbCD45nA)
			accCorrectedInbCDFT_VGG = divideHist(histBHDVCSInbCDFT*histVGGGenInbCDFT45nA , histVGGInbCDFT45nA)
			accCorrectedInbCR_VGG = divideHist(histBHDVCSInbCR*histVGGGenInbCR45nA , histVGGInbCR45nA)
			accCorrectedInb_VGG = accCorrectedInbFD_VGG + accCorrectedInbCD_VGG + accCorrectedInbCDFT_VGG + accCorrectedInbCR_VGG
			#FD inbending threshold [:, :, 0, :]
			accCorrectedInb_VGG[:, :, 0, :] = accCorrectedInbCD_VGG[:, :, 0, :] + accCorrectedInbCDFT_VGG[:, :, 0, :] + accCorrectedInbCR_VGG[:, :, 0, :]

			uncStatInbFD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD**2+histExpUncInbFD**2), histBHDVCSInbFD)**2 + inverseHist(histVGGInbFD45nA) + inverseHist(histVGGGenInbFD45nA))
			uncStatInbCD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD**2+histExpUncInbCD**2), histBHDVCSInbCD)**2 + inverseHist(histVGGInbCD45nA) + inverseHist(histVGGGenInbCD45nA))
			uncStatInbCDFT_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT**2+histExpUncInbCDFT**2), histBHDVCSInbCDFT)**2 + inverseHist(histVGGInbCDFT45nA) + inverseHist(histVGGGenInbCDFT45nA))
			uncStatInbCR_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR**2+histExpUncInbCR**2), histBHDVCSInbCR)**2 + inverseHist(histVGGInbCR45nA) + inverseHist(histVGGGenInbCR45nA))
			uncStatInb_VGG = np.sqrt(accCorrectedInbFD_VGG**2 *uncStatInbFD_VGG**2 + accCorrectedInbCD_VGG**2 * uncStatInbCD_VGG**2 + accCorrectedInbCDFT_VGG**2 * uncStatInbCDFT_VGG**2 + accCorrectedInbCR_VGG**2 * uncStatInbCR_VGG**2)
			uncStatInb_VGG = divideHist(uncStatInb_VGG, accCorrectedInb_VGG)
			#FD inbending threshold [:, :, 0, :]
			uncStatInb_VGG[:, :, 0, :] = np.sqrt(accCorrectedInbCD_VGG**2 * uncStatInbCD_VGG**2 + accCorrectedInbCDFT_VGG**2 * uncStatInbCDFT_VGG**2 + accCorrectedInbCR_VGG**2 * uncStatInbCR_VGG**2)[:, :, 0, :]
			uncStatInb_VGG[:, :, 0, :] = divideHist(uncStatInb_VGG[:, :, 0, :], accCorrectedInb_VGG[:, :, 0, :])
			#Correct by the ratio of covered/uncovered
			accCorrectedInb_VGG = divideHist(accCorrectedInb_VGG*histVGGGenInb45nA, histVGGGenInbFD45nA + histVGGGenInbCD45nA + histVGGGenInbCDFT45nA + histVGGGenInbCR45nA)
			accCorrectedInb_VGG[:, :, 0, :] = divideHist(accCorrectedInb_VGG[:, :, 0, :]*histVGGGenInb45nA[:, :, 0, :], histVGGGenInbCD45nA[:, :, 0, :] + histVGGGenInbCDFT45nA[:, :, 0, :] + histVGGGenInbCR45nA[:, :, 0, :])

			# plus
			accCorrectedInbFD_VGG_plus = divideHist(histBHDVCSInbFD_plus*histVGGGenInbFD45nA_plus , histVGGInbFD45nA_plus)
			accCorrectedInbCD_VGG_plus = divideHist(histBHDVCSInbCD_plus*histVGGGenInbCD45nA_plus , histVGGInbCD45nA_plus)
			accCorrectedInbCDFT_VGG_plus = divideHist(histBHDVCSInbCDFT_plus*histVGGGenInbCDFT45nA_plus , histVGGInbCDFT45nA_plus)
			accCorrectedInbCR_VGG_plus = divideHist(histBHDVCSInbCR_plus*histVGGGenInbCR45nA_plus , histVGGInbCR45nA_plus)
			accCorrectedInb_VGG_plus = accCorrectedInbFD_VGG_plus + accCorrectedInbCD_VGG_plus + accCorrectedInbCDFT_VGG_plus + accCorrectedInbCR_VGG_plus
			#FD inbending threshold [:, :, 0, :]
			accCorrectedInb_VGG_plus[:, :, 0, :] = accCorrectedInbCD_VGG_plus[:, :, 0, :] + accCorrectedInbCDFT_VGG_plus[:, :, 0, :] + accCorrectedInbCR_VGG_plus[:, :, 0, :]

			uncStatInbFD_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD_plus**2+histExpUncInbFD_plus**2), histBHDVCSInbFD_plus)**2 + inverseHist(histVGGInbFD45nA_plus) + inverseHist(histVGGGenInbFD45nA_plus))
			uncStatInbCD_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD_plus**2+histExpUncInbCD_plus**2), histBHDVCSInbCD_plus)**2 + inverseHist(histVGGInbCD45nA_plus) + inverseHist(histVGGGenInbCD45nA_plus))
			uncStatInbCDFT_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT_plus**2+histExpUncInbCDFT_plus**2), histBHDVCSInbCDFT_plus)**2 + inverseHist(histVGGInbCDFT45nA_plus) + inverseHist(histVGGGenInbCDFT45nA_plus))
			uncStatInbCR_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR_plus**2+histExpUncInbCR_plus**2), histBHDVCSInbCR_plus)**2 + inverseHist(histVGGInbCR45nA_plus) + inverseHist(histVGGGenInbCR45nA_plus))
			uncStatInb_VGG_plus = np.sqrt(accCorrectedInbFD_VGG_plus**2 *uncStatInbFD_VGG_plus**2 + accCorrectedInbCD_VGG_plus**2 * uncStatInbCD_VGG_plus**2 + accCorrectedInbCDFT_VGG_plus**2 * uncStatInbCDFT_VGG_plus**2 + accCorrectedInbCR_VGG_plus**2 * uncStatInbCR_VGG_plus**2)
			uncStatInb_VGG_plus = divideHist(uncStatInb_VGG_plus, accCorrectedInb_VGG_plus)
			#FD inbending threshold [:, :, 0, :]
			uncStatInb_VGG_plus[:, :, 0, :] = np.sqrt(accCorrectedInbCD_VGG_plus**2 * uncStatInbCD_VGG_plus**2 + accCorrectedInbCDFT_VGG_plus**2 * uncStatInbCDFT_VGG_plus**2 + accCorrectedInbCR_VGG_plus**2 * uncStatInbCR_VGG_plus**2)[:, :, 0, :]
			uncStatInb_VGG_plus[:, :, 0, :] = divideHist(uncStatInb_VGG_plus[:, :, 0, :], accCorrectedInb_VGG_plus[:, :, 0, :])
			#Correct by the ratio of covered/uncovered
			accCorrectedInb_VGG_plus = divideHist(accCorrectedInb_VGG_plus*histVGGGenInb45nA_plus, histVGGGenInbFD45nA_plus + histVGGGenInbCD45nA_plus + histVGGGenInbCDFT45nA_plus + histVGGGenInbCR45nA_plus)
			accCorrectedInb_VGG_plus[:, :, 0, :] = divideHist(accCorrectedInb_VGG_plus[:, :, 0, :]*histVGGGenInb45nA_plus[:, :, 0, :], histVGGGenInbCD45nA_plus[:, :, 0, :] + histVGGGenInbCDFT45nA_plus[:, :, 0, :] + histVGGGenInbCR45nA_plus[:, :, 0, :])

			# minus
			accCorrectedInbFD_VGG_minus = divideHist(histBHDVCSInbFD_minus*histVGGGenInbFD45nA_minus , histVGGInbFD45nA_minus)
			accCorrectedInbCD_VGG_minus = divideHist(histBHDVCSInbCD_minus*histVGGGenInbCD45nA_minus , histVGGInbCD45nA_minus)
			accCorrectedInbCDFT_VGG_minus = divideHist(histBHDVCSInbCDFT_minus*histVGGGenInbCDFT45nA_minus , histVGGInbCDFT45nA_minus)
			accCorrectedInbCR_VGG_minus = divideHist(histBHDVCSInbCR_minus*histVGGGenInbCR45nA_minus , histVGGInbCR45nA_minus)
			accCorrectedInb_VGG_minus = accCorrectedInbFD_VGG_minus + accCorrectedInbCD_VGG_minus + accCorrectedInbCDFT_VGG_minus + accCorrectedInbCR_VGG_minus
			#FD inbending threshold [:, :, 0, :]
			accCorrectedInb_VGG_minus[:, :, 0, :] = accCorrectedInbCD_VGG_minus[:, :, 0, :] + accCorrectedInbCDFT_VGG_minus[:, :, 0, :] + accCorrectedInbCR_VGG_minus[:, :, 0, :]

			uncStatInbFD_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD_minus**2+histExpUncInbFD_minus**2), histBHDVCSInbFD_minus)**2 + inverseHist(histVGGInbFD45nA_minus) + inverseHist(histVGGGenInbFD45nA_minus))
			uncStatInbCD_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD_minus**2+histExpUncInbCD_minus**2), histBHDVCSInbCD_minus)**2 + inverseHist(histVGGInbCD45nA_minus) + inverseHist(histVGGGenInbCD45nA_minus))
			uncStatInbCDFT_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT_minus**2+histExpUncInbCDFT_minus**2), histBHDVCSInbCDFT_minus)**2 + inverseHist(histVGGInbCDFT45nA_minus) + inverseHist(histVGGGenInbCDFT45nA_minus))
			uncStatInbCR_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR_minus**2+histExpUncInbCR_minus**2), histBHDVCSInbCR_minus)**2 + inverseHist(histVGGInbCR45nA_minus) + inverseHist(histVGGGenInbCR45nA_minus))
			uncStatInb_VGG_minus = np.sqrt(accCorrectedInbFD_VGG_minus**2 *uncStatInbFD_VGG_minus**2 + accCorrectedInbCD_VGG_minus**2 * uncStatInbCD_VGG_minus**2 + accCorrectedInbCDFT_VGG_minus**2 * uncStatInbCDFT_VGG_minus**2 + accCorrectedInbCR_VGG_minus**2 * uncStatInbCR_VGG_minus**2)
			uncStatInb_VGG_minus = divideHist(uncStatInb_VGG_minus, accCorrectedInb_VGG_minus)
			#FD inbending threshold [:, :, 0, :]
			uncStatInb_VGG_minus[:, :, 0, :] = np.sqrt(accCorrectedInbCD_VGG_minus**2 * uncStatInbCD_VGG_minus**2 + accCorrectedInbCDFT_VGG_minus**2 * uncStatInbCDFT_VGG_minus**2 + accCorrectedInbCR_VGG_minus**2 * uncStatInbCR_VGG_minus**2)[:, :, 0, :]
			uncStatInb_VGG_minus[:, :, 0, :] = divideHist(uncStatInb_VGG_minus[:, :, 0, :], accCorrectedInb_VGG_minus[:, :, 0, :])
			#Correct by the ratio of covered/uncovered
			accCorrectedInb_VGG_minus = divideHist(accCorrectedInb_VGG_minus*histVGGGenInb45nA_minus, histVGGGenInbFD45nA_minus + histVGGGenInbCD45nA_minus + histVGGGenInbCDFT45nA_minus + histVGGGenInbCR45nA_minus)
			accCorrectedInb_VGG_minus[:, :, 0, :] = divideHist(accCorrectedInb_VGG_minus[:, :, 0, :]*histVGGGenInb45nA_minus[:, :, 0, :], histVGGGenInbCD45nA_minus[:, :, 0, :] + histVGGGenInbCDFT45nA_minus[:, :, 0, :] + histVGGGenInbCR45nA_minus[:, :, 0, :])

			#config specific
			accCorrectedInbFDonly_VGG = divideHist(histBHDVCSInbFD*histVGGGenInb45nA , histVGGInbFD45nA)
			accCorrectedInbCDonly_VGG = divideHist(histBHDVCSInbCD*histVGGGenInb45nA , histVGGInbCD45nA)
			accCorrectedInbCDFTonly_VGG = divideHist(histBHDVCSInbCDFT*histVGGGenInb45nA , histVGGInbCDFT45nA)
			accCorrectedInbCRonly_VGG = divideHist(histBHDVCSInbCR*histVGGGenInb45nA , histVGGInbCR45nA)

			uncStatInbFDonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD**2+histExpUncInbFD**2), histBHDVCSInbFD)**2 + inverseHist(histVGGInbFD45nA) + inverseHist(histVGGGenInb45nA))
			uncStatInbCDonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD**2+histExpUncInbCD**2), histBHDVCSInbCD)**2 + inverseHist(histVGGInbCD45nA) + inverseHist(histVGGGenInb45nA))
			uncStatInbCDFTonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT**2+histExpUncInbCDFT**2), histBHDVCSInbCDFT)**2 + inverseHist(histVGGInbCDFT45nA) + inverseHist(histVGGGenInb45nA))
			uncStatInbCRonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR**2+histExpUncInbCR**2), histBHDVCSInbCR)**2 + inverseHist(histVGGInbCR45nA) + inverseHist(histVGGGenInb45nA))

			#don't use the low Q^2
			accCorrectedInb_VGG[~ActiveInb] = 0
			uncStatInb_VGG[~ActiveInb] = 0
			accCorrectedInb_VGG_plus[~ActiveInb] = 0
			uncStatInb_VGG_plus[~ActiveInb] = 0
			accCorrectedInb_VGG_minus[~ActiveInb] = 0
			uncStatInb_VGG_minus[~ActiveInb] = 0

			# acc correction inbending, with BH
			accCorrectedInbFD_BH = divideHist(histBHDVCSInbFD*histBHGenInbFD45nA , histBHInbFD45nA)
			accCorrectedInbCD_BH = divideHist(histBHDVCSInbCD*histBHGenInbCD45nA , histBHInbCD45nA)
			accCorrectedInbCDFT_BH = divideHist(histBHDVCSInbCDFT*histBHGenInbCDFT45nA , histBHInbCDFT45nA)
			accCorrectedInbCR_BH = divideHist(histBHDVCSInbCR*histBHGenInbCR45nA , histBHInbCR45nA)
			accCorrectedInb_BH = accCorrectedInbFD_BH + accCorrectedInbCD_BH + accCorrectedInbCDFT_BH + accCorrectedInbCR_BH
			#FD inbending threshold [:, :, 0, :]
			accCorrectedInb_BH[:, :, 0, :] = accCorrectedInbCD_BH[:, :, 0, :] + accCorrectedInbCDFT_BH[:, :, 0, :] + accCorrectedInbCR_BH[:, :, 0, :]

			uncStatInbFD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD**2+histExpUncInbFD**2), histBHDVCSInbFD)**2 + inverseHist(histBHInbFD45nA) + inverseHist(histBHGenInbFD45nA))
			uncStatInbCD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD**2+histExpUncInbCD**2), histBHDVCSInbCD)**2 + inverseHist(histBHInbCD45nA) + inverseHist(histBHGenInbCD45nA))
			uncStatInbCDFT_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT**2+histExpUncInbCDFT**2), histBHDVCSInbCDFT)**2 + inverseHist(histBHInbCDFT45nA) + inverseHist(histBHGenInbCDFT45nA))
			uncStatInbCR_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR**2+histExpUncInbCR**2), histBHDVCSInbCR)**2 + inverseHist(histBHInbCR45nA) + inverseHist(histBHGenInbCR45nA))
			uncStatInb_BH = np.sqrt(accCorrectedInbFD_BH**2 *uncStatInbFD_BH**2 + accCorrectedInbCD_BH**2 * uncStatInbCD_BH**2 + accCorrectedInbCDFT_BH**2 * uncStatInbCDFT_BH**2 + accCorrectedInbCR_BH**2 * uncStatInbCR_BH**2)
			uncStatInb_BH = divideHist(uncStatInb_BH, accCorrectedInb_BH)
			#FD inbending threshold [:, :, 0, :]
			uncStatInb_BH[:, :, 0, :] = np.sqrt(accCorrectedInbCD_BH**2 * uncStatInbCD_BH**2 + accCorrectedInbCDFT_BH**2 * uncStatInbCDFT_BH**2 + accCorrectedInbCR_BH**2 * uncStatInbCR_BH**2)[:, :, 0, :]
			uncStatInb_BH[:, :, 0, :] = divideHist(uncStatInb_BH[:, :, 0, :], accCorrectedInb_BH[:, :, 0, :])
			#Correct by the ratio of covered/uncovered
			accCorrectedInb_BH = divideHist(accCorrectedInb_BH*histBHGenInb45nA, histBHGenInbFD45nA + histBHGenInbCD45nA + histBHGenInbCDFT45nA + histBHGenInbCR45nA)
			accCorrectedInb_BH[:, :, 0, :] = divideHist(accCorrectedInb_BH[:, :, 0, :]*histBHGenInb45nA[:, :, 0, :], histBHGenInbCD45nA[:, :, 0, :] + histBHGenInbCDFT45nA[:, :, 0, :] + histBHGenInbCR45nA[:, :, 0, :])

			# plus
			accCorrectedInbFD_BH_plus = divideHist(histBHDVCSInbFD_plus*histBHGenInbFD45nA_plus , histBHInbFD45nA_plus)
			accCorrectedInbCD_BH_plus = divideHist(histBHDVCSInbCD_plus*histBHGenInbCD45nA_plus , histBHInbCD45nA_plus)
			accCorrectedInbCDFT_BH_plus = divideHist(histBHDVCSInbCDFT_plus*histBHGenInbCDFT45nA_plus , histBHInbCDFT45nA_plus)
			accCorrectedInbCR_BH_plus = divideHist(histBHDVCSInbCR_plus*histBHGenInbCR45nA_plus , histBHInbCR45nA_plus)
			accCorrectedInb_BH_plus = accCorrectedInbFD_BH_plus + accCorrectedInbCD_BH_plus + accCorrectedInbCDFT_BH_plus + accCorrectedInbCR_BH_plus
			#FD inbending threshold [:, :, 0, :]
			accCorrectedInb_BH_plus[:, :, 0, :] = accCorrectedInbCD_BH_plus[:, :, 0, :] + accCorrectedInbCDFT_BH_plus[:, :, 0, :] + accCorrectedInbCR_BH_plus[:, :, 0, :]

			uncStatInbFD_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD_plus**2+histExpUncInbFD_plus**2), histBHDVCSInbFD_plus)**2 + inverseHist(histBHInbFD45nA_plus) + inverseHist(histBHGenInbFD45nA_plus))
			uncStatInbCD_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD_plus**2+histExpUncInbCD_plus**2), histBHDVCSInbCD_plus)**2 + inverseHist(histBHInbCD45nA_plus) + inverseHist(histBHGenInbCD45nA_plus))
			uncStatInbCDFT_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT_plus**2+histExpUncInbCDFT_plus**2), histBHDVCSInbCDFT_plus)**2 + inverseHist(histBHInbCDFT45nA_plus) + inverseHist(histBHGenInbCDFT45nA_plus))
			uncStatInbCR_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR_plus**2+histExpUncInbCR_plus**2), histBHDVCSInbCR_plus)**2 + inverseHist(histBHInbCR45nA_plus) + inverseHist(histBHGenInbCR45nA_plus))
			uncStatInb_BH_plus = np.sqrt(accCorrectedInbFD_BH_plus**2 *uncStatInbFD_BH_plus**2 + accCorrectedInbCD_BH_plus**2 * uncStatInbCD_BH_plus**2 + accCorrectedInbCDFT_BH_plus**2 * uncStatInbCDFT_BH_plus**2 + accCorrectedInbCR_BH_plus**2 * uncStatInbCR_BH_plus**2)
			uncStatInb_BH_plus = divideHist(uncStatInb_BH_plus, accCorrectedInb_BH_plus)
			#FD inbending threshold [:, :, 0, :]
			uncStatInb_BH_plus[:, :, 0, :] = np.sqrt(accCorrectedInbCD_BH_plus**2 * uncStatInbCD_BH_plus**2 + accCorrectedInbCDFT_BH_plus**2 * uncStatInbCDFT_BH_plus**2 + accCorrectedInbCR_BH_plus**2 * uncStatInbCR_BH_plus**2)[:, :, 0, :]
			uncStatInb_BH_plus[:, :, 0, :] = divideHist(uncStatInb_BH_plus[:, :, 0, :], accCorrectedInb_BH_plus[:, :, 0, :])
			#Correct by the ratio of covered/uncovered
			accCorrectedInb_BH_plus = divideHist(accCorrectedInb_BH_plus*histBHGenInb45nA_plus, histBHGenInbFD45nA_plus + histBHGenInbCD45nA_plus + histBHGenInbCDFT45nA_plus + histBHGenInbCR45nA_plus)
			accCorrectedInb_BH_plus[:, :, 0, :] = divideHist(accCorrectedInb_BH_plus[:, :, 0, :]*histBHGenInb45nA_plus[:, :, 0, :], histBHGenInbCD45nA_plus[:, :, 0, :] + histBHGenInbCDFT45nA_plus[:, :, 0, :] + histBHGenInbCR45nA_plus[:, :, 0, :])

			# minus
			accCorrectedInbFD_BH_minus = divideHist(histBHDVCSInbFD_minus*histBHGenInbFD45nA_minus , histBHInbFD45nA_minus)
			accCorrectedInbCD_BH_minus = divideHist(histBHDVCSInbCD_minus*histBHGenInbCD45nA_minus , histBHInbCD45nA_minus)
			accCorrectedInbCDFT_BH_minus = divideHist(histBHDVCSInbCDFT_minus*histBHGenInbCDFT45nA_minus , histBHInbCDFT45nA_minus)
			accCorrectedInbCR_BH_minus = divideHist(histBHDVCSInbCR_minus*histBHGenInbCR45nA_minus , histBHInbCR45nA_minus)
			accCorrectedInb_BH_minus = accCorrectedInbFD_BH_minus + accCorrectedInbCD_BH_minus + accCorrectedInbCDFT_BH_minus + accCorrectedInbCR_BH_minus
			#FD inbending threshold [:, :, 0, :]
			accCorrectedInb_BH_minus[:, :, 0, :] = accCorrectedInbCD_BH_minus[:, :, 0, :] + accCorrectedInbCDFT_BH_minus[:, :, 0, :] + accCorrectedInbCR_BH_minus[:, :, 0, :]

			uncStatInbFD_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD_minus**2+histExpUncInbFD_minus**2), histBHDVCSInbFD_minus)**2 + inverseHist(histBHInbFD45nA_minus) + inverseHist(histBHGenInbFD45nA_minus))
			uncStatInbCD_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD_minus**2+histExpUncInbCD_minus**2), histBHDVCSInbCD_minus)**2 + inverseHist(histBHInbCD45nA_minus) + inverseHist(histBHGenInbCD45nA_minus))
			uncStatInbCDFT_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT_minus**2+histExpUncInbCDFT_minus**2), histBHDVCSInbCDFT_minus)**2 + inverseHist(histBHInbCDFT45nA_minus) + inverseHist(histBHGenInbCDFT45nA_minus))
			uncStatInbCR_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR_minus**2+histExpUncInbCR_minus**2), histBHDVCSInbCR_minus)**2 + inverseHist(histBHInbCR45nA_minus) + inverseHist(histBHGenInbCR45nA_minus))
			uncStatInb_BH_minus = np.sqrt(accCorrectedInbFD_BH_minus**2 *uncStatInbFD_BH_minus**2 + accCorrectedInbCD_BH_minus**2 * uncStatInbCD_BH_minus**2 + accCorrectedInbCDFT_BH_minus**2 * uncStatInbCDFT_BH_minus**2 + accCorrectedInbCR_BH_minus**2 * uncStatInbCR_BH_minus**2)
			uncStatInb_BH_minus = divideHist(uncStatInb_BH_minus, accCorrectedInb_BH_minus)
			#FD inbending threshold [:, :, 0, :]
			uncStatInb_BH_minus[:, :, 0, :] = np.sqrt(accCorrectedInbCD_BH_minus**2 * uncStatInbCD_BH_minus**2 + accCorrectedInbCDFT_BH_minus**2 * uncStatInbCDFT_BH_minus**2 + accCorrectedInbCR_BH_minus**2 * uncStatInbCR_BH_minus**2)[:, :, 0, :]
			uncStatInb_BH_minus[:, :, 0, :] = divideHist(uncStatInb_BH_minus[:, :, 0, :], accCorrectedInb_BH_minus[:, :, 0, :])
			#Correct by the ratio of covered/uncovered
			accCorrectedInb_BH_minus = divideHist(accCorrectedInb_BH_minus*histBHGenInb45nA_minus, histBHGenInbFD45nA_minus + histBHGenInbCD45nA_minus + histBHGenInbCDFT45nA_minus + histBHGenInbCR45nA_minus)
			accCorrectedInb_BH_minus[:, :, 0, :] = divideHist(accCorrectedInb_BH_minus[:, :, 0, :]*histBHGenInb45nA_minus[:, :, 0, :], histBHGenInbCD45nA_minus[:, :, 0, :] + histBHGenInbCDFT45nA_minus[:, :, 0, :] + histBHGenInbCR45nA_minus[:, :, 0, :])

			#config specific
			accCorrectedInbFDonly_BH = divideHist(histBHDVCSInbFD*histBHGenInb45nA , histBHInbFD45nA)
			accCorrectedInbCDonly_BH = divideHist(histBHDVCSInbCD*histBHGenInb45nA , histBHInbCD45nA)
			accCorrectedInbCDFTonly_BH = divideHist(histBHDVCSInbCDFT*histBHGenInb45nA , histBHInbCDFT45nA)
			accCorrectedInbCRonly_BH = divideHist(histBHDVCSInbCR*histBHGenInb45nA , histBHInbCR45nA)

			uncStatInbFDonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD**2+histExpUncInbFD**2), histBHDVCSInbFD)**2 + inverseHist(histBHInbFD45nA) + inverseHist(histBHGenInb45nA))
			uncStatInbCDonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD**2+histExpUncInbCD**2), histBHDVCSInbCD)**2 + inverseHist(histBHInbCD45nA) + inverseHist(histBHGenInb45nA))
			uncStatInbCDFTonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT**2+histExpUncInbCDFT**2), histBHDVCSInbCDFT)**2 + inverseHist(histBHInbCDFT45nA) + inverseHist(histBHGenInb45nA))
			uncStatInbCRonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCR**2+histExpUncInbCR**2), histBHDVCSInbCR)**2 + inverseHist(histBHInbCR45nA) + inverseHist(histBHGenInb45nA))

			accCorrectedInb_BH[~ActiveInb] = 0
			uncStatInb_BH[~ActiveInb] = 0
			accCorrectedInb_BH_plus[~ActiveInb] = 0
			uncStatInb_BH_plus[~ActiveInb] = 0
			accCorrectedInb_BH_minus[~ActiveInb] = 0
			uncStatInb_BH_minus[~ActiveInb] = 0

			# acceptance correction with stat error - outbending
			# with VGG
			accCorrectedOutbFD_VGG = divideHist(histBHDVCSOutbFD*histVGGGenOutbFD50nA , histVGGOutbFD50nA)
			accCorrectedOutbCD_VGG = divideHist(histBHDVCSOutbCD*histVGGGenOutbCD50nA , histVGGOutbCD50nA)
			accCorrectedOutbCDFT_VGG = divideHist(histBHDVCSOutbCDFT*histVGGGenOutbCDFT50nA , histVGGOutbCDFT50nA)
			accCorrectedOutbCR_VGG = divideHist(histBHDVCSOutbCR*histVGGGenOutbCR50nA , histVGGOutbCR50nA)
			accCorrectedOutb_VGG = accCorrectedOutbFD_VGG + accCorrectedOutbCD_VGG + accCorrectedOutbCDFT_VGG + accCorrectedOutbCR_VGG
			#FD inbending threshold [:, :, [0, 1], :]
			accCorrectedOutb_VGG[:, :, [0, 1], :] = accCorrectedOutbCD_VGG[:, :, [0, 1], :] + accCorrectedOutbCDFT_VGG[:, :, [0, 1], :] + accCorrectedOutbCR_VGG[:, :, [0, 1], :]

			uncStatOutbFD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD**2+histExpUncOutbFD**2), histBHDVCSOutbFD)**2 + inverseHist(histVGGOutbFD50nA) + inverseHist(histVGGGenOutbFD50nA))
			uncStatOutbCD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD**2+histExpUncOutbCD**2), histBHDVCSOutbCD)**2 + inverseHist(histVGGOutbCD50nA) + inverseHist(histVGGGenOutbCD50nA))
			uncStatOutbCDFT_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT**2+histExpUncOutbCDFT**2), histBHDVCSOutbCDFT)**2 + inverseHist(histVGGOutbCDFT50nA) + inverseHist(histVGGGenOutbCDFT50nA))
			uncStatOutbCR_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR**2+histExpUncOutbCR**2), histBHDVCSOutbCR)**2 + inverseHist(histVGGOutbCR50nA) + inverseHist(histVGGGenOutbCR50nA))
			uncStatOutb_VGG = np.sqrt(accCorrectedOutbFD_VGG**2 *uncStatOutbFD_VGG**2 + accCorrectedOutbCD_VGG**2 * uncStatOutbCD_VGG**2 + accCorrectedOutbCDFT_VGG**2 * uncStatOutbCDFT_VGG**2 + accCorrectedOutbCR_VGG**2 * uncStatOutbCR_VGG**2)
			uncStatOutb_VGG = divideHist(uncStatOutb_VGG, accCorrectedOutb_VGG)
			#FD inbending threshold [:, :, [0, 1], :]
			uncStatOutb_VGG[:, :, [0, 1], :] = np.sqrt(accCorrectedOutbCD_VGG**2 * uncStatOutbCD_VGG**2 + accCorrectedOutbCDFT_VGG**2 * uncStatOutbCDFT_VGG**2 + accCorrectedOutbCR_VGG**2 * uncStatOutbCR_VGG**2)[:, :, [0, 1], :]
			uncStatOutb_VGG[:, :, [0, 1], :] = divideHist(uncStatOutb_VGG[:, :, [0, 1], :], accCorrectedOutb_VGG[:, :, [0, 1], :])
			#Correct by the ratio of covered/uncovered
			accCorrectedOutb_VGG = divideHist(accCorrectedOutb_VGG*histVGGGenOutb50nA, histVGGGenOutbFD50nA + histVGGGenOutbCD50nA + histVGGGenOutbCDFT50nA + histVGGGenOutbCR50nA)
			accCorrectedOutb_VGG[:, :, [0, 1], :] = divideHist(accCorrectedOutb_VGG[:, :, [0, 1], :]*histVGGGenOutb50nA[:, :, [0, 1], :], histVGGGenOutbCD50nA[:, :, [0, 1], :] + histVGGGenOutbCDFT50nA[:, :, [0, 1], :] + histVGGGenOutbCR50nA[:, :, [0, 1], :])

			# plus
			accCorrectedOutbFD_VGG_plus = divideHist(histBHDVCSOutbFD_plus*histVGGGenOutbFD50nA_plus , histVGGOutbFD50nA_plus)
			accCorrectedOutbCD_VGG_plus = divideHist(histBHDVCSOutbCD_plus*histVGGGenOutbCD50nA_plus , histVGGOutbCD50nA_plus)
			accCorrectedOutbCDFT_VGG_plus = divideHist(histBHDVCSOutbCDFT_plus*histVGGGenOutbCDFT50nA_plus , histVGGOutbCDFT50nA_plus)
			accCorrectedOutbCR_VGG_plus = divideHist(histBHDVCSOutbCR_plus*histVGGGenOutbCR50nA_plus , histVGGOutbCR50nA_plus)
			accCorrectedOutb_VGG_plus = accCorrectedOutbFD_VGG_plus + accCorrectedOutbCD_VGG_plus + accCorrectedOutbCDFT_VGG_plus + accCorrectedOutbCR_VGG_plus
			#FD inbending threshold [:, :, [0, 1], :]
			accCorrectedOutb_VGG_plus[:, :, [0, 1], :] = accCorrectedOutbCD_VGG_plus[:, :, [0, 1], :] + accCorrectedOutbCDFT_VGG_plus[:, :, [0, 1], :] + accCorrectedOutbCR_VGG_plus[:, :, [0, 1], :]

			uncStatOutbFD_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD_plus**2+histExpUncOutbFD_plus**2), histBHDVCSOutbFD_plus)**2 + inverseHist(histVGGOutbFD50nA_plus) + inverseHist(histVGGGenOutbFD50nA_plus))
			uncStatOutbCD_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD_plus**2+histExpUncOutbCD_plus**2), histBHDVCSOutbCD_plus)**2 + inverseHist(histVGGOutbCD50nA_plus) + inverseHist(histVGGGenOutbCD50nA_plus))
			uncStatOutbCDFT_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT_plus**2+histExpUncOutbCDFT_plus**2), histBHDVCSOutbCDFT_plus)**2 + inverseHist(histVGGOutbCDFT50nA_plus) + inverseHist(histVGGGenOutbCDFT50nA_plus))
			uncStatOutbCR_VGG_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR_plus**2+histExpUncOutbCR_plus**2), histBHDVCSOutbCR_plus)**2 + inverseHist(histVGGOutbCR50nA_plus) + inverseHist(histVGGGenOutbCR50nA_plus))
			uncStatOutb_VGG_plus = np.sqrt(accCorrectedOutbFD_VGG_plus**2 *uncStatOutbFD_VGG_plus**2 + accCorrectedOutbCD_VGG_plus**2 * uncStatOutbCD_VGG_plus**2 + accCorrectedOutbCDFT_VGG_plus**2 * uncStatOutbCDFT_VGG_plus**2 + accCorrectedOutbCR_VGG_plus**2 * uncStatOutbCR_VGG_plus**2)
			uncStatOutb_VGG_plus = divideHist(uncStatOutb_VGG_plus, accCorrectedOutb_VGG_plus)
			#FD inbending threshold [:, :, [0, 1], :]
			uncStatOutb_VGG_plus[:, :, [0, 1], :] = np.sqrt(accCorrectedOutbCD_VGG_plus**2 * uncStatOutbCD_VGG_plus**2 + accCorrectedOutbCDFT_VGG_plus**2 * uncStatOutbCDFT_VGG_plus**2 + accCorrectedOutbCR_VGG_plus**2 * uncStatOutbCR_VGG_plus**2)[:, :, [0, 1], :]
			uncStatOutb_VGG_plus[:, :, [0, 1], :] = divideHist(uncStatOutb_VGG_plus[:, :, [0, 1], :], accCorrectedOutb_VGG_plus[:, :, [0, 1], :])

			accCorrectedOutb_VGG_plus = divideHist(accCorrectedOutb_VGG_plus*histVGGGenOutb50nA_plus, histVGGGenOutbFD50nA_plus + histVGGGenOutbCD50nA_plus + histVGGGenOutbCDFT50nA_plus + histVGGGenOutbCR50nA_plus)
			accCorrectedOutb_VGG_plus[:, :, [0, 1], :] = divideHist(accCorrectedOutb_VGG_plus[:, :, [0, 1], :]*histVGGGenOutb50nA_plus[:, :, [0, 1], :], histVGGGenOutbCD50nA_plus[:, :, [0, 1], :] + histVGGGenOutbCDFT50nA_plus[:, :, [0, 1], :] + histVGGGenOutbCR50nA_plus[:, :, [0, 1], :])

			# minus
			accCorrectedOutbFD_VGG_minus = divideHist(histBHDVCSOutbFD_minus*histVGGGenOutbFD50nA_minus , histVGGOutbFD50nA_minus)
			accCorrectedOutbCD_VGG_minus = divideHist(histBHDVCSOutbCD_minus*histVGGGenOutbCD50nA_minus , histVGGOutbCD50nA_minus)
			accCorrectedOutbCDFT_VGG_minus = divideHist(histBHDVCSOutbCDFT_minus*histVGGGenOutbCDFT50nA_minus , histVGGOutbCDFT50nA_minus)
			accCorrectedOutbCR_VGG_minus = divideHist(histBHDVCSOutbCR_minus*histVGGGenOutbCR50nA_minus , histVGGOutbCR50nA_minus)
			accCorrectedOutb_VGG_minus = accCorrectedOutbFD_VGG_minus + accCorrectedOutbCD_VGG_minus + accCorrectedOutbCDFT_VGG_minus + accCorrectedOutbCR_VGG_minus
			#FD inbending threshold [:, :, [0, 1], :]
			accCorrectedOutb_VGG_minus[:, :, [0, 1], :] = accCorrectedOutbCD_VGG_minus[:, :, [0, 1], :] + accCorrectedOutbCDFT_VGG_minus[:, :, [0, 1], :] + accCorrectedOutbCR_VGG_minus[:, :, [0, 1], :]

			uncStatOutbFD_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD_minus**2+histExpUncOutbFD_minus**2), histBHDVCSOutbFD_minus)**2 + inverseHist(histVGGOutbFD50nA_minus) + inverseHist(histVGGGenOutbFD50nA_minus))
			uncStatOutbCD_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD_minus**2+histExpUncOutbCD_minus**2), histBHDVCSOutbCD_minus)**2 + inverseHist(histVGGOutbCD50nA_minus) + inverseHist(histVGGGenOutbCD50nA_minus))
			uncStatOutbCDFT_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT_minus**2+histExpUncOutbCDFT_minus**2), histBHDVCSOutbCDFT_minus)**2 + inverseHist(histVGGOutbCDFT50nA_minus) + inverseHist(histVGGGenOutbCDFT50nA_minus))
			uncStatOutbCR_VGG_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR_minus**2+histExpUncOutbCR_minus**2), histBHDVCSOutbCR_minus)**2 + inverseHist(histVGGOutbCR50nA_minus) + inverseHist(histVGGGenOutbCR50nA_minus))
			uncStatOutb_VGG_minus = np.sqrt(accCorrectedOutbFD_VGG_minus**2 *uncStatOutbFD_VGG_minus**2 + accCorrectedOutbCD_VGG_minus**2 * uncStatOutbCD_VGG_minus**2 + accCorrectedOutbCDFT_VGG_minus**2 * uncStatOutbCDFT_VGG_minus**2 + accCorrectedOutbCR_VGG_minus**2 * uncStatOutbCR_VGG_minus**2)
			uncStatOutb_VGG_minus = divideHist(uncStatOutb_VGG_minus, accCorrectedOutb_VGG_minus)
			#FD inbending threshold [:, :, [0, 1], :]
			uncStatOutb_VGG_minus[:, :, [0, 1], :] = np.sqrt(accCorrectedOutbCD_VGG_minus**2 * uncStatOutbCD_VGG_minus**2 + accCorrectedOutbCDFT_VGG_minus**2 * uncStatOutbCDFT_VGG_minus**2 + accCorrectedOutbCR_VGG_minus**2 * uncStatOutbCR_VGG_minus**2)[:, :, [0, 1], :]
			uncStatOutb_VGG_minus[:, :, [0, 1], :] = divideHist(uncStatOutb_VGG_minus[:, :, [0, 1], :], accCorrectedOutb_VGG_minus[:, :, [0, 1], :])

			accCorrectedOutb_VGG_minus = divideHist(accCorrectedOutb_VGG_minus*histVGGGenOutb50nA_minus, histVGGGenOutbFD50nA_minus + histVGGGenOutbCD50nA_minus + histVGGGenOutbCDFT50nA_minus + histVGGGenOutbCR50nA_minus)
			accCorrectedOutb_VGG_minus[:, :, [0, 1], :] = divideHist(accCorrectedOutb_VGG_minus[:, :, [0, 1], :]*histVGGGenOutb50nA_minus[:, :, [0, 1], :], histVGGGenOutbCD50nA_minus[:, :, [0, 1], :] + histVGGGenOutbCDFT50nA_minus[:, :, [0, 1], :] + histVGGGenOutbCR50nA_minus[:, :, [0, 1], :])

			#config specific
			accCorrectedOutbFDonly_VGG = divideHist(histBHDVCSOutbFD*histVGGGenOutb50nA , histVGGOutbFD50nA)
			accCorrectedOutbCDonly_VGG = divideHist(histBHDVCSOutbCD*histVGGGenOutb50nA , histVGGOutbCD50nA)
			accCorrectedOutbCDFTonly_VGG = divideHist(histBHDVCSOutbCDFT*histVGGGenOutb50nA , histVGGOutbCDFT50nA)
			accCorrectedOutbCRonly_VGG = divideHist(histBHDVCSOutbCR*histVGGGenOutb50nA , histVGGOutbCR50nA)

			uncStatOutbFDonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD**2+histExpUncOutbFD**2), histBHDVCSOutbFD)**2 + inverseHist(histVGGOutbFD50nA) + inverseHist(histVGGGenOutb50nA))
			uncStatOutbCDonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD**2+histExpUncOutbCD**2), histBHDVCSOutbCD)**2 + inverseHist(histVGGOutbCD50nA) + inverseHist(histVGGGenOutb50nA))
			uncStatOutbCDFTonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT**2+histExpUncOutbCDFT**2), histBHDVCSOutbCDFT)**2 + inverseHist(histVGGOutbCDFT50nA) + inverseHist(histVGGGenOutb50nA))
			uncStatOutbCRonly_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR**2+histExpUncOutbCR**2), histBHDVCSOutbCR)**2 + inverseHist(histVGGOutbCR50nA) + inverseHist(histVGGGenOutb50nA))

			#remove inactive
			accCorrectedOutb_VGG[~ActiveOutb] = 0
			uncStatOutb_VGG[~ActiveOutb] = 0
			accCorrectedOutb_VGG_plus[~ActiveOutb] = 0
			uncStatOutb_VGG_plus[~ActiveOutb] = 0
			accCorrectedOutb_VGG_minus[~ActiveOutb] = 0
			uncStatOutb_VGG_minus[~ActiveInb] = 0

			# acc correction with BH, outbending
			accCorrectedOutbFD_BH = divideHist(histBHDVCSOutbFD*histBHGenOutbFD50nA , histBHOutbFD50nA)
			accCorrectedOutbCD_BH = divideHist(histBHDVCSOutbCD*histBHGenOutbCD50nA , histBHOutbCD50nA)
			accCorrectedOutbCDFT_BH = divideHist(histBHDVCSOutbCDFT*histBHGenOutbCDFT50nA , histBHOutbCDFT50nA)
			accCorrectedOutbCR_BH = divideHist(histBHDVCSOutbCR*histBHGenOutbCR50nA , histBHOutbCR50nA)
			accCorrectedOutb_BH = accCorrectedOutbFD_BH + accCorrectedOutbCD_BH + accCorrectedOutbCDFT_BH + accCorrectedOutbCR_BH
			#FD inbending threshold [:, :, [0, 1], :]
			accCorrectedOutb_BH[:, :, [0, 1], :] = accCorrectedOutbCD_BH[:, :, [0, 1], :] + accCorrectedOutbCDFT_BH[:, :, [0, 1], :] + accCorrectedOutbCR_BH[:, :, [0, 1], :]

			uncStatOutbFD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD**2+histExpUncOutbFD**2), histBHDVCSOutbFD)**2 + inverseHist(histBHOutbFD50nA) + inverseHist(histBHGenOutbFD50nA))
			uncStatOutbCD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD**2+histExpUncOutbCD**2), histBHDVCSOutbCD)**2 + inverseHist(histBHOutbCD50nA) + inverseHist(histBHGenOutbCD50nA))
			uncStatOutbCDFT_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT**2+histExpUncOutbCDFT**2), histBHDVCSOutbCDFT)**2 + inverseHist(histBHOutbCDFT50nA) + inverseHist(histBHGenOutbCDFT50nA))
			uncStatOutbCR_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR**2+histExpUncOutbCR**2), histBHDVCSOutbCR)**2 + inverseHist(histBHOutbCR50nA) + inverseHist(histBHGenOutbCR50nA))
			uncStatOutb_BH = np.sqrt(accCorrectedOutbFD_BH**2 *uncStatOutbFD_BH**2 + accCorrectedOutbCD_BH**2 * uncStatOutbCD_BH**2 + accCorrectedOutbCDFT_BH**2 * uncStatOutbCDFT_BH**2 + accCorrectedOutbCR_BH**2 * uncStatOutbCR_BH**2)
			uncStatOutb_BH = divideHist(uncStatOutb_BH, accCorrectedOutb_BH)
			#FD inbending threshold [:, :, [0, 1], :]
			uncStatOutb_BH[:, :, [0, 1], :] = np.sqrt(accCorrectedOutbCD_BH**2 * uncStatOutbCD_BH**2 + accCorrectedOutbCDFT_BH**2 * uncStatOutbCDFT_BH**2 + accCorrectedOutbCR_BH**2 * uncStatOutbCR_BH**2)[:, :, [0, 1], :]
			uncStatOutb_BH[:, :, [0, 1], :] = divideHist(uncStatOutb_BH[:, :, [0, 1], :], accCorrectedOutb_BH[:, :, [0, 1], :])
			#Correct by the ratio of covered/uncovered
			accCorrectedOutb_BH = divideHist(accCorrectedOutb_BH*histBHGenOutb50nA, histBHGenOutbFD50nA + histBHGenOutbCD50nA + histBHGenOutbCDFT50nA + histBHGenOutbCR50nA)
			accCorrectedOutb_BH[:, :, [0, 1], :] = divideHist(accCorrectedOutb_BH[:, :, [0, 1], :]*histBHGenOutb50nA[:, :, [0, 1], :], histBHGenOutbCD50nA[:, :, [0, 1], :] + histBHGenOutbCDFT50nA[:, :, [0, 1], :] + histBHGenOutbCR50nA[:, :, [0, 1], :])

			# plus
			accCorrectedOutbFD_BH_plus = divideHist(histBHDVCSOutbFD_plus*histBHGenOutbFD50nA_plus , histBHOutbFD50nA_plus)
			accCorrectedOutbCD_BH_plus = divideHist(histBHDVCSOutbCD_plus*histBHGenOutbCD50nA_plus , histBHOutbCD50nA_plus)
			accCorrectedOutbCDFT_BH_plus = divideHist(histBHDVCSOutbCDFT_plus*histBHGenOutbCDFT50nA_plus , histBHOutbCDFT50nA_plus)
			accCorrectedOutbCR_BH_plus = divideHist(histBHDVCSOutbCR_plus*histBHGenOutbCR50nA_plus , histBHOutbCR50nA_plus)
			accCorrectedOutb_BH_plus = accCorrectedOutbFD_BH_plus + accCorrectedOutbCD_BH_plus + accCorrectedOutbCDFT_BH_plus + accCorrectedOutbCR_BH_plus
			#FD inbending threshold [:, :, [0, 1], :]
			accCorrectedOutb_BH_plus[:, :, [0, 1], :] = accCorrectedOutbCD_BH_plus[:, :, [0, 1], :] + accCorrectedOutbCDFT_BH_plus[:, :, [0, 1], :] + accCorrectedOutbCR_BH_plus[:, :, [0, 1], :]

			uncStatOutbFD_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD_plus**2+histExpUncOutbFD_plus**2), histBHDVCSOutbFD_plus)**2 + inverseHist(histBHOutbFD50nA_plus) + inverseHist(histBHGenOutbFD50nA_plus))
			uncStatOutbCD_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD_plus**2+histExpUncOutbCD_plus**2), histBHDVCSOutbCD_plus)**2 + inverseHist(histBHOutbCD50nA_plus) + inverseHist(histBHGenOutbCD50nA_plus))
			uncStatOutbCDFT_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT_plus**2+histExpUncOutbCDFT_plus**2), histBHDVCSOutbCDFT_plus)**2 + inverseHist(histBHOutbCDFT50nA_plus) + inverseHist(histBHGenOutbCDFT50nA_plus))
			uncStatOutbCR_BH_plus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR_plus**2+histExpUncOutbCR_plus**2), histBHDVCSOutbCR_plus)**2 + inverseHist(histBHOutbCR50nA_plus) + inverseHist(histBHGenOutbCR50nA_plus))
			uncStatOutb_BH_plus = np.sqrt(accCorrectedOutbFD_BH_plus**2 *uncStatOutbFD_BH_plus**2 + accCorrectedOutbCD_BH_plus**2 * uncStatOutbCD_BH_plus**2 + accCorrectedOutbCDFT_BH_plus**2 * uncStatOutbCDFT_BH_plus**2 + accCorrectedOutbCR_BH_plus**2 * uncStatOutbCR_BH_plus**2)
			uncStatOutb_BH_plus = divideHist(uncStatOutb_BH_plus, accCorrectedOutb_BH_plus)
			#FD inbending threshold [:, :, [0, 1], :]
			uncStatOutb_BH_plus[:, :, [0, 1], :] = np.sqrt(accCorrectedOutbCD_BH_plus**2 * uncStatOutbCD_BH_plus**2 + accCorrectedOutbCDFT_BH_plus**2 * uncStatOutbCDFT_BH_plus**2 + accCorrectedOutbCR_BH_plus**2 * uncStatOutbCR_BH_plus**2)[:, :, [0, 1], :]
			uncStatOutb_BH_plus[:, :, [0, 1], :] = divideHist(uncStatOutb_BH_plus[:, :, [0, 1], :], accCorrectedOutb_BH_plus[:, :, [0, 1], :])

			accCorrectedOutb_BH_plus = divideHist(accCorrectedOutb_BH_plus*histBHGenOutb50nA_plus, histBHGenOutbFD50nA_plus + histBHGenOutbCD50nA_plus + histBHGenOutbCDFT50nA_plus + histBHGenOutbCR50nA_plus)
			accCorrectedOutb_BH_plus[:, :, [0, 1], :] = divideHist(accCorrectedOutb_BH_plus[:, :, [0, 1], :]*histBHGenOutb50nA_plus[:, :, [0, 1], :], histBHGenOutbCD50nA_plus[:, :, [0, 1], :] + histBHGenOutbCDFT50nA_plus[:, :, [0, 1], :] + histBHGenOutbCR50nA_plus[:, :, [0, 1], :])

			# minus
			accCorrectedOutbFD_BH_minus = divideHist(histBHDVCSOutbFD_minus*histBHGenOutbFD50nA_minus , histBHOutbFD50nA_minus)
			accCorrectedOutbCD_BH_minus = divideHist(histBHDVCSOutbCD_minus*histBHGenOutbCD50nA_minus , histBHOutbCD50nA_minus)
			accCorrectedOutbCDFT_BH_minus = divideHist(histBHDVCSOutbCDFT_minus*histBHGenOutbCDFT50nA_minus , histBHOutbCDFT50nA_minus)
			accCorrectedOutbCR_BH_minus = divideHist(histBHDVCSOutbCR_minus*histBHGenOutbCR50nA_minus , histBHOutbCR50nA_minus)
			accCorrectedOutb_BH_minus = accCorrectedOutbFD_BH_minus + accCorrectedOutbCD_BH_minus + accCorrectedOutbCDFT_BH_minus + accCorrectedOutbCR_BH_minus
			#FD inbending threshold [:, :, [0, 1], :]
			accCorrectedOutb_BH_minus[:, :, [0, 1], :] = accCorrectedOutbCD_BH_minus[:, :, [0, 1], :] + accCorrectedOutbCDFT_BH_minus[:, :, [0, 1], :] + accCorrectedOutbCR_BH_minus[:, :, [0, 1], :]

			uncStatOutbFD_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD_minus**2+histExpUncOutbFD_minus**2), histBHDVCSOutbFD_minus)**2 + inverseHist(histBHOutbFD50nA_minus) + inverseHist(histBHGenOutbFD50nA_minus))
			uncStatOutbCD_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD_minus**2+histExpUncOutbCD_minus**2), histBHDVCSOutbCD_minus)**2 + inverseHist(histBHOutbCD50nA_minus) + inverseHist(histBHGenOutbCD50nA_minus))
			uncStatOutbCDFT_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT_minus**2+histExpUncOutbCDFT_minus**2), histBHDVCSOutbCDFT_minus)**2 + inverseHist(histBHOutbCDFT50nA_minus) + inverseHist(histBHGenOutbCDFT50nA_minus))
			uncStatOutbCR_BH_minus = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR_minus**2+histExpUncOutbCR_minus**2), histBHDVCSOutbCR_minus)**2 + inverseHist(histBHOutbCR50nA_minus) + inverseHist(histBHGenOutbCR50nA_minus))
			uncStatOutb_BH_minus = np.sqrt(accCorrectedOutbFD_BH_minus**2 *uncStatOutbFD_BH_minus**2 + accCorrectedOutbCD_BH_minus**2 * uncStatOutbCD_BH_minus**2 + accCorrectedOutbCDFT_BH_minus**2 * uncStatOutbCDFT_BH_minus**2 + accCorrectedOutbCR_BH_minus**2 * uncStatOutbCR_BH_minus**2)
			uncStatOutb_BH_minus = divideHist(uncStatOutb_BH_minus, accCorrectedOutb_BH_minus)
			#FD inbending threshold [:, :, [0, 1], :]
			uncStatOutb_BH_minus[:, :, [0, 1], :] = np.sqrt(accCorrectedOutbCD_BH_minus**2 * uncStatOutbCD_BH_minus**2 + accCorrectedOutbCDFT_BH_minus**2 * uncStatOutbCDFT_BH_minus**2 + accCorrectedOutbCR_BH_minus**2 * uncStatOutbCR_BH_minus**2)[:, :, [0, 1], :]
			uncStatOutb_BH_minus[:, :, [0, 1], :] = divideHist(uncStatOutb_BH_minus[:, :, [0, 1], :], accCorrectedOutb_BH_minus[:, :, [0, 1], :])

			accCorrectedOutb_BH_minus = divideHist(accCorrectedOutb_BH_minus*histBHGenOutb50nA_minus, histBHGenOutbFD50nA_minus + histBHGenOutbCD50nA_minus + histBHGenOutbCDFT50nA_minus + histBHGenOutbCR50nA_minus)
			accCorrectedOutb_BH_minus[:, :, [0, 1], :] = divideHist(accCorrectedOutb_BH_minus[:, :, [0, 1], :]*histBHGenOutb50nA_minus[:, :, [0, 1], :], histBHGenOutbCD50nA_minus[:, :, [0, 1], :] + histBHGenOutbCDFT50nA_minus[:, :, [0, 1], :] + histBHGenOutbCR50nA_minus[:, :, [0, 1], :])

			#config specific
			accCorrectedOutbFDonly_BH = divideHist(histBHDVCSOutbFD*histBHGenOutb50nA , histBHOutbFD50nA)
			accCorrectedOutbCDonly_BH = divideHist(histBHDVCSOutbCD*histBHGenOutb50nA , histBHOutbCD50nA)
			accCorrectedOutbCDFTonly_BH = divideHist(histBHDVCSOutbCDFT*histBHGenOutb50nA , histBHOutbCDFT50nA)
			accCorrectedOutbCRonly_BH = divideHist(histBHDVCSOutbCR*histBHGenOutb50nA , histBHOutbCR50nA)

			uncStatOutbFDonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD**2+histExpUncOutbFD**2), histBHDVCSOutbFD)**2 + inverseHist(histBHOutbFD50nA) + inverseHist(histBHGenOutb50nA))
			uncStatOutbCDonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD**2+histExpUncOutbCD**2), histBHDVCSOutbCD)**2 + inverseHist(histBHOutbCD50nA) + inverseHist(histBHGenOutb50nA))
			uncStatOutbCDFTonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT**2+histExpUncOutbCDFT**2), histBHDVCSOutbCDFT)**2 + inverseHist(histBHOutbCDFT50nA) + inverseHist(histBHGenOutb50nA))
			uncStatOutbCRonly_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCR**2+histExpUncOutbCR**2), histBHDVCSOutbCR)**2 + inverseHist(histBHOutbCR50nA) + inverseHist(histBHGenOutb50nA))
			#remove inactive
			accCorrectedOutb_BH[~ActiveOutb] = 0
			uncStatOutb_BH[~ActiveOutb] = 0
			accCorrectedOutb_BH_plus[~ActiveOutb] = 0
			uncStatOutb_BH_plus[~ActiveOutb] = 0
			accCorrectedOutb_BH_minus[~ActiveOutb] = 0
			uncStatOutb_BH_minus[~ActiveInb] = 0

			#stat error - all
			accCorrected_VGG = accCorrectedInb_VGG + accCorrectedOutb_VGG
			uncStat_VGG = np.sqrt((accCorrectedInb_VGG*uncStatInb_VGG)**2 + (accCorrectedOutb_VGG*uncStatOutb_VGG)**2)
			uncStat_VGG = divideHist(uncStat_VGG, accCorrected_VGG)

			accCorrected_BH = accCorrectedInb_BH + accCorrectedOutb_BH
			uncStat_BH = np.sqrt((accCorrectedInb_BH*uncStatInb_BH)**2 + (accCorrectedOutb_BH*uncStatOutb_BH)**2)
			uncStat_BH = divideHist(uncStat_BH, accCorrected_BH)

			accCorrected_VGG_plus = accCorrectedInb_VGG_plus + accCorrectedOutb_VGG_plus
			uncStat_VGG_plus = np.sqrt((accCorrectedInb_VGG_plus*uncStatInb_VGG_plus)**2 + (accCorrectedOutb_VGG_plus*uncStatOutb_VGG_plus)**2)
			uncStat_VGG_plus = divideHist(uncStat_VGG_plus, accCorrected_VGG_plus)

			accCorrected_BH_plus = accCorrectedInb_BH_plus + accCorrectedOutb_BH_plus
			uncStat_BH_plus = np.sqrt((accCorrectedInb_BH_plus*uncStatInb_BH_plus)**2 + (accCorrectedOutb_BH_plus*uncStatOutb_BH_plus)**2)
			uncStat_BH_plus = divideHist(uncStat_BH_plus, accCorrected_BH_plus)

			accCorrected_VGG_minus = accCorrectedInb_VGG_minus + accCorrectedOutb_VGG_minus
			uncStat_VGG_minus = np.sqrt((accCorrectedInb_VGG_minus*uncStatInb_VGG_minus)**2 + (accCorrectedOutb_VGG_minus*uncStatOutb_VGG_minus)**2)
			uncStat_VGG_minus = divideHist(uncStat_VGG_minus, accCorrected_VGG_minus)

			accCorrected_BH_minus = accCorrectedInb_BH_minus + accCorrectedOutb_BH_minus
			uncStat_BH_minus = np.sqrt((accCorrectedInb_BH_minus*uncStatInb_BH_minus)**2 + (accCorrectedOutb_BH_minus*uncStatOutb_BH_minus)**2)
			uncStat_BH_minus = divideHist(uncStat_BH_minus, accCorrected_BH_minus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_VGG.npz".format(optionaltag, k, i), hist = accCorrectedInb_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_BH.npz".format(optionaltag, k, i), hist = accCorrectedInb_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_VGG.npz".format(optionaltag, k, i), hist = accCorrectedOutb_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_BH.npz".format(optionaltag, k, i), hist = accCorrectedOutb_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_VGG.npz".format(optionaltag, k, i), hist = accCorrected_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_BH.npz".format(optionaltag, k, i), hist = accCorrected_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_VGG_plus.npz".format(optionaltag, k, i), hist = accCorrectedInb_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_BH_plus.npz".format(optionaltag, k, i), hist = accCorrectedInb_BH_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_VGG_plus.npz".format(optionaltag, k, i), hist = accCorrectedOutb_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_BH_plus.npz".format(optionaltag, k, i), hist = accCorrectedOutb_BH_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_VGG_plus.npz".format(optionaltag, k, i), hist = accCorrected_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_BH_plus.npz".format(optionaltag, k, i), hist = accCorrected_BH_plus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_VGG_minus.npz".format(optionaltag, k, i), hist = accCorrectedInb_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_BH_minus.npz".format(optionaltag, k, i), hist = accCorrectedInb_BH_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_VGG_minus.npz".format(optionaltag, k, i), hist = accCorrectedOutb_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_BH_minus.npz".format(optionaltag, k, i), hist = accCorrectedOutb_BH_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_VGG_minus.npz".format(optionaltag, k, i), hist = accCorrected_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_BH_minus.npz".format(optionaltag, k, i), hist = accCorrected_BH_minus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbFDonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedInbFDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbFDonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedInbFDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedInbCDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedInbCDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDFTonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedInbCDFTonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDFTonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedInbCDFTonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCRonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedInbCRonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCRonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedInbCRonly_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbFDonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedOutbFDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbFDonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedOutbFDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedOutbCDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedOutbCDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDFTonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedOutbCDFTonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDFTonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedOutbCDFTonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCRonly_VGG.npz".format(optionaltag, k, i), hist = accCorrectedOutbCRonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCRonly_BH.npz".format(optionaltag, k, i), hist = accCorrectedOutbCRonly_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(optionaltag, k, i), hist = uncStatInb_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(optionaltag, k, i), hist = uncStatInb_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(optionaltag, k, i), hist = uncStatOutb_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(optionaltag, k, i), hist = uncStatOutb_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(optionaltag, k, i), hist = uncStat_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i), hist = uncStat_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_VGG_plus.npz".format(optionaltag, k, i), hist = uncStatInb_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_BH_plus.npz".format(optionaltag, k, i), hist = uncStatInb_BH_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_VGG_plus.npz".format(optionaltag, k, i), hist = uncStatOutb_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_BH_plus.npz".format(optionaltag, k, i), hist = uncStatOutb_BH_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG_plus.npz".format(optionaltag, k, i), hist = uncStat_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH_plus.npz".format(optionaltag, k, i), hist = uncStat_BH_plus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_VGG_minus.npz".format(optionaltag, k, i), hist = uncStatInb_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_BH_minus.npz".format(optionaltag, k, i), hist = uncStatInb_BH_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_VGG_minus.npz".format(optionaltag, k, i), hist = uncStatOutb_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_BH_minus.npz".format(optionaltag, k, i), hist = uncStatOutb_BH_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG_minus.npz".format(optionaltag, k, i), hist = uncStat_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH_minus.npz".format(optionaltag, k, i), hist = uncStat_BH_minus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbFDonly_VGG.npz".format(optionaltag, k, i), hist = uncStatInbFDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbFDonly_BH.npz".format(optionaltag, k, i), hist = uncStatInbFDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDonly_VGG.npz".format(optionaltag, k, i), hist = uncStatInbCDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDonly_BH.npz".format(optionaltag, k, i), hist = uncStatInbCDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDFTonly_VGG.npz".format(optionaltag, k, i), hist = uncStatInbCDFTonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDFTonly_BH.npz".format(optionaltag, k, i), hist = uncStatInbCDFTonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCRonly_VGG.npz".format(optionaltag, k, i), hist = uncStatInbCRonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCRonly_BH.npz".format(optionaltag, k, i), hist = uncStatInbCRonly_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbFDonly_VGG.npz".format(optionaltag, k, i), hist = uncStatOutbFDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbFDonly_BH.npz".format(optionaltag, k, i), hist = uncStatOutbFDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDonly_VGG.npz".format(optionaltag, k, i), hist = uncStatOutbCDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDonly_BH.npz".format(optionaltag, k, i), hist = uncStatOutbCDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDFTonly_VGG.npz".format(optionaltag, k, i), hist = uncStatOutbCDFTonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDFTonly_BH.npz".format(optionaltag, k, i), hist = uncStatOutbCDFTonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCRonly_VGG.npz".format(optionaltag, k, i), hist = uncStatOutbCRonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCRonly_BH.npz".format(optionaltag, k, i), hist = uncStatOutbCRonly_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll.npz".format(optionaltag, k, i), hist = ActiveAll)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i), hist = ActiveAny)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb.npz".format(optionaltag, k, i), hist = ActiveInb)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb.npz".format(optionaltag, k, i), hist = ActiveOutb)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(optionaltag, k, i), hist = ActiveAll_int)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(optionaltag, k, i), hist = ActiveAny_int)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(optionaltag, k, i), hist = ActiveInb_int)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(optionaltag, k, i), hist = ActiveOutb_int)

if args.savekine:
	# nominal only
	for k in range(kstart, kend):
		print("bin scheme {}".format(k))
		ActiveAll       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAll.npz".format(k, 0))["hist"]
		ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAny.npz".format(k, 0))["hist"]
		ActiveInb          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveInb.npz".format(k, 0))["hist"]
		ActiveOutb         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveOutb.npz".format(k, 0))["hist"]

		ActiveAll_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(k, 0))["hist"]
		ActiveAny_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(k, 0))["hist"]
		ActiveInb_int          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(k, 0))["hist"]
		ActiveOutb_int         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(k, 0))["hist"]

		#read the <xB>, <Q2>, <t>, <phi>, <sigma>
		histVGGGenInbxB45nA, histVGGGenInbQ245nA, histVGGGenInbt145nA, histVGGGenInbphi45nA = [], [], [], []
		histVGGGenInbrad45nA, histVGGGenInbborn45nA = [], []
		histVGGGenInbrad45nA_plus, histVGGGenInbborn45nA_plus = [], []
		histVGGGenInbrad45nA_minus, histVGGGenInbborn45nA_minus = [], []
		for jobNum in runs_inb_vgg45nA:
			histVGGGenInbxB45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"])
			histVGGGenInbQ245nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"])
			histVGGGenInbt145nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"])
			histVGGGenInbphi45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"])
			histVGGGenInbrad45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
			histVGGGenInbborn45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
			histVGGGenInbrad45nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradplus.npz".format(k, jobNum))["hist"])
			histVGGGenInbborn45nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornplus.npz".format(k, jobNum))["hist"])
			histVGGGenInbrad45nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradminus.npz".format(k, jobNum))["hist"])
			histVGGGenInbborn45nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornminus.npz".format(k, jobNum))["hist"])

		histBHGenInbxB45nA, histBHGenInbQ245nA, histBHGenInbt145nA, histBHGenInbphi45nA = [], [], [], []
		histBHGenInbrad45nA, histBHGenInbborn45nA = [], []
		histBHGenInbrad45nA_plus, histBHGenInbborn45nA_plus = [], []
		histBHGenInbrad45nA_minus, histBHGenInbborn45nA_minus = [], []
		for jobNum in runs_inb_bh45nA:
			histBHGenInbxB45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"])
			histBHGenInbQ245nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"])
			histBHGenInbt145nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"])
			histBHGenInbphi45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"])
			histBHGenInbrad45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
			histBHGenInbborn45nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
			histBHGenInbrad45nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradplus.npz".format(k, jobNum))["hist"])
			histBHGenInbborn45nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornplus.npz".format(k, jobNum))["hist"])
			histBHGenInbrad45nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradminus.npz".format(k, jobNum))["hist"])
			histBHGenInbborn45nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornminus.npz".format(k, jobNum))["hist"])

		histVGGGenOutbxB50nA, histVGGGenOutbQ250nA, histVGGGenOutbt150nA, histVGGGenOutbphi50nA = [], [], [], []
		histVGGGenOutbrad50nA, histVGGGenOutbborn50nA = [], []
		histVGGGenOutbrad50nA_plus, histVGGGenOutbborn50nA_plus = [], []
		histVGGGenOutbrad50nA_minus, histVGGGenOutbborn50nA_minus = [], []
		for jobNum in runs_outb_vgg50nA:
			histVGGGenOutbxB50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"])
			histVGGGenOutbQ250nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"])
			histVGGGenOutbt150nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"])
			histVGGGenOutbphi50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"])
			histVGGGenOutbrad50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
			histVGGGenOutbborn50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
			histVGGGenOutbrad50nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradplus.npz".format(k, jobNum))["hist"])
			histVGGGenOutbborn50nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornplus.npz".format(k, jobNum))["hist"])
			histVGGGenOutbrad50nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradminus.npz".format(k, jobNum))["hist"])
			histVGGGenOutbborn50nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornminus.npz".format(k, jobNum))["hist"])

		histBHGenOutbxB50nA, histBHGenOutbQ250nA, histBHGenOutbt150nA, histBHGenOutbphi50nA = [], [], [], []
		histBHGenOutbrad50nA, histBHGenOutbborn50nA = [], []
		histBHGenOutbrad50nA_plus, histBHGenOutbborn50nA_plus = [], []
		histBHGenOutbrad50nA_minus, histBHGenOutbborn50nA_minus = [], []
		for jobNum in runs_outb_bh50nA:
			histBHGenOutbxB50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"])
			histBHGenOutbQ250nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"])
			histBHGenOutbt150nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"])
			histBHGenOutbphi50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"])
			histBHGenOutbrad50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
			histBHGenOutbborn50nA.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
			histBHGenOutbrad50nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradplus.npz".format(k, jobNum))["hist"])
			histBHGenOutbborn50nA_plus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornplus.npz".format(k, jobNum))["hist"])
			histBHGenOutbrad50nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genradminus.npz".format(k, jobNum))["hist"])
			histBHGenOutbborn50nA_minus.append(np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/{}Genbornminus.npz".format(k, jobNum))["hist"])

		phi1avg_VGG = np.mean(histVGGGenOutbphi50nA, axis = 0)
		xBavg_VGG = np.stack([np.mean(histVGGGenOutbxB50nA, axis = 0)]*(len(phibins)-1), axis = -1)
		Q2avg_VGG = np.stack([np.mean(histVGGGenOutbQ250nA, axis = 0)]*(len(phibins)-1), axis = -1)
		t1avg_VGG = np.stack([np.mean(histVGGGenOutbt150nA, axis = 0)]*(len(phibins)-1), axis = -1)
		phi1avg_VGG[~ActiveAny_int] = 0
		xBavg_VGG[~ActiveAny_int] = 0
		Q2avg_VGG[~ActiveAny_int] = 0
		t1avg_VGG[~ActiveAny_int] = 0

		phi1avg_BH = np.mean(histBHGenOutbphi50nA, axis = 0)
		xBavg_BH = np.stack([np.mean(histBHGenOutbxB50nA, axis = 0)]*(len(phibins)-1), axis = -1)
		Q2avg_BH = np.stack([np.mean(histBHGenOutbQ250nA, axis = 0)]*(len(phibins)-1), axis = -1)
		t1avg_BH = np.stack([np.mean(histBHGenOutbt150nA, axis = 0)]*(len(phibins)-1), axis = -1)
		phi1avg_BH[~ActiveAny_int] = 0
		xBavg_BH[~ActiveAny_int] = 0
		Q2avg_BH[~ActiveAny_int] = 0
		t1avg_BH[~ActiveAny_int] = 0

		xsecTh_KM = np.zeros(phi1avg_BH.shape)
		xsecTh_BH = np.zeros(phi1avg_BH.shape)
		xsecTh_VGG = np.zeros(phi1avg_BH.shape)
		xsecTh_KM_plus = np.zeros(phi1avg_BH.shape)
		xsecTh_VGG_plus = np.zeros(phi1avg_BH.shape)
		xsecTh_KM_minus = np.zeros(phi1avg_BH.shape)
		xsecTh_VGG_minus = np.zeros(phi1avg_BH.shape)
		xsecTh_VGG[ActiveAny_int] = np.array(printVGGarray(xBavg_VGG[ActiveAny_int], Q2avg_VGG[ActiveAny_int], t1avg_VGG[ActiveAny_int], np.radians(phi1avg_VGG[ActiveAny_int]), globalfit = True))
		xsecTh_BH[ActiveAny_int] = np.array(printBHarray(xBavg_BH[ActiveAny_int], Q2avg_BH[ActiveAny_int], t1avg_BH[ActiveAny_int], np.radians(phi1avg_BH[ActiveAny_int]), globalfit = True))
		xsecTh_KM[ActiveAny_int] = np.array(printKMarray(xBavg_BH[ActiveAny_int], Q2avg_BH[ActiveAny_int], t1avg_BH[ActiveAny_int], np.radians(phi1avg_BH[ActiveAny_int])))
		xsecTh_VGG_plus[ActiveAny_int] = np.array(printVGGarray(xBavg_VGG[ActiveAny_int], Q2avg_VGG[ActiveAny_int], t1avg_VGG[ActiveAny_int], np.radians(phi1avg_VGG[ActiveAny_int]), globalfit = True, pol = 1))
		xsecTh_KM_plus[ActiveAny_int] = np.array(printKMarray(xBavg_BH[ActiveAny_int], Q2avg_BH[ActiveAny_int], t1avg_BH[ActiveAny_int], np.radians(phi1avg_BH[ActiveAny_int]), pol = 1))
		xsecTh_VGG_minus[ActiveAny_int] = np.array(printVGGarray(xBavg_VGG[ActiveAny_int], Q2avg_VGG[ActiveAny_int], t1avg_VGG[ActiveAny_int], np.radians(phi1avg_VGG[ActiveAny_int]), globalfit = True, pol = -1))
		xsecTh_KM_minus[ActiveAny_int] = np.array(printKMarray(xBavg_BH[ActiveAny_int], Q2avg_BH[ActiveAny_int], t1avg_BH[ActiveAny_int], np.radians(phi1avg_BH[ActiveAny_int]), pol = -1))
		xsecTh_BH_plus = xsecTh_BH
		xsecTh_BH_minus = xsecTh_BH

		#use the outbending only for now.
		integratedRad_VGG = np.mean(histVGGGenOutbrad50nA, axis = 0)
		integratedBorn_VGG = np.mean(histVGGGenOutbborn50nA, axis = 0)
		rcfactors_VGG = divideHist(integratedRad_VGG, xsecTh_VGG)
		integratedRad_BH = np.mean(histBHGenOutbrad50nA, axis = 0)
		integratedBorn_BH = np.mean(histBHGenOutbborn50nA, axis = 0)
		rcfactors_BH = divideHist(integratedRad_BH, xsecTh_BH)

		integratedRad_VGG_plus = np.mean(histVGGGenOutbrad50nA_plus, axis = 0)
		integratedBorn_VGG_plus = np.mean(histVGGGenOutbborn50nA_plus, axis = 0)
		rcfactors_VGG_plus = divideHist(integratedRad_VGG_plus, xsecTh_VGG_plus)
		integratedRad_BH_plus = np.mean(histBHGenOutbrad50nA_plus, axis = 0)
		integratedBorn_BH_plus = np.mean(histBHGenOutbborn50nA_plus, axis = 0)
		rcfactors_BH_plus = divideHist(integratedRad_BH_plus, xsecTh_BH_plus)

		integratedRad_VGG_minus = np.mean(histVGGGenOutbrad50nA_minus, axis = 0)
		integratedBorn_VGG_minus = np.mean(histVGGGenOutbborn50nA_minus, axis = 0)
		rcfactors_VGG_minus = divideHist(integratedRad_VGG_minus, xsecTh_VGG_minus)
		integratedRad_BH_minus = np.mean(histBHGenOutbrad50nA_minus, axis = 0)
		integratedBorn_BH_minus = np.mean(histBHGenOutbborn50nA_minus, axis = 0)
		rcfactors_BH_minus = divideHist(integratedRad_BH_minus, xsecTh_BH_minus)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_VGG.npz".format(k), hist = phi1avg_VGG)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_VGG.npz".format(k), hist = xBavg_VGG)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_VGG.npz".format(k), hist = Q2avg_VGG)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_VGG.npz".format(k), hist = t1avg_VGG)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_BH.npz".format(k), hist = phi1avg_BH)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_BH.npz".format(k), hist = xBavg_BH)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_BH.npz".format(k), hist = Q2avg_BH)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_BH.npz".format(k), hist = t1avg_BH)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM.npz".format(k), hist = xsecTh_KM)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH.npz".format(k), hist = xsecTh_BH)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG.npz".format(k), hist = xsecTh_VGG)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM_plus.npz".format(k), hist = xsecTh_KM_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH_plus.npz".format(k), hist = xsecTh_BH_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG_plus.npz".format(k), hist = xsecTh_VGG_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM_minus.npz".format(k), hist = xsecTh_KM_minus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH_minus.npz".format(k), hist = xsecTh_BH_minus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG_minus.npz".format(k), hist = xsecTh_VGG_minus)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_VGG.npz".format(k), hist = integratedRad_VGG)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_VGG.npz".format(k), hist = integratedBorn_VGG)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG.npz".format(k), hist = rcfactors_VGG)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_BH.npz".format(k), hist = integratedRad_BH)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_BH.npz".format(k), hist = integratedBorn_BH)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH.npz".format(k), hist = rcfactors_BH)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_VGG_plus.npz".format(k), hist = integratedRad_VGG_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_VGG_plus.npz".format(k), hist = integratedBorn_VGG_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG_plus.npz".format(k), hist = rcfactors_VGG_plus)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_BH_plus.npz".format(k), hist = integratedRad_BH_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_BH_plus.npz".format(k), hist = integratedBorn_BH_plus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH_plus.npz".format(k), hist = rcfactors_BH_plus)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_VGG_minus.npz".format(k), hist = integratedRad_VGG_minus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_VGG_minus.npz".format(k), hist = integratedBorn_VGG_minus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG_minus.npz".format(k), hist = rcfactors_VGG_minus)

		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_BH_minus.npz".format(k), hist = integratedRad_BH_minus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_BH_minus.npz".format(k), hist = integratedBorn_BH_minus)
		np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH_minus.npz".format(k), hist = rcfactors_BH_minus)

if args.savexsec:
	print("read exp...")
	epgExp = pd.read_pickle("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/epgExp.pkl".format(optionaltag))

	for k in range(kstart, kend):

		print("bin scheme {}".format(k))
		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		binVolume = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/binVolume.npz".format(k))["hist"]
		rcfactors_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH.npz".format(k))["hist"]
		rcfactors_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG.npz".format(k))["hist"]
		rcfactors_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH_plus.npz".format(k))["hist"]
		rcfactors_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG_plus.npz".format(k))["hist"]
		rcfactors_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH_minus.npz".format(k))["hist"]
		rcfactors_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG_minus.npz".format(k))["hist"]

		for i in range(istart, len(collection_cont_xBbins)):
			#Inbending cross sections 
			# i = 0 #selected background estimation
			accCorrectedInb_VGG	 = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInb_BH	 = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutb_BH	 = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrected_VGG	 = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrected_BH	 = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_BH.npz".format(optionaltag, k, i))["hist"]

			accCorrectedInb_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInb_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_BH_plus.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutb_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutb_BH_plus= np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_BH_plus.npz".format(optionaltag, k, i))["hist"]
			accCorrected_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			accCorrected_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_BH_plus.npz".format(optionaltag, k, i))["hist"]

			accCorrectedInb_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInb_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInb_BH_minus.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutb_VGG_minus =  np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutb_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutb_BH_minus.npz".format(optionaltag, k, i))["hist"]
			accCorrected_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			accCorrected_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrected_BH_minus.npz".format(optionaltag, k, i))["hist"]

			accCorrectedInbFDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbFDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbFDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbFDonly_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbCDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbCDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDonly_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbCDFTonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDFTonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbCDFTonly_BH =  np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCDFTonly_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbCRonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCRonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedInbCRonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedInbCRonly_BH.npz".format(optionaltag, k, i))["hist"]

			accCorrectedOutbFDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbFDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbFDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbFDonly_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbCDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbCDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDonly_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbCDFTonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDFTonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbCDFTonly_BH =  np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCDFTonly_BH.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbCRonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCRonly_VGG.npz".format(optionaltag, k, i))["hist"]
			accCorrectedOutbCRonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}accCorrectedOutbCRonly_BH.npz".format(optionaltag, k, i))["hist"]

			ActiveAll       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll.npz".format(optionaltag, k, i))["hist"]
			ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
			ActiveInb          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb.npz".format(optionaltag, k, i))["hist"]
			ActiveOutb         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb.npz".format(optionaltag, k, i))["hist"]

			ActiveAll_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(optionaltag, k, i))["hist"]
			ActiveAny_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(optionaltag, k, i))["hist"]
			ActiveInb_int          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(optionaltag, k, i))["hist"]
			ActiveOutb_int         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(optionaltag, k, i))["hist"]

			print("bkg scheme {}".format(i))
			#calculate cross sections
			charges = charge_epg*np.zeros(ActiveAll.shape)
			charges[ActiveAll] = charge_epg
			charges[(ActiveInb)&(~ActiveOutb)] = inbcharge_epg
			charges[(ActiveOutb)&(~ActiveInb)] = outbcharge_epg
			xsecInb_VGG = divideHist(accCorrectedInb_VGG, binVolume*rcfactors_VGG)/(1.324*inbcharge_epg)
			xsecInb_BH = divideHist(accCorrectedInb_BH, binVolume*rcfactors_BH)/(1.324*inbcharge_epg)
			xsecInb_VGG_plus = divideHist(accCorrectedInb_VGG_plus, binVolume*rcfactors_VGG_plus)/(1.324*inbcharge_epg/2) # charge/2 for plus/minus
			xsecInb_BH_plus = divideHist(accCorrectedInb_BH_plus, binVolume*rcfactors_BH_plus)/(1.324*inbcharge_epg/2)
			xsecInb_VGG_minus = divideHist(accCorrectedInb_VGG_minus, binVolume*rcfactors_VGG_minus)/(1.324*inbcharge_epg/2)
			xsecInb_BH_minus = divideHist(accCorrectedInb_BH_minus, binVolume*rcfactors_BH_minus)/(1.324*inbcharge_epg/2)
			xsecInbFDonly_VGG = divideHist(accCorrectedInbFDonly_VGG, binVolume*rcfactors_VGG)/(1.324*inbcharge_epg)
			xsecInbFDonly_BH = divideHist(accCorrectedInbFDonly_BH, binVolume*rcfactors_BH)/(1.324*inbcharge_epg)
			xsecInbCDonly_VGG = divideHist(accCorrectedInbCDonly_VGG, binVolume*rcfactors_VGG)/(1.324*inbcharge_epg)
			xsecInbCDonly_BH = divideHist(accCorrectedInbCDonly_BH, binVolume*rcfactors_BH)/(1.324*inbcharge_epg)
			xsecInbCDFTonly_VGG = divideHist(accCorrectedInbCDFTonly_VGG, binVolume*rcfactors_VGG)/(1.324*inbcharge_epg)
			xsecInbCDFTonly_BH = divideHist(accCorrectedInbCDFTonly_BH, binVolume*rcfactors_BH)/(1.324*inbcharge_epg)
			xsecInbCRonly_VGG = divideHist(accCorrectedInbCRonly_VGG, binVolume*rcfactors_VGG)/(1.324*inbcharge_epg)
			xsecInbCRonly_BH = divideHist(accCorrectedInbCRonly_BH, binVolume*rcfactors_BH)/(1.324*inbcharge_epg)
			xsecOutb_VGG = divideHist(accCorrectedOutb_VGG, binVolume*rcfactors_VGG)/(1.324*outbcharge_epg)
			xsecOutb_BH = divideHist(accCorrectedOutb_BH, binVolume*rcfactors_BH)/(1.324*outbcharge_epg)
			xsecOutb_VGG_plus = divideHist(accCorrectedOutb_VGG_plus, binVolume*rcfactors_VGG_plus)/(1.324*outbcharge_epg/2)
			xsecOutb_BH_plus = divideHist(accCorrectedOutb_BH_plus, binVolume*rcfactors_BH_plus)/(1.324*outbcharge_epg/2)
			xsecOutb_VGG_minus = divideHist(accCorrectedOutb_VGG_minus, binVolume*rcfactors_VGG_minus)/(1.324*outbcharge_epg/2)
			xsecOutb_BH_minus = divideHist(accCorrectedOutb_BH_minus, binVolume*rcfactors_BH_minus)/(1.324*outbcharge_epg/2)
			xsecOutbFDonly_VGG = divideHist(accCorrectedOutbFDonly_VGG, binVolume*rcfactors_VGG)/(1.324*outbcharge_epg)
			xsecOutbFDonly_BH = divideHist(accCorrectedOutbFDonly_BH, binVolume*rcfactors_BH)/(1.324*outbcharge_epg)
			xsecOutbCDonly_VGG = divideHist(accCorrectedOutbCDonly_VGG, binVolume*rcfactors_VGG)/(1.324*outbcharge_epg)
			xsecOutbCDonly_BH = divideHist(accCorrectedOutbCDonly_BH, binVolume*rcfactors_BH)/(1.324*outbcharge_epg)
			xsecOutbCDFTonly_VGG = divideHist(accCorrectedOutbCDFTonly_VGG, binVolume*rcfactors_VGG)/(1.324*outbcharge_epg)
			xsecOutbCDFTonly_BH = divideHist(accCorrectedOutbCDFTonly_BH, binVolume*rcfactors_BH)/(1.324*outbcharge_epg)
			xsecOutbCRonly_VGG = divideHist(accCorrectedOutbCRonly_VGG, binVolume*rcfactors_VGG)/(1.324*outbcharge_epg)
			xsecOutbCRonly_BH = divideHist(accCorrectedOutbCRonly_BH, binVolume*rcfactors_BH)/(1.324*outbcharge_epg)
			xsec_VGG = divideHist(accCorrected_VGG, binVolume*rcfactors_VGG*(1.324*charges))
			xsec_BH = divideHist(accCorrected_BH, binVolume*rcfactors_BH*(1.324*charges))
			xsec_VGG_plus = divideHist(accCorrected_VGG_plus, binVolume*rcfactors_VGG_plus*(1.324*charges/2))
			xsec_BH_plus = divideHist(accCorrected_BH_plus, binVolume*rcfactors_BH_plus*(1.324*charges/2))
			xsec_VGG_minus = divideHist(accCorrected_VGG_minus, binVolume*rcfactors_VGG_minus*(1.324*charges/2))
			xsec_BH_minus = divideHist(accCorrected_BH_minus, binVolume*rcfactors_BH_minus*(1.324*charges/2))


			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(optionaltag, k, i), hist = xsecInb_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(optionaltag, k, i), hist = xsecInb_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(optionaltag, k, i), hist = xsecOutb_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(optionaltag, k, i), hist = xsecOutb_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG.npz".format(optionaltag, k, i), hist = xsec_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i), hist = xsec_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_VGG_plus.npz".format(optionaltag, k, i), hist = xsecInb_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_BH_plus.npz".format(optionaltag, k, i), hist = xsecInb_BH_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_VGG_plus.npz".format(optionaltag, k, i), hist = xsecOutb_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_BH_plus.npz".format(optionaltag, k, i), hist = xsecOutb_BH_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG_plus.npz".format(optionaltag, k, i), hist = xsec_VGG_plus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH_plus.npz".format(optionaltag, k, i), hist = xsec_BH_plus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_VGG_minus.npz".format(optionaltag, k, i), hist = xsecInb_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_BH_minus.npz".format(optionaltag, k, i), hist = xsecInb_BH_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_VGG_minus.npz".format(optionaltag, k, i), hist = xsecOutb_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_BH_minus.npz".format(optionaltag, k, i), hist = xsecOutb_BH_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG_minus.npz".format(optionaltag, k, i), hist = xsec_VGG_minus)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH_minus.npz".format(optionaltag, k, i), hist = xsec_BH_minus)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbFDonly_VGG.npz".format(optionaltag, k, i), hist = xsecInbFDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbFDonly_BH.npz".format(optionaltag, k, i), hist = xsecInbFDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDonly_VGG.npz".format(optionaltag, k, i), hist = xsecInbCDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDonly_BH.npz".format(optionaltag, k, i), hist = xsecInbCDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDFTonly_VGG.npz".format(optionaltag, k, i), hist = xsecInbCDFTonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDFTonly_BH.npz".format(optionaltag, k, i), hist = xsecInbCDFTonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCRonly_VGG.npz".format(optionaltag, k, i), hist = xsecInbCRonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCRonly_BH.npz".format(optionaltag, k, i), hist = xsecInbCRonly_BH)

			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbFDonly_VGG.npz".format(optionaltag, k, i), hist = xsecOutbFDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbFDonly_BH.npz".format(optionaltag, k, i), hist = xsecOutbFDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDonly_VGG.npz".format(optionaltag, k, i), hist = xsecOutbCDonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDonly_BH.npz".format(optionaltag, k, i), hist = xsecOutbCDonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDFTonly_VGG.npz".format(optionaltag, k, i), hist = xsecOutbCDFTonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDFTonly_BH.npz".format(optionaltag, k, i), hist = xsecOutbCDFTonly_BH)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCRonly_VGG.npz".format(optionaltag, k, i), hist = xsecOutbCRonly_VGG)
			np.savez("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCRonly_BH.npz".format(optionaltag, k, i), hist = xsecOutbCRonly_BH)

if args.savesyst:
	k = 3

	xBbins  = collection_xBbins[k]
	Q2bins  = collection_Q2bins[k]
	tbins   = collection_tbins [k]
	phibins = collection_phibins[k]

	# phi1avg_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_VGG.npz".format(k))["hist"]
	# xBavg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_VGG.npz".format(k))["hist"]
	# Q2avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_VGG.npz".format(k))["hist"]
	# t1avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_VGG.npz".format(k))["hist"]

	phi1avg_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_BH.npz".format(k))["hist"]
	xBavg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_BH.npz".format(k))["hist"]
	Q2avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_BH.npz".format(k))["hist"]
	t1avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_BH.npz".format(k))["hist"]

	xsecTh_KM          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM.npz".format(k))["hist"]
	xsecTh_BH          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH.npz".format(k))["hist"]
	xsecTh_VGG         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG.npz".format(k))["hist"]

	xBbin = 3
	Q2bin = 2
	tbin = 2

	fig, axs = plt.subplots(1, 1, figsize = (10, 6))
	colorscheme = ['k', 'purple', 'b', 'g', 'brown', 'orange', 'pink']
	scheme = 0
	for optionaltag in ['', '_bkg', '_2sigma', '_4sigma', '_sm09', '_sm11', '_tightfid']:
		if optionaltag == '_bkg':
			i = 0
			ActiveAll       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll.npz".format('', k, i))["hist"]
			ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format('', k, i))["hist"]
			ActiveInb       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb.npz".format('', k, i))["hist"]
			ActiveOutb      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb.npz".format('', k, i))["hist"]

			xsec_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG.npz".format('', k, i))["hist"]
			xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format('', k, i))["hist"]
			uncStat_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG.npz".format('', k, i))["hist"]
			uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format('', k, i))["hist"]
		else:
			i = 3
			ActiveAll       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll.npz".format(optionaltag, k, i))["hist"]
			ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
			ActiveInb       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb.npz".format(optionaltag, k, i))["hist"]
			ActiveOutb      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb.npz".format(optionaltag, k, i))["hist"]

			xsec_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG.npz".format(optionaltag, k, i))["hist"]
			xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]
			uncStat_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]

		phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
		P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
		P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])

		# axs.errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental data')
		res_lsq = curve_fit(FourierSeries, [0, 0, 0], args=(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin]))
		N = 40
		phi1s = np.linspace(0, 360, N)[:-1]
		phi1s = phi1s + np.diff(phi1s)[0]
		xBs = np.ones(len(phi1s))*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
		Q2s = np.ones(len(phi1s))*Q2avg_BH[xBbin, Q2bin, tbin, phibin][0]
		t1s = np.ones(len(phi1s))*t1avg_BH[xBbin, Q2bin, tbin, phibin][0]

		if optionaltag == '':
			axs.plot(phi1s, FourierSeries(res_lsq.x, phi1s), label = 'Fitting results', color = 'k', linestyle = '--')
			axs.errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH)[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = P1b*P2b*(xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental Data')
			axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsecTh_KM[xBbin, Q2bin, tbin, phibin], color = 'cyan', label = 'Theory (KM15)')
			axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsecTh_BH[xBbin, Q2bin, tbin, phibin], color = 'r', label = 'Theory (BH)')
			P1_zero = P1(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0)
			P2_zero = P2(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0)
			BHDVCS_zero = getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0, mode = 1)
			Normalization = FourierSeries(res_lsq.x, 0)/(P1_zero*P2_zero*BHDVCS_zero)
		else:
			axs.plot(phi1s, FourierSeries(res_lsq.x, phi1s), label = 'Fitting results {}'.format(optionaltag[1:]), color = colorscheme[scheme], linestyle = '--')
		scheme = scheme+1

	handles, labels = axs.get_legend_handles_labels()
	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, bbox_to_anchor = (1.0, 0.8))
	axs.set_xlim([0, 360])
	axs.set_xticks([0, 90, 180, 270, 360])
	axs.set_xlabel(r"$\phi$" + " ["+degree+"]")
	axs.set_ylabel(r"$\mathcal{P}_1(\phi)\mathcal{P}_2(\phi)\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")

	xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}, ".format(xBbins[xBbin+1]) +r"$~<x_B>=$"+ "{:.3f}\n".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
	Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f}, ".format(Q2bins[Q2bin+1])+ r"$~<Q^2>=$"+"{:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"$~(\mathrm{GeV/c})^2$"+ "\n"
	theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f}, ".format(tbins[tbin+1]) + r"$~<|t|>=$"+"{:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"$~\mathrm{GeV}^2$"
	header = xBheader + Q2header + theader
	axs.set_title(header, loc = 'left')

	# fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
	plt.savefig("plots/systematics.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')

	# P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])
	# P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])

		# tbin = 2
		# xBbin = 3

		# fig, axs = plt.subplots(4, 1, figsize = (12, 20))
		# for Q2bin in range(4):

		# 	print(Q2bin)

		# 	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()		
		# 	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
		# 	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])

		# 	axs[4-Q2bin-1].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], (xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental data')
		# 	res_lsq = curve_fit(FourierSeries, [0, 0, 0], args=(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin]))
		# 	N = 40
		# 	phi1s = np.linspace(0, 360, N)[:-1]
		# 	phi1s = phi1s + np.diff(phi1s)[0]
		# 	xBs = np.ones(len(phi1s))*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
		# 	Q2s = np.ones(len(phi1s))*Q2avg_BH[xBbin, Q2bin, tbin, phibin][0]
		# 	t1s = np.ones(len(phi1s))*t1avg_BH[xBbin, Q2bin, tbin, phibin][0]
		# 	P1b = P1(xBs, Q2s, t1s, phi1s)
		# 	P2b = P2(xBs, Q2s, t1s, phi1s)

		# 	axs[4-Q2bin-1].plot(phi1s, 1/P1b*1/P2b*FourierSeries(res_lsq.x, phi1s), label = 'Fitting results', color = 'k', linestyle = '--')
		# 	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])
		# 	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])

		# 	axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 5), color = 'cyan', label = 'Theory (BH+Int.+DVCS'+r"${}^{2}$"+")")
		# 	axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 1), color = 'r', label = 'Theory (Pure BH)')
		# 	# axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 2), color = 'g', label = 'Interference')
		# 	# axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 3), color = 'purple', label = 'DVCS'+r"${}^{2}$")
		# 	xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
		# 	Q2header = r"$,~<Q^2>=$"+" {:.3f} (GeV/c)".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"${}^{2}$"
		# 	theader = r"$,~<|t|>=$"+" {:.3f} GeV".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"${}^{2}$"
		# 	header = xBheader +Q2header + theader
		# 	axs[4-Q2bin-1].set_xlim([0, 360])
		# 	axs[4-Q2bin-1].set_ylim([0.005, 2])
		# 	axs[4-Q2bin-1].set_title(header, fontsize = 20)
		# 	axs[4-Q2bin-1].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")
		# 	axs[4-Q2bin-1].set_yscale('log')
		# 	axs[4-Q2bin-1].set_xticks([0, 90, 180, 270, 360])
		# 	axs[4-Q2bin-1].set_xlabel(r"$\phi$" + " [" + degree + "]")
		# 	if Q2bin == 0:
		# 		handles, labels = axs[4-Q2bin-1].get_legend_handles_labels()
		# 		order = [3, 0, 1, 2]

		# lgd = plt.figlegend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'upper left', bbox_to_anchor =(1.02, 0.5), title = "")
		# fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
		# plt.savefig("plots/richard_rolf_Q2.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
		# # plt.savefig("plots/richard_rolf{}.pdf".format(Q2bin))
		# plt.clf()

if args.savesyst2:
	k = 3
	i = 3

	xBbins  = collection_xBbins[k]
	Q2bins  = collection_Q2bins[k]
	tbins   = collection_tbins [k]
	phibins = collection_phibins[k]

	phi1avg_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_VGG.npz".format(k))["hist"]
	xBavg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_VGG.npz".format(k))["hist"]
	Q2avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_VGG.npz".format(k))["hist"]
	t1avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_VGG.npz".format(k))["hist"]

	phi1avg_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_BH.npz".format(k))["hist"]
	xBavg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_BH.npz".format(k))["hist"]
	Q2avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_BH.npz".format(k))["hist"]
	t1avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_BH.npz".format(k))["hist"]

	xsecTh_KM          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM.npz".format(k))["hist"]
	xsecTh_BH          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH.npz".format(k))["hist"]
	xsecTh_VGG         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG.npz".format(k))["hist"]


	fig, axs = plt.subplots(1, 1, figsize = (10, 6))
	colorscheme = ['k', 'purple', 'b', 'g', 'brown', 'orange', 'pink']
	scheme = 0
	# for optionaltag in ['', '_bkg', '_2sigma', '_4sigma', '_sm09', '_sm11', '_tightfid']:
	# 	if optionaltag == '_bkg':
	# 		i = 0
	xBbin = 3
	Q2bin = 2
	tbin = 2

	# P1_zero = P1(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0)
	# P2_zero = P2(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0)
	# BHDVCS_zero = getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0, mode = 1)
	reduced_zero = getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0, mode = 0)

	optionaltag = '_2sigma'
	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]
	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	tightexcl = popt

	optionaltag = '_4sigma'
	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]
	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	looseexcl = popt

	optionaltag = '_sm09'
	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]
	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	sm09 = popt

	optionaltag = '_sm11'
	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]
	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	sm11 = popt

	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format('', k, 0))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format('', k, 0))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format('', k, 0))["hist"]
	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	bkg = popt

	optionaltag = '_tightfid'
	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]
	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	tightfid = popt

	UncNorm = (getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0, mode = 5)/getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, 0], Q2avg_BH[xBbin, Q2bin, tbin, 0], t1avg_BH[xBbin, Q2bin, tbin, 0], 0, mode = 1)-1)/2

	#nominal
	ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format('', k, i))["hist"]

	xsec_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG.npz".format('', k, i))["hist"]
	xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format('', k, i))["hist"]
	uncStat_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG.npz".format('', k, i))["hist"]
	uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format('', k, i))["hist"]

	phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
	P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
	P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])

	popt, pcov = curve_fit(FourierSeries, phi1avg_VGG[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_VGG[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_VGG[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	VGG = popt
	popt, pcov = curve_fit(FourierSeries, phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*xsec_BH[xBbin, Q2bin, tbin, phibin], p0 =[0, 0, 0], sigma = P1b*P2b*uncStat_BH[xBbin, Q2bin, tbin, phibin], absolute_sigma = True)
	nominal = popt

	Normalization = FourierSeries(0, *nominal)/reduced_zero
	Normalization_2sigma = FourierSeries(0, *tightexcl)/reduced_zero
	Normalization_4sigma = FourierSeries(0, *looseexcl)/reduced_zero
	Normalization_sm09 = FourierSeries(0, *sm09)/reduced_zero
	Normalization_sm11 = FourierSeries(0, *sm11)/reduced_zero
	Normalization_tightfid = FourierSeries(0, *tightfid)/reduced_zero

	nominal = np.divide(nominal, Normalization)
	tightexcl = np.divide(tightexcl,Normalization_2sigma)
	looseexcl = np.divide(looseexcl,Normalization_4sigma)
	sm09 = np.divide(sm09,Normalization_sm09)
	sm11 = np.divide(sm11,Normalization_sm11)
	tightfid = np.divide(tightfid,Normalization_tightfid)

	UncModel = np.abs(divideHist(FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(VGG-nominal)),FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(nominal)), threshold=-np.inf))
	UncExcl = 0.5*np.abs(divideHist(FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(tightexcl-looseexcl)),FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(nominal)), threshold=-np.inf))
	UncSmear = 0.5*np.abs(divideHist(FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(sm11-sm09)),FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(nominal)), threshold=-np.inf))
	UncFid = np.abs(divideHist(FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(tightfid-nominal)),FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(nominal)), threshold=-np.inf))
	UncBkg = np.abs(divideHist(FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(bkg-nominal)),FourierSeries(phi1avg_BH[xBbin, Q2bin, tbin, phibin],*(nominal)), threshold=-np.inf))

	SystUnc = np.sqrt(UncNorm**2+ UncModel**2 + UncExcl**2 + UncSmear**2 + UncFid**2 + UncBkg**2 + 0.04**2)
	# Unc = np.sqrt(uncStat_BH[xBbin, Q2bin, tbin, phibin]**2 + SystUnc**2)

	axs.errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], (xsec_BH/Normalization)[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (uncStat_BH*xsec_BH/Normalization)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental data')
	for i in range(-500, 500):
		axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, phibin], (1+i*SystUnc/500)*(xsec_BH/Normalization)[xBbin, Q2bin, tbin, phibin], color = 'g', alpha = 1/100)

	axs.plot(np.linspace(0, 360, 40), 1/(P1b*P2b)*FourierSeries(np.linspace(0, 360, 40),*(nominal)), label = 'Fitting results', color = 'k', linestyle = '--')
	axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin], color = 'cyan', label = 'Theory (KM15)')
	Nplot = 40
	xBs = np.ones(Nplot)*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
	xBs = np.ones(Nplot)*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
	xBs = np.ones(Nplot)*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
	xBs = np.ones(Nplot)*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]

	axs.plot(np.linspace(0, 360, Nplot), getBHDVCS(np.linspace(0, 360, Nplot), mode  =1), color = 'r', label = 'Theory (BH)')
	axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin], color = 'purple', label = 'Theory (BH)')

	handles, labels = axs.get_legend_handles_labels()
	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, bbox_to_anchor = (1.0, 0.8))
	axs.set_xlim([0, 360])
	axs.set_xticks([0, 90, 180, 270, 360])
	axs.set_xlabel(r"$\phi$" + " ["+degree+"]")
	axs.set_ylabel(r"$\mathcal{P}_1(\phi)\mathcal{P}_2(\phi)\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")

	xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}, ".format(xBbins[xBbin+1]) +r"$~<x_B>=$"+ "{:.3f}\n".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
	Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f}, ".format(Q2bins[Q2bin+1])+ r"$~<Q^2>=$"+"{:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"$~(\mathrm{GeV/c})^2$"+ "\n"
	theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f}, ".format(tbins[tbin+1]) + r"$~<|t|>=$"+"{:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"$~\mathrm{GeV}^2$"
	header = xBheader + Q2header + theader
	axs.set_title(header, loc = 'left')

	# fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
	plt.savefig("plots/systematics2.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')

if args.saveplot:

	for k in range(kstart, kend):

		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		os.makedirs("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}".format(optionaltag, k), exist_ok = True)

		phi1avg_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_VGG.npz".format(k))["hist"]
		xBavg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_VGG.npz".format(k))["hist"]
		Q2avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_VGG.npz".format(k))["hist"]
		t1avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_VGG.npz".format(k))["hist"]

		phi1avg_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_BH.npz".format(k))["hist"]
		xBavg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_BH.npz".format(k))["hist"]
		Q2avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_BH.npz".format(k))["hist"]
		t1avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_BH.npz".format(k))["hist"]

		xsecTh_KM          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM.npz".format(k))["hist"]
		xsecTh_BH          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH.npz".format(k))["hist"]
		xsecTh_VGG         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG.npz".format(k))["hist"]
		xsecTh_KM_plus     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM_plus.npz".format(k))["hist"]
		xsecTh_BH_plus     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH_plus.npz".format(k))["hist"]
		xsecTh_VGG_plus    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG_plus.npz".format(k))["hist"]
		xsecTh_KM_minus     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_KM_minus.npz".format(k))["hist"]
		xsecTh_BH_minus     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH_minus.npz".format(k))["hist"]
		xsecTh_VGG_minus    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG_minus.npz".format(k))["hist"]
		binVolume          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/binVolume.npz".format(k))["hist"]

		for i in range(istart, len(collection_cont_xBbins)):

			print("reading the xsec vars")

			ActiveAll       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll.npz".format(optionaltag, k, i))["hist"]
			ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny.npz".format(optionaltag, k, i))["hist"]
			ActiveInb       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb.npz".format(optionaltag, k, i))["hist"]
			ActiveOutb      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb.npz".format(optionaltag, k, i))["hist"]

			ActiveAll_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(optionaltag, k, i))["hist"]
			ActiveAny_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(optionaltag, k, i))["hist"]
			ActiveInb_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(optionaltag, k, i))["hist"]
			ActiveOutb_int      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(optionaltag, k, i))["hist"]

			xsecInb_VGG  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecInb_BH   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(optionaltag, k, i))["hist"]
			xsecOutb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecOutb_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(optionaltag, k, i))["hist"]
			xsec_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG.npz".format(optionaltag, k, i))["hist"]
			xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH.npz".format(optionaltag, k, i))["hist"]

			uncStatInb_VGG  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatInb_BH   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatOutb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatOutb_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(optionaltag, k, i))["hist"]
			uncStat_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH.npz".format(optionaltag, k, i))["hist"]

			xsecInb_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			xsecInb_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_BH_plus.npz".format(optionaltag, k, i))["hist"]
			xsecOutb_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			xsecOutb_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_BH_plus.npz".format(optionaltag, k, i))["hist"]
			xsec_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			xsec_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH_plus.npz".format(optionaltag, k, i))["hist"]

			xsecInb_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			xsecInb_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInb_BH_minus.npz".format(optionaltag, k, i))["hist"]
			xsecOutb_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			xsecOutb_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutb_BH_minus.npz".format(optionaltag, k, i))["hist"]
			xsec_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			xsec_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsec_BH_minus.npz".format(optionaltag, k, i))["hist"]

			xsecInbFDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbFDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecInbFDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbFDonly_BH.npz".format(optionaltag, k, i))["hist"]
			xsecInbCDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecInbCDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDonly_BH.npz".format(optionaltag, k, i))["hist"]
			xsecInbCDFTonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDFTonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecInbCDFTonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCDFTonly_BH.npz".format(optionaltag, k, i))["hist"]
			xsecInbCRonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCRonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecInbCRonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecInbCRonly_BH.npz".format(optionaltag, k, i))["hist"]

			xsecOutbFDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbFDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecOutbFDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbFDonly_BH.npz".format(optionaltag, k, i))["hist"]
			xsecOutbCDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecOutbCDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDonly_BH.npz".format(optionaltag, k, i))["hist"]
			xsecOutbCDFTonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDFTonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecOutbCDFTonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCDFTonly_BH.npz".format(optionaltag, k, i))["hist"]
			xsecOutbCRonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCRonly_VGG.npz".format(optionaltag, k, i))["hist"]
			xsecOutbCRonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}xsecOutbCRonly_BH.npz".format(optionaltag, k, i))["hist"]

			uncStatInb_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			uncStatInb_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_BH_plus.npz".format(optionaltag, k, i))["hist"]
			uncStatOutb_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			uncStatOutb_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_BH_plus.npz".format(optionaltag, k, i))["hist"]
			uncStat_VGG_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG_plus.npz".format(optionaltag, k, i))["hist"]
			uncStat_BH_plus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH_plus.npz".format(optionaltag, k, i))["hist"]

			uncStatInb_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			uncStatInb_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInb_BH_minus.npz".format(optionaltag, k, i))["hist"]
			uncStatOutb_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			uncStatOutb_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutb_BH_minus.npz".format(optionaltag, k, i))["hist"]
			uncStat_VGG_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_VGG_minus.npz".format(optionaltag, k, i))["hist"]
			uncStat_BH_minus = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStat_BH_minus.npz".format(optionaltag, k, i))["hist"]

			uncStatInbFDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbFDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatInbFDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbFDonly_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatInbCDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatInbCDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDonly_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatInbCDFTonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDFTonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatInbCDFTonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCDFTonly_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatInbCRonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCRonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatInbCRonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatInbCRonly_BH.npz".format(optionaltag, k, i))["hist"]

			uncStatOutbFDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbFDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbFDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbFDonly_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbCDonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbCDonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDonly_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbCDFTonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDFTonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbCDFTonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCDFTonly_BH.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbCRonly_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCRonly_VGG.npz".format(optionaltag, k, i))["hist"]
			uncStatOutbCRonly_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms{}/binscheme{}/bkgscheme{}uncStatOutbCRonly_BH.npz".format(optionaltag, k, i))["hist"]

			print("plotting...")

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


			# #all config bh
			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInbFDonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbFDonly_BH*uncStatInbFDonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Inb. FD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInbCDonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbCDonly_BH*uncStatInbCDonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'r', label = 'Inb. CD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInbCDFTonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbCDFTonly_BH*uncStatInbCDFTonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb. CDFT only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInbCRonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbCRonly_BH*uncStatInbCRonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'b', label = 'Inb. CR only')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutbFDonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbFDonly_BH*uncStatOutbFDonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'purple', label = 'Outb. FD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutbCDonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbCDonly_BH*uncStatOutbCDonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'brown', label = 'Outb. CD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutbCDFTonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbCDFTonly_BH*uncStatOutbCDFTonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'orange', label = 'Outb. CDFT only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutbCRonly_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbCRonly_BH*uncStatOutbCRonly_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb. CR only')

			# 			# phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			# axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'BH')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			# if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 			# 	handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 			# 	active = 1
			# 	# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	# plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/all_configs_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/all_configs_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin))
			# 	plt.clf()

			# #all config vgg
			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecInbFDonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbFDonly_VGG*uncStatInbFDonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Inb. FD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecInbCDonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbCDonly_VGG*uncStatInbCDonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'r', label = 'Inb. CD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecInbCDFTonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbCDFTonly_VGG*uncStatInbCDFTonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb. CDFT only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecInbCRonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInbCRonly_VGG*uncStatInbCRonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'b', label = 'Inb. CR only')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecOutbFDonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbFDonly_VGG*uncStatOutbFDonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'purple', label = 'Outb. FD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecOutbCDonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbCDonly_VGG*uncStatOutbCDonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'brown', label = 'Outb. CD only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecOutbCDFTonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbCDFTonly_VGG*uncStatOutbCDFTonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'orange', label = 'Outb. CDFT only')
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsecOutbCRonly_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutbCRonly_VGG*uncStatOutbCRonly_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb. CR only')

			# 			# phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			# axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_VGG[xBbin, Q2bin, tbin, phibin], xsec_VGG[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_VGG[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_VGG[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_VGG*uncStat_VGG)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_VGG[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_VGG[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'BH')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_VGG[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_VGG[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_VGG[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			# if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 			# 	handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 			# 	active = 1
			# 	# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	# plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/all_configs_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/all_configs_VGG_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin))
			# 	plt.clf()


			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH_plus[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH_plus*uncStatInb_BH_plus)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH_plus[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH_plus*uncStatOutb_BH_plus)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb.')

			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH_plus[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH_plus*uncStat_BH_plus)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], ((1+polrate)/2*xsecTh_KM_plus+(1-polrate)/2*xsecTh_KM_minus)[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], ((1+polrate)/2*xsecTh_VGG_plus+(1-polrate)/2*xsecTh_VGG_minus)[xBbin, Q2bin, tbin, :], color = 'orange', label = 'VGG')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH_plus[xBbin, Q2bin, tbin, :], color = 'r', label = 'BH')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma_{+}}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/plus_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()

			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH_minus[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH_minus*uncStatInb_BH_minus)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH_minus[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH_minus*uncStatOutb_BH_minus)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb.')

			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH_minus[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH_minus*uncStat_BH_minus)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], ((1-polrate)/2*xsecTh_KM_plus+(1+polrate)/2*xsecTh_KM_minus)[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], ((1-polrate)/2*xsecTh_VGG_plus+(1+polrate)/2*xsecTh_VGG_minus)[xBbin, Q2bin, tbin, :], color = 'orange', label = 'VGG')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma_{-}}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/minus_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()

			for tbin in range(num_plott):
				active = 0
				ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
				fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
				for xBbin in range(num_plotxB):
					for Q2bin in range(num_plotQ2):
						#skip inactive bins
						if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
							axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
							axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
							continue
						# if ActiveInb[xBbin, Q2bin, tbin, :].any():
						# 	phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						# 	axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], 0.5*(1/polrate)*(xsecInb_BH_plus - xsecInb_BH_minus)[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = np.sqrt((xsecInb_BH_plus*uncStatInb_BH_plus)**2+(xsecInb_BH_minus*uncStatInb_BH_minus)**2)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb.')
						# if ActiveOutb[xBbin, Q2bin, tbin, :].any():
						# 	phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						# 	axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], 0.5*(1/polrate)*(xsecOutb_BH_plus - xsecOutb_BH_minus)[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = np.sqrt((xsecOutb_BH_plus*uncStatOutb_BH_plus)**2+(xsecOutb_BH_minus*uncStatOutb_BH_minus)**2)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb.')

						phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], 0.5*(xsec_BH_plus-xsec_BH_minus)[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = np.sqrt((xsec_BH_plus*uncStat_BH_plus)**2+(xsec_BH_minus*uncStat_BH_minus)**2)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental Data')
						axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], 0.5*(xsecTh_KM_plus - xsecTh_KM_minus)[xBbin, Q2bin, tbin, :], color = 'cyan', label = 'Theory (KM15)')
						# axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], 0.5*(xsecTh_VGG_plus - xsecTh_VGG_minus)[xBbin, Q2bin, tbin, :], color = 'orange', label = 'VGG')

						xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
						Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
						theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
						header = xBheader +Q2header + theader
						axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
						axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma_{pol.}}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
						# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
						axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
						axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
						if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
							handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
							active = 1
				lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
				fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
				plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/pol_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
				plt.clf()

			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				BSAInb_BH = (1/polrate)*divideHist(xsecInb_BH_plus - xsecInb_BH_minus, xsecInb_BH_plus + xsecInb_BH_minus, threshold = -np.inf)
			# 				uncBSAInb_BH = BSAInb_BH*np.sqrt((xsecInb_BH_plus*uncStatInb_BH_plus)**2+(xsecInb_BH_minus*uncStatInb_BH_minus)**2)*np.sqrt(divideHist(1, (xsecInb_BH_plus-xsecInb_BH_minus)**2) + divideHist(1, (xsecInb_BH_plus+xsecInb_BH_minus)**2))
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], BSAInb_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = uncBSAInb_BH[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				BSAOutb_BH = (1/polrate)*divideHist(xsecOutb_BH_plus - xsecOutb_BH_minus, xsecOutb_BH_plus + xsecOutb_BH_minus, threshold = -np.inf)
			# 				uncBSAOutb_BH = BSAOutb_BH*np.sqrt((xsecOutb_BH_plus*uncStatOutb_BH_plus)**2+(xsecOutb_BH_minus*uncStatOutb_BH_minus)**2)*np.sqrt(divideHist(1, (xsecOutb_BH_plus-xsecOutb_BH_minus)**2) + divideHist(1, (xsecOutb_BH_plus+xsecOutb_BH_minus)**2))
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], BSAOutb_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = uncBSAOutb_BH[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb.')

			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			BSA_BH = (1/polrate)*divideHist(xsec_BH_plus - xsec_BH_minus, xsec_BH_plus + xsec_BH_minus, threshold = -np.inf)
			# 			uncBSA_BH = BSA_BH*np.sqrt((xsec_BH_plus*uncStat_BH_plus)**2+(xsec_BH_minus*uncStat_BH_minus)**2)*np.sqrt(divideHist(1, (xsec_BH_plus-xsec_BH_minus)**2) + divideHist(1, (xsec_BH_plus+xsec_BH_minus)**2))
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], BSA_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = uncBSA_BH[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Merged.')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_KM_plus - xsecTh_KM_minus, xsecTh_KM_plus + xsecTh_KM_minus, threshold = -np.inf)[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_VGG_plus - xsecTh_VGG_minus, xsecTh_VGG_plus + xsecTh_VGG_minus, threshold = -np.inf)[xBbin, Q2bin, tbin, :], color = 'orange', label = 'VGG')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel("BSA")
			# 			# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/bsa_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()


			for tbin in range(num_plott):
				active = 0
				ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
				fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
				for xBbin in range(num_plotxB):
					for Q2bin in range(num_plotQ2):
						#skip inactive bins
						if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
							axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
							axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
							continue
						# if ActiveInb[xBbin, Q2bin, tbin, :].any():
						# 	phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						# 	axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb.')
						# if ActiveOutb[xBbin, Q2bin, tbin, :].any():
						# 	phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						# 	axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb.')

						phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental Data')
						axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'cyan', label = 'Theory (KM15)')
						axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'Theory (BH)')

						xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
						Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
						theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
						header = xBheader +Q2header + theader
						axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
						axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
						axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
						axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
						axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
						if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
							handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
							active = 1
				lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
				fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
				plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/unpol_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
				plt.clf()

			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecInb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecOutb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'cyan', label = 'Outb.')
			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsec_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_KM[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :]), color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :]), color = 'r', label = 'BH')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylim([0.01, 100])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/KM_ratio_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()

			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue
			# 			if ActiveInb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecInb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb[xBbin, Q2bin, tbin, :].any():
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecOutb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'cyan', label = 'Outb.')
			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsec_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_KM[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :]), color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :]), color = 'r', label = 'BH')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylim([0.01, 100])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/BH_ratio_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()

			# Normalization = np.zeros(xsecTh_BH.shape[:-1])
			# Normalization_Inb = np.zeros(xsecTh_BH.shape[:-1])
			# Normalization_Outb = np.zeros(xsecTh_BH.shape[:-1])
			# Normalization_KM = np.zeros(xsecTh_BH.shape[:-1])

			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue

			# 			if ActiveInb_int[xBbin, Q2bin, tbin, :].any():
			# 				wings_Inb = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsecInb_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
			# 				Normalization_Inb[xBbin, Q2bin, tbin] = np.mean(xsecInb_BH[xBbin, Q2bin, tbin, wings_Inb], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings_Inb], axis = -1)
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin], linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb_int[xBbin, Q2bin, tbin, :].any():
			# 				wings_Outb = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsecOutb_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
			# 				Normalization_Outb[xBbin, Q2bin, tbin] = np.mean(xsecOutb_BH[xBbin, Q2bin, tbin, wings_Outb], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings_Outb], axis = -1)
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin], linestyle ='', color = 'cyan', label = 'Outb.')
			# 			wings = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsec_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
			# 			Normalization[xBbin, Q2bin, tbin] = np.mean(xsec_BH[xBbin, Q2bin, tbin, wings], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings], axis = -1)
			# 			Normalization_KM[xBbin, Q2bin, tbin] = np.mean(xsecTh_KM[xBbin, Q2bin, tbin, wings], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings], axis = -1)
			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], linestyle ='', color = 'k', label = 'Merged')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'BH')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll_int[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/Normalized_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()

			# IntegratedDiff = np.zeros(xsecTh_BH.shape[:-1])
			# IntegratedDiff_Inb = np.zeros(xsecTh_BH.shape[:-1])
			# IntegratedDiff_Outb = np.zeros(xsecTh_BH.shape[:-1])
			# IntegratedDiff_KM = np.zeros(xsecTh_BH.shape[:-1])

			# for tbin in range(num_plott):
			# 	active = 0
			# 	ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			# 	fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# 	for xBbin in range(num_plotxB):
			# 		for Q2bin in range(num_plotQ2):
			# 			#skip inactive bins
			# 			if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 				continue

			# 			centers = [8, 9, 10, 11, 12, 13, 14, 15]
			# 			if ActiveInb_int[xBbin, Q2bin, tbin, :].any():
			# 				IntegratedDiff_Inb[xBbin, Q2bin, tbin] = np.sum(xsecInb_BH[xBbin, Q2bin, tbin, centers]/Normalization_Inb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
			# 				phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin], linestyle ='', color = 'g', label = 'Inb.')
			# 			if ActiveOutb_int[xBbin, Q2bin, tbin, :].any():
			# 				IntegratedDiff_Outb[xBbin, Q2bin, tbin] = np.sum(xsecOutb_BH[xBbin, Q2bin, tbin, centers]/Normalization_Outb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
			# 				phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
			# 				axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin], linestyle ='', color = 'cyan', label = 'Outb.')
			# 			IntegratedDiff[xBbin, Q2bin, tbin] = np.sum(xsec_BH[xBbin, Q2bin, tbin, centers]/Normalization[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
			# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], linestyle ='', color = 'k', label = 'Merged')
			# 			IntegratedDiff_KM[xBbin, Q2bin, tbin] = np.sum(xsecTh_KM[xBbin, Q2bin, tbin, centers] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], (xsecTh_KM-xsecTh_BH)[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')

			# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
			# 			Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 			theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 			header = xBheader +Q2header + theader
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
			# 			# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
			# 			axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
			# 			if (active == 0) and ActiveAll_int[xBbin, Q2bin, tbin, :].any():
			# 				handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 				active = 1
			# 	lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			# 	fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# 	plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/Differences_bkgscheme{}tbin{}.pdf".format(optionaltag, k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# 	plt.clf()

			# active = 0
			# fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# for xBbin in range(num_plotxB):
			# 	for Q2bin in range(num_plotQ2):
			# 		#skip inactive bins
			# 		if badBinCondxBQ2(xBbin, Q2bin, k):
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 			continue
			# 		tbin = np.argwhere(ActiveInb_int[xBbin, Q2bin, :, 0]).flatten()
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Inb[xBbin, Q2bin, tbin], color = 'g', label = "Inb", marker = 'o', markersize = 20)
			# 		tbin = np.argwhere(ActiveOutb_int[xBbin, Q2bin, :, 0]).flatten()
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Outb[xBbin, Q2bin, tbin], color = 'cyan', label = "Outb.", marker = 'o', markersize = 20)
			# 		tbin = np.argwhere(ActiveAny_int[xBbin, Q2bin, :, 0]).flatten()
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff[xBbin, Q2bin, tbin], color = 'k', label = "Merged", marker = 'o', markersize = 20)
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_KM[xBbin, Q2bin, tbin], color = 'b', label = "KM", marker = 'o', markersize = 20)

			# 		xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0].mean())
			# 		Q2header = r"$<Q^2>=$"+" {:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0].mean())
			# 		# theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 		header = xBheader +Q2header
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 		# axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$"))
			# 		# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_xlim([0, 2])
			# 		# axs[num_plotQ2-Q2bin-1, xBbin].set_ylim([0.2, 1.])
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 0.5, 1, 1.5, 2])
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$|t|$" + " [" + GeV2 + "]")
			# 		if (active == 0) and ActiveAll_int[xBbin, Q2bin, :, 0].any():
			# 			handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 			active = 1
			# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = r"$\int (d\sigma-d\sigma_{BH})$", bbox_to_anchor = (1.0, 0.6))
			# fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/Integrals_bkgscheme{}.pdf".format(optionaltag, k, i), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# # plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/NormScale_bkgscheme{}.pdf".format(optionaltag, k, i, tbin), bbox_inches = 'tight')
			# plt.clf()

			# active = 0
			# fig, axs = plt.subplots(num_plott, num_plotxB, figsize = (7.5*(num_plotxB), 7.5*(num_plott)))
			# for xBbin in range(num_plotxB):
			# 	for tbin in range(num_plott):
			# 		#skip inactive bins
			# 		if badBinCondxBt(xBbin, tbin, k):
			# 			axs[num_plott-tbin-1 , xBbin].yaxis.set_visible(False)
			# 			axs[num_plott-tbin-1 , xBbin].xaxis.set_visible(False)
			# 			continue

			# 		Q2bin = np.argwhere(ActiveInb_int[xBbin, :, tbin, 0]).flatten()
			# 		axs[num_plott-tbin-1 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Inb[xBbin, Q2bin, tbin], color = 'g', label = "Inb.", marker = 'o', markersize = 20)
			# 		Q2bin = np.argwhere(ActiveOutb_int[xBbin, :, tbin, 0]).flatten()
			# 		axs[num_plott-tbin-1 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Outb[xBbin, Q2bin, tbin], color = 'cyan', label = "Outb.", marker = 'o', markersize = 20)
			# 		Q2bin = np.argwhere(ActiveAny_int[xBbin, :, tbin, 0]).flatten()
			# 		axs[num_plott-tbin-1 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff[xBbin, Q2bin, tbin], color = 'k', label = "Merged", marker = 'o', markersize = 20)
			# 		axs[num_plott-tbin-1 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_KM[xBbin, Q2bin, tbin], color = 'b', label = "KM", marker = 'o', markersize = 20)

			# 		xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0].mean())
			# 		# Q2header = r"$<Q^2>=$"+" {:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
			# 		theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0].mean())
			# 		header = xBheader +theader
			# 		axs[num_plott-tbin-1, xBbin].set_title(header, fontsize = 20)
			# 		# axs[num_plott-tbin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$"))
			# 		# axs[num_plott-tbin-1, xBbin].set_yscale('log')
			# 		# axs[num_plott-tbin-1, xBbin].set_xlim([0, 2])
			# 		# axs[num_plott-tbin-1, xBbin].set_ylim([0.2, 1.])
			# 		# axs[num_plott-tbin-1, xBbin].set_xticks([0, 0.5, 1, 1.5, 2])
			# 		axs[num_plott-tbin-1, xBbin].set_xlabel(r"$Q^2$" + " [" + GeVc2 + "]")
			# 		axs[num_plott-tbin-1, xBbin].set_xlim([1, 5])
			# 		axs[num_plott-tbin-1, xBbin].set_xticks([1, 2, 3, 4, 5])
			# 		if (active == 0) and ActiveAll_int[xBbin, Q2bin, :, 0].any():
			# 			handles, labels = axs[num_plott-tbin-1, xBbin].get_legend_handles_labels()
			# 			active = 1
			# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = r"$\int (d\sigma-d\sigma_{BH})$", bbox_to_anchor = (1.0, 0.6))
			# fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/Integrals_bkgscheme{}_inQ2.pdf".format(optionaltag, k, i), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# # plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/NormScale_bkgscheme{}.pdf".format(optionaltag, k, i, tbin), bbox_inches = 'tight')
			# plt.clf()

			# active = 0
			# fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
			# for xBbin in range(num_plotxB):
			# 	for Q2bin in range(num_plotQ2):
			# 		#skip inactive bins
			# 		if badBinCondxBQ2(xBbin, Q2bin, k):
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
			# 			axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
			# 			continue
			# 		tbin = np.argwhere(ActiveInb_int[xBbin, Q2bin, :, 0]).flatten()
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization_Inb[xBbin, Q2bin, tbin], color = 'g', label = "Inb.")
			# 		tbin = np.argwhere(ActiveOutb_int[xBbin, Q2bin, :, 0]).flatten()
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization_Outb[xBbin, Q2bin, tbin], color = 'cyan', label = "Outb.")
			# 		tbin = np.argwhere(ActiveAny_int[xBbin, Q2bin, :, 0]).flatten()
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization[xBbin, Q2bin, tbin], color = 'r', label = "Inb.")
			# 		axs[num_plotQ2-Q2bin-1 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization_KM[xBbin, Q2bin, tbin], color = 'b', label = "KM")

			# 		# axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, :], Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, :]-phibins[:-1], phibins[1:]-phi1avg_BH[xBbin, Q2bin, tbin, :]], yerr = 0, linestyle ='', color = 'cyan', label = 'Outb.')

			# 		# axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, :], Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, :]-phibins[:-1], phibins[1:]-phi1avg_BH[xBbin, Q2bin, tbin, :]], yerr = 0, linestyle ='', color = 'k', label = 'Merged')

			# 		xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0].mean())
			# 		Q2header = r"$<Q^2>=$"+" {:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0].mean())
			# 		# theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
			# 		header = xBheader +Q2header
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
			# 		# axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$"))
			# 		# axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_xlim([0, 2])
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_ylim([0.2, 1.])
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 0.5, 1, 1.5, 2])
			# 		axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$|t|$" + " [" + GeV2 + "]")
			# 		if (active == 0) and ActiveAll_int[xBbin, Q2bin, :, 0].any():
			# 			handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
			# 			active = 1
			# lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = "Norm. to BH", bbox_to_anchor = (1.0, 0.6))
			# fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			# plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/NormScale_bkgscheme{}.pdf".format(optionaltag, k, i), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			# # plt.savefig("/volatile/clas12/sangbaek/clas12DVCS/plots{}/binscheme{}/NormScale_bkgscheme{}.pdf".format(optionaltag, k, i, tbin), bbox_inches = 'tight')
			# plt.clf()


#not actively used
if args.contplot:
	print("read exp...")
	epgExp = pd.read_pickle("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/epgExp.pkl")

	k = 3
	xBbins  = collection_xBbins[k]
	Q2bins  = collection_Q2bins[k]
	tbins   = collection_tbins [k]
	phibins = collection_phibins[k]

	histExpInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	histExpInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	histExpInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

	histExpOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	histExpOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	histExpOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

	histExpInb = histExpInbFD + histExpInbCD + histExpInbCDFT
	histExpOutb = histExpOutbFD + histExpOutbCD + histExpOutbCDFT

	histBHDVCSInb, histBHDVCSInbFD, histBHDVCSInbCD, histBHDVCSInbCDFT = {}, {}, {}, {}
	histBHDVCSOutb, histBHDVCSOutbFD, histBHDVCSOutbCD, histBHDVCSOutbCDFT = {}, {}, {}, {}
	for i in range(3):

	    histBHDVCSInbFD[i], bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)])
	    histBHDVCSInbCD[i], bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)])
	    histBHDVCSInbCDFT[i], bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)])

	    histBHDVCSInb[i] = histBHDVCSInbFD[i] + histBHDVCSInbCD[i] + histBHDVCSInbCDFT[i]
	    
	    histBHDVCSOutbFD[i], bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)])
	    histBHDVCSOutbCD[i], bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)])
	    histBHDVCSOutbCDFT[i], bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)])

	    histBHDVCSOutb[i] = histBHDVCSOutbFD[i] + histBHDVCSOutbCD[i] + histBHDVCSOutbCDFT[i]

	plt.rcParams["figure.figsize"] = (10,6)
	plt.rcParams['legend.title_fontsize'] = 'small'

	xBbin = 5
	Q2bin = 2
	tbin = 2

	plt.hist(phibins[:-1], phibins, weights = histExpInb[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'k', label = "Sig+Bkg")
	plt.hist(phibins[:-1], phibins, weights = histBHDVCSInb[1][xBbin, Q2bin, tbin,:], histtype = 'step', color = 'r', label = "Sig")
	plt.hist(phibins[:-1], phibins, weights = histExpInb[xBbin, Q2bin, tbin, :] - histBHDVCSInb[1][xBbin, Q2bin, tbin,:], histtype = 'stepfilled', color = 'brown', alpha = 0.5, label = "Bkg(1)")
	plt.hist(phibins[:-1], phibins, weights = histExpInb[xBbin, Q2bin, tbin, :] - histBHDVCSInb[2][xBbin, Q2bin, tbin,:], histtype = 'stepfilled', color ='b', alpha = 0.5, label = "Bkg(2)")


	plt.xlim([90, 270])
	plt.ylim([0, 70])
	# plt.yscale('log')
	plt.xticks([90, 180, 270])
	plt.xlabel(r"$\phi$" + " ["+degree+"]")

	xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}".format(xBbins[xBbin+1]) + "\n"
	Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f} ".format(Q2bins[Q2bin+1])+ "\n"
	theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f} ".format(tbins[tbin+1])
	header = xBheader + Q2header + theader
	leg = plt.legend(loc = 'upper right', title = header, ncol = 2)
	plt.savefig("plots/contamination{}{}{}.pdf".format(xBbin, Q2bin, tbin))

	xBbin = 3
	Q2bin = 2
	tbin = 2

	plt.hist(phibins[:-1], phibins, weights = histExpInb[xBbin, Q2bin, tbin,:], histtype = 'step', color = 'k', label = "Sig+Bkg")
	plt.hist(phibins[:-1], phibins, weights = histBHDVCSInb[1][xBbin, Q2bin, tbin,:], histtype = 'step', color = 'r', label = "Sig")
	plt.hist(phibins[:-1], phibins, weights = histExpInb[xBbin, Q2bin, tbin, :] - histBHDVCSInb[1][xBbin, Q2bin, tbin,:], histtype = 'stepfilled', color = 'brown', alpha = 0.5, label = "Bkg(1)")
	plt.hist(phibins[:-1], phibins, weights = histExpInb[xBbin, Q2bin, tbin, :] - histBHDVCSInb[2][xBbin, Q2bin, tbin,:], histtype = 'stepfilled', color ='b', alpha = 0.5, label = "Bkg(2)")


	plt.xlim([90, 270])
	plt.ylim([0, 70])
	# plt.yscale('log')
	plt.xticks([90, 180, 270])
	plt.xlabel(r"$\phi$" + " ["+degree+"]")

	xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}".format(xBbins[xBbin+1]) + "\n"
	Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f} ".format(Q2bins[Q2bin+1])+ "\n"
	theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f} ".format(tbins[tbin+1])
	header = xBheader + Q2header + theader
	leg = plt.legend(loc = 'upper right', title = header, ncol = 2)
	plt.savefig("plots/contamination{}{}{}.pdf".format(xBbin, Q2bin, tbin))

#not actively used
if args.radplot:
	print("read rad factors...")

	k = 3

	xBbins  = collection_xBbins[k]
	Q2bins  = collection_Q2bins[k]
	tbins   = collection_tbins [k]
	phibins = collection_phibins[k]

	xsecTh_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_BH.npz".format(k))["hist"]
	xsecTh_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xsecTh_VGG.npz".format(k))["hist"]

	integratedRad_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_VGG.npz".format(k))["hist"]
	integratedBorn_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_VGG.npz".format(k))["hist"]
	rcfactors_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_VGG.npz".format(k))["hist"]

	integratedRad_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedRad_BH.npz".format(k))["hist"]
	integratedBorn_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/integratedBorn_BH.npz".format(k))["hist"]
	rcfactors_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/rcfactors_BH.npz".format(k))["hist"]

	phi1avg_BH = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_BH.npz".format(k))["hist"]
	xBavg_BH   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_BH.npz".format(k))["hist"]
	Q2avg_BH   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_BH.npz".format(k))["hist"]
	t1avg_BH   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_BH.npz".format(k))["hist"]

	phi1avg_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/phi1avg_VGG.npz".format(k))["hist"]
	xBavg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/xBavg_VGG.npz".format(k))["hist"]
	Q2avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/Q2avg_VGG.npz".format(k))["hist"]
	t1avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/t1avg_VGG.npz".format(k))["hist"]

	xBbin = 3
	Q2bin = 2
	tbin = 2

	plt.scatter(phi1avg_VGG[xBbin, Q2bin, tbin, :], divideHist(integratedBorn_BH,xsecTh_BH)[xBbin, Q2bin, tbin, :] , color = 'k')
	plt.xlim([0, 360])
	plt.xticks([0, 90, 180, 270, 360])
	plt.xlabel(r"$\phi$" + " ["+degree+"]")
	plt.ylim([0.9, 1.1])
	plt.ylabel(r"$\frac{<d\sigma>}{d\sigma}$", rotation = 0, labelpad=40, fontsize = 40)

	xBheader = "{:.3f} ".format(xBbins[xBbin])+r"$<~~~~~~~~~~x_B~~~~~~~~~<$"+ " {:.3f}, ".format(xBbins[xBbin+1]) +r"$~<x_B>=$"+ "{:.3f}\n".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
	Q2header = "{:.3f} ".format(Q2bins[Q2bin])+ r"$<Q^2/(1~(\mathrm{GeV/c})^2<$"+ " {:.3f}, ".format(Q2bins[Q2bin+1])+ r"$~<Q^2>=$"+"{:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"$~(\mathrm{GeV/c})^2$"+ "\n"
	theader = "{:.3f} ".format(tbins[tbin])+ r"$<~~|t|/(1~\mathrm{GeV}^2)~~~<$"+ " {:.3f}, ".format(tbins[tbin+1]) + r"$~<|t|>=$"+"{:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"$~\mathrm{GeV}^2$"
	header = xBheader + Q2header + theader
	plt.title(header, loc = 'left')
	plt.savefig("plots/finitebin{}{}{}.pdf".format(xBbin, Q2bin, tbin))


# if args.readonly: 
# 	phi1avg_VGG, xBavg_VGG, Q2avg_VGG, t1avg_VGG = {}, {}, {}, {}  
# 	phi1avg_BH , xBavg_BH , Q2avg_BH,  t1avg_BH  = {}, {}, {}, {} 
# 	xsecInb_VGG, xsecInb_BH, xsecOutb_VGG, xsecOutb_BH, xsec_VGG, xsec_BH = {}, {}, {}, {}, {}, {}
# 	uncStatInb_VGG, uncStatInb_BH, uncStatOutb_VGG, uncStatOutb_BH, uncStat_VGG, uncStat_BH = {}, {}, {}, {}, {}, {}
# 	xsecTh_KM, xsecTh_BH, xsecTh_VGG = {}, {}, {}     
# binVolume      
# ActiveAll      
# ActiveAny      
# ActiveInb      
# ActiveOutb     
# ActiveAll_int  
# ActiveAny_int  
# ActiveInb_int  
# ActiveOutb_int 
# 	for k in range(kstart, kend):
# 		for i in range(istart, len(collection_cont_xBbins)):
# 			phi1avg_VGG [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}phi1avg_VGG.npz".format(k, i))["hist"]
# 			xBavg_VGG   [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xBavg_VGG.npz".format(k, i))["hist"]
# 			Q2avg_VGG   [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}Q2avg_VGG.npz".format(k, i))["hist"]
# 			t1avg_VGG   [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}t1avg_VGG.npz".format(k, i))["hist"]

# 			phi1avg_BH  [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}phi1avg_BH.npz".format(k, i))["hist"]
# 			xBavg_BH    [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xBavg_BH.npz".format(k, i))["hist"]
# 			Q2avg_BH    [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}Q2avg_BH.npz".format(k, i))["hist"]
# 			t1avg_BH    [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}t1avg_BH.npz".format(k, i))["hist"]

# 			xsecInb_VGG [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(k, i))["hist"]
# 			xsecInb_BH  [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(k, i))["hist"]
# 			xsecOutb_VGG [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(k, i))["hist"]
# 			xsecOutb_BH  [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(k, i))["hist"]
# 			xsec_VGG     [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsec_VGG.npz".format(k, i))["hist"]
# 			xsec_BH      [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsec_BH.npz".format(k, i))["hist"]

# 			uncStatInb_VGG [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(k, i))["hist"]
# 			uncStatInb_BH  [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(k, i))["hist"]
# 			uncStatOutb_VGG [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(k, i))["hist"]
# 			uncStatOutb_BH  [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(k, i))["hist"]
# 			uncStat_VGG     [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(k, i))["hist"]
# 			uncStat_BH      [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStat_BH.npz".format(k, i))["hist"]

# 			xsecTh_KM          [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecTh_KM.npz".format(k, i))["hist"]
# 			xsecTh_BH          [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecTh_BH.npz".format(k, i))["hist"]
# 			xsecTh_VGG         [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecTh_VGG.npz".format(k, i))["hist"]
# 			binVolume          [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}binVolume.npz".format(k, i))["hist"]

# 			ActiveAll       [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAll.npz".format(k, i))["hist"]
# 			ActiveAny       [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAny.npz".format(k, i))["hist"]
# 			ActiveInb          [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveInb.npz".format(k, i))["hist"]
# 			ActiveOutb         [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveOutb.npz".format(k, i))["hist"]

# 			ActiveAll_int       [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(k, i))["hist"]
# 			ActiveAny_int       [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(k, i))["hist"]
# 			ActiveInb_int          [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(k, i))["hist"]
# 			ActiveOutb_int         [k, i] = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(k, i))["hist"]

# for k in range(2, kend):

# 	xBbins  = collection_xBbins[k]
# 	Q2bins  = collection_Q2bins[k]
# 	tbins   = collection_tbins [k]
# 	phibins = collection_phibins[k]

# 	os.makedirs("plots/binscheme{}".format(k), exist_ok = True)

# 	for i in range(2, len(collection_cont_xBbins)):

# 		print("reading the xsec vars")
# 		phi1avg_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}phi1avg_VGG.npz".format(k, i))["hist"]
# 		xBavg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xBavg_VGG.npz".format(k, i))["hist"]
# 		Q2avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}Q2avg_VGG.npz".format(k, i))["hist"]
# 		t1avg_VGG   = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}t1avg_VGG.npz".format(k, i))["hist"]

# 		phi1avg_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}phi1avg_BH.npz".format(k, i))["hist"]
# 		xBavg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xBavg_BH.npz".format(k, i))["hist"]
# 		Q2avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}Q2avg_BH.npz".format(k, i))["hist"]
# 		t1avg_BH    = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}t1avg_BH.npz".format(k, i))["hist"]

# 		xsecInb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(k, i))["hist"]
# 		xsecInb_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(k, i))["hist"]
# 		xsecOutb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(k, i))["hist"]
# 		xsecOutb_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(k, i))["hist"]
# 		xsec_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsec_VGG.npz".format(k, i))["hist"]
# 		xsec_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsec_BH.npz".format(k, i))["hist"]

# 		uncStatInb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(k, i))["hist"]
# 		uncStatInb_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(k, i))["hist"]
# 		uncStatOutb_VGG = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(k, i))["hist"]
# 		uncStatOutb_BH  = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(k, i))["hist"]
# 		uncStat_VGG     = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(k, i))["hist"]
# 		uncStat_BH      = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}uncStat_BH.npz".format(k, i))["hist"]

# 		xsecTh_KM          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecTh_KM.npz".format(k, i))["hist"]
# 		xsecTh_BH          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecTh_BH.npz".format(k, i))["hist"]
# 		xsecTh_VGG         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}xsecTh_VGG.npz".format(k, i))["hist"]
# 		binVolume          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}binVolume.npz".format(k, i))["hist"]

# 		ActiveAll       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAll.npz".format(k, i))["hist"]
# 		ActiveAny       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAny.npz".format(k, i))["hist"]
# 		ActiveInb          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveInb.npz".format(k, i))["hist"]
# 		ActiveOutb         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveOutb.npz".format(k, i))["hist"]

# 		ActiveAll_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(k, i))["hist"]
# 		ActiveAny_int       = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(k, i))["hist"]
# 		ActiveInb_int          = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(k, i))["hist"]
# 		ActiveOutb_int         = np.load("/volatile/clas12/sangbaek/clas12DVCS/nphistograms/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(k, i))["hist"]

# 		print("plotting...")

# 		def badBinCondxBQ2t(xBbin, Q2bin, tbin, k = 0):
# 			if k ==0:
# 				return (xBbin==1 and Q2bin == 0) or (xBbin==0 and Q2bin==4) or (tbin==0 and xBbin==1)
# 			else:
# 				return ~ActiveAny_int[xBbin, Q2bin, tbin, :].any()

# 		def badBinCondxBQ2(xBbin, Q2bin, k = 0):
# 			if k ==0:
# 				return (xBbin==1 and Q2bin == 0) or (xBbin==0 and Q2bin==4)
# 			else:
# 				return ~ActiveAny_int[xBbin, Q2bin, :, :].any()

# 		def badBinCondxBt(xBbin, tbin, k = 0):
# 			if k ==0:
# 				return (xBbin==1 and tbin == 0)
# 			else:
# 				return ~ActiveAny_int[xBbin, :, tbin, :].any()

# 		num_plotQ2 = len(Q2bins) - 1
# 		num_plotxB = len(xBbins) - 1
# 		num_plott = len(tbins) - 1

# 		# if k == 2:
# 		# 	num_plotQ2 = 4
# 		# 	num_plotxB = 4

# 		Normalization = np.zeros(xsecTh_BH.shape[:-1])
# 		Normalization_Inb = np.zeros(xsecTh_BH.shape[:-1])
# 		Normalization_Outb = np.zeros(xsecTh_BH.shape[:-1])
# 		Normalization_KM = np.zeros(xsecTh_BH.shape[:-1])
# 		order = [2, 0, 1]
# 		for tbin in range(num_plott):
# 			active = 0
# 			ttitle = "{:.3f} GeV".format(tbins[tbin])+r"${}^{2}<|t|<$"+"{:.3f} GeV".format(tbins[tbin+1])+r"${}^{2}$"
# 			fig, axs = plt.subplots(num_plotQ2, num_plotxB, figsize = (7.5*(num_plotxB), 6*(num_plotQ2)))
# 			for xBbin in range(num_plotxB):
# 				for Q2bin in range(num_plotQ2):
# 					#skip inactive bins
# 					if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
# 						axs[num_plotQ2-Q2bin-1 , xBbin].yaxis.set_visible(False)
# 						axs[num_plotQ2-Q2bin-1 , xBbin].xaxis.set_visible(False)
# 						continue
# 					wings = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsec_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
# 					Normalization[xBbin, Q2bin, tbin] = np.mean(xsec_BH[xBbin, Q2bin, tbin, wings], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings], axis = -1)
# 					Normalization_KM[xBbin, Q2bin, tbin] = np.mean(xsecTh_KM[xBbin, Q2bin, tbin, wings], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings], axis = -1)
# 					phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
# 					axs[num_plotQ2-Q2bin-1 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], linestyle ='', color = 'k', label = 'Experimental data')
# 					axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'cyan', label = 'Theory (BH+Int.+DVCS'+r"${}^{2}$"+")")
# 					axs[num_plotQ2-Q2bin-1 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'Theory (Pure BH)')

# 					xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
# 					Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
# 					theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
# 					header = xBheader +Q2header + theader
# 					axs[num_plotQ2-Q2bin-1, xBbin].set_title(header, fontsize = 20)
# 					axs[num_plotQ2-Q2bin-1, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")
# 					axs[num_plotQ2-Q2bin-1, xBbin].set_yscale('log')
# 					axs[num_plotQ2-Q2bin-1, xBbin].set_xticks([0, 90, 180, 270, 360])
# 					axs[num_plotQ2-Q2bin-1, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
# 					if (active == 0) and ActiveAll_int[xBbin, Q2bin, tbin, :].any():
# 						handles, labels = axs[num_plotQ2-Q2bin-1, xBbin].get_legend_handles_labels()
# 						active = 1
# 			lgd = plt.figlegend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
# 			fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
# 			plt.savefig("plots/richard_rolf_tbin{}.pdf".format(tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
# 			plt.clf()

# 		# xBbin, Q2bin, tbin = (3, 2, 2)

# 		# fig, axs = plt.subplots(1, 1, figsize = (10, 6))
# 		# phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()		
# 		# P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
# 		# P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])

# 		# axs.errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental data')
# 		# res_lsq = curve_fit(FourierSeries, [0, 0, 0], args=(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin]))
# 		# N = 40
# 		# phi1s = np.linspace(0, 360, N)[:-1]
# 		# phi1s = phi1s + np.diff(phi1s)[0]
# 		# xBs = np.ones(len(phi1s))*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
# 		# Q2s = np.ones(len(phi1s))*Q2avg_BH[xBbin, Q2bin, tbin, phibin][0]
# 		# t1s = np.ones(len(phi1s))*t1avg_BH[xBbin, Q2bin, tbin, phibin][0]

# 		# axs.plot(phi1s, FourierSeries(res_lsq.x, phi1s), label = 'Fitting results', color = 'k', linestyle = '--')
# 		# P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])
# 		# P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])

# 		# axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, :], P1b*P2b*getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 1), color = 'r', label = 'Theory (Pure BH)')
# 		# axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, :], P1b*P2b*getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 5), color = 'cyan', label = 'Theory (BH+Int.+DVCS'+r"${}^{2}$"+")")
# 		# axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, :], P1b*P2b*getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 2), color = 'g', label = 'Interference')
# 		# axs.plot(phi1avg_BH[xBbin, Q2bin, tbin, :], P1b*P2b*getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 3), color = 'orange', label = 'DVCS'+r"${}^{2}$")
# 		# xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
# 		# Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
# 		# theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
# 		# header = xBheader +Q2header + theader
# 		# axs.set_xlim([0, 360])
# 		# axs.set_title(header, fontsize = 20)
# 		# axs.set_ylabel(r"$P_1(\phi)P_2(\phi)\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")
# 		# # axs.set_yscale('log')
# 		# axs.set_xticks([0, 90, 180, 270, 360])
# 		# axs.set_xlabel(r"$\phi$" + " [" + degree + "]")
# 		# handles, labels = axs.get_legend_handles_labels()
# 		# order = [5, 0, 1, 2, 3, 4]

# 		# lgd = plt.figlegend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'upper left', bbox_to_anchor =(1.1, 0.9), title = "Reduced Cross Sections")
# 		# plt.savefig("plots/richard_rolf.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')

# 		# xBbin, Q2bin, tbin = (3, 2, 2)
# 		tbin = 2
# 		xBbin = 3

# 		fig, axs = plt.subplots(4, 1, figsize = (12, 20))
# 		for Q2bin in range(4):

# 			print(Q2bin)

# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()		
# 			P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
# 			P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])

# 			axs[4-Q2bin-1].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], (xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental data')
# 			res_lsq = curve_fit(FourierSeries, [0, 0, 0], args=(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin]))
# 			N = 40
# 			phi1s = np.linspace(0, 360, N)[:-1]
# 			phi1s = phi1s + np.diff(phi1s)[0]
# 			xBs = np.ones(len(phi1s))*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
# 			Q2s = np.ones(len(phi1s))*Q2avg_BH[xBbin, Q2bin, tbin, phibin][0]
# 			t1s = np.ones(len(phi1s))*t1avg_BH[xBbin, Q2bin, tbin, phibin][0]
# 			P1b = P1(xBs, Q2s, t1s, phi1s)
# 			P2b = P2(xBs, Q2s, t1s, phi1s)

# 			axs[4-Q2bin-1].plot(phi1s, 1/P1b*1/P2b*FourierSeries(res_lsq.x, phi1s), label = 'Fitting results', color = 'k', linestyle = '--')
# 			P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])
# 			P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])

# 			axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 5), color = 'cyan', label = 'Theory (BH+Int.+DVCS'+r"${}^{2}$"+")")
# 			axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 1), color = 'r', label = 'Theory (Pure BH)')
# 			# axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 2), color = 'g', label = 'Interference')
# 			# axs[4-Q2bin-1].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 3), color = 'purple', label = 'DVCS'+r"${}^{2}$")
# 			xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
# 			Q2header = r"$,~<Q^2>=$"+" {:.3f} (GeV/c)".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"${}^{2}$"
# 			theader = r"$,~<|t|>=$"+" {:.3f} GeV".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"${}^{2}$"
# 			header = xBheader +Q2header + theader
# 			axs[4-Q2bin-1].set_xlim([0, 360])
# 			axs[4-Q2bin-1].set_ylim([0.005, 2])
# 			axs[4-Q2bin-1].set_title(header, fontsize = 20)
# 			axs[4-Q2bin-1].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")
# 			axs[4-Q2bin-1].set_yscale('log')
# 			axs[4-Q2bin-1].set_xticks([0, 90, 180, 270, 360])
# 			axs[4-Q2bin-1].set_xlabel(r"$\phi$" + " [" + degree + "]")
# 			if Q2bin == 0:
# 				handles, labels = axs[4-Q2bin-1].get_legend_handles_labels()
# 				order = [3, 0, 1, 2]

# 		lgd = plt.figlegend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'upper left', bbox_to_anchor =(1.02, 0.5), title = "")
# 		fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
# 		plt.savefig("plots/richard_rolf_Q2.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
# 		# plt.savefig("plots/richard_rolf{}.pdf".format(Q2bin))
# 		plt.clf()

# 		Q2bin = 2
# 		xBbin = 3

# 		fig, axs = plt.subplots(4, 1, figsize = (12, 20))
# 		for tbin in range(1, 5):

# 			print(tbin)

# 			phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()		
# 			P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])
# 			P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, phibin], Q2avg_BH[xBbin, Q2bin, tbin, phibin], t1avg_BH[xBbin, Q2bin, tbin, phibin], phi1avg_BH[xBbin, Q2bin, tbin, phibin])

# 			axs[4-tbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], (xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Experimental data')
# 			res_lsq = curve_fit(FourierSeries, [0, 0, 0], args=(phi1avg_BH[xBbin, Q2bin, tbin, phibin], P1b*P2b*(xsec_BH/Normalization[xBbin, Q2bin, tbin])[xBbin, Q2bin, tbin, phibin]))
# 			N = 40
# 			phi1s = np.linspace(0, 360, N)[:-1]
# 			phi1s = phi1s + np.diff(phi1s)[0]
# 			xBs = np.ones(len(phi1s))*xBavg_BH[xBbin, Q2bin, tbin, phibin][0]
# 			Q2s = np.ones(len(phi1s))*Q2avg_BH[xBbin, Q2bin, tbin, phibin][0]
# 			t1s = np.ones(len(phi1s))*t1avg_BH[xBbin, Q2bin, tbin, phibin][0]
# 			P1b = P1(xBs, Q2s, t1s, phi1s)
# 			P2b = P2(xBs, Q2s, t1s, phi1s)

# 			axs[4-tbin].plot(phi1s, 1/P1b*1/P2b*FourierSeries(res_lsq.x, phi1s), label = 'Fitting results', color = 'k', linestyle = '--')
# 			P1b = P1(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])
# 			P2b = P2(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :])

# 			axs[4-tbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 5), color = 'cyan', label = 'Theory (BH+Int.+DVCS'+r"${}^{2}$"+")")
# 			axs[4-tbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 1), color = 'r', label = 'Theory (Pure BH)')
# 			# axs[4-tbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 2), color = 'g', label = 'Interference')
# 			# axs[4-tbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], getBHDVCS(xBavg_BH[xBbin, Q2bin, tbin, :], Q2avg_BH[xBbin, Q2bin, tbin, :], t1avg_BH[xBbin, Q2bin, tbin, :], phi1avg_BH[xBbin, Q2bin, tbin, :], mode = 3), color = 'purple', label = 'DVCS'+r"${}^{2}$")
# 			xBheader = r"$<x_B>=$"+" {:.3f}".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
# 			Q2header = r"$,~<Q^2>=$"+" {:.3f} (GeV/c)".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])+r"${}^{2}$"
# 			theader = r"$,~<|t|>=$"+" {:.3f} GeV".format(t1avg_BH[xBbin, Q2bin, tbin, 0])+r"${}^{2}$"
# 			header = xBheader +Q2header + theader
# 			axs[4-tbin].set_xlim([0, 360])
# 			axs[4-tbin].set_ylim([0.005, 2])
# 			axs[4-tbin].set_title(header, fontsize = 20)
# 			axs[4-tbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + " [nb/GeV"+r"$^4$"+"]")
# 			axs[4-tbin].set_yscale('log')
# 			axs[4-tbin].set_xticks([0, 90, 180, 270, 360])
# 			axs[4-tbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
# 			if tbin == 0:
# 				handles, labels = axs[4-tbin].get_legend_handles_labels()
# 				order = [3, 0, 1, 2]

# 		lgd = plt.figlegend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'upper left', bbox_to_anchor =(1.02, 0.5), title = "")
# 		fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
# 		plt.savefig("plots/richard_rolf_t.pdf", bbox_extra_artists=[lgd], bbox_inches = 'tight')
# 		# plt.savefig("plots/richard_rolf{}.pdf".format(Q2bin))
# 		plt.clf()
