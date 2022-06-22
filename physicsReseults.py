#!/usr/bin/env python3
"""
Main Script to save the cross sections
"""
import pandas as pd
import numpy as np
import gc
import matplotlib.pyplot as plt
from copy import copy
cmap = copy(plt.cm.get_cmap("jet"))
from scipy.optimize import least_squares
from utils.const import *
from matplotlib.colors import LogNorm
import argparse
from gepard.fits import th_KM15
import gepard as g
import os, subprocess

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

def printKMarray(xBarray, Q2array, tarray, phiarray, **kwargs):
    BHarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi in zip(xBarray, Q2array, tarray, phiarray):
        BHarray.append(printKM(xB, Q2, t, phi, **kwargs))
    return np.array(BHarray)

def printKM(xB, Q2, t, phi, frame = 'trento', pol = 0 ):
	phi = np.pi - phi
	pt1 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phi,
					process='ep2epgamma', exptype='fixed target', frame =frame,
					in1energy=10.604, in1charge=-1, in1polarization=pol)
	#     print(pt1.frame, pol)
	#     pt2 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phi,
	#                    process='ep2epgamma', exptype='fixed target', frame = 'trento',
	#                    in1energy=10.604, in1charge=-1, in1polarization=+1)
	try:
		return th_KM15.XS(pt1)
	except:
		print(xB, Q2, t, phi)
		return 0

def printVGGarray(xBarray, Q2array, tarray, phiarray, **kwargs):
    VGGarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi in zip(xBarray, Q2array, tarray, phiarray):
        VGGarray.append(printVGG(xB, Q2, t, phi, **kwargs))
    return VGGarray

def printVGG(xB, Q2, t, phi, globalfit = True):
	my_env = os.environ.copy()
	my_env["PATH"] = "/Users/sangbaek/CLAS12/dvcs/print:" + my_env["PATH"]
	my_env["CLASDVCS_PDF"] = "/Users/sangbaek/CLAS12/dvcs/print"
	if globalfit:
		dstot = subprocess.check_output(['/home/sangbaek/printDVCSBH/dvcsgen', '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '3', '--phi', str(phi), '--gpd', '101', '--globalfit'], env = my_env)
	else:
		dstot = subprocess.check_output(['/home/sangbaek/printDVCSBH/dvcsgen', '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '3', '--phi', str(phi), '--gpd', '101'], env = my_env)
	try:
		dstot = float(dstot.splitlines()[-1].decode("utf-8"))
		return dstot
	except:
		print(xB, Q2, t, phi)
		return 0

def printBHarray(xBarray, Q2array, tarray, phiarray, **kwargs):
    BHarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi in zip(xBarray, Q2array, tarray, phiarray):
        BHarray.append(printBHonly(xB, Q2, t, phi, **kwargs))
    return np.array(BHarray)

def printBHonly(xB, Q2, t, phi, globalfit = True):
	if globalfit:
		dstot = subprocess.check_output(['/home/sangbaek/printDVCSBH/dvcsgen', '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '1', '--phi', str(phi), '--globalfit'])
	else:
		dstot = subprocess.check_output(['/home/sangbaek/printDVCSBH/dvcsgen', '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '1', '--phi', str(phi)])
	try:
		dstot = float(dstot.splitlines()[-1].decode("utf-8"))
		return dstot
	except:
		print(xB, Q2, t, phi)
		return 0

def nphistmean(hist, bins):
    s=0
    for i in range(len(hist)):
        s += hist[i] * ((bins[i] + bins[i+1]) / 2) 
    mean = s / np.sum(hist)
    return mean

def createBinEdges(binCenters):
    start = binCenters[0] - np.diff(binCenters)[0]/2
    end = binCenters[-1] + np.diff(binCenters)[-1]/2
    middle = binCenters[:-1] + np.diff(binCenters)/2
    return np.array([start, *middle, end])

def makeReduced(df):
    columns_needed = ["polarity", "config", "beamCurrent", "xB", "Q2", "t1", "phi1"]
    return df.loc[:, columns_needed]

def readReduced(parent, jobNum, polarity, beamCurrent):
    df = pd.read_pickle(parent + "{}.pkl".format(jobNum))
    df.loc[:, "polarity"] = polarity
    df.loc[:, "beamCurrent"] = beamCurrent
    columns_needed = ["polarity", "config", "beamCurrent", "xB", "Q2", "t1", "phi1"]
    return df.loc[:, columns_needed]

def divideHist(df1, df2, threshold = 0):
	return np.divide(df1, df2, where = (df2!=0) & (df1>threshold), out = np.zeros(df2.shape, dtype = float))

def inverseHist(df1):
	return np.divide(np.ones(df1.shape), df1, where = df1!=0, out = np.zeros_like(df1))

def binVolumes(xBbin, Q2bin, tbin, finehist, k=0):
	xBbins  = collection_xBbins[k]
	Q2bins  = collection_Q2bins[k]
	tbins   = collection_tbins [k]
	phibins = collection_phibins[k]
	fineVols = []
	for phibin in range(len(phibins)-1):
		fineVol = finehist[6*xBbin:6*(xBbin+1), 6*Q2bin:6*(Q2bin+1), 6*tbin:6*(tbin+1), 6*phibin:6*(phibin+1)].flatten()
		fineVols.append(fineVol)
	return (np.sum(fineVols, axis = 1)/6**4)*np.diff(np.radians(phibins))*np.diff(xBbins)[xBbin]*np.diff(Q2bins)[Q2bin]*np.diff(tbins)[tbin]

parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-sc","--savecont", help="save cont", action = "store_true")
parser.add_argument("-sx","--savexsec", help="save xsec", action = "store_true")
parser.add_argument("-ct","--contplot", help="save cont plots", action = "store_true")
parser.add_argument("-sp","--saveplot", help="save plots", action = "store_true")
parser.add_argument("-pr","--parent", help="parent", default = "/volatile/clas12/sangbaek/nov2021/")
parser.add_argument("-k","--binscheme", help="binning scheme number", default=None)
parser.add_argument("-kstart","--kstart", help="binning scheme start", default=None)

args = parser.parse_args()

if args.kstart:
	kstart = int(args.kstart)
else:
	kstart = len(collection_xBbins) - 1 

inbcharge_epg, outbcharge_epg = 30398709.057119943, 32085430.046131916
charge_epg = inbcharge_epg + outbcharge_epg

if args.savecont:
	# read exp
	parent = args.parent

	parent_MC_inb = parent + "convPkl_full/inb/dvcs/"
	parent_MC_BH_inb = parent + "convPkl_full/inb/bh/"
	parent_Gen_inb = parent + "convPkl_Gen/inb/dvcs/"
	parent_Gen_BH_inb = parent + "convPkl_Gen/inb/bh/"
	parent_MC_bkg1g_inb = parent + "convPkl_full/inb/bkg_1g/"
	parent_MC_bkg2g_inb = parent + "convPkl_full/inb/bkg_2g/"
	parent_exp_inb = parent + "convPkl_full/inb/exp/"

	parent_MC_outb = parent + "convPkl_full/outb/dvcs/"
	parent_MC_BH_outb = parent + "convPkl_full/outb/bh/"
	parent_Gen_outb = parent + "convPkl_Gen/outb/dvcs/"
	parent_Gen_BH_outb = parent + "convPkl_Gen/outb/bh/"
	parent_MC_bkg1g_outb = parent + "convPkl_full/outb/bkg_1g/"
	parent_MC_bkg2g_outb = parent + "convPkl_full/outb/bkg_2g/"
	parent_exp_outb = parent + "convPkl_full/outb/exp/"

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
	for i in range(len(collection_cont_xBbins)):

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

		epgExp.loc[:, "newxBbin{}".format(i)] = (len(phibins)-1)*(len(tbins)-1)*(len(Q2bins)-1)*(np.digitize(epgExp.xB, xBbins)-1)
		epgExp.loc[:, "newQ2bin{}".format(i)] = (len(phibins)-1)*(len(tbins)-1)*(np.digitize(epgExp.Q2, Q2bins)-1)
		epgExp.loc[:, "newtbin{}".format(i)] = (len(phibins)-1)*(np.digitize(epgExp.t1, tbins)-1)
		epgExp.loc[:, "newphibin{}".format(i)] = np.digitize(epgExp.phi1, phibins)-1
		epgExp.loc[:, "newbin{}".format(i)] = np.sum(epgExp.loc[:, ["newxBbin{}".format(i), "newQ2bin{}".format(i), "newtbin{}".format(i), "newphibin{}".format(i)]], axis = 1)

		contInbFD = divideHist(histBkgInbFD*histPi0InbFD, histRefInbFD*histExpInbFD)
		contInbCD = divideHist(histBkgInbCD*histPi0InbCD, histRefInbCD*histExpInbCD)
		contInbCDFT = divideHist(histBkgInbCDFT*histPi0InbCDFT, histRefInbCDFT*histExpInbCDFT)
		contOutbFD = divideHist(histBkgOutbFD*histPi0OutbFD, histRefOutbFD*histExpOutbFD)
		contOutbCD = divideHist(histBkgOutbCD*histPi0OutbCD, histRefOutbCD*histExpOutbCD)
		contOutbCDFT = divideHist(histBkgOutbCDFT*histPi0OutbCDFT, histRefOutbCDFT*histExpOutbCDFT)

		unccontInbFD = contInbFD*np.sqrt(inverseHist(histBkgInbFD)+inverseHist(histPi0InbFD)+inverseHist(histRefInbFD)+inverseHist(histExpInbFD))
		unccontInbCD = contInbCD*np.sqrt(inverseHist(histBkgInbCD)+inverseHist(histPi0InbCD)+inverseHist(histRefInbCD)+inverseHist(histExpInbCD))
		unccontInbCDFT = contInbCDFT*np.sqrt(inverseHist(histBkgInbCDFT)+inverseHist(histPi0InbCDFT)+inverseHist(histRefInbCDFT)+inverseHist(histExpInbCDFT))
		unccontOutbFD = contOutbFD*np.sqrt(inverseHist(histBkgOutbFD)+inverseHist(histPi0OutbFD)+inverseHist(histRefOutbFD)+inverseHist(histExpOutbFD))
		unccontOutbCD = contOutbCD*np.sqrt(inverseHist(histBkgOutbCD)+inverseHist(histPi0OutbCD)+inverseHist(histRefOutbCD)+inverseHist(histExpOutbCD))
		unccontOutbCDFT = contOutbCDFT*np.sqrt(inverseHist(histBkgOutbCDFT)+inverseHist(histPi0OutbCDFT)+inverseHist(histRefOutbCDFT)+inverseHist(histExpOutbCDFT))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT.flatten())))

		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontInbFD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontInbCD.flatten())))
		epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontInbCDFT.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(unccontOutbFD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCD.flatten())))
		epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "unccont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(unccontOutbCDFT.flatten())))

		epgExp.loc[epgExp.loc[:, "cont{}".format(i)] > 1, "cont{}".format(i)] = 1

		print("saved contaminations {}".format(i))

	print("clear memory...")

	del df_bkg1gs_inb
	del df_bkg1gs_outb
	del df_bkg2gs_inb
	del df_bkg2gs_outb
	gc.collect()

	epgExp.to_pickle("nphistograms/epgExp.pkl")

if args.savexsec:
	print("read exp...")
	epgExp = pd.read_pickle("nphistograms/epgExp.pkl")

	for k in range(kstart, len(collection_xBbins)):

		print("bin scheme {}".format(k))
		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		for i in range(len(collection_cont_xBbins)):
			#Inbending cross sections 
			# i = 0 #selected background estimation

			print("bkg scheme {}".format(i))
			#inbending
			histExpInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)])
			histBHDVCSInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)])
			histBHDVCSInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)])

			histExpInb = histExpInbFD + histExpInbCD + histExpInbCDFT
			histBHDVCSInb = histBHDVCSInbFD + histBHDVCSInbCD + histBHDVCSInbCDFT

			ActiveInb = histBHDVCSInb>20
			if k == 2:
				ActiveInb[:, 0, :, :] = False
			ActiveInb_int = np.stack([np.sum(ActiveInb, axis=-1)>8]*(len(phibins)-1), axis = -1)

			print("reading bhs - inbending ")

			histBHInb50nA, histBHInbFD50nA, histBHInbCD50nA, histBHInbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_inb_bh50nA:
				histBHInb50nA = histBHInb50nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histBHInbFD50nA = histBHInbFD50nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histBHInbCD50nA = histBHInbCD50nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histBHInbCDFT50nA = histBHInbCDFT50nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histBHInb45nA, histBHInbFD45nA, histBHInbCD45nA, histBHInbCDFT45nA = 0, 0, 0, 0
			for jobNum in runs_inb_bh45nA:
				histBHInb45nA = histBHInb45nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histBHInbFD45nA = histBHInbFD45nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histBHInbCD45nA = histBHInbCD45nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histBHInbCDFT45nA = histBHInbCDFT45nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histBHGenInb50nA, histBHGenInbFD50nA, histBHGenInbCD50nA, histBHGenInbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_inb_bh50nA:
				histBHGenInb50nA = histBHGenInb50nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histBHGenInbFD50nA = histBHGenInbFD50nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histBHGenInbCD50nA = histBHGenInbCD50nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histBHGenInbCDFT50nA = histBHGenInbCDFT50nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histBHGenInb45nA, histBHGenInbFD45nA, histBHGenInbCD45nA, histBHGenInbCDFT45nA = 0, 0, 0, 0
			for jobNum in runs_inb_bh45nA:
				histBHGenInb45nA = histBHGenInb45nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histBHGenInbFD45nA = histBHGenInbFD45nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histBHGenInbCD45nA = histBHGenInbCD45nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histBHGenInbCDFT45nA = histBHGenInbCDFT45nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			print("reading vggs  - inbending")

			histVGGInb50nA, histVGGInbFD50nA, histVGGInbCD50nA, histVGGInbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg50nA:
				histVGGInb50nA = histVGGInb50nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGInbFD50nA = histVGGInbFD50nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGInbCD50nA = histVGGInbCD50nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGInbCDFT50nA = histVGGInbCDFT50nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGInb45nA, histVGGInbFD45nA, histVGGInbCD45nA, histVGGInbCDFT45nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg45nA:
				histVGGInb45nA = histVGGInb45nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGInbFD45nA = histVGGInbFD45nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGInbCD45nA = histVGGInbCD45nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGInbCDFT45nA = histVGGInbCDFT45nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGInb55nA, histVGGInbFD55nA, histVGGInbCD55nA, histVGGInbCDFT55nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg55nA:
				histVGGInb55nA = histVGGInb55nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGInbFD55nA = histVGGInbFD55nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGInbCD55nA = histVGGInbCD55nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGInbCDFT55nA = histVGGInbCDFT55nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGInb0nA, histVGGInbFD0nA, histVGGInbCD0nA, histVGGInbCDFT0nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg0nA:
				histVGGInb0nA = histVGGInb0nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGInbFD0nA = histVGGInbFD0nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGInbCD0nA = histVGGInbCD0nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGInbCDFT0nA = histVGGInbCDFT0nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGGenInb50nA, histVGGGenInbFD50nA, histVGGGenInbCD50nA, histVGGGenInbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg50nA:
				histVGGGenInb50nA = histVGGGenInb50nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD50nA = histVGGGenInbFD50nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD50nA = histVGGGenInbCD50nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT50nA = histVGGGenInbCDFT50nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histVGGGenInb45nA, histVGGGenInbFD45nA, histVGGGenInbCD45nA, histVGGGenInbCDFT45nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg45nA:
				histVGGGenInb45nA = histVGGGenInb45nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD45nA = histVGGGenInbFD45nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD45nA = histVGGGenInbCD45nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT45nA = histVGGGenInbCDFT45nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histVGGGenInb55nA, histVGGGenInbFD55nA, histVGGGenInbCD55nA, histVGGGenInbCDFT55nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg55nA:
				histVGGGenInb55nA = histVGGGenInb55nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD55nA = histVGGGenInbFD55nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD55nA = histVGGGenInbCD55nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT55nA = histVGGGenInbCDFT55nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histVGGGenInb0nA, histVGGGenInbFD0nA, histVGGGenInbCD0nA, histVGGGenInbCDFT0nA = 0, 0, 0, 0
			for jobNum in runs_inb_vgg0nA:
				histVGGGenInb0nA = histVGGGenInb0nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenInbFD0nA = histVGGGenInbFD0nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenInbCD0nA = histVGGGenInbCD0nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenInbCDFT0nA = histVGGGenInbCDFT0nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histBkgUncInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "unccont{}".format(i)])
			histBkgUncInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "unccont{}".format(i)])
			histBkgUncInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "unccont{}".format(i)])
			histExpUncInbFD = np.sqrt(histExpInbFD)
			histExpUncInbCD = np.sqrt(histExpInbCD)
			histExpUncInbCDFT = np.sqrt(histExpInbCDFT)

			histVGGGenInbInt50nA, histVGGGenInbxB50nA, histVGGGenInbQ250nA, histVGGGenInbt150nA = 0, 0, 0, 0
			histVGGGenInbphi50nA, histVGGGenInbbinVol50nA = 0, 0
			histVGGGenInbrad50nA, histVGGGenInbborn50nA = [], []
			for jobNum in runs_inb_vgg50nA:
				histVGGGenInbInt50nA = histVGGGenInbInt50nA + np.load("nphistograms/binscheme{}/{}GenInt.npz".format(k, jobNum))["hist"]
				histVGGGenInbxB50nA = histVGGGenInbxB50nA + np.load("nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"]
				histVGGGenInbQ250nA = histVGGGenInbQ250nA + np.load("nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"]
				histVGGGenInbt150nA = histVGGGenInbt150nA + np.load("nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"]
				histVGGGenInbphi50nA = histVGGGenInbphi50nA + np.load("nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"]
				histVGGGenInbrad50nA.append(np.load("nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
				histVGGGenInbborn50nA.append(np.load("nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
				histVGGGenInbbinVol50nA = histVGGGenInbbinVol50nA | np.load("nphistograms/binscheme{}/{}GenbinVolume.npz".format(k, jobNum))["hist"].astype(int)

			histBHGenInbInt45nA, histBHGenInbxB45nA, histBHGenInbQ245nA, histBHGenInbt145nA = 0, 0, 0, 0
			histBHGenInbphi45nA, histBHGenInbbinVol45nA = 0, 0
			histBHGenInbrad45nA, histBHGenInbborn45nA = [], []
			for jobNum in runs_inb_bh45nA:
				histBHGenInbInt45nA = histBHGenInbInt45nA + np.load("nphistograms/binscheme{}/{}GenInt.npz".format(k, jobNum))["hist"]
				histBHGenInbxB45nA = histBHGenInbxB45nA + np.load("nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"]
				histBHGenInbQ245nA = histBHGenInbQ245nA + np.load("nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"]
				histBHGenInbt145nA = histBHGenInbt145nA + np.load("nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"]
				histBHGenInbphi45nA = histBHGenInbphi45nA + np.load("nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"]
				histBHGenInbrad45nA.append(np.load("nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
				histBHGenInbborn45nA.append(np.load("nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
				histBHGenInbbinVol45nA = histBHGenInbbinVol45nA | np.load("nphistograms/binscheme{}/{}GenbinVolume.npz".format(k, jobNum))["hist"].astype(int)

			#outbending
			histExpOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histExpOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			histBHDVCSOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)])
			histBHDVCSOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)])
			histBHDVCSOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)])

			histExpOutb = histExpOutbFD + histExpOutbCD + histExpOutbCDFT
			histBHDVCSOutb = histBHDVCSOutbFD + histBHDVCSOutbCD + histBHDVCSOutbCDFT

			ActiveOutb = histBHDVCSOutb>20
			ActiveOutb_int = np.stack([np.sum(ActiveOutb, axis=-1)>8]*(len(phibins)-1), axis = -1)
			ActiveAll = ActiveInb & ActiveOutb
			ActiveAny = ActiveInb | ActiveOutb
			ActiveAll_int = ActiveInb_int & ActiveOutb_int
			ActiveAny_int = ActiveInb_int | ActiveOutb_int

			print("reading bhs  - outbending")

			histBHOutb50nA, histBHOutbFD50nA, histBHOutbCD50nA, histBHOutbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_outb_bh50nA:
				histBHOutb50nA = histBHOutb50nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histBHOutbFD50nA = histBHOutbFD50nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histBHOutbCD50nA = histBHOutbCD50nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histBHOutbCDFT50nA = histBHOutbCDFT50nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histBHGenOutb50nA, histBHGenOutbFD50nA, histBHGenOutbCD50nA, histBHGenOutbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_outb_bh50nA:
				histBHGenOutb50nA = histBHGenOutb50nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histBHGenOutbFD50nA = histBHGenOutbFD50nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histBHGenOutbCD50nA = histBHGenOutbCD50nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histBHGenOutbCDFT50nA = histBHGenOutbCDFT50nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			print("reading vggs  - outbending")

			histVGGOutb50nA, histVGGOutbFD50nA, histVGGOutbCD50nA, histVGGOutbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_outb_vgg50nA:
				histVGGOutb50nA = histVGGOutb50nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGOutbFD50nA = histVGGOutbFD50nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGOutbCD50nA = histVGGOutbCD50nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGOutbCDFT50nA = histVGGOutbCDFT50nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGOutb40nA, histVGGOutbFD40nA, histVGGOutbCD40nA, histVGGOutbCDFT40nA = 0, 0, 0, 0
			for jobNum in runs_outb_vgg40nA:
				histVGGOutb40nA = histVGGOutb40nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGOutbFD40nA = histVGGOutbFD40nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGOutbCD40nA = histVGGOutbCD40nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGOutbCDFT40nA = histVGGOutbCDFT40nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGOutb0nA, histVGGOutbFD0nA, histVGGOutbCD0nA, histVGGOutbCDFT0nA = 0, 0, 0, 0
			for jobNum in runs_outb_vgg0nA:
				histVGGOutb0nA = histVGGOutb0nA + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGOutbFD0nA = histVGGOutbFD0nA + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGOutbCD0nA = histVGGOutbCD0nA + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGOutbCDFT0nA = histVGGOutbCDFT0nA + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGOutb40naT, histVGGOutbFD40naT, histVGGOutbCD40naT, histVGGOutbCDFT40naT = 0, 0, 0, 0
			for jobNum in runs_outb_vgg40nAT:
				histVGGOutb40naT = histVGGOutb40naT + np.load("nphistograms/binscheme{}/{}Rec.npz".format(k, jobNum))["hist"]
				histVGGOutbFD40naT = histVGGOutbFD40naT + np.load("nphistograms/binscheme{}/{}Rec1.npz".format(k, jobNum))["hist"]
				histVGGOutbCD40naT = histVGGOutbCD40naT + np.load("nphistograms/binscheme{}/{}Rec2.npz".format(k, jobNum))["hist"]
				histVGGOutbCDFT40naT = histVGGOutbCDFT40naT + np.load("nphistograms/binscheme{}/{}Rec3.npz".format(k, jobNum))["hist"]

			histVGGGenOutb50nA, histVGGGenOutbFD50nA, histVGGGenOutbCD50nA, histVGGGenOutbCDFT50nA = 0, 0, 0, 0
			for jobNum in runs_outb_vgg50nA:
				histVGGGenOutb50nA = histVGGGenOutb50nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD50nA = histVGGGenOutbFD50nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD50nA = histVGGGenOutbCD50nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT50nA = histVGGGenOutbCDFT50nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histVGGGenOutb40nA, histVGGGenOutbFD40nA, histVGGGenOutbCD40nA, histVGGGenOutbCDFT40nA = 0, 0, 0, 0
			for jobNum in runs_outb_vgg40nA:
				histVGGGenOutb40nA = histVGGGenOutb40nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD40nA = histVGGGenOutbFD40nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD40nA = histVGGGenOutbCD40nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT40nA = histVGGGenOutbCDFT40nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histVGGGenOutb0nA, histVGGGenOutbFD0nA, histVGGGenOutbCD0nA, histVGGGenOutbCDFT0nA = 0, 0, 0, 0
			for jobNum in runs_outb_vgg0nA:
				histVGGGenOutb0nA = histVGGGenOutb0nA + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD0nA = histVGGGenOutbFD0nA + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD0nA = histVGGGenOutbCD0nA + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT0nA = histVGGGenOutbCDFT0nA + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histVGGGenOutb40nAT, histVGGGenOutbFD40nAT, histVGGGenOutbCD40nAT, histVGGGenOutbCDFT40nAT = 0, 0, 0, 0
			for jobNum in runs_outb_vgg40nAT:
				histVGGGenOutb40nAT = histVGGGenOutb40nAT + np.load("nphistograms/binscheme{}/{}Gen.npz".format(k, jobNum))["hist"]
				histVGGGenOutbFD40nAT = histVGGGenOutbFD40nAT + np.load("nphistograms/binscheme{}/{}Gen1.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCD40nAT = histVGGGenOutbCD40nAT + np.load("nphistograms/binscheme{}/{}Gen2.npz".format(k, jobNum))["hist"]
				histVGGGenOutbCDFT40nAT = histVGGGenOutbCDFT40nAT + np.load("nphistograms/binscheme{}/{}Gen3.npz".format(k, jobNum))["hist"]

			histBkgUncOutbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "unccont{}".format(i)])
			histBkgUncOutbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "unccont{}".format(i)])
			histBkgUncOutbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "unccont{}".format(i)])
			histExpUncOutbFD = np.sqrt(histExpOutbFD)
			histExpUncOutbCD = np.sqrt(histExpOutbCD)
			histExpUncOutbCDFT = np.sqrt(histExpOutbCDFT)

			histVGGGenOutbInt50nA, histVGGGenOutbxB50nA, histVGGGenOutbQ250nA, histVGGGenOutbt150nA = 0, 0, 0, 0
			histVGGGenOutbphi50nA, histVGGGenOutbbinVol50nA = 0, 0
			histVGGGenOutbrad50nA, histVGGGenOutbborn50nA = [], []
			for jobNum in runs_outb_vgg50nA:
				histVGGGenOutbInt50nA = histVGGGenOutbInt50nA + np.load("nphistograms/binscheme{}/{}GenInt.npz".format(k, jobNum))["hist"]
				histVGGGenOutbxB50nA = histVGGGenOutbxB50nA + np.load("nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"]
				histVGGGenOutbQ250nA = histVGGGenOutbQ250nA + np.load("nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"]
				histVGGGenOutbt150nA = histVGGGenOutbt150nA + np.load("nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"]
				histVGGGenOutbphi50nA = histVGGGenOutbphi50nA + np.load("nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"]
				histVGGGenOutbrad50nA.append(np.load("nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
				histVGGGenOutbborn50nA.append(np.load("nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
				histVGGGenOutbbinVol50nA = histVGGGenOutbbinVol50nA | np.load("nphistograms/binscheme{}/{}GenbinVolume.npz".format(k, jobNum))["hist"].astype(int)

			histBHGenOutbInt50nA, histBHGenOutbxB50nA, histBHGenOutbQ250nA, histBHGenOutbt150nA = 0, 0, 0, 0
			histBHGenOutbphi50nA, histBHGenOutbbinVol50nA = 0, 0
			histBHGenOutbrad50nA, histBHGenOutbborn50nA = [], []
			for jobNum in runs_outb_bh50nA:
				histBHGenOutbInt50nA = histBHGenOutbInt50nA + np.load("nphistograms/binscheme{}/{}GenInt.npz".format(k, jobNum))["hist"]
				histBHGenOutbxB50nA = histBHGenOutbxB50nA + np.load("nphistograms/binscheme{}/{}GenxB.npz".format(k, jobNum))["hist"]
				histBHGenOutbQ250nA = histBHGenOutbQ250nA + np.load("nphistograms/binscheme{}/{}GenQ2.npz".format(k, jobNum))["hist"]
				histBHGenOutbt150nA = histBHGenOutbt150nA + np.load("nphistograms/binscheme{}/{}Gent1.npz".format(k, jobNum))["hist"]
				histBHGenOutbphi50nA = histBHGenOutbphi50nA + np.load("nphistograms/binscheme{}/{}Genphi.npz".format(k, jobNum))["hist"]
				histBHGenOutbrad50nA.append(np.load("nphistograms/binscheme{}/{}Genrad.npz".format(k, jobNum))["hist"])
				histBHGenOutbborn50nA.append(np.load("nphistograms/binscheme{}/{}Genborn.npz".format(k, jobNum))["hist"])
				histBHGenOutbbinVol50nA = histBHGenOutbbinVol50nA | np.load("nphistograms/binscheme{}/{}GenbinVolume.npz".format(k, jobNum))["hist"].astype(int)

			#stat error - inbending
			accCorrectedInbFD_VGG = divideHist(histBHDVCSInbFD*histVGGGenInbFD50nA , histVGGInbFD50nA)
			accCorrectedInbCD_VGG = divideHist(histBHDVCSInbCD*histVGGGenInbCD50nA , histVGGInbCD50nA)
			accCorrectedInbCDFT_VGG = divideHist(histBHDVCSInbCDFT*histVGGGenInbCDFT50nA , histVGGInbCDFT50nA)
			uncStatInbFD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD**2+histExpUncInbFD**2), histBHDVCSInbFD)**2 + inverseHist(histVGGInbFD50nA) + inverseHist(histVGGGenInbFD50nA))
			uncStatInbCD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD**2+histExpUncInbCD**2), histBHDVCSInbCD)**2 + inverseHist(histVGGInbCD50nA) + inverseHist(histVGGGenInbCD50nA))
			uncStatInbCDFT_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT**2+histExpUncInbCDFT**2), histBHDVCSInbCDFT)**2 + inverseHist(histVGGInbCDFT50nA) + inverseHist(histVGGGenInbCDFT50nA))
			uncStatInb_VGG = np.sqrt(accCorrectedInbFD_VGG**2 *uncStatInbFD_VGG**2 + accCorrectedInbCD_VGG**2 * uncStatInbCD_VGG**2 + accCorrectedInbCDFT_VGG**2 * uncStatInbCDFT_VGG**2)

			accCorrectedInb_VGG = accCorrectedInbFD_VGG + accCorrectedInbCD_VGG + accCorrectedInbCDFT_VGG
			accCorrectedInb_VGG = divideHist(accCorrectedInb_VGG*histVGGGenInb50nA, histVGGGenInbFD50nA + histVGGGenInbCD50nA+ histVGGGenInbCDFT50nA)
			uncStatInb_VGG = divideHist(uncStatInb_VGG, accCorrectedInb_VGG)
			accCorrectedInb_VGG[~ActiveInb] = 0
			uncStatInb_VGG[~ActiveInb] = 0

			accCorrectedInbFD_BH = divideHist(histBHDVCSInbFD*histBHGenInbFD45nA , histBHInbFD45nA)
			accCorrectedInbCD_BH = divideHist(histBHDVCSInbCD*histBHGenInbCD45nA , histBHInbCD45nA)
			accCorrectedInbCDFT_BH = divideHist(histBHDVCSInbCDFT*histBHGenInbCDFT45nA , histBHInbCDFT45nA)
			uncStatInbFD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbFD**2+histExpUncInbFD**2), histBHDVCSInbFD)**2 + inverseHist(histBHInbFD45nA) + inverseHist(histBHGenInbFD45nA))
			uncStatInbCD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCD**2+histExpUncInbCD**2), histBHDVCSInbCD)**2 + inverseHist(histBHInbCD45nA) + inverseHist(histBHGenInbCD45nA))
			uncStatInbCDFT_BH = np.sqrt(divideHist(np.sqrt(histBkgUncInbCDFT**2+histExpUncInbCDFT**2), histBHDVCSInbCDFT)**2 + inverseHist(histBHInbCDFT45nA) + inverseHist(histBHGenInbCDFT45nA))
			uncStatInb_BH = np.sqrt(accCorrectedInbFD_BH**2 *uncStatInbFD_BH**2 + accCorrectedInbCD_BH**2 * uncStatInbCD_BH**2 + accCorrectedInbCDFT_BH**2 * uncStatInbCDFT_BH**2)

			accCorrectedInb_BH = accCorrectedInbFD_BH + accCorrectedInbCD_BH + accCorrectedInbCDFT_BH
			accCorrectedInb_BH = divideHist(accCorrectedInb_BH*histBHGenInb45nA, histBHGenInbFD45nA + histBHGenInbCD45nA+ histBHGenInbCDFT45nA)
			uncStatInb_BH = divideHist(uncStatInb_BH, accCorrectedInb_BH)
			accCorrectedInb_BH[~ActiveInb] = 0
			uncStatInb_BH[~ActiveInb] = 0

			#stat error - outbending
			accCorrectedOutbFD_VGG = divideHist(histBHDVCSOutbFD*histVGGGenOutbFD50nA , histVGGOutbFD50nA)
			accCorrectedOutbCD_VGG = divideHist(histBHDVCSOutbCD*histVGGGenOutbCD50nA , histVGGOutbCD50nA)
			accCorrectedOutbCDFT_VGG = divideHist(histBHDVCSOutbCDFT*histVGGGenOutbCDFT50nA , histVGGOutbCDFT50nA)
			uncStatOutbFD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD**2+histExpUncOutbFD**2), histBHDVCSOutbFD)**2 + inverseHist(histVGGOutbFD50nA) + inverseHist(histVGGGenOutbFD50nA))
			uncStatOutbCD_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD**2+histExpUncOutbCD**2), histBHDVCSOutbCD)**2 + inverseHist(histVGGOutbCD50nA) + inverseHist(histVGGGenOutbCD50nA))
			uncStatOutbCDFT_VGG = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT**2+histExpUncOutbCDFT**2), histBHDVCSOutbCDFT)**2 + inverseHist(histVGGOutbCDFT50nA) + inverseHist(histVGGGenOutbCDFT50nA))
			uncStatOutb_VGG = np.sqrt(accCorrectedOutbFD_VGG**2 *uncStatOutbFD_VGG**2 + accCorrectedOutbCD_VGG**2 * uncStatOutbCD_VGG**2 + accCorrectedOutbCDFT_VGG**2 * uncStatOutbCDFT_VGG**2)

			accCorrectedOutb_VGG = accCorrectedOutbFD_VGG + accCorrectedOutbCD_VGG + accCorrectedOutbCDFT_VGG
			accCorrectedOutb_VGG = divideHist(accCorrectedOutb_VGG*histVGGGenOutb50nA, histVGGGenOutbFD50nA + histVGGGenOutbCD50nA+ histVGGGenOutbCDFT50nA)
			uncStatOutb_VGG = divideHist(uncStatOutb_VGG, accCorrectedOutb_VGG)
			accCorrectedOutb_VGG[~ActiveOutb] = 0
			uncStatOutb_VGG[~ActiveOutb] = 0

			accCorrectedOutbFD_BH = divideHist(histBHDVCSOutbFD*histBHGenOutbFD50nA , histBHOutbFD50nA)
			accCorrectedOutbCD_BH = divideHist(histBHDVCSOutbCD*histBHGenOutbCD50nA , histBHOutbCD50nA)
			accCorrectedOutbCDFT_BH = divideHist(histBHDVCSOutbCDFT*histBHGenOutbCDFT50nA , histBHOutbCDFT50nA)
			uncStatOutbFD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbFD**2+histExpUncOutbFD**2), histBHDVCSOutbFD)**2 + inverseHist(histBHOutbFD50nA) + inverseHist(histBHGenOutbFD50nA))
			uncStatOutbCD_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCD**2+histExpUncOutbCD**2), histBHDVCSOutbCD)**2 + inverseHist(histBHOutbCD50nA) + inverseHist(histBHGenOutbCD50nA))
			uncStatOutbCDFT_BH = np.sqrt(divideHist(np.sqrt(histBkgUncOutbCDFT**2+histExpUncOutbCDFT**2), histBHDVCSOutbCDFT)**2 + inverseHist(histBHOutbCDFT50nA) + inverseHist(histBHGenOutbCDFT50nA))
			uncStatOutb_BH = np.sqrt(accCorrectedOutbFD_BH**2 *uncStatOutbFD_BH**2 + accCorrectedOutbCD_BH**2 * uncStatOutbCD_BH**2 + accCorrectedOutbCDFT_BH**2 * uncStatOutbCDFT_BH**2)

			accCorrectedOutb_BH = accCorrectedOutbCDFT_BH + accCorrectedOutbCD_BH + accCorrectedOutbFD_BH
			accCorrectedOutb_BH = divideHist(accCorrectedOutb_BH*histBHGenOutb50nA, histBHGenOutbFD50nA + histBHGenOutbCD50nA+ histBHGenOutbCDFT50nA)
			uncStatOutb_BH = divideHist(uncStatOutb_BH, accCorrectedOutb_BH)
			accCorrectedOutb_BH[~ActiveOutb] = 0
			uncStatOutb_BH[~ActiveOutb] = 0

			#stat error - all
			accCorrected_VGG = accCorrectedInb_VGG + accCorrectedOutb_VGG
			uncStat_VGG = np.sqrt((accCorrectedInb_VGG*uncStatInb_VGG)**2 + (accCorrectedOutb_VGG*uncStatOutb_VGG)**2)
			uncStat_VGG = divideHist(uncStat_VGG, accCorrected_VGG)

			accCorrected_BH = accCorrectedInb_BH + accCorrectedOutb_BH
			uncStat_BH = np.sqrt((accCorrectedInb_BH*uncStatInb_BH)**2 + (accCorrectedOutb_BH*uncStatOutb_BH)**2)
			uncStat_BH = divideHist(uncStat_BH, accCorrected_BH)

			phi1avg_VGG = np.zeros(accCorrected_VGG.shape)
			xBavg_VGG = np.zeros(accCorrected_VGG.shape)
			Q2avg_VGG = np.zeros(accCorrected_VGG.shape)
			t1avg_VGG = np.zeros(accCorrected_VGG.shape)
			integratedRad_VGG = np.zeros(accCorrected_VGG.shape)
			pointBorn_VGG = np.zeros(accCorrected_VGG.shape)
			rcfactors_VGG = np.zeros(accCorrected_VGG.shape)

			phi1avg_BH = np.zeros(accCorrected_BH.shape)
			xBavg_BH = np.zeros(accCorrected_BH.shape)
			Q2avg_BH = np.zeros(accCorrected_BH.shape)
			t1avg_BH = np.zeros(accCorrected_BH.shape)
			integratedRad_BH = np.zeros(accCorrected_BH.shape)
			pointBorn_BH = np.zeros(accCorrected_BH.shape)
			rcfactors_BH = np.zeros(accCorrected_BH.shape)

			xsecTh_KM = np.zeros(accCorrected_VGG.shape)
			xsecTh_BH = np.zeros(accCorrected_VGG.shape)
			xsecTh_VGG = np.zeros(accCorrected_VGG.shape)

			binVolume = np.zeros(accCorrected_VGG.shape)

			phi1avg_VGG = divideHist(histVGGGenInbphi50nA+histVGGGenOutbphi50nA, histVGGGenInb50nA+histVGGGenOutb50nA)
			xBavg_VGG = np.stack([divideHist(histVGGGenInbxB50nA+histVGGGenOutbxB50nA, histVGGGenInbInt50nA+histVGGGenOutbInt50nA)]*(len(phibins)-1), axis = -1)
			Q2avg_VGG = np.stack([divideHist(histVGGGenInbQ250nA+histVGGGenOutbQ250nA, histVGGGenInbInt50nA+histVGGGenOutbInt50nA)]*(len(phibins)-1), axis = -1)
			t1avg_VGG = np.stack([divideHist(histVGGGenInbt150nA+histVGGGenOutbt150nA, histVGGGenInbInt50nA+histVGGGenOutbInt50nA)]*(len(phibins)-1), axis = -1)
			phi1avg_VGG[~ActiveAny_int] = 0
			xBavg_VGG[~ActiveAny_int] = 0
			Q2avg_VGG[~ActiveAny_int] = 0
			t1avg_VGG[~ActiveAny_int] = 0

			xsecTh_VGG[ActiveAny_int] = np.array(printVGGarray(xBavg_VGG[ActiveAny_int], Q2avg_VGG[ActiveAny_int], t1avg_VGG[ActiveAny_int], np.radians(phi1avg_VGG[ActiveAny_int]), globalfit = True))

			phi1avg_BH = divideHist(histBHGenInbphi45nA+histBHGenOutbphi50nA, histBHGenInb45nA+histBHGenOutb50nA)
			xBavg_BH = np.stack([divideHist(histBHGenInbxB45nA+histBHGenOutbxB50nA, histBHGenInbInt45nA+histBHGenOutbInt50nA)]*(len(phibins)-1), axis = -1)
			Q2avg_BH = np.stack([divideHist(histBHGenInbQ245nA+histBHGenOutbQ250nA, histBHGenInbInt45nA+histBHGenOutbInt50nA)]*(len(phibins)-1), axis = -1)
			t1avg_BH = np.stack([divideHist(histBHGenInbt145nA+histBHGenOutbt150nA, histBHGenInbInt45nA+histBHGenOutbInt50nA)]*(len(phibins)-1), axis = -1)
			phi1avg_BH[~ActiveAny_int] = 0
			xBavg_BH[~ActiveAny_int] = 0
			Q2avg_BH[~ActiveAny_int] = 0
			t1avg_BH[~ActiveAny_int] = 0

			binVolume = np.zeros(phi1avg_BH.shape)
			xsecTh_BH[ActiveAny_int] = np.array(printBHarray(xBavg_BH[ActiveAny_int], Q2avg_BH[ActiveAny_int], t1avg_BH[ActiveAny_int], np.radians(phi1avg_BH[ActiveAny_int]), globalfit = True))
			xsecTh_KM[ActiveAny_int] = np.array(printKMarray(xBavg_BH[ActiveAny_int], Q2avg_BH[ActiveAny_int], t1avg_BH[ActiveAny_int], np.radians(phi1avg_BH[ActiveAny_int])))
			
			for xBbin in range(len(xBbins)-1):
				for Q2bin in range(len(Q2bins)-1):
					for tbin in range(len(tbins)-1):
						binVolume[xBbin, Q2bin, tbin, :] = binVolumes(xBbin, Q2bin, tbin, histBHGenInbbinVol45nA|histVGGGenInbbinVol50nA|histBHGenOutbbinVol50nA|histVGGGenOutbbinVol50nA, k = k)

			integratedRad_VGG = np.mean([*histVGGGenInbrad50nA, *histVGGGenOutbrad50nA], axis = 0)
			rcfactors_VGG = divideHist(integratedRad_VGG, xsecTh_VGG)
			integratedRad_BH = np.mean([*histBHGenInbrad45nA, *histBHGenOutbrad50nA], axis = 0)
			rcfactors_BH = divideHist(integratedRad_BH, xsecTh_BH)

			charges = charge_epg*np.zeros(xsecTh_BH.shape)
			charges[ActiveAll] = charge_epg
			charges[(ActiveInb)&(~ActiveOutb)] = inbcharge_epg
			charges[(ActiveOutb)&(~ActiveInb)] = outbcharge_epg
			xsecInb_VGG = divideHist(accCorrectedInb_VGG, binVolume*rcfactors_VGG)/(1.324*inbcharge_epg)
			xsecInb_BH = divideHist(accCorrectedInb_BH, binVolume*rcfactors_BH)/(1.324*inbcharge_epg)
			xsecOutb_VGG = divideHist(accCorrectedOutb_VGG, binVolume*rcfactors_VGG)/(1.324*outbcharge_epg)
			xsecOutb_BH = divideHist(accCorrectedOutb_BH, binVolume*rcfactors_BH)/(1.324*outbcharge_epg)
			xsec_VGG = divideHist(accCorrected_VGG, binVolume*rcfactors_VGG*(1.324*charges))
			xsec_BH = divideHist(accCorrected_BH, binVolume*rcfactors_BH*(1.324*charges))

			np.savez("nphistograms/binscheme{}/bkgscheme{}phi1avg_VGG.npz".format(k, i), hist = phi1avg_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xBavg_VGG.npz".format(k, i), hist = xBavg_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}Q2avg_VGG.npz".format(k, i), hist = Q2avg_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}t1avg_VGG.npz".format(k, i), hist = t1avg_VGG)

			np.savez("nphistograms/binscheme{}/bkgscheme{}phi1avg_BH.npz".format(k, i), hist = phi1avg_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xBavg_BH.npz".format(k, i), hist = xBavg_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}Q2avg_BH.npz".format(k, i), hist = Q2avg_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}t1avg_BH.npz".format(k, i), hist = t1avg_BH)

			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(k, i), hist = xsecInb_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(k, i), hist = xsecInb_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(k, i), hist = xsecOutb_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(k, i), hist = xsecOutb_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsec_VGG.npz".format(k, i), hist = xsec_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsec_BH.npz".format(k, i), hist = xsec_BH)

			np.savez("nphistograms/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(k, i), hist = uncStatInb_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(k, i), hist = uncStatInb_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(k, i), hist = uncStatOutb_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(k, i), hist = uncStatOutb_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(k, i), hist = uncStat_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}uncStat_BH.npz".format(k, i), hist = uncStat_BH)

			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecTh_KM.npz".format(k, i), hist = xsecTh_KM)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecTh_BH.npz".format(k, i), hist = xsecTh_BH)
			np.savez("nphistograms/binscheme{}/bkgscheme{}xsecTh_VGG.npz".format(k, i), hist = xsecTh_VGG)
			np.savez("nphistograms/binscheme{}/bkgscheme{}binVolume.npz".format(k, i), hist = binVolume)

			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveAll.npz".format(k, i), hist = ActiveAll)
			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveAny.npz".format(k, i), hist = ActiveAny)
			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveInb.npz".format(k, i), hist = ActiveInb)
			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveOutb.npz".format(k, i), hist = ActiveOutb)

			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(k, i), hist = ActiveAll_int)
			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(k, i), hist = ActiveAny_int)
			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(k, i), hist = ActiveInb_int)
			np.savez("nphistograms/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(k, i), hist = ActiveOutb_int)

if args.contplot:
	print("read exp...")
	epgExp = pd.read_pickle("nphistograms/epgExp.pkl")

	k = 2
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

	xBbin = 4
	Q2bin = 2
	tbin = 1

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
	tbin = 1

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
# 	for k in range(kstart, len(collection_xBbins)):
# 		for i in range(len(collection_cont_xBbins)):
# 			phi1avg_VGG [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}phi1avg_VGG.npz".format(k, i))["hist"]
# 			xBavg_VGG   [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xBavg_VGG.npz".format(k, i))["hist"]
# 			Q2avg_VGG   [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}Q2avg_VGG.npz".format(k, i))["hist"]
# 			t1avg_VGG   [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}t1avg_VGG.npz".format(k, i))["hist"]

# 			phi1avg_BH  [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}phi1avg_BH.npz".format(k, i))["hist"]
# 			xBavg_BH    [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xBavg_BH.npz".format(k, i))["hist"]
# 			Q2avg_BH    [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}Q2avg_BH.npz".format(k, i))["hist"]
# 			t1avg_BH    [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}t1avg_BH.npz".format(k, i))["hist"]

# 			xsecInb_VGG [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(k, i))["hist"]
# 			xsecInb_BH  [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(k, i))["hist"]
# 			xsecOutb_VGG [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(k, i))["hist"]
# 			xsecOutb_BH  [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(k, i))["hist"]
# 			xsec_VGG     [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsec_VGG.npz".format(k, i))["hist"]
# 			xsec_BH      [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsec_BH.npz".format(k, i))["hist"]

# 			uncStatInb_VGG [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(k, i))["hist"]
# 			uncStatInb_BH  [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(k, i))["hist"]
# 			uncStatOutb_VGG [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(k, i))["hist"]
# 			uncStatOutb_BH  [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(k, i))["hist"]
# 			uncStat_VGG     [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(k, i))["hist"]
# 			uncStat_BH      [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}uncStat_BH.npz".format(k, i))["hist"]

# 			xsecTh_KM          [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecTh_KM.npz".format(k, i))["hist"]
# 			xsecTh_BH          [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecTh_BH.npz".format(k, i))["hist"]
# 			xsecTh_VGG         [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}xsecTh_VGG.npz".format(k, i))["hist"]
# 			binVolume          [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}binVolume.npz".format(k, i))["hist"]

# 			ActiveAll       [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAll.npz".format(k, i))["hist"]
# 			ActiveAny       [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAny.npz".format(k, i))["hist"]
# 			ActiveInb          [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveInb.npz".format(k, i))["hist"]
# 			ActiveOutb         [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveOutb.npz".format(k, i))["hist"]

# 			ActiveAll_int       [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(k, i))["hist"]
# 			ActiveAny_int       [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(k, i))["hist"]
# 			ActiveInb_int          [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(k, i))["hist"]
# 			ActiveOutb_int         [k, i] = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(k, i))["hist"]

if args.saveplot:

	for k in range(kstart, len(collection_xBbins)):

		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		i = 0

		os.makedirs("plots/binscheme{}".format(k), exist_ok = True)

		print("reading the xsec vars")
		phi1avg_VGG = np.load("nphistograms/binscheme{}/bkgscheme{}phi1avg_VGG.npz".format(k, i))["hist"]
		xBavg_VGG   = np.load("nphistograms/binscheme{}/bkgscheme{}xBavg_VGG.npz".format(k, i))["hist"]
		Q2avg_VGG   = np.load("nphistograms/binscheme{}/bkgscheme{}Q2avg_VGG.npz".format(k, i))["hist"]
		t1avg_VGG   = np.load("nphistograms/binscheme{}/bkgscheme{}t1avg_VGG.npz".format(k, i))["hist"]

		phi1avg_BH  = np.load("nphistograms/binscheme{}/bkgscheme{}phi1avg_BH.npz".format(k, i))["hist"]
		xBavg_BH    = np.load("nphistograms/binscheme{}/bkgscheme{}xBavg_BH.npz".format(k, i))["hist"]
		Q2avg_BH    = np.load("nphistograms/binscheme{}/bkgscheme{}Q2avg_BH.npz".format(k, i))["hist"]
		t1avg_BH    = np.load("nphistograms/binscheme{}/bkgscheme{}t1avg_BH.npz".format(k, i))["hist"]

		xsecInb_VGG = np.load("nphistograms/binscheme{}/bkgscheme{}xsecInb_VGG.npz".format(k, i))["hist"]
		xsecInb_BH  = np.load("nphistograms/binscheme{}/bkgscheme{}xsecInb_BH.npz".format(k, i))["hist"]
		xsecOutb_VGG = np.load("nphistograms/binscheme{}/bkgscheme{}xsecOutb_VGG.npz".format(k, i))["hist"]
		xsecOutb_BH  = np.load("nphistograms/binscheme{}/bkgscheme{}xsecOutb_BH.npz".format(k, i))["hist"]
		xsec_VGG     = np.load("nphistograms/binscheme{}/bkgscheme{}xsec_VGG.npz".format(k, i))["hist"]
		xsec_BH      = np.load("nphistograms/binscheme{}/bkgscheme{}xsec_BH.npz".format(k, i))["hist"]

		uncStatInb_VGG = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatInb_VGG.npz".format(k, i))["hist"]
		uncStatInb_BH  = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatInb_BH.npz".format(k, i))["hist"]
		uncStatOutb_VGG = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatOutb_VGG.npz".format(k, i))["hist"]
		uncStatOutb_BH  = np.load("nphistograms/binscheme{}/bkgscheme{}uncStatOutb_BH.npz".format(k, i))["hist"]
		uncStat_VGG     = np.load("nphistograms/binscheme{}/bkgscheme{}uncStat_VGG.npz".format(k, i))["hist"]
		uncStat_BH      = np.load("nphistograms/binscheme{}/bkgscheme{}uncStat_BH.npz".format(k, i))["hist"]

		xsecTh_KM          = np.load("nphistograms/binscheme{}/bkgscheme{}xsecTh_KM.npz".format(k, i))["hist"]
		xsecTh_BH          = np.load("nphistograms/binscheme{}/bkgscheme{}xsecTh_BH.npz".format(k, i))["hist"]
		xsecTh_VGG         = np.load("nphistograms/binscheme{}/bkgscheme{}xsecTh_VGG.npz".format(k, i))["hist"]
		binVolume          = np.load("nphistograms/binscheme{}/bkgscheme{}binVolume.npz".format(k, i))["hist"]

		ActiveAll       = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAll.npz".format(k, i))["hist"]
		ActiveAny       = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAny.npz".format(k, i))["hist"]
		ActiveInb          = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveInb.npz".format(k, i))["hist"]
		ActiveOutb         = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveOutb.npz".format(k, i))["hist"]

		ActiveAll_int       = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAll_int.npz".format(k, i))["hist"]
		ActiveAny_int       = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveAny_int.npz".format(k, i))["hist"]
		ActiveInb_int          = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveInb_int.npz".format(k, i))["hist"]
		ActiveOutb_int         = np.load("nphistograms/binscheme{}/bkgscheme{}ActiveOutb_int.npz".format(k, i))["hist"]

		print("plotting...")

		def badBinCondxBQ2t(xBbin, Q2bin, tbin, k = 0):
			if k ==0:
				return (xBbin==1 and Q2bin == 0) or (xBbin==0 and Q2bin==4) or (tbin==0 and xBbin==1)
			else:
				return ~ActiveAny_int[xBbin, Q2bin, tbin, :].any()

		def badBinCondxBQ2(xBbin, Q2bin, k = 0):
			if k ==0:
				return (xBbin==1 and Q2bin == 0) or (xBbin==0 and Q2bin==4)
			else:
				return ~ActiveAny_int[xBbin, Q2bin, :, :].any()

		def badBinCondxBt(xBbin, tbin, k = 0):
			if k ==0:
				return (xBbin==1 and tbin == 0)
			else:
				return ~ActiveAny_int[xBbin, :, tbin, :].any()

		num_plotQ2 = len(Q2bins)
		num_plotxB = len(xBbins)
		num_plott = len(tbins)

		if k == 0:
			num_plotQ2 = 6
			num_plotxB = 3


		for tbin in range(num_plott -1):
			active = 0
			ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
			for xBbin in range(num_plotxB - 1):
				for Q2bin in range(num_plotQ2 - 1):
					#skip inactive bins
					if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
						axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
						axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
						continue
					if ActiveInb[xBbin, Q2bin, tbin, :].any():
						phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'g', label = 'Inb.')
					if ActiveOutb[xBbin, Q2bin, tbin, :].any():
						phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'cyan', label = 'Outb.')

					phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
					axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], linestyle ='', color = 'k', label = 'Merged')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'BH')

					xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
					Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
					theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
					header = xBheader +Q2header + theader
					axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
					axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
					axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 90, 180, 270, 360])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
					if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
						handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
						active = 1
			lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			plt.savefig("plots/binscheme{}/bkgscheme{}tbin{}.pdf".format(k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			plt.clf()

		for tbin in range(num_plott -1):
			active = 0
			ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
			for xBbin in range(num_plotxB - 1):
				for Q2bin in range(num_plotQ2 - 1):
					#skip inactive bins
					if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
						axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
						axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
						continue
					if ActiveInb[xBbin, Q2bin, tbin, :].any():
						phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecInb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'g', label = 'Inb.')
					if ActiveOutb[xBbin, Q2bin, tbin, :].any():
						phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecOutb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'cyan', label = 'Outb.')
					phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
					axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsec_BH[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_KM[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'k', label = 'Merged')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_KM[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :]), color = 'b', label = 'KM15')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :]), color = 'r', label = 'BH')

					xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
					Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
					theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
					header = xBheader +Q2header + theader
					axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
					axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylim([0.01, 100])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 90, 180, 270, 360])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
					if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
						handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
						active = 1
			lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			plt.savefig("plots/binscheme{}/KM_ratio_bkgscheme{}tbin{}.pdf".format(k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			plt.clf()

		for tbin in range(num_plott -1):
			active = 0
			ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
			for xBbin in range(num_plotxB - 1):
				for Q2bin in range(num_plotQ2 - 1):
					#skip inactive bins
					if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
						axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
						axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
						continue
					if ActiveInb[xBbin, Q2bin, tbin, :].any():
						phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecInb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'g', label = 'Inb.')
					if ActiveOutb[xBbin, Q2bin, tbin, :].any():
						phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsecOutb_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'cyan', label = 'Outb.')
					phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
					axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], divideHist(xsec_BH[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = divideHist((xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin], xsecTh_BH[xBbin, Q2bin, tbin, phibin]), linestyle ='', color = 'k', label = 'Merged')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_KM[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :]), color = 'b', label = 'KM15')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], divideHist(xsecTh_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :]), color = 'r', label = 'BH')

					xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
					Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
					theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
					header = xBheader +Q2header + theader
					axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
					axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylim([0.01, 100])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 90, 180, 270, 360])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
					if (active == 0) and ActiveAll[xBbin, Q2bin, tbin, :].any():
						handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
						active = 1
			lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			plt.savefig("plots/binscheme{}/BH_ratio_bkgscheme{}tbin{}.pdf".format(k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			plt.clf()

		Normalization = np.zeros(xsecTh_BH.shape[:-1])
		Normalization_Inb = np.zeros(xsecTh_BH.shape[:-1])
		Normalization_Outb = np.zeros(xsecTh_BH.shape[:-1])
		Normalization_KM = np.zeros(xsecTh_BH.shape[:-1])

		for tbin in range(num_plott -1):
			active = 0
			ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
			for xBbin in range(num_plotxB - 1):
				for Q2bin in range(num_plotQ2 - 1):
					#skip inactive bins
					if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
						axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
						axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
						continue

					if ActiveInb_int[xBbin, Q2bin, tbin, :].any():
						wings_Inb = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsecInb_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
						Normalization_Inb[xBbin, Q2bin, tbin] = np.mean(xsecInb_BH[xBbin, Q2bin, tbin, wings_Inb], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings_Inb], axis = -1)
						phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin], linestyle ='', color = 'g', label = 'Inb.')
					if ActiveOutb_int[xBbin, Q2bin, tbin, :].any():
						wings_Outb = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsecOutb_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
						Normalization_Outb[xBbin, Q2bin, tbin] = np.mean(xsecOutb_BH[xBbin, Q2bin, tbin, wings_Outb], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings_Outb], axis = -1)
						phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin], linestyle ='', color = 'cyan', label = 'Outb.')
					wings = np.argwhere((xsecTh_BH[xBbin, Q2bin, tbin, :]>0)&(xsec_BH[xBbin, Q2bin, tbin, :]>0))[[0,1,-2,-1]].flatten()
					Normalization[xBbin, Q2bin, tbin] = np.mean(xsec_BH[xBbin, Q2bin, tbin, wings], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings], axis = -1)
					Normalization_KM[xBbin, Q2bin, tbin] = np.mean(xsecTh_KM[xBbin, Q2bin, tbin, wings], axis = -1)/np.mean(xsecTh_BH[xBbin, Q2bin, tbin, wings], axis = -1)
					phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
					axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], linestyle ='', color = 'k', label = 'Merged')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_KM[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], xsecTh_BH[xBbin, Q2bin, tbin, :], color = 'r', label = 'BH')

					xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
					Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
					theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
					header = xBheader +Q2header + theader
					axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
					axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
					axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 90, 180, 270, 360])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
					if (active == 0) and ActiveAll_int[xBbin, Q2bin, tbin, :].any():
						handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
						active = 1
			lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			plt.savefig("plots/binscheme{}/Normalized_bkgscheme{}tbin{}.pdf".format(k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			plt.clf()

		IntegratedDiff = np.zeros(xsecTh_BH.shape[:-1])
		IntegratedDiff_Inb = np.zeros(xsecTh_BH.shape[:-1])
		IntegratedDiff_Outb = np.zeros(xsecTh_BH.shape[:-1])
		IntegratedDiff_KM = np.zeros(xsecTh_BH.shape[:-1])

		for tbin in range(num_plott -1):
			active = 0
			ttitle = "{:.3f}".format(tbins[tbin])+r"$<|t|<$"+"{:.3f}".format(tbins[tbin+1])
			fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
			for xBbin in range(num_plotxB - 1):
				for Q2bin in range(num_plotQ2 - 1):
					#skip inactive bins
					if badBinCondxBQ2t(xBbin, Q2bin, tbin, k):
						axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
						axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
						continue

					centers = [8, 9, 10, 11, 12, 13, 14, 15]
					if ActiveInb_int[xBbin, Q2bin, tbin, :].any():
						IntegratedDiff_Inb[xBbin, Q2bin, tbin] = np.sum(xsecInb_BH[xBbin, Q2bin, tbin, centers]/Normalization_Inb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
						phibin = np.argwhere(ActiveInb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecInb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecInb_BH*uncStatInb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Inb[xBbin, Q2bin, tbin], linestyle ='', color = 'g', label = 'Inb.')
					if ActiveOutb_int[xBbin, Q2bin, tbin, :].any():
						IntegratedDiff_Outb[xBbin, Q2bin, tbin] = np.sum(xsecOutb_BH[xBbin, Q2bin, tbin, centers]/Normalization_Outb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
						phibin = np.argwhere(ActiveOutb[xBbin, Q2bin, tbin, :]).flatten()
						axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsecOutb_BH[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsecOutb_BH*uncStatOutb_BH)[xBbin, Q2bin, tbin, phibin]/Normalization_Outb[xBbin, Q2bin, tbin], linestyle ='', color = 'cyan', label = 'Outb.')
					IntegratedDiff[xBbin, Q2bin, tbin] = np.sum(xsec_BH[xBbin, Q2bin, tbin, centers]/Normalization[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
					phibin = np.argwhere(ActiveAny[xBbin, Q2bin, tbin, :]).flatten()
					axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, phibin], xsec_BH[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin] - xsecTh_BH[xBbin, Q2bin, tbin, phibin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, phibin]-phibins[:-1][phibin], phibins[1:][phibin]-phi1avg_BH[xBbin, Q2bin, tbin, phibin]], yerr = (xsec_BH*uncStat_BH)[xBbin, Q2bin, tbin, phibin]/Normalization[xBbin, Q2bin, tbin], linestyle ='', color = 'k', label = 'Merged')
					IntegratedDiff_KM[xBbin, Q2bin, tbin] = np.sum(xsecTh_KM[xBbin, Q2bin, tbin, centers] - xsecTh_BH[xBbin, Q2bin, tbin, centers], axis = -1)
					axs[num_plotQ2-Q2bin-2 , xBbin].plot(phi1avg_BH[xBbin, Q2bin, tbin, :], (xsecTh_KM-xsecTh_BH)[xBbin, Q2bin, tbin, :], color = 'b', label = 'KM15')

					xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0])
					Q2header = r"$<Q^2>=$"+" {:.3f}, ".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
					theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
					header = xBheader +Q2header + theader
					axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
					axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$")
					# axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
					axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 90, 180, 270, 360])
					axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$\phi$" + " [" + degree + "]")
					if (active == 0) and ActiveAll_int[xBbin, Q2bin, tbin, :].any():
						handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
						active = 1
			lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = ttitle, bbox_to_anchor = (1.0, 0.6))
			fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
			plt.savefig("plots/binscheme{}/Differences_bkgscheme{}tbin{}.pdf".format(k, i, tbin), bbox_extra_artists=[lgd], bbox_inches = 'tight')
			plt.clf()

		active = 0
		fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
		for xBbin in range(num_plotxB - 1):
			for Q2bin in range(num_plotQ2 - 1):
				#skip inactive bins
				if badBinCondxBQ2(xBbin, Q2bin, k):
					axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
					axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
					continue
				tbin = np.argwhere(ActiveInb_int[xBbin, Q2bin, :, 0]).flatten()
				axs[num_plotQ2-Q2bin-2 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Inb[xBbin, Q2bin, tbin], color = 'g', label = "Inb", marker = 'o', markersize = 20)
				tbin = np.argwhere(ActiveOutb_int[xBbin, Q2bin, :, 0]).flatten()
				axs[num_plotQ2-Q2bin-2 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Outb[xBbin, Q2bin, tbin], color = 'cyan', label = "Outb.", marker = 'o', markersize = 20)
				tbin = np.argwhere(ActiveAny_int[xBbin, Q2bin, :, 0]).flatten()
				axs[num_plotQ2-Q2bin-2 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff[xBbin, Q2bin, tbin], color = 'k', label = "Merged", marker = 'o', markersize = 20)
				axs[num_plotQ2-Q2bin-2 , xBbin].plot(t1avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_KM[xBbin, Q2bin, tbin], color = 'b', label = "KM", marker = 'o', markersize = 20)

				xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0].mean())
				Q2header = r"$<Q^2>=$"+" {:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0].mean())
				# theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
				header = xBheader +Q2header
				axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
				# axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$"))
				# axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
				axs[num_plotQ2-Q2bin-2, xBbin].set_xlim([0, 2])
				# axs[num_plotQ2-Q2bin-2, xBbin].set_ylim([0.2, 1.])
				axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 0.5, 1, 1.5, 2])
				axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$|t|$" + " [" + GeV2 + "]")
				if (active == 0) and ActiveAll_int[xBbin, Q2bin, :, 0].any():
					handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
					active = 1
		lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = r"$\int (d\sigma-d\sigma_{BH})$", bbox_to_anchor = (1.0, 0.6))
		fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
		plt.savefig("plots/binscheme{}/Integrals_bkgscheme{}.pdf".format(k, i), bbox_extra_artists=[lgd], bbox_inches = 'tight')
		# plt.savefig("plots/binscheme{}/NormScale_bkgscheme{}.pdf".format(k, i, tbin), bbox_inches = 'tight')
		plt.clf()

		active = 0
		fig, axs = plt.subplots(num_plott-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 7.5*(num_plott-1)))
		for xBbin in range(num_plotxB - 1):
			for tbin in range(num_plott - 1):
				#skip inactive bins
				if badBinCondxBt(xBbin, tbin, k):
					axs[num_plott-tbin-2 , xBbin].yaxis.set_visible(False)
					axs[num_plott-tbin-2 , xBbin].xaxis.set_visible(False)
					continue

				Q2bin = np.argwhere(ActiveInb_int[xBbin, :, tbin, 0]).flatten()
				axs[num_plott-tbin-2 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Inb[xBbin, Q2bin, tbin], color = 'g', label = "Inb.", marker = 'o', markersize = 20)
				Q2bin = np.argwhere(ActiveOutb_int[xBbin, :, tbin, 0]).flatten()
				axs[num_plott-tbin-2 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_Outb[xBbin, Q2bin, tbin], color = 'cyan', label = "Outb.", marker = 'o', markersize = 20)
				Q2bin = np.argwhere(ActiveAny_int[xBbin, :, tbin, 0]).flatten()
				axs[num_plott-tbin-2 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff[xBbin, Q2bin, tbin], color = 'k', label = "Merged", marker = 'o', markersize = 20)
				axs[num_plott-tbin-2 , xBbin].plot(Q2avg_BH[xBbin, Q2bin, tbin, 0], IntegratedDiff_KM[xBbin, Q2bin, tbin], color = 'b', label = "KM", marker = 'o', markersize = 20)

				xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0].mean())
				# Q2header = r"$<Q^2>=$"+" {:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0])
				theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0].mean())
				header = xBheader +theader
				axs[num_plott-tbin-2, xBbin].set_title(header, fontsize = 20)
				# axs[num_plott-tbin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$"))
				# axs[num_plott-tbin-2, xBbin].set_yscale('log')
				# axs[num_plott-tbin-2, xBbin].set_xlim([0, 2])
				# axs[num_plott-tbin-2, xBbin].set_ylim([0.2, 1.])
				# axs[num_plott-tbin-2, xBbin].set_xticks([0, 0.5, 1, 1.5, 2])
				axs[num_plott-tbin-2, xBbin].set_xlabel(r"$Q^2$" + " [" + GeVc2 + "]")
				axs[num_plott-tbin-2, xBbin].set_xlim([1, 5])
				axs[num_plott-tbin-2, xBbin].set_xticks([1, 2, 3, 4, 5])
				if (active == 0) and ActiveAll_int[xBbin, Q2bin, :, 0].any():
					handles, labels = axs[num_plott-tbin-2, xBbin].get_legend_handles_labels()
					active = 1
		lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = r"$\int (d\sigma-d\sigma_{BH})$", bbox_to_anchor = (1.0, 0.6))
		fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
		plt.savefig("plots/binscheme{}/Integrals_bkgscheme{}_inQ2.pdf".format(k, i), bbox_extra_artists=[lgd], bbox_inches = 'tight')
		# plt.savefig("plots/binscheme{}/NormScale_bkgscheme{}.pdf".format(k, i, tbin), bbox_inches = 'tight')
		plt.clf()

		active = 0
		fig, axs = plt.subplots(num_plotQ2-1, num_plotxB-1, figsize = (7.5*(num_plotxB-1), 6*(num_plotQ2-1)))
		for xBbin in range(num_plotxB - 1):
			for Q2bin in range(num_plotQ2 - 1):
				#skip inactive bins
				if badBinCondxBQ2(xBbin, Q2bin, k):
					axs[num_plotQ2-Q2bin-2 , xBbin].yaxis.set_visible(False)
					axs[num_plotQ2-Q2bin-2 , xBbin].xaxis.set_visible(False)
					continue
				tbin = np.argwhere(ActiveInb_int[xBbin, Q2bin, :, 0]).flatten()
				axs[num_plotQ2-Q2bin-2 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization_Inb[xBbin, Q2bin, tbin], color = 'g', label = "Inb.")
				tbin = np.argwhere(ActiveOutb_int[xBbin, Q2bin, :, 0]).flatten()
				axs[num_plotQ2-Q2bin-2 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization_Outb[xBbin, Q2bin, tbin], color = 'cyan', label = "Outb.")
				tbin = np.argwhere(ActiveAny_int[xBbin, Q2bin, :, 0]).flatten()
				axs[num_plotQ2-Q2bin-2 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization[xBbin, Q2bin, tbin], color = 'r', label = "Inb.")
				axs[num_plotQ2-Q2bin-2 , xBbin].scatter(t1avg_BH[xBbin, Q2bin, tbin, 0], Normalization_KM[xBbin, Q2bin, tbin], color = 'b', label = "KM")

				# axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, :], Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, :]-phibins[:-1], phibins[1:]-phi1avg_BH[xBbin, Q2bin, tbin, :]], yerr = 0, linestyle ='', color = 'cyan', label = 'Outb.')

				# axs[num_plotQ2-Q2bin-2 , xBbin].errorbar(phi1avg_BH[xBbin, Q2bin, tbin, :], Normalization[xBbin, Q2bin, tbin], xerr = [phi1avg_BH[xBbin, Q2bin, tbin, :]-phibins[:-1], phibins[1:]-phi1avg_BH[xBbin, Q2bin, tbin, :]], yerr = 0, linestyle ='', color = 'k', label = 'Merged')

				xBheader = r"$<x_B>=$"+" {:.3f}, ".format(xBavg_BH[xBbin, Q2bin, tbin, 0].mean())
				Q2header = r"$<Q^2>=$"+" {:.3f}".format(Q2avg_BH[xBbin, Q2bin, tbin, 0].mean())
				# theader = r"$<|t|>=$"+" {:.3f}".format(t1avg_BH[xBbin, Q2bin, tbin, 0])
				header = xBheader +Q2header
				axs[num_plotQ2-Q2bin-2, xBbin].set_title(header, fontsize = 20)
				# axs[num_plotQ2-Q2bin-2, xBbin].set_ylabel(r"$\frac{d\sigma}{dx_B dQ^2 d|t|d\phi}$" + "nb/GeV"+r"$^4$"))
				# axs[num_plotQ2-Q2bin-2, xBbin].set_yscale('log')
				axs[num_plotQ2-Q2bin-2, xBbin].set_xlim([0, 2])
				axs[num_plotQ2-Q2bin-2, xBbin].set_ylim([0.2, 1.])
				axs[num_plotQ2-Q2bin-2, xBbin].set_xticks([0, 0.5, 1, 1.5, 2])
				axs[num_plotQ2-Q2bin-2, xBbin].set_xlabel(r"$|t|$" + " [" + GeV2 + "]")
				if (active == 0) and ActiveAll_int[xBbin, Q2bin, :, 0].any():
					handles, labels = axs[num_plotQ2-Q2bin-2, xBbin].get_legend_handles_labels()
					active = 1
		lgd = plt.figlegend(handles, labels, loc='upper left', fontsize= 20, title = "Norm. to BH", bbox_to_anchor = (1.0, 0.6))
		fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
		plt.savefig("plots/binscheme{}/NormScale_bkgscheme{}.pdf".format(k, i), bbox_extra_artists=[lgd], bbox_inches = 'tight')
		# plt.savefig("plots/binscheme{}/NormScale_bkgscheme{}.pdf".format(k, i, tbin), bbox_inches = 'tight')
		plt.clf()