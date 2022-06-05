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
		"figure.autolayout": True
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
    return th_KM15.XS(pt1)

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
        
    dstot = float(dstot.splitlines()[-1].decode("utf-8"))
    return dstot

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

    dstot = float(dstot.splitlines()[-1].decode("utf-8"))
    return dstot

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

def divideHist(df1, df2):
	return np.divide(df1, df2, where = df2>0, out = np.zeros(df2.shape))

def inverseHist(df1):
	return np.divide(np.ones(df1.shape), df1, where = df1>0, out = np.zeros(df1))

def binVolumes(xBbin, Q2bin, tbin, i=0):
	xBbins  = collection_xBbins[i]
	Q2bins  = collection_Q2bins[i]
	tbins   = collection_tbins [i]
	phibins = collection_phibins[i]
	return np.diff(np.radians(phibins))*np.diff(xBbins)[xBbin]*np.diff(Q2bins)[Q2bin]*np.diff(tbins)[tbin]

parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-sc","--skipcont", help="skip cont", action = "store_true")

args = parser.parse_args()

# read exp
parent = "/volatile/clas12/sangbaek/nov2021/"

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

inbcharge_epg, outbcharge_epg = 30398709.057119943, 32085430.046131916

epgExpInb = makeReduced(epgExpInb)
pi0ExpInb = makeReduced(pi0ExpInb)
epgExpOutb = makeReduced(epgExpOutb)
pi0ExpOutb = makeReduced(pi0ExpOutb)

epgExp = pd.concat([epgExpInb, epgExpOutb])
pi0Exp = pd.concat([pi0ExpInb, pi0ExpOutb])

if args.skipcont:
	pass

else:
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
	for xBbins in collection_cont_xBbins:
		for Q2bins in collection_cont_Q2bins:
			for tbins in collection_cont_tbins:
				for phibins in collection_cont_phibins:

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

					epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD.flatten())))
					epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD.flatten())))
					epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT.flatten())))
					epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD.flatten())))
					epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD.flatten())))
					epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT.flatten())))

					epgExp.loc[epgExp.loc[:, "cont{}".format(i)] > 1, "cont{}".format(i)] = 1

					i = i+1

					print("saved contaminations {}".format(i))

	print("clear memory...")

	del df_bkg1gs_inb
	del df_bkg1gs_outb
	del df_bkg2gs_inb
	del df_bkg2gs_outb
	gc.collect()

	epgExp.to_pickle("nphistograms/epgExp.pkl")


epgExp = pd.read_pickle("nphistograms/epgExp.pkl")
# trial 0
i = 0

xBbins  = collection_xBbins[i]
Q2bins  = collection_Q2bins[i]
tbins   = collection_tbins [i]
phibins = collection_phibins[i]

#Inbending cross sections
histBHDVCSInbFD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)])
histBHDVCSInbCD, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)])
histBHDVCSInbCDFT, bins = np.histogramdd(epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1 - epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)])

histBHDVCSInb = histBHDVCSInbFD + histBHDVCSInbCD + histBHDVCSInbCDFT

print("reading bhs")

histBHInb50nA, histBHInbFD50nA, histBHInbCD50nA, histBHInbCDFT50nA = 0, 0, 0, 0
for jobNum in runs_inb_bh50nA:
	histBHInb50nA = histBHInb50nA + np.load("nphistograms/{}Rec.npz".format(jobNum))["hist"]
	histBHInbFD50nA = histBHInbFD50nA + np.load("nphistograms/{}Rec1.npz".format(jobNum))["hist"]
	histBHInbCD50nA = histBHInbCD50nA + np.load("nphistograms/{}Rec2.npz".format(jobNum))["hist"]
	histBHInbCDFT50nA = histBHInbCDFT50nA + np.load("nphistograms/{}Rec3.npz".format(jobNum))["hist"]

histBHInb45nA, histBHInbFD45nA, histBHInbCD45nA, histBHInbCDFT45nA = 0, 0, 0, 0
for jobNum in runs_inb_bh45nA:
	histBHInb45nA = histBHInb45nA + np.load("nphistograms/{}Rec.npz".format(jobNum))["hist"]
	histBHInbFD45nA = histBHInbFD45nA + np.load("nphistograms/{}Rec1.npz".format(jobNum))["hist"]
	histBHInbCD45nA = histBHInbCD45nA + np.load("nphistograms/{}Rec2.npz".format(jobNum))["hist"]
	histBHInbCDFT45nA = histBHInbCDFT45nA + np.load("nphistograms/{}Rec3.npz".format(jobNum))["hist"]

histBHGenInb50nA, histBHGenInbFD50nA, histBHGenInbCD50nA, histBHGenInbCDFT50nA = 0, 0, 0, 0
for jobNum in runs_inb_bh50nA:
	histBHGenInb50nA = histBHGenInb50nA + np.load("nphistograms/{}Gen.npz".format(jobNum))["hist"]
	histBHGenInbFD50nA = histBHGenInbFD50nA + np.load("nphistograms/{}Gen1.npz".format(jobNum))["hist"]
	histBHGenInbCD50nA = histBHGenInbCD50nA + np.load("nphistograms/{}Gen2.npz".format(jobNum))["hist"]
	histBHGenInbCDFT50nA = histBHGenInbCDFT50nA + np.load("nphistograms/{}Gen3.npz".format(jobNum))["hist"]

histBHGenInb45nA, histBHGenInbFD45nA, histBHGenInbCD45nA, histBHGenInbCDFT45nA = 0, 0, 0, 0
for jobNum in runs_inb_bh45nA:
	histBHGenInb45nA = histBHGenInb45nA + np.load("nphistograms/{}Gen.npz".format(jobNum))["hist"]
	histBHGenInbFD45nA = histBHGenInbFD45nA + np.load("nphistograms/{}Gen1.npz".format(jobNum))["hist"]
	histBHGenInbCD45nA = histBHGenInbCD45nA + np.load("nphistograms/{}Gen2.npz".format(jobNum))["hist"]
	histBHGenInbCDFT45nA = histBHGenInbCDFT45nA + np.load("nphistograms/{}Gen3.npz".format(jobNum))["hist"]

print("reading vggs")

histVGGInb50nA, histVGGInbFD50nA, histVGGInbCD50nA, histVGGInbCDFT50nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg50nA:
	histVGGInb50nA = histVGGInb50nA + np.load("nphistograms/{}Rec.npz".format(jobNum))["hist"]
	histVGGInbFD50nA = histVGGInbFD50nA + np.load("nphistograms/{}Rec1.npz".format(jobNum))["hist"]
	histVGGInbCD50nA = histVGGInbCD50nA + np.load("nphistograms/{}Rec2.npz".format(jobNum))["hist"]
	histVGGInbCDFT50nA = histVGGInbCDFT50nA + np.load("nphistograms/{}Rec3.npz".format(jobNum))["hist"]

histVGGInb45nA, histVGGInbFD45nA, histVGGInbCD45nA, histVGGInbCDFT45nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg45nA:
	histVGGInb45nA = histVGGInb45nA + np.load("nphistograms/{}Rec.npz".format(jobNum))["hist"]
	histVGGInbFD45nA = histVGGInbFD45nA + np.load("nphistograms/{}Rec1.npz".format(jobNum))["hist"]
	histVGGInbCD45nA = histVGGInbCD45nA + np.load("nphistograms/{}Rec2.npz".format(jobNum))["hist"]
	histVGGInbCDFT45nA = histVGGInbCDFT45nA + np.load("nphistograms/{}Rec3.npz".format(jobNum))["hist"]

histVGGInb55nA, histVGGInbFD55nA, histVGGInbCD55nA, histVGGInbCDFT55nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg55nA:
	histVGGInb55nA = histVGGInb55nA + np.load("nphistograms/{}Rec.npz".format(jobNum))["hist"]
	histVGGInbFD55nA = histVGGInbFD55nA + np.load("nphistograms/{}Rec1.npz".format(jobNum))["hist"]
	histVGGInbCD55nA = histVGGInbCD55nA + np.load("nphistograms/{}Rec2.npz".format(jobNum))["hist"]
	histVGGInbCDFT55nA = histVGGInbCDFT55nA + np.load("nphistograms/{}Rec3.npz".format(jobNum))["hist"]

histVGGInb0nA, histVGGInbFD0nA, histVGGInbCD0nA, histVGGInbCDFT0nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg0nA:
	histVGGInb0nA = histVGGInb0nA + np.load("nphistograms/{}Rec.npz".format(jobNum))["hist"]
	histVGGInbFD0nA = histVGGInbFD0nA + np.load("nphistograms/{}Rec1.npz".format(jobNum))["hist"]
	histVGGInbCD0nA = histVGGInbCD0nA + np.load("nphistograms/{}Rec2.npz".format(jobNum))["hist"]
	histVGGInbCDFT0nA = histVGGInbCDFT0nA + np.load("nphistograms/{}Rec3.npz".format(jobNum))["hist"]

	
histVGGGenInb50nA, histVGGGenInbFD50nA, histVGGGenInbCD50nA, histVGGGenInbCDFT50nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg50nA:
	histVGGGenInb50nA = histVGGGenInb50nA + np.load("nphistograms/{}Gen.npz".format(jobNum))["hist"]
	histVGGGenInbFD50nA = histVGGGenInbFD50nA + np.load("nphistograms/{}Gen1.npz".format(jobNum))["hist"]
	histVGGGenInbCD50nA = histVGGGenInbCD50nA + np.load("nphistograms/{}Gen2.npz".format(jobNum))["hist"]
	histVGGGenInbCDFT50nA = histVGGGenInbCDFT50nA + np.load("nphistograms/{}Gen3.npz".format(jobNum))["hist"]

histVGGGenInb45nA, histVGGGenInbFD45nA, histVGGGenInbCD45nA, histVGGGenInbCDFT45nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg45nA:
	histVGGGenInb45nA = histVGGGenInb45nA + np.load("nphistograms/{}Gen.npz".format(jobNum))["hist"]
	histVGGGenInbFD45nA = histVGGGenInbFD45nA + np.load("nphistograms/{}Gen1.npz".format(jobNum))["hist"]
	histVGGGenInbCD45nA = histVGGGenInbCD45nA + np.load("nphistograms/{}Gen2.npz".format(jobNum))["hist"]
	histVGGGenInbCDFT45nA = histVGGGenInbCDFT45nA + np.load("nphistograms/{}Gen3.npz".format(jobNum))["hist"]

histVGGGenInb55nA, histVGGGenInbFD55nA, histVGGGenInbCD55nA, histVGGGenInbCDFT55nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg55nA:
	histVGGGenInb55nA = histVGGGenInb55nA + np.load("nphistograms/{}Gen.npz".format(jobNum))["hist"]
	histVGGGenInbFD55nA = histVGGGenInbFD55nA + np.load("nphistograms/{}Gen1.npz".format(jobNum))["hist"]
	histVGGGenInbCD55nA = histVGGGenInbCD55nA + np.load("nphistograms/{}Gen2.npz".format(jobNum))["hist"]
	histVGGGenInbCDFT55nA = histVGGGenInbCDFT55nA + np.load("nphistograms/{}Gen3.npz".format(jobNum))["hist"]

histVGGGenInb0nA, histVGGGenInbFD0nA, histVGGGenInbCD0nA, histVGGGenInbCDFT0nA = 0, 0, 0, 0
for jobNum in runs_inb_vgg0nA:
	histVGGGenInb0nA = histVGGGenInb0nA + np.load("nphistograms/{}Gen.npz".format(jobNum))["hist"]
	histVGGGenInbFD0nA = histVGGGenInbFD0nA + np.load("nphistograms/{}Gen1.npz".format(jobNum))["hist"]
	histVGGGenInbCD0nA = histVGGGenInbCD0nA + np.load("nphistograms/{}Gen2.npz".format(jobNum))["hist"]
	histVGGGenInbCDFT0nA = histVGGGenInbCDFT0nA + np.load("nphistograms/{}Gen3.npz".format(jobNum))["hist"]

CDFTcontribution = divideHist(histBHDVCSInbCDFT*histVGGGenInbCDFT50nA , histVGGInbCDFT50nA)
CDcontribution = divideHist(histBHDVCSInbCD*histVGGGenInbCD50nA , histVGGInbCD50nA)
FDcontribution = divideHist(histBHDVCSInbFD*histVGGGenInbFD50nA , histVGGInbFD50nA)

accCorrected_VGG = FDcontribution + CDcontribution + CDFTcontribution

CDFTcontribution = divideHist(histBHDVCSInbCDFT*histBHGenInbCDFT50nA , histBHInbCDFT50nA)
CDcontribution = divideHist(histBHDVCSInbCD*histBHGenInbCD50nA , histBHInbCD50nA)
FDcontribution = divideHist(histBHDVCSInbFD*histBHGenInbFD50nA , histBHInbFD50nA)

accCorrected_BH = FDcontribution + CDcontribution + CDFTcontribution

xBbin = 1
Q2bin = 2
tbin = 1
binVolume = binVolumes(xBbin, Q2bin, tbin)

histVGGGenInbInt50nA, histVGGGenInbxB50nA, histVGGGenInbQ250nA, histVGGGenInbt150nA = 0, 0, 0, 0
histVGGGenInbphi50nA, histVGGGenInbrad50nA, histVGGGenInbborn50nA = 0, 0, 0

for jobNum in runs_inb_vgg50nA:
	histVGGGenInbInt50nA = histVGGGenInbInt50nA + np.load("nphistograms/{}GenInt.npz".format(jobNum))["hist"]
	histVGGGenInbxB50nA = histVGGGenInbxB50nA + np.load("nphistograms/{}GenxB.npz".format(jobNum))["hist"]
	histVGGGenInbQ250nA = histVGGGenInbQ250nA + np.load("nphistograms/{}GenQ2.npz".format(jobNum))["hist"]
	histVGGGenInbt150nA = histVGGGenInbt150nA + np.load("nphistograms/{}Gent1.npz".format(jobNum))["hist"]
	histVGGGenInbphi50nA = histVGGGenInbphi50nA + np.load("nphistograms/{}Genphi.npz".format(jobNum))["hist"]
	histVGGGenInbrad50nA = histVGGGenInbrad50nA + np.load("nphistograms/{}Genrad.npz".format(jobNum))["hist"]
	histVGGGenInbborn50nA = histVGGGenInbborn50nA + np.load("nphistograms/{}Genborn.npz".format(jobNum))["hist"]

for jobNum in runs_inb_bh45nA:
	histBHGenInbInt45nA = histBHGenInbInt45nA + np.load("nphistograms/{}GenInt.npz".format(jobNum))["hist"]
	histBHGenInbxB45nA = histBHGenInbxB45nA + np.load("nphistograms/{}GenxB.npz".format(jobNum))["hist"]
	histBHGenInbQ245nA = histBHGenInbQ245nA + np.load("nphistograms/{}GenQ2.npz".format(jobNum))["hist"]
	histBHGenInbt145nA = histBHGenInbt145nA + np.load("nphistograms/{}Gent1.npz".format(jobNum))["hist"]
	histBHGenInbphi45nA = histBHGenInbphi45nA + np.load("nphistograms/{}Genphi.npz".format(jobNum))["hist"]
	histBHGenInbrad45nA = histBHGenInbrad45nA + np.load("nphistograms/{}Genrad.npz".format(jobNum))["hist"]
	histBHGenInbborn45nA = histBHGenInbborn45nA + np.load("nphistograms/{}Genborn.npz".format(jobNum))["hist"]


phi1avg_VGG = divideHist(histVGGGenInbphi50nA, histVGGGenInb50nA)[xBbin, Q2bin, tbin, :]
xBavg_VGG = divideHist(histVGGGenInbxB50nA, histVGGGenInbInt50nA)[xBbin, Q2bin, tbin]*np.ones(phi1avg_VGG.shape)
Q2avg_VGG = divideHist(histVGGGenInbQ250nA, histVGGGenInbInt50nA)[xBbin, Q2bin, tbin]*np.ones(phi1avg_VGG.shape)
t1avg_VGG = divideHist(histVGGGenInbt150nA, histVGGGenInbInt50nA)[xBbin, Q2bin, tbin]*np.ones(phi1avg_VGG.shape)
integratedRad_VGG = divideHist(histVGGGenInb45nA, histVGGGenInbrad45nA)[xBbin, Q2bin, tbin, :]
pointBorn_VGG = printVGGarray(xBavg_VGG, Q2avg_VGG, t1avg_VGG, np.radians(phi1avg_VGG), globalfit = True)
rcfactors_VGG = divideHist(integratedRad_VGG, pointBorn_VGG, where = pointBorn_VGG>0, out = np.zeros(pointBorn_VGG.shape))


phi1avg_BH = divideHist(histBHGenInbphi50nA, histBHGenInb50nA)[xBbin, Q2bin, tbin, :]
xBavg_BH = divideHist(histBHGenInbxB50nA, histBHGenInbInt50nA)[xBbin, Q2bin, tbin]*np.ones(phi1avg_BH.shape)
Q2avg_BH = divideHist(histBHGenInbQ250nA, histBHGenInbInt50nA)[xBbin, Q2bin, tbin]*np.ones(phi1avg_BH.shape)
t1avg_BH = divideHist(histBHGenInbt150nA, histBHGenInbInt50nA)[xBbin, Q2bin, tbin]*np.ones(phi1avg_BH.shape)
integratedRad_BH = divideHist(histBHGenInb45nA, histBHGenInbrad45nA)[xBbin, Q2bin, tbin, :]
pointBorn_BH = printBHarray(xBavg_BH, Q2avg_BH, t1avg_BH, np.radians(phi1avg_BH), globalfit = True)
rcfactors_BH = divideHist(integratedRad_BH, pointBorn_BH, where = pointBorn_BH>0, out = np.zeros(pointBorn_BH.shape))

plt.hist(phi1avg_VGG[:-1], phi1avg_VGG, weights = accCorrected_VGG[xBbin, Q2bin, tbin, :]/binVolume/inbcharge_epg/rcfactors_VGG)
plt.hist(phi1avg_BH[:-1], phi1avg_BH, weights = accCorrected_BH[xBbin, Q2bin, tbin, :]/binVolume/inbcharge_epg/rcfactors_BH)

plt.plot(phi1avg, printKMarray(xBavg, Q2avg, t1avg, np.radians(phi1avg)), color = 'b')
plt.plot(phi1avg, printBHarray(xBavg, Q2avg, t1avg, np.radians(phi1avg)), color = 'r')
plt.plot(phi1avg, printVGGarray(xBavg, Q2avg, t1avg, np.radians(phi1avg)), color = 'g')

plt.show()