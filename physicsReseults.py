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
from matplotlib.colors import LogNorm
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
	return np.divide(df1, df2, where = df1>0, out = np.zeros(df1.shape))

def inverseHist(df1):
	return np.divide(np.ones(df1.shape), df1, where = df1>0, out = np.zeros(df1))

# simulation run numbers
runs_inb_vgg50nA = [3987, 4124, 4139, 4181, 4182, 4397, 4528, 4529, 4535, 4539]
runs_inb_vgg55nA = [4186, 4545]
runs_inb_vgg45nA = [4188, 4547]
runs_inb_vgg0nA = [4192, 4561]
runs_inb_bh50nA = [4238, 4542]
runs_inb_bh45nA = [4740, 4742, 4745, 4751, 4760]
runs_inb_bkg50nA = [4076, 4202, 4209]
runs_inb_bkg55nA = [4212]
runs_inb_bkg45nA = [4217]
runs_inb_bkg0nA = [4231]

runs_outb_vgg50nA = [4240, 4250, 4251, 4252, 4255, 4398, 4532, 4534, 4540, 4541, 4717]
runs_outb_vgg40nA = [4263, 4546]
runs_outb_vgg0nA = [4262, 4554]
runs_outb_vgg40nAT = [4266, 4562]
runs_outb_bh50nA = [4249, 4544, 4780, 4808, 4812, 4818]
runs_outb_bkg50nA = [4243, 4271, 4290]
runs_outb_bkg40nA = [4293]
runs_outb_bkg0nA = [4304]
runs_outb_bkg40nAT = [4306]

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

#define the collecion of bin edges
collection_cont_xBbins = [np.linspace(0.05, 0.85, 6)]
collection_cont_Q2bins = [np.array([1, 1.5, 2, 2.5, 3.5, 4.5, 6, 7.5, 10.3])]
collection_cont_tbins = [np.array([0.09, 0.2, 0.4, 0.8, 1.8])]
collection_cont_phibins = [np.linspace(0, 360, 6)]


# trial 0
i = 0

xBbins  = collection_cont_xBbins[i]
Q2bins  = collection_cont_Q2bins[i]
tbins   = collection_cont_tbins [i]
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

epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD.flatten())))
epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD.flatten())))
epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT.flatten())))
epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD.flatten())))
epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD.flatten())))
epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT.flatten())))

epgExp.loc[epgExp.loc[:, "cont{}".format(i)] > 1, "cont{}".format(i)] = 1

print("saved contaminations...")

print("clear memory...")

del df_bkg1gs_inb
del df_bkg1gs_outb
del df_bkg2gs_inb
del df_bkg2gs_outb
gc.collect()


#define the collecion of bin edges
collection_xBbins = [np.linspace(0.05, 0.85, 6)]
collection_Q2bins = [np.array([1, 1.5, 2, 2.5, 3.5, 4.5, 6, 7.5, 10.3])]
collection_tbins = [np.array([0.09, 0.2, 0.4, 0.8, 1.8])]
collection_phibins = [np.linspace(0, 360, 25)]


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


print("reading bhs")

df_bhs_inb = []
for jobNum in runs_inb_bh50nA:
    df_bhs_inb.append(readReduced(parent_MC_BH_inb, jobNum, -1, 50))
for jobNum in runs_inb_bh45nA:
    df_bhs_inb.append(readReduced(parent_MC_BH_inb, jobNum, -1, 45))
    
df_bhs_inb = pd.concat(df_bhs_inb)

print("reading bh Gens")

df_bhs_Gen_inb = []
for jobNum in runs_inb_bh50nA:
    df_bhs_Gen_inb.append(readReduced(parent_Gen_BH_inb, jobNum, -1, 50))
for jobNum in runs_inb_bh45nA:
    df_bhs_Gen_inb.append(readReduced(parent_Gen_BH_inb, jobNum, -1, 45))
    
df_bhs_Gen_inb = pd.concat(df_bhs_Gen_inb)

histBHInb, bins = np.histogramdd(df_bhs_inb.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHGenInb, bins = np.histogramdd(df_bhs_Gen_inb.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHInbFD, bins = np.histogramdd(df_bhs_inb.loc[df_bhs_inb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHGenInbFD, bins = np.histogramdd(df_bhs_Gen_inb.loc[df_bhs_Gen_inb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHInbCD, bins = np.histogramdd(df_bhs_inb.loc[df_bhs_inb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHGenInbCD, bins = np.histogramdd(df_bhs_Gen_inb.loc[df_bhs_Gen_inb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHInbCDFT, bins = np.histogramdd(df_bhs_inb.loc[df_bhs_inb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histBHGenInbCDFT, bins = np.histogramdd(df_bhs_Gen_inb.loc[df_bhs_Gen_inb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

df_vggs_inb = []
for jobNum in runs_inb_vgg50nA:
    df_vggs_inb.append(readReduced(parent_MC_inb, jobNum, -1, 50))
for jobNum in runs_inb_vgg55nA:
    df_vggs_inb.append(readReduced(parent_MC_inb, jobNum, -1, 55))
for jobNum in runs_inb_vgg45nA:
    df_vggs_inb.append(readReduced(parent_MC_inb, jobNum, -1, 45))
for jobNum in runs_inb_vgg0nA:
    df_vggs_inb.append(readReduced(parent_MC_inb, jobNum, -1, 0))
    
df_vggs_inb = pd.concat(df_vggs_inb)

print("reading vgg Gens")

df_vggs_Gen_inb = []
for jobNum in runs_inb_vgg50nA:
    df_vggs_Gen_inb.append(readReduced(parent_Gen_inb, jobNum, -1, 50))
for jobNum in runs_inb_vgg55nA:
    df_vggs_Gen_inb.append(readReduced(parent_Gen_inb, jobNum, -1, 55))
for jobNum in runs_inb_vgg45nA:
    df_vggs_Gen_inb.append(readReduced(parent_Gen_inb, jobNum, -1, 45))
for jobNum in runs_inb_vgg0nA:
    df_vggs_Gen_inb.append(readReduced(parent_Gen_inb, jobNum, -1, 0))
    
df_vggs_Gen_inb = pd.concat(df_vggs_Gen_inb)


histVGGInb, bins = np.histogramdd(df_vggs_inb.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGGenInb, bins = np.histogramdd(df_vggs_Gen_inb.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGInbFD, bins = np.histogramdd(df_vggs_inb.loc[df_vggs_inb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGGenInbFD, bins = np.histogramdd(df_vggs_Gen_inb.loc[df_vggs_Gen_inb.config == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGInbCD, bins = np.histogramdd(df_vggs_inb.loc[df_vggs_inb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGGenInbCD, bins = np.histogramdd(df_vggs_Gen_inb.loc[df_vggs_Gen_inb.config == 2 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGInbCDFT, bins = np.histogramdd(df_vggs_inb.loc[df_vggs_inb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
histVGGGenInbCDFT, bins = np.histogramdd(df_vggs_Gen_inb.loc[df_vggs_Gen_inb.config == 3 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])

#Method 1
accVGGInbFD = divideHist(histVGGInbFD, histVGGGenInbFD)
accVGGInbCD = divideHist(histVGGInbCD, histVGGGenInbCD)
accVGGInbCDFT = divideHist(histVGGInbCDFT, histVGGGenInbCDFT)

