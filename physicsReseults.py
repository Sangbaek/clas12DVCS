#!/usr/bin/env python3
"""
Main Script to save the cross sections
"""
import pandas as pd
import numpy as np


def makeReduced(df):
    columns_needed = ["polarity", "config", "beamCurrent", "xB", "Q2", "t1", "phi1"]
    return df.loc[:, columns_needed]

def readReduced(parent, jobNum, polarity, beamCurrent):
    df = pd.read_pickle(parent + "{}.pkl".format(jobNum))
    df.loc[:, "polarity"] = polarity
    df.loc[:, "beamCurrent"] = beamCurrent
    columns_needed = ["polarity", "config", "beamCurrent", "xB", "Q2", "t1", "phi1"]
    return df.loc[:, columns_needed]


# simulation run numbers
runs_inb_dvcs50nA = [3987, 4124, 4139, 4181, 4182, 4397, 4528, 4529, 4535, 4539]
runs_inb_dvcs55nA = [4186, 4545]
runs_inb_dvcs45nA = [4188, 4547]
runs_inb_dvcs0nA = [4192, 4561]
runs_inb_bh50nA = [4238, 4542]
runs_inb_bh45nA = [4740, 4742, 4745, 4751, 4760]
runs_inb_bkg50nA = [4076, 4202, 4209]
runs_inb_bkg55nA = [4212]
runs_inb_bkg45nA = [4217]
runs_inb_bkg0nA = [4231]

runs_outb_dvcs50nA = [4240, 4250, 4251, 4252, 4255, 4398, 4532, 4534, 4540, 4541, 4717]
runs_outb_dvcs40nA = [4263, 4546]
runs_outb_dvcs0nA = [4262, 4554]
runs_outb_dvcs40nAT = [4266, 4562]
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
    df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 50))
for jobNum in runs_inb_bkg55nA:
    df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 55))
for jobNum in runs_inb_bkg45nA:
    df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 45))
for jobNum in runs_inb_bkg0nA:
    df_bkg1gs_inb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 0))

df_bkg2gs_inb = []
for jobNum in runs_inb_bkg50nA:
    df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 50))
for jobNum in runs_inb_bkg55nA:
    df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 55))
for jobNum in runs_inb_bkg45nA:
    df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 45))
for jobNum in runs_inb_bkg0nA:
    df_bkg2gs_inb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 0))
    
df_bkg1gs_inb = pd.concat(df_bkg1gs_inb)
df_bkg2gs_inb = pd.concat(df_bkg2gs_inb)

df_bkg1gs_outb = []
for jobNum in runs_outb_bkg50nA:
    df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 50))
for jobNum in runs_outb_bkg40nA:
    df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 40))
for jobNum in runs_outb_bkg0nA:
    df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 0))
for jobNum in runs_outb_bkg40nAT:
    df_bkg1gs_outb.append(readReduced(parent_MC_bkg1g, jobNum, -1, 40))

df_bkg2gs_outb = []
for jobNum in runs_outb_bkg50nA:
    df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 50))
for jobNum in runs_outb_bkg40nA:
    df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 40))
for jobNum in runs_outb_bkg0nA:
    df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 0))
for jobNum in runs_outb_bkg40nAT:
    df_bkg2gs_outb.append(readReduced(parent_MC_bkg2g, jobNum, -1, 40))
    
df_bkg1gs_outb = pd.concat(df_bkg1gs_outb)
df_bkg2gs_outb = pd.concat(df_bkg2gs_outb)

#define the collecion of bin edges
collection_xBbins = [np.linspace(0.05, 0.85, 6)]
collection_Q2bins = [np.array([1, 1.5, 2, 2.5, 3.5, 4.5, 6, 7.5, 10.3])]
collection_tbins = [np.array([0.09, 0.2, 0.4, 0.8, 1.8])]
collection_phibins = [np.linspace(0, 360, 6)]


# trial 0
i = 0

xBbins  = collection_xBbins[i]
Q2bins  = collection_Q2bins[i]
tbins   = collection_tbins [i]
phibins = collection_phibin[i]

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

contInbFD = np.divide(histBkgInbFD*histPi0InbFD, histRefInbFD*histExpInbFD, where = histRefInbFD*histExpInbFD>0, out = np.zeros(histRefInbFD.shape))
contInbCD = np.divide(histBkgInbCD*histPi0InbCD, histRefInbCD*histExpInbCD, where = histRefInbCD*histExpInbCD>0, out = np.zeros(histRefInbCD.shape))
contInbCDFT = np.divide(histBkgInbCDFT*histPi0InbCDFT, histRefInbCDFT*histExpInbCDFT, where = histRefInbCDFT*histExpInbCDFT>0, out = np.zeros(histRefInbCDFT.shape))
contOutbFD = np.divide(histBkgOutbFD*histPi0OutbFD, histRefOutbFD*histExpOutbFD, where = histRefOutbFD*histExpOutbFD>0, out = np.zeros(histRefOutbFD.shape))
contOutbCD = np.divide(histBkgOutbCD*histPi0OutbCD, histRefOutbCD*histExpOutbCD, where = histRefOutbCD*histExpOutbCD>0, out = np.zeros(histRefOutbCD.shape))
contOutbCDFT = np.divide(histBkgOutbCDFT*histPi0OutbCDFT, histRefOutbCDFT*histExpOutbCDFT, where = histRefOutbCDFT*histExpOutbCDFT>0, out = np.zeros(histRefOutbCDFT.shape))

epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contInbFD.flatten())))
epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contInbCD.flatten())))
epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == -1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contInbCDFT.flatten())))
epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 1), "newbin{}".format(i)].map(dict(enumerate(contOutbFD.flatten())))
epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 2), "newbin{}".format(i)].map(dict(enumerate(contOutbCD.flatten())))
epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "cont{}".format(i)] = epgExp.loc[(epgExp.polarity == 1) & (epgExp.config == 3), "newbin{}".format(i)].map(dict(enumerate(contOutbCDFT.flatten())))

epgExp.loc[epgExp.loc[:, "cont{}".format(i)] > 1, "cont{}".format(i)] = 1

print("saved contaminations...")