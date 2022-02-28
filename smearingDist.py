#!/usr/bin/env python3
"""
A script to compare sim-to-dat
"""

import uproot
import pandas as pd
import numpy as np
import argparse
from copy import copy
from utils.const import *
from utils.physics import *
import matplotlib.pyplot as plt
from copy import copy
cmap = copy(plt.cm.get_cmap("jet"))
from scipy.optimize import least_squares
from scipy.stats import entropy

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

import warnings
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

def nphistmean(hist, bins):
    s=0
    for i in range(len(hist)):
        s += hist[i] * ((bins[i] + bins[i+1]) / 2) 
    mean = s / np.sum(hist)
    return mean

class smearingDist():

	def __init__(self, version, outDir = "SimtoDat/v0/", pol = "inbending"):
		self.GetVersion(version)
		if version == "v0":
			self.MakeV0(outDir = outDir)
		if version == "v1":
			self.MakeV1()
		if version == "v2":
			self.MakeV2()
		if version == "v3":
			self.MakeV3()
		if version == "v4":
			self.MakeV4(pol = pol)

	def GetVersion(self, version):
		self.version = version
		self.outDir = "SimtoDat/{}".format(version)

	def MakeV0(self, outDir = "SimtoDat/v0/"):

		if self.version == "v0":
			inDirInb = "../nov2021/convPkl_full/inb/"
			inDirOutb = "../nov2021/convPkl_full/outb/"

		dvcsInbJobs = [	3987]#, 4124, 4139, 4181, 4182, 4186, 4188, 4192 ]#, 4397 ]
		dvcsOutbJobs = [ 4240]#, 4250, 4251, 4252, 4255, 4263, 4262, 4266 ]#, 4398 ]
		pi0InbJobs = [ 4076]#, 4202, 4209, 4212, 4217, 4231 ]
		pi0OutbJobs = [ 4243]#, 4271, 4290, 4293, 4304, 4306 ]

		epgExpInb = pd.read_pickle(inDirInb+"/exp/dvcs.pkl")
		pi0ExpInb = pd.read_pickle(inDirInb+"/exp/pi0.pkl")
		df = []
		for job in dvcsInbJobs:
			df.append(pd.read_pickle(inDirInb+"/dvcs/{}.pkl".format(job)))
		dvcsSimInb = pd.concat(df)
		
		df = []
		for job in pi0InbJobs:
			df.append(pd.read_pickle(inDirInb+"/bkg_1g/{}.pkl".format(job)))
		bkgSimInb = pd.concat(df)
		df = []
		for job in pi0InbJobs:
			df.append(pd.read_pickle(inDirInb+"/bkg_2g/{}.pkl".format(job)))
		pi0SimInb = pd.concat(df)

		epgExpOutb = pd.read_pickle(inDirOutb+"/exp/dvcs.pkl")
		pi0ExpOutb = pd.read_pickle(inDirOutb+"/exp/pi0.pkl")
		
		df = []
		for job in dvcsOutbJobs:
			df.append(pd.read_pickle(inDirOutb+"/dvcs/{}.pkl".format(job)))
		dvcsSimOutb = pd.concat(df)
		
		df = []
		for job in pi0OutbJobs:
			df.append(pd.read_pickle(inDirOutb+"/bkg_1g/{}.pkl".format(job)))
		bkgSimOutb = pd.concat(df)

		df = []
		for job in pi0OutbJobs:
			df.append(pd.read_pickle(inDirOutb+"/bkg_2g/{}.pkl".format(job)))
		pi0SimOutb = pd.concat(df)

		ParticleVars1 = ["Epx", "Epy", "Epz", "Ep", "Etheta", "Ephi", "Esector"]
		ParticleVars2 = ["Ppx", "Ppy", "Ppz", "Pp", "Ptheta", "Pphi", "Psector", "Pstat"]
		ParticleVars3 = ["Gpx", "Gpy", "Gpz", "Gp", "Gtheta", "Gphi", "Gsector", "GFid"]
		ParticleVars4 = ["Gpx2", "Gpy2", "Gpz2", "Gp2", "Gtheta2", "Gphi2", "Gsector2", "GFid2"]
		Exclvars1 = ['event', 'Mpx', 'Mpy', 'Mpz','Q2','nu','y','xB','t1','t2','W','phi1','phi2','MM2_epg','ME_epg','MM2_ep','MM2_eg','MPt','coneAngle','reconGam','coplanarity', 'config']
		Exclvars2 = ['event', 'Mpx', 'Mpy', 'Mpz','Q2','nu','xB','t1','t2','W','phi1','phi2','MPt', 'MM2_ep', 'MM2_egg', 'MM2_epgg', 'ME_epgg', 'Mpi0', 'reconPi', "Pie", 'coplanarity', 'coneAngle1', 'coneAngle2', 'closeness', 'config']
		DVCSvars = ParticleVars1 + ParticleVars2+ ParticleVars3 + Exclvars1
		DVpi0Pvars = ParticleVars1 + ParticleVars2 + ParticleVars3 + ParticleVars4 + Exclvars2

		epgExpInb = epgExpInb.loc[:, DVCSvars]
		dvcsSimInb = dvcsSimInb.loc[:, DVCSvars]
		bkgSimInb = bkgSimInb.loc[:, DVCSvars]
		pi0ExpInb = pi0ExpInb.loc[:, DVpi0Pvars]
		pi0SimInb = pi0SimInb.loc[:, DVpi0Pvars]

		epgExpOutb = epgExpOutb.loc[:, DVCSvars]
		dvcsSimOutb = dvcsSimOutb.loc[:, DVCSvars]
		bkgSimOutb = bkgSimOutb.loc[:, DVCSvars]
		pi0ExpOutb = pi0ExpOutb.loc[:, DVpi0Pvars]
		pi0SimOutb = pi0SimOutb.loc[:, DVpi0Pvars]

		#CDFT
		epgExpInbCDFT = epgExpInb.loc[epgExpInb.config == 3]
		pi0ExpInbCDFT = pi0ExpInb.loc[pi0ExpInb.config == 3]
		dvcsSimInbCDFT = dvcsSimInb.loc[dvcsSimInb.config == 3]
		bkgSimInbCDFT = bkgSimInb.loc[bkgSimInb.config == 3]
		pi0SimInbCDFT = pi0SimInb.loc[pi0SimInb.config == 3]

		epgExpOutbCDFT = epgExpOutb.loc[epgExpOutb.config == 3]
		pi0ExpOutbCDFT = pi0ExpOutb.loc[pi0ExpOutb.config == 3]
		dvcsSimOutbCDFT = dvcsSimOutb.loc[dvcsSimOutb.config == 3]
		bkgSimOutbCDFT = bkgSimOutb.loc[bkgSimOutb.config == 3]
		pi0SimOutbCDFT = pi0SimOutb.loc[pi0SimOutb.config == 3]
		#CD
		epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
		pi0ExpInbCD = pi0ExpInb.loc[pi0ExpInb.config == 2]
		dvcsSimInbCD = dvcsSimInb.loc[dvcsSimInb.config == 2]
		bkgSimInbCD = bkgSimInb.loc[bkgSimInb.config == 2]
		pi0SimInbCD = pi0SimInb.loc[pi0SimInb.config == 2]

		epgExpOutbCD = epgExpOutb.loc[epgExpOutb.config == 2]
		pi0ExpOutbCD = pi0ExpOutb.loc[pi0ExpOutb.config == 2]
		dvcsSimOutbCD = dvcsSimOutb.loc[dvcsSimOutb.config == 2]
		bkgSimOutbCD = bkgSimOutb.loc[bkgSimOutb.config == 2]
		pi0SimOutbCD = pi0SimOutb.loc[pi0SimOutb.config == 2]
		#CDFT
		epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]
		pi0ExpInbFD = pi0ExpInb.loc[pi0ExpInb.config == 1]
		dvcsSimInbFD = dvcsSimInb.loc[dvcsSimInb.config == 1]
		bkgSimInbFD = bkgSimInb.loc[bkgSimInb.config == 1]
		pi0SimInbFD = pi0SimInb.loc[pi0SimInb.config == 1]

		epgExpOutbFD = epgExpOutb.loc[epgExpOutb.config == 1]
		pi0ExpOutbFD = pi0ExpOutb.loc[pi0ExpOutb.config == 1]
		dvcsSimOutbFD = dvcsSimOutb.loc[dvcsSimOutb.config == 1]
		bkgSimOutbFD = bkgSimOutb.loc[bkgSimOutb.config == 1]
		pi0SimOutbFD = pi0SimOutb.loc[pi0SimOutb.config == 1]

		epgExpInbCDFT.to_pickle(outDir+ "epgExpInbCDFT")
		dvcsSimInbCDFT.to_pickle(outDir+ "dvcsSimInbCDFT")
		pi0ExpInbCDFT.to_pickle(outDir+ "pi0ExpInbCDFT")
		pi0SimInbCDFT.to_pickle(outDir+ "pi0SimInbCDFT")
		bkgSimInbCDFT.to_pickle(outDir+ "bkgSimInbCDFT")

		epgExpInbCD.to_pickle(outDir+ "epgExpInbCD")
		dvcsSimInbCD.to_pickle(outDir+ "dvcsSimInbCD")
		pi0ExpInbCD.to_pickle(outDir+ "pi0ExpInbCD")
		pi0SimInbCD.to_pickle(outDir+ "pi0SimInbCD")
		bkgSimInbCD.to_pickle(outDir+ "bkgSimInbCD")

		epgExpInbFD.to_pickle(outDir+ "epgExpInbFD")
		dvcsSimInbFD.to_pickle(outDir+ "dvcsSimInbFD")
		pi0ExpInbFD.to_pickle(outDir+ "pi0ExpInbFD")
		pi0SimInbFD.to_pickle(outDir+ "pi0SimInbFD")
		bkgSimInbFD.to_pickle(outDir+ "bkgSimInbFD")

		epgExpOutbCDFT.to_pickle(outDir+ "epgExpOutbCDFT")
		dvcsSimOutbCDFT.to_pickle(outDir+ "dvcsSimOutbCDFT")
		pi0ExpOutbCDFT.to_pickle(outDir+ "pi0ExpOutbCDFT")
		pi0SimOutbCDFT.to_pickle(outDir+ "pi0SimOutbCDFT")
		bkgSimOutbCDFT.to_pickle(outDir+ "bkgSimOutbCDFT")

		epgExpOutbCD.to_pickle(outDir+ "epgExpOutbCD")
		dvcsSimOutbCD.to_pickle(outDir+ "dvcsSimOutbCD")
		pi0ExpOutbCD.to_pickle(outDir+ "pi0ExpOutbCD")
		pi0SimOutbCD.to_pickle(outDir+ "pi0SimOutbCD")
		bkgSimOutbCD.to_pickle(outDir+ "bkgSimOutbCD")

		epgExpOutbFD.to_pickle(outDir+ "epgExpOutbFD")
		dvcsSimOutbFD.to_pickle(outDir+ "dvcsSimOutbFD")
		pi0ExpOutbFD.to_pickle(outDir+ "pi0ExpOutbFD")
		pi0SimOutbFD.to_pickle(outDir+ "pi0SimOutbFD")
		bkgSimOutbFD.to_pickle(outDir+ "bkgSimOutbFD")

	def MakeV1Exp(self, inDir = "SimtoDat/v0", outDir = "SimtoDat/v1/"):
		epgExpInbCDFT = pd.read_pickle(inDir+"/epgExpInbCDFT")
		epgExpOutbCDFT = pd.read_pickle(inDir+"/epgExpOutbCDFT")
		#divide into fd photon energy range
		epgExpInbCDFT0 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>2)&(epgExpInbCDFT.Ge<3)]
		epgExpInbCDFT1 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>3)&(epgExpInbCDFT.Ge<3.5)]
		epgExpInbCDFT2 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>3.5)&(epgExpInbCDFT.Ge<4)]
		epgExpInbCDFT3 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>4)&(epgExpInbCDFT.Ge<4.5)]
		epgExpInbCDFT4 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>4.5)&(epgExpInbCDFT.Ge<5)]
		epgExpInbCDFT5 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>5)&(epgExpInbCDFT.Ge<5.5)]
		epgExpInbCDFT6 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>5.5)&(epgExpInbCDFT.Ge<6)]
		epgExpInbCDFT7 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>6)&(epgExpInbCDFT.Ge<6.5)]
		epgExpInbCDFT8 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>6.5)&(epgExpInbCDFT.Ge<7)]
		epgExpInbCDFT9 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>7)&(epgExpInbCDFT.Ge<7.5)]
		epgExpInbCDFT10 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>7.5)&(epgExpInbCDFT.Ge<8)]
		epgExpInbCDFT11 = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>8)&(epgExpInbCDFT.Ge<9)]

		epgExpOutbCDFT0 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>2)&(epgExpOutbCDFT.Ge<3)]
		epgExpOutbCDFT1 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>3)&(epgExpOutbCDFT.Ge<3.5)]
		epgExpOutbCDFT2 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>3.5)&(epgExpOutbCDFT.Ge<4)]
		epgExpOutbCDFT3 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>4)&(epgExpOutbCDFT.Ge<4.5)]
		epgExpOutbCDFT4 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>4.5)&(epgExpOutbCDFT.Ge<5)]
		epgExpOutbCDFT5 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>5)&(epgExpOutbCDFT.Ge<5.5)]
		epgExpOutbCDFT6 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>5.5)&(epgExpOutbCDFT.Ge<6)]
		epgExpOutbCDFT7 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>6)&(epgExpOutbCDFT.Ge<6.5)]
		epgExpOutbCDFT8 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>6.5)&(epgExpOutbCDFT.Ge<7)]
		epgExpOutbCDFT9 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>7)&(epgExpOutbCDFT.Ge<7.5)]
		epgExpOutbCDFT10 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>7.5)&(epgExpOutbCDFT.Ge<8)]
		epgExpOutbCDFT11 = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>8)&(epgExpOutbCDFT.Ge<9)]

		epgExpInbCDFT0.to_pickle(outDir+ "epgExpInbCDFT0")
		epgExpInbCDFT1.to_pickle(outDir+ "epgExpInbCDFT1")
		epgExpInbCDFT2.to_pickle(outDir+ "epgExpInbCDFT2")
		epgExpInbCDFT3.to_pickle(outDir+ "epgExpInbCDFT3")
		epgExpInbCDFT4.to_pickle(outDir+ "epgExpInbCDFT4")
		epgExpInbCDFT5.to_pickle(outDir+ "epgExpInbCDFT5")
		epgExpInbCDFT6.to_pickle(outDir+ "epgExpInbCDFT6")
		epgExpInbCDFT7.to_pickle(outDir+ "epgExpInbCDFT7")
		epgExpInbCDFT8.to_pickle(outDir+ "epgExpInbCDFT8")
		epgExpInbCDFT9.to_pickle(outDir+ "epgExpInbCDFT9")
		epgExpInbCDFT10.to_pickle(outDir+ "epgExpInbCDFT10")
		epgExpInbCDFT11.to_pickle(outDir+ "epgExpInbCDFT11")

		epgExpOutbCDFT0.to_pickle(outDir+ "epgExpOutbCDFT0")
		epgExpOutbCDFT1.to_pickle(outDir+ "epgExpOutbCDFT1")
		epgExpOutbCDFT2.to_pickle(outDir+ "epgExpOutbCDFT2")
		epgExpOutbCDFT3.to_pickle(outDir+ "epgExpOutbCDFT3")
		epgExpOutbCDFT4.to_pickle(outDir+ "epgExpOutbCDFT4")
		epgExpOutbCDFT5.to_pickle(outDir+ "epgExpOutbCDFT5")
		epgExpOutbCDFT6.to_pickle(outDir+ "epgExpOutbCDFT6")
		epgExpOutbCDFT7.to_pickle(outDir+ "epgExpOutbCDFT7")
		epgExpOutbCDFT8.to_pickle(outDir+ "epgExpOutbCDFT8")
		epgExpOutbCDFT9.to_pickle(outDir+ "epgExpOutbCDFT9")
		epgExpOutbCDFT10.to_pickle(outDir+ "epgExpOutbCDFT10")
		epgExpOutbCDFT11.to_pickle(outDir+ "epgExpOutbCDFT11")

	def MakeV1(self, inDir = "SimtoDat/v0/", outDir = "SimtoDat/v1/", exp = None):

		if exp:
			self.MakeV1Exp(inDir, outDir)
			exit()

		# binsMEepgInb = np.linspace(-0.439, 0.484, 101)
		# binsMM2egInb = np.linspace(0.246, 1.569, 101)
		# binsMEepgOutb = np.linspace(-0.436, 0.481, 101) 
		# binsMM2egOutb = np.linspace(0.223, 1.602, 101) 

		def distance(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1, _ = np.histogram(df1.loc[:, var], bins = bins)
			unchist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			unchist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1 = np.sqrt(unchist1)
			unchist2 = np.sqrt(unchist2)
			unchist_exp = np.sqrt(unchist_exp)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			uncdist1 = unchist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			uncdist2 = unchist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist_exp = unchist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist = np.sqrt((1-cont)**2 * uncdist1**2 + cont**2 * uncdist2 **2 + unchist_exp**2)
			uncdist = np.where(uncdist>0, uncdist, np.inf)
			distance = np.sum(((1-cont)*dist1 + cont*dist2 -dist_exp)**2)
			return distance


		def corr(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			simDist = (1-cont)*dist1 + cont*dist2
			expDist = dist_exp
			simPeak = bincenters[np.argmax(simDist)]
			expPeak = bincenters[np.argmax(expDist)]
			return expPeak - simPeak

		# if isinstance(sigma, str):
		# 	sigma = float(sigma)

		GeEdges = [3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 9]
		sigmas = np.linspace(0.013, 0.018, 6)
		corrections = []
		sigmas_opt = []

		epgExpInbCDFT = pd.read_pickle(inDir + "/epgExpInbCDFT")
		epgExpOutbCDFT = pd.read_pickle(inDir + "/epgExpOutbCDFT")

		pi0ExpInbCDFT = pd.read_pickle(inDir+"/pi0ExpInbCDFT")
		pi0ExpOutbCDFT = pd.read_pickle(inDir+"/pi0ExpOutbCDFT")

		dvcsSimInbCDFT = pd.read_pickle(inDir+"/dvcsSimInbCDFT")
		dvcsSimOutbCDFT = pd.read_pickle(inDir+"/dvcsSimOutbCDFT")

		bkgSimInbCDFT = pd.read_pickle(inDir+"/bkgSimInbCDFT")
		bkgSimOutbCDFT = pd.read_pickle(inDir+"/bkgSimOutbCDFT")

		pi0SimInbCDFT = pd.read_pickle(inDir+"/pi0SimInbCDFT")
		pi0SimOutbCDFT = pd.read_pickle(inDir+"/pi0SimOutbCDFT")


		# #performing correcting
		# correction = 0.186/(1+np.exp(1.254*(epgExpInbCDFT.Gp-4.5)))
		# self.CorrectionV0(epgExpInbCDFT, correction, mode = "epg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "inbending")
		# epgExpInbCDFT = self.df_epg
		# correction = 0.186/(1+np.exp(1.254*(pi0ExpInbCDFT.Gp-4.5)))
		# self.CorrectionV0(pi0ExpInbCDFT, correction, mode = "epgg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "inbending")
		# pi0ExpInbCDFT = self.df_epg
		# #performing correcting
		# correction = 0.186/(1+np.exp(1.254*(epgExpOutbCDFT.Gp-4.5)))
		# self.CorrectionV0(epgExpOutbCDFT, correction, mode = "epg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "outbending")
		# epgExpOutbCDFT = self.df_epg
		# correction = 0.186/(1+np.exp(1.254*(pi0ExpOutbCDFT.Gp-4.5)))
		# self.CorrectionV0(pi0ExpOutbCDFT, correction, mode = "epgg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "outbending")
		# pi0ExpOutbCDFT = self.df_epg

		# epgExpInbCDFT.to_pickle(outDir + "/epgExpInbCDFT")
		# epgExpOutbCDFT.to_pickle(outDir + "/epgExpOutbCDFT")

		# pi0ExpInbCDFT.to_pickle(outDir+"/pi0ExpInbCDFT")
		# pi0ExpOutbCDFT.to_pickle(outDir+"/pi0ExpOutbCDFT")


		for i in range(len(GeEdges)-1):

			distances = []
			GeMin = GeEdges[i]
			GeMax = GeEdges[i+1]

			epgExpInbCDFT_selected = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>GeMin) & (epgExpInbCDFT.Ge<GeMax)]
			pi0ExpInbCDFT_selected = pi0ExpInbCDFT.loc[(pi0ExpInbCDFT.Ge>GeMin) & (pi0ExpInbCDFT.Ge<GeMax)]
			dvcsSimInbCDFT_selected = dvcsSimInbCDFT.loc[(dvcsSimInbCDFT.Ge>GeMin) & (dvcsSimInbCDFT.Ge<GeMax)]
			pi0SimInbCDFT_selected = pi0SimInbCDFT.loc[(pi0SimInbCDFT.Ge>GeMin) & (pi0SimInbCDFT.Ge<GeMax)]
			bkgSimInbCDFT_selected = bkgSimInbCDFT.loc[(bkgSimInbCDFT.Ge>GeMin) & (bkgSimInbCDFT.Ge<GeMax)]

			epgExpOutbCDFT_selected = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>GeMin) & (epgExpOutbCDFT.Ge<GeMax)]
			pi0ExpOutbCDFT_selected = pi0ExpOutbCDFT.loc[(pi0ExpOutbCDFT.Ge>GeMin) & (pi0ExpOutbCDFT.Ge<GeMax)]
			dvcsSimOutbCDFT_selected = dvcsSimOutbCDFT.loc[(dvcsSimOutbCDFT.Ge>GeMin) & (dvcsSimOutbCDFT.Ge<GeMax)]
			pi0SimOutbCDFT_selected = pi0SimOutbCDFT.loc[(pi0SimOutbCDFT.Ge>GeMin) & (pi0SimOutbCDFT.Ge<GeMax)]
			bkgSimOutbCDFT_selected = bkgSimOutbCDFT.loc[(bkgSimOutbCDFT.Ge>GeMin) & (bkgSimOutbCDFT.Ge<GeMax)]

			contInb = 0
			contOutb = 0
			if len(epgExpInbCDFT_selected)*len(pi0SimInbCDFT_selected) > 0:
				contInb = len(bkgSimInbCDFT_selected)/len(pi0SimInbCDFT_selected)*len(pi0ExpInbCDFT_selected)/len(epgExpInbCDFT_selected)
			if len(epgExpInbCDFT_selected)*len(pi0SimOutbCDFT_selected) > 0:
				contOutb = len(bkgSimOutbCDFT_selected)/len(pi0SimOutbCDFT_selected)*len(pi0ExpOutbCDFT_selected)/len(epgExpOutbCDFT_selected)

			correction1 = epgExpInbCDFT_selected.ME_epg.mean() - (1-contInb)*dvcsSimInbCDFT_selected.ME_epg.mean() - contInb*bkgSimInbCDFT_selected.ME_epg.mean()
			correction2 = epgExpOutbCDFT_selected.ME_epg.mean() - (1-contOutb)*dvcsSimOutbCDFT_selected.ME_epg.mean() - contOutb*bkgSimOutbCDFT_selected.ME_epg.mean()
			# correction1 = corr(dvcsSimInbCDFT_selected, bkgSimInbCDFT_selected, epgExpInbCDFT_selected, cont = contInb, var = "ME_epg")
			# correction2 = corr(dvcsSimOutbCDFT_selected, bkgSimOutbCDFT_selected, epgExpOutbCDFT_selected, cont = contOutb, var = "ME_epg")
			correction = (correction1+correction2)/2
			print(correction1, contInb, correction2, contOutb)
			corrections.append(correction)

			#performing correcting
			self.CorrectionV0(epgExpInbCDFT, correction, mode = "epg")
			self.saveDVCSvars()
			self.makeDVCS()
			epgExpInbCDFT_corrected = self.df_epg
			self.CorrectionV0(pi0ExpInbCDFT, correction, mode = "epgg")
			self.saveDVCSvars()
			self.makeDVCS()
			pi0ExpInbCDFT_corrected = self.df_epg

			#performing correcting
			self.CorrectionV0(epgExpOutbCDFT, correction, mode = "epg")
			self.saveDVCSvars()
			self.makeDVCS(pol = "outbending")
			epgExpOutbCDFT_corrected = self.df_epg
			self.CorrectionV0(pi0ExpOutbCDFT, correction, mode = "epgg")
			self.saveDVCSvars()
			self.makeDVCS(pol = "outbending")
			pi0ExpOutbCDFT_corrected = self.df_epg

			epgExpInbCDFT_corrected = epgExpInbCDFT_corrected.loc[(epgExpInbCDFT_corrected.Ge>GeMin) & (epgExpInbCDFT_corrected.Ge<GeMax)]
			epgExpOutbCDFT_corrected = epgExpOutbCDFT_corrected.loc[(epgExpOutbCDFT_corrected.Ge>GeMin) & (epgExpOutbCDFT_corrected.Ge<GeMax)]

			# epgExpInbCDFT_corrected = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>GeMin) & (epgExpInbCDFT.Ge<GeMax)]
			# epgExpOutbCDFT_corrected = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>GeMin) & (epgExpOutbCDFT.Ge<GeMax)]


			for sigma in sigmas:

				print("smearing with {:.3f}".format(sigma))

				#performing smearing
				self.SmearingV0(dvcsSimInbCDFT, sigma, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "inbending")
				dvcsSimInbCDFT_smeared = self.df_epg

				self.SmearingV0(bkgSimInbCDFT, sigma, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "inbending")
				bkgSimInbCDFT_smeared = self.df_epg

				self.SmearingV0(dvcsSimOutbCDFT, sigma, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "outbending")
				dvcsSimOutbCDFT_smeared = self.df_epg

				self.SmearingV0(bkgSimOutbCDFT, sigma, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "outbending")
				bkgSimOutbCDFT_smeared = self.df_epg

				dvcsSimInbCDFT_smeared = dvcsSimInbCDFT_smeared.loc[(dvcsSimInbCDFT_smeared.Ge>GeMin) & (dvcsSimInbCDFT_smeared.Ge<GeMax)]
				dvcsSimOutbCDFT_smeared = dvcsSimOutbCDFT_smeared.loc[(dvcsSimOutbCDFT_smeared.Ge>GeMin) & (dvcsSimOutbCDFT_smeared.Ge<GeMax)]

				distance1 = distance(dvcsSimInbCDFT_smeared, bkgSimInbCDFT_smeared, epgExpInbCDFT_corrected, cont = contInb, var = "ME_epg")
				distance2 = distance(dvcsSimOutbCDFT_smeared, bkgSimOutbCDFT_smeared, epgExpOutbCDFT_corrected, cont = contOutb, var = "ME_epg")

				distances.append(( distance1+ distance2 )/2) 

			sigma_opt = sigmas[np.argmin(distances)]
			sigmas_opt.append(sigma_opt)
			self.SmearingV0(dvcsSimInbCDFT, sigma_opt, mode = "epg")
			self.saveDVCSvars()
			self.makeDVCS(pol = "inbending")
			dvcsSimInbCDFT_opt = self.df_epg

			self.SmearingV0(dvcsSimOutbCDFT, sigma_opt, mode = "epg")
			self.saveDVCSvars()
			self.makeDVCS(pol = "outbending")
			dvcsSimOutbCDFT_opt = self.df_epg
			dvcsSimInbCDFT_opt = dvcsSimInbCDFT_opt.loc[(dvcsSimInbCDFT_opt.Ge>GeMin) & (dvcsSimInbCDFT_opt.Ge<GeMax)]
			dvcsSimOutbCDFT_opt = dvcsSimOutbCDFT_opt.loc[(dvcsSimOutbCDFT_opt.Ge>GeMin) & (dvcsSimOutbCDFT_opt.Ge<GeMax)]
			print(len(dvcsSimInbCDFT_opt), len(dvcsSimOutbCDFT_opt))
			print(len(epgExpInbCDFT_corrected), len(epgExpOutbCDFT_corrected))
			print(GeMin, GeMax, distances, sigma_opt, correction) 

			varstoplot = ["Gp", "Gtheta", "Gphi", "coneAngle",  "reconGam", "MPt", "ME_epg", "MM2_epg", "MM2_eg", "coplanarity"]
			title = [r"$p_{\gamma}$", r"$\theta_{\gamma}$", r"$\phi_{\gamma}$", r"$\theta_{e'\gamma}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", "MPt"+r"${}_{epg}$", "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{eg}$", r"$\Delta \phi$"]
			unit = [GeV, degree, degree, degree, degree, GeV, GeV, GeV2, GeV2, degree]
			# binstarts = [GeMin, 0, -180, 0, 0, 0, -0.5, -0.01, 0.1, 0]
			# binends = [GeMax, 7, 180, 30, 0.8, .1, 1.2, 0.01, 1.7, 10]

			fig, axs = plt.subplots(2, 5, figsize = (15,10))
			for yind in range(0, 2):
			    for xind in range(0,5):
			        ind = 5*yind + xind
			        # start = binstarts[ind]
			        # end = binends[ind]
			        # bins = np.linspace(start, end, 101)
			        simDist_dvcs, bins = np.histogram(dvcsSimInbCDFT_opt[varstoplot[ind]], 100, density = True)
			        simDist_dvpi0, bins = np.histogram(bkgSimInbCDFT[varstoplot[ind]], bins, density = True)
			        simDist = (1-contInb)*simDist_dvcs + contInb*simDist_dvpi0
			        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
			        axs[yind, xind].hist(epgExpInbCDFT_corrected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
			        axs[yind, xind].set_title(title[ind])
			        # axs[yind, xind].set_xlim([start, end])
			        if (unit[ind]):
			            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
			        else:
			            axs[yind, xind].set_xlabel(title[ind])
			plt.tight_layout()
			plt.savefig(outDir+"InbCDFT{}_{:.3f}.pdf".format(i, sigma_opt))
			plt.clf()

			fig, axs = plt.subplots(2, 5, figsize = (15,10))
			for yind in range(0, 2):
			    for xind in range(0,5):
			        ind = 5*yind + xind
			        # start = binstarts[ind]
			        # end = binends[ind]
			        # bins = np.linspace(start, end, 101)
			        simDist_dvcs, bins = np.histogram(dvcsSimOutbCDFT_opt[varstoplot[ind]], 100, density = True)
			        simDist_dvpi0, bins = np.histogram(bkgSimOutbCDFT[varstoplot[ind]], bins, density = True)
			        simDist = (1-contOutb)*simDist_dvcs + contOutb*simDist_dvpi0
			        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
			        axs[yind, xind].hist(epgExpOutbCDFT_corrected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
			        axs[yind, xind].set_title(title[ind])
			        # axs[yind, xind].set_xlim([start, end])
			        if (unit[ind]):
			            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
			        else:
			            axs[yind, xind].set_xlabel(title[ind])
			plt.tight_layout()
			plt.savefig(outDir+"OutbCDFT{}_{:.3f}.pdf".format(i, sigma_opt))

		print(sigmas_opt, corrections)

	def CorrectionV0(self, df, correction, mode = "epg"):
		df_epg = copy(df)
		if mode == "epg":
			df_epg.loc[df_epg.Gsector>7, 'Gp'] = df_epg.Gp + correction
			df_epg.loc[df_epg.Gsector>7, 'Ge'] = df_epg.loc[df_epg.Gsector>7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))
		elif mode == "epgg":
			df_epg.loc[df_epg.Gsector>7, 'Gp'] = df_epg.Gp + correction
			df_epg.loc[df_epg.Gsector>7, 'Ge'] = df_epg.loc[df_epg.Gsector>7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))			
		self.df_epg = df_epg

	def SmearingV0(self, df, sigma, mode = "epg"):
		df_epg = copy(df)
		if mode == "epg":
			df_epg.loc[df_epg.Gsector>7, 'Gp'] = np.random.normal(1, sigma, len(df_epg.loc[df_epg.Gsector>7]))*df_epg.loc[df_epg.Gsector>7, 'Gp']
			df_epg.loc[df_epg.Gsector>7, 'Ge'] = df_epg.loc[df_epg.Gsector>7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))
		if mode == "epgg":
			df_epg.loc[df_epg.Gsector>7, 'Gp'] = np.random.normal(1, sigma, len(df_epg.loc[df_epg.Gsector>7]))*df_epg.loc[df_epg.Gsector>7, 'Gp']
			df_epg.loc[df_epg.Gsector>7, 'Ge'] = df_epg.loc[df_epg.Gsector>7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))

			df_epg.loc[df_epg.Gsector>7, 'Gp2'] = np.random.normal(1, 0.014, len(df_epg.loc[df_epg.Gsector>7]))*df_epg.loc[df_epg.Gsector>7, 'Gp2']
			df_epg.loc[df_epg.Gsector>7, 'Ge2'] = df_epg.loc[df_epg.Gsector>7, 'Gp2']
			df_epg.loc[:, "Gpx2"] = df_epg.loc[:, "Gp2"]*np.sin(np.radians(df_epg.loc[:, "Gtheta2"]))*np.cos(np.radians(df_epg.loc[:, "Gphi2"]))
			df_epg.loc[:, "Gpy2"] = df_epg.loc[:, "Gp2"]*np.sin(np.radians(df_epg.loc[:, "Gtheta2"]))*np.sin(np.radians(df_epg.loc[:, "Gphi2"]))
			df_epg.loc[:, "Gpz2"] = df_epg.loc[:, "Gp2"]*np.cos(np.radians(df_epg.loc[:, "Gtheta2"]))

		self.df_epg = df_epg

	def MakeV2(self, inDir = "SimtoDat/v1/", outDir = "SimtoDat/v2/"):

		def distance(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1, _ = np.histogram(df1.loc[:, var], bins = bins)
			unchist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			unchist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1 = np.sqrt(unchist1)
			unchist2 = np.sqrt(unchist2)
			unchist_exp = np.sqrt(unchist_exp)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			uncdist1 = unchist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			uncdist2 = unchist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist_exp = unchist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist = np.sqrt((1-cont)**2 * uncdist1**2 + cont**2 * uncdist2 **2 + unchist_exp**2)
			uncdist = np.where(uncdist>0, uncdist, np.inf)
			distance = np.sum(((1-cont)*dist1 + cont*dist2 -dist_exp)**2)
			return distance


		PpEdges = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
		sigma1s = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08]
		sigma2s = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
		# PthetaEdges = [40, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65]
		# sigma3s = np.linspace(0, 3, 31)
		# corrections = []
		sigma1s_temp = []
		sigma2s_temp = []
		# sigma3s_temp = []

		sigmas_opt = []

		epgExpInbCDFT = pd.read_pickle(inDir + "/epgExpInbCDFT")
		epgExpOutbCDFT = pd.read_pickle(inDir + "/epgExpOutbCDFT")

		pi0ExpInbCDFT = pd.read_pickle(inDir+"/pi0ExpInbCDFT")
		pi0ExpOutbCDFT = pd.read_pickle(inDir+"/pi0ExpOutbCDFT")

		dvcsSimInbCDFT = pd.read_pickle(inDir+"/dvcsSimInbCDFT")
		dvcsSimOutbCDFT = pd.read_pickle(inDir+"/dvcsSimOutbCDFT")

		bkgSimInbCDFT = pd.read_pickle(inDir+"/bkgSimInbCDFT")
		bkgSimOutbCDFT = pd.read_pickle(inDir+"/bkgSimOutbCDFT")

		pi0SimInbCDFT = pd.read_pickle(inDir+"/pi0SimInbCDFT")
		pi0SimOutbCDFT = pd.read_pickle(inDir+"/pi0SimOutbCDFT")

		# corrections = []

		for i in range(len(PpEdges)-1):
		# for i in range(len(PthetaEdges)-1):

			distances1 = []
			distances2 = []
			# distances3 = []
			PpMin = PpEdges[i]
			PpMax = PpEdges[i+1]
			# PthetaMin = PthetaEdges[i]
			# PthetaMax = PthetaEdges[i+1]

			epgExpInbCDFT_selected = epgExpInbCDFT.loc[(epgExpInbCDFT.Pp>PpMin) & (epgExpInbCDFT.Pp<PpMax)]
			pi0ExpInbCDFT_selected = pi0ExpInbCDFT.loc[(pi0ExpInbCDFT.Pp>PpMin) & (pi0ExpInbCDFT.Pp<PpMax)]
			dvcsSimInbCDFT_selected = dvcsSimInbCDFT.loc[(dvcsSimInbCDFT.Pp>PpMin) & (dvcsSimInbCDFT.Pp<PpMax)]
			pi0SimInbCDFT_selected = pi0SimInbCDFT.loc[(pi0SimInbCDFT.Pp>PpMin) & (pi0SimInbCDFT.Pp<PpMax)]
			bkgSimInbCDFT_selected = bkgSimInbCDFT.loc[(bkgSimInbCDFT.Pp>PpMin) & (bkgSimInbCDFT.Pp<PpMax)]

			epgExpOutbCDFT_selected = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Pp>PpMin) & (epgExpOutbCDFT.Pp<PpMax)]
			pi0ExpOutbCDFT_selected = pi0ExpOutbCDFT.loc[(pi0ExpOutbCDFT.Pp>PpMin) & (pi0ExpOutbCDFT.Pp<PpMax)]
			dvcsSimOutbCDFT_selected = dvcsSimOutbCDFT.loc[(dvcsSimOutbCDFT.Pp>PpMin) & (dvcsSimOutbCDFT.Pp<PpMax)]
			pi0SimOutbCDFT_selected = pi0SimOutbCDFT.loc[(pi0SimOutbCDFT.Pp>PpMin) & (pi0SimOutbCDFT.Pp<PpMax)]
			bkgSimOutbCDFT_selected = bkgSimOutbCDFT.loc[(bkgSimOutbCDFT.Pp>PpMin) & (bkgSimOutbCDFT.Pp<PpMax)]

			# epgExpInbCDFT_selected = epgExpInbCDFT.loc[(epgExpInbCDFT.Ptheta>PthetaMin) & (epgExpInbCDFT.Ptheta<PthetaMax)]
			# pi0ExpInbCDFT_selected = pi0ExpInbCDFT.loc[(pi0ExpInbCDFT.Ptheta>PthetaMin) & (pi0ExpInbCDFT.Ptheta<PthetaMax)]
			# dvcsSimInbCDFT_selected = dvcsSimInbCDFT.loc[(dvcsSimInbCDFT.Ptheta>PthetaMin) & (dvcsSimInbCDFT.Ptheta<PthetaMax)]
			# pi0SimInbCDFT_selected = pi0SimInbCDFT.loc[(pi0SimInbCDFT.Ptheta>PthetaMin) & (pi0SimInbCDFT.Ptheta<PthetaMax)]
			# bkgSimInbCDFT_selected = bkgSimInbCDFT.loc[(bkgSimInbCDFT.Ptheta>PthetaMin) & (bkgSimInbCDFT.Ptheta<PthetaMax)]

			# epgExpOutbCDFT_selected = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ptheta>PthetaMin) & (epgExpOutbCDFT.Ptheta<PthetaMax)]
			# pi0ExpOutbCDFT_selected = pi0ExpOutbCDFT.loc[(pi0ExpOutbCDFT.Ptheta>PthetaMin) & (pi0ExpOutbCDFT.Ptheta<PthetaMax)]
			# dvcsSimOutbCDFT_selected = dvcsSimOutbCDFT.loc[(dvcsSimOutbCDFT.Ptheta>PthetaMin) & (dvcsSimOutbCDFT.Ptheta<PthetaMax)]
			# pi0SimOutbCDFT_selected = pi0SimOutbCDFT.loc[(pi0SimOutbCDFT.Ptheta>PthetaMin) & (pi0SimOutbCDFT.Ptheta<PthetaMax)]
			# bkgSimOutbCDFT_selected = bkgSimOutbCDFT.loc[(bkgSimOutbCDFT.Ptheta>PthetaMin) & (bkgSimOutbCDFT.Ptheta<PthetaMax)]

			contInb = 0
			contOutb = 0
			if len(epgExpInbCDFT_selected)*len(pi0SimInbCDFT_selected) > 0:
				contInb = len(bkgSimInbCDFT_selected)/len(pi0SimInbCDFT_selected)*len(pi0ExpInbCDFT_selected)/len(epgExpInbCDFT_selected)
			if len(epgExpInbCDFT_selected)*len(pi0SimOutbCDFT_selected) > 0:
				contOutb = len(bkgSimOutbCDFT_selected)/len(pi0SimOutbCDFT_selected)*len(pi0ExpOutbCDFT_selected)/len(epgExpOutbCDFT_selected)

			# scores = []
			# correction_all = np.linspace(-0.8, 0, 17)
			# for correction_temp in correction_all:
			# 	#performing correcting
			# 	self.CorrectionV1(epgExpInbCDFT, correction_temp)
			# 	self.saveDVCSvars()
			# 	self.makeDVCS()
			# 	epgExpInbCDFT_corrected = self.df_epg

			# 	self.CorrectionV1(epgExpOutbCDFT, correction_temp)
			# 	self.saveDVCSvars()
			# 	self.makeDVCS()
			# 	epgExpOutbCDFT_corrected = self.df_epg

			# 	# epgExpInbCDFT_corrected = epgExpInbCDFT_corrected.loc[(epgExpInbCDFT_corrected.Pp>PpMin) & (epgExpInbCDFT_corrected.Pp<PpMax)]
			# 	# epgExpOutbCDFT_corrected = epgExpOutbCDFT_corrected.loc[(epgExpOutbCDFT_corrected.Pp>PpMin) & (epgExpOutbCDFT_corrected.Pp<PpMax)]

			# 	epgExpInbCDFT_corrected = epgExpInbCDFT_corrected.loc[(epgExpInbCDFT_corrected.Ptheta>PthetaMin) & (epgExpInbCDFT_corrected.Ptheta<PthetaMax)]
			# 	epgExpOutbCDFT_corrected = epgExpOutbCDFT_corrected.loc[(epgExpOutbCDFT_corrected.Ptheta>PthetaMin) & (epgExpOutbCDFT_corrected.Ptheta<PthetaMax)]

			# 	score1 = epgExpInbCDFT_corrected.MM2_ep.mean() - (1-contInb)*dvcsSimInbCDFT_selected.MM2_ep.mean() - contInb*bkgSimInbCDFT_selected.MM2_ep.mean()
			# 	score2 = epgExpOutbCDFT_corrected.MM2_ep.mean() - (1-contOutb)*dvcsSimOutbCDFT_selected.MM2_ep.mean() - contOutb*bkgSimOutbCDFT_selected.MM2_ep.mean()
			# 	score = score1**2 +score2**2
			# 	scores.append(score)
			# print(scores)
			# correction_opt = correction_all[np.argmin(scores)]
			# corrections.append(correction_opt)
			# print(correction_opt)

			# self.CorrectionV1(epgExpInbCDFT, 0.8/(1+np.exp(58.657*(epgExpInbCDFT.Pp-0.41))))
			# self.saveDVCSvars()
			# self.makeDVCS()
			# epgExpInbCDFT_corrected = self.df_epg

			# self.CorrectionV1(epgExpOutbCDFT, 0.8/(1+np.exp(58.657*(epgExpOutbCDFT.Pp-0.41))))
			# self.saveDVCSvars()
			# self.makeDVCS()
			# epgExpOutbCDFT_corrected = self.df_epg

			# self.CorrectionV1(epgExpInbCDFT, correction_opt)
			# self.saveDVCSvars()
			# self.makeDVCS()
			# epgExpInbCDFT_corrected = self.df_epg

			# self.CorrectionV1(epgExpOutbCDFT, correction_opt)
			# self.saveDVCSvars()
			# self.makeDVCS()
			# epgExpOutbCDFT_corrected = self.df_epg

			# epgExpInbCDFT_selected = epgExpInbCDFT_corrected.loc[(epgExpInbCDFT_corrected.Pp>PpMin) & (epgExpInbCDFT_corrected.Pp<PpMax)]
			# epgExpOutbCDFT_selected = epgExpOutbCDFT_corrected.loc[(epgExpOutbCDFT_corrected.Pp>PpMin) & (epgExpOutbCDFT_corrected.Pp<PpMax)]

			# epgExpInbCDFT_selected = epgExpInbCDFT_corrected.loc[(epgExpInbCDFT_corrected.Ptheta>PthetaMin) & (epgExpInbCDFT_corrected.Ptheta<PthetaMax)]
			# epgExpOutbCDFT_selected = epgExpOutbCDFT_corrected.loc[(epgExpOutbCDFT_corrected.Ptheta>PthetaMin) & (epgExpOutbCDFT_corrected.Ptheta<PthetaMax)]

			for sigma1 in sigma1s:
				for sigma2 in sigma2s:
						# for sigma3 in sigma3s:
							# print("smearing with {:.3f}, {:.3f}, {:.3f}".format(sigma1, sigma2, sigma3))
					print("smearing with {:.3f}, {:.3f}".format(sigma1, sigma2))
					# print("smearing with {:.3f}".format(sigma1))
					# print("smearing with {:.3f}".format(sigma2))

					#performing smearing
					self.SmearingV1(dvcsSimInbCDFT, sigma1, sigma2, 0)
					self.saveDVCSvars()
					self.makeDVCS(pol = "inbending")
					dvcsSimInbCDFT_smeared = self.df_epg
					self.SmearingV1(bkgSimInbCDFT, sigma1, sigma2, 0)
					self.saveDVCSvars()
					self.makeDVCS(pol = "inbending")
					bkgSimInbCDFT_smeared = self.df_epg

					self.SmearingV1(pi0SimInbCDFT, sigma1, sigma2, 0)
					self.saveDVpi0Pvars()
					self.makeDVpi0P(pol = "inbending")
					pi0SimInbCDFT_smeared = self.df_epg

					self.SmearingV1(dvcsSimOutbCDFT, sigma1, sigma2, 0)
					self.saveDVCSvars()
					self.makeDVCS(pol = "outbending")
					dvcsSimOutbCDFT_smeared = self.df_epg

					self.SmearingV1(bkgSimOutbCDFT, sigma1, sigma2, 0)
					self.saveDVCSvars()
					self.makeDVCS(pol = "outbending")
					bkgSimOutbCDFT_smeared = self.df_epg

					self.SmearingV1(pi0SimOutbCDFT, sigma1, sigma2, 0)
					self.saveDVpi0Pvars()
					self.makeDVpi0P(pol = "outbending")
					pi0SimOutbCDFT_smeared = self.df_epg

					print(len(dvcsSimInbCDFT_smeared), len(bkgSimInbCDFT_smeared), len(pi0SimInbCDFT_smeared), len(epgExpInbCDFT_selected))
					print(len(dvcsSimOutbCDFT_smeared), len(bkgSimOutbCDFT_smeared), len(pi0SimOutbCDFT_smeared), len(epgExpOutbCDFT_selected))
					# dvcsSimInbCDFT_smeared = dvcsSimInbCDFT_smeared.loc[(dvcsSimInbCDFT_smeared.Ptheta>PthetaMin) & (dvcsSimInbCDFT_smeared.Ptheta<PthetaMax)]
					# dvcsSimOutbCDFT_smeared = dvcsSimOutbCDFT_smeared.loc[(dvcsSimOutbCDFT_smeared.Ptheta>PthetaMin) & (dvcsSimOutbCDFT_smeared.Ptheta<PthetaMax)]

					dvcsSimInbCDFT_smeared = dvcsSimInbCDFT_smeared.loc[(dvcsSimInbCDFT_smeared.Pp>PpMin) & (dvcsSimInbCDFT_smeared.Pp<PpMax)]
					dvcsSimOutbCDFT_smeared = dvcsSimOutbCDFT_smeared.loc[(dvcsSimOutbCDFT_smeared.Pp>PpMin) & (dvcsSimOutbCDFT_smeared.Pp<PpMax)]

					bkgSimInbCDFT_smeared = bkgSimInbCDFT_smeared.loc[(bkgSimInbCDFT_smeared.Pp>PpMin) & (bkgSimInbCDFT_smeared.Pp<PpMax)]
					bkgSimOutbCDFT_smeared = bkgSimOutbCDFT_smeared.loc[(bkgSimOutbCDFT_smeared.Pp>PpMin) & (bkgSimOutbCDFT_smeared.Pp<PpMax)]

					pi0SimInbCDFT_smeared = pi0SimInbCDFT_smeared.loc[(pi0SimInbCDFT_smeared.Pp>PpMin) & (pi0SimInbCDFT_smeared.Pp<PpMax)]
					pi0SimOutbCDFT_smeared = pi0SimOutbCDFT_smeared.loc[(pi0SimOutbCDFT_smeared.Pp>PpMin) & (pi0SimOutbCDFT_smeared.Pp<PpMax)]
					
					contInb = 0
					contOutb = 0
					if len(epgExpInbCDFT_selected)*len(pi0SimInbCDFT_smeared) > 0:
						contInb = len(bkgSimInbCDFT_smeared)/len(pi0SimInbCDFT_smeared)*len(pi0ExpInbCDFT_selected)/len(epgExpInbCDFT_selected)
					if len(epgExpInbCDFT_selected)*len(pi0SimOutbCDFT_smeared) > 0:
						contOutb = len(bkgSimOutbCDFT_smeared)/len(pi0SimOutbCDFT_smeared)*len(pi0ExpOutbCDFT_selected)/len(epgExpOutbCDFT_selected)

					distance1 = distance(dvcsSimInbCDFT_smeared, bkgSimInbCDFT_smeared, epgExpInbCDFT_selected, cont = contInb, var = "MM2_ep")
					distance2 = distance(dvcsSimOutbCDFT_smeared, bkgSimOutbCDFT_smeared, epgExpOutbCDFT_selected, cont = contOutb, var = "MM2_ep")
					distances1.append(( distance1+ distance2 )/2) 

					distance1 = distance(dvcsSimInbCDFT_smeared, bkgSimInbCDFT_smeared, epgExpInbCDFT_selected, cont = contInb, var = "reconGam")
					distance2 = distance(dvcsSimOutbCDFT_smeared, bkgSimOutbCDFT_smeared, epgExpOutbCDFT_selected, cont = contOutb, var = "reconGam")
					distances2.append(( distance1+ distance2 )/2) 

					# distance1 = distance(dvcsSimInbCDFT_smeared, bkgSimInbCDFT_smeared, epgExpInbCDFT_selected, cont = contInb, var = "coplanarity")
					# distance2 = distance(dvcsSimOutbCDFT_smeared, bkgSimOutbCDFT_smeared, epgExpOutbCDFT_selected, cont = contOutb, var = "coplanarity")
					# distances3.append(( distance1+ distance2 )/2) 

					sigma1s_temp.append(sigma1) 
					sigma2s_temp.append(sigma2) 
					# sigma3s_temp.append(sigma3) 

			distances = ((np.array(distances1) - np.min(distances1))/np.std(distances1))**2 + ((np.array(distances2) - np.min(distances2))/np.std(distances2))**2
			sigma1_opt = sigma1s_temp[np.argmin(distance1)]
			sigma2_opt = sigma2s_temp[np.argmin(distances)]
			# sigma3_opt = sigma3s_temp[np.argmin(distances3)]
			sigmas_opt.append([sigma1_opt, sigma2_opt])#, sigma3_opt])
			# sigmas_opt.append(sigma2_opt)#, sigma3_opt])
			# sigmas_opt.append(sigma1_opt)#, sigma3_opt])
			self.SmearingV1(dvcsSimInbCDFT, sigma1_opt, sigma2_opt, 0)#, sigma3_opt)
			self.saveDVCSvars()
			self.makeDVCS(pol = "inbending")
			dvcsSimInbCDFT_opt = self.df_epg

			self.SmearingV1(bkgSimInbCDFT, sigma1_opt, sigma2_opt, 0)#, sigma3_opt)
			self.saveDVCSvars()
			self.makeDVCS(pol = "inbending")
			bkgSimInbCDFT_opt = self.df_epg

			self.SmearingV1(pi0SimInbCDFT, sigma1_opt, sigma2_opt, 0)#, sigma3_opt)
			self.saveDVpi0Pvars()
			self.makeDVpi0P(pol = "inbending")
			pi0SimInbCDFT_opt = self.df_epg

			self.SmearingV1(dvcsSimOutbCDFT, sigma1_opt, sigma2_opt, 0)#sigma3_opt)
			self.saveDVCSvars()
			self.makeDVCS(pol = "outbending")
			dvcsSimOutbCDFT_opt = self.df_epg

			self.SmearingV1(bkgSimOutbCDFT, sigma1_opt, sigma2_opt, 0)#sigma3_opt)
			self.saveDVCSvars()
			self.makeDVCS(pol = "outbending")
			bkgSimOutbCDFT_opt = self.df_epg

			self.SmearingV1(pi0SimOutbCDFT, sigma1_opt, sigma2_opt, 0)#, sigma3_opt)
			self.saveDVpi0Pvars()
			self.makeDVpi0P(pol = "outbending")
			pi0SimOutbCDFT_opt = self.df_epg

			dvcsSimInbCDFT_opt = dvcsSimInbCDFT_opt.loc[(dvcsSimInbCDFT_opt.Pp>PpMin) & (dvcsSimInbCDFT_opt.Pp<PpMax)]
			bkgSimInbCDFT_opt = bkgSimInbCDFT_opt.loc[(bkgSimInbCDFT_opt.Pp>PpMin) & (bkgSimInbCDFT_opt.Pp<PpMax)]
			pi0SimInbCDFT_opt = pi0SimInbCDFT_opt.loc[(pi0SimInbCDFT_opt.Pp>PpMin) & (pi0SimInbCDFT_opt.Pp<PpMax)]
			dvcsSimOutbCDFT_opt = dvcsSimOutbCDFT_opt.loc[(dvcsSimOutbCDFT_opt.Pp>PpMin) & (dvcsSimOutbCDFT_opt.Pp<PpMax)]
			bkgSimOutbCDFT_opt = bkgSimOutbCDFT_opt.loc[(bkgSimOutbCDFT_opt.Pp>PpMin) & (bkgSimOutbCDFT_opt.Pp<PpMax)]
			pi0SimOutbCDFT_opt = pi0SimOutbCDFT_opt.loc[(pi0SimOutbCDFT_opt.Pp>PpMin) & (pi0SimOutbCDFT_opt.Pp<PpMax)]

			contInb = 0
			contOutb = 0
			if len(epgExpInbCDFT_selected)*len(pi0SimInbCDFT_opt) > 0:
				contInb = len(bkgSimInbCDFT_opt)/len(pi0SimInbCDFT_opt)*len(pi0ExpInbCDFT_selected)/len(epgExpInbCDFT_selected)
			if len(epgExpInbCDFT_selected)*len(pi0SimOutbCDFT_opt) > 0:
				contOutb = len(bkgSimOutbCDFT_opt)/len(pi0SimOutbCDFT_opt)*len(pi0ExpOutbCDFT_selected)/len(epgExpOutbCDFT_selected)

			print(len(dvcsSimInbCDFT_opt), len(dvcsSimOutbCDFT_opt))
			print(len(epgExpInbCDFT_selected), len(epgExpOutbCDFT_selected))
			# print(np.array([sigma1s_temp, sigma2s_temp, sigma3s_temp, distances1, distances2, distances3]).T)
			# print(np.array([sigma1s_temp, sigma2s_temp, distances1, distances2, distances]).T)
			# print(np.array([sigma2s_temp, distances2]).T)
			# print(np.array([sigma1s_temp, distances1]).T)
			# print(PpMin, PpMax, sigma1_opt, sigma2_opt, correction_opt) #, sigma3_opt)#, correction_opt) 
			print(PpMin, PpMax, sigma1_opt, sigma2_opt)#, correction_opt) #, sigma3_opt)#, correction_opt) 
			# print(PthetaMin, PthetaMax, sigma2_opt)#, sigma3_opt)#, correction_opt) 
			# print(PthetaMin, PthetaMax, sigma1_opt)#, sigma3_opt)#, correction_opt) 
			# print(PpMin, PpMax, sigma2_opt)#, sigma3_opt)#, correction_opt) 
##

			varstoplot = ["Pp", "Ptheta", "Pphi","Gp", "Gtheta", "Gphi",  "coneAngle", "MM2_eg", "", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
			title = [r"$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", r"$p_{\gamma}$", r"$\theta_{\gamma}$", r"$\phi_{\gamma}$", r"$\theta_{e'\gamma}$", r"$MM^2_{e'\gamma}$", "", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
			unit = [GeV, degree, degree, GeV, degree, degree, degree, GeV2, "", GeV2, degree, degree, degree, degree, GeV, GeV2, GeV2, GeV2, GeVc]
			# binstarts = [GeMin, 0, -180, 0, 0, 0, -0.5, -0.01, 0.1, 0]
			# binends = [GeMax, 7, 180, 30, 0.8, .1, 1.2, 0.01, 1.7, 10]

			fig, axs = plt.subplots(5, 3, figsize = (15,25))
			for yind in range(0, 5):
			    for xind in range(0, 3):
			        ind = 3*yind + xind
			        if varstoplot[ind]:
			            pass
			        else:
			            continue
			        simDist_dvcs, bins = np.histogram(dvcsSimInbCDFT_opt[varstoplot[ind]], 100, density = True)
			        simDist_dvpi0, _ = np.histogram(bkgSimInbCDFT_opt[varstoplot[ind]], bins, density = True)
			        simDist = (1-contInb)*simDist_dvcs + contInb*simDist_dvpi0
			        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
			        axs[yind, xind].hist(epgExpInbCDFT_selected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
			        axs[yind, xind].set_title(title[ind])
			        # axs[yind, xind].set_xlim([start, end])
			        if (unit[ind]):
			            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
			        else:
			            axs[yind, xind].set_xlabel(title[ind])
			plt.tight_layout()
			# plt.savefig(outDir+"InbCDFT{}_{:.3f}.pdf".format(i, sigma2_opt))
			plt.savefig(outDir+"InbCDFT{}_{:.3f}_{:.3f}.pdf".format(i, sigma1_opt, sigma2_opt))
			plt.clf()

			fig, axs = plt.subplots(5, 3, figsize = (15,25))
			for yind in range(0, 5):
			    for xind in range(0, 3):
			        ind = 3*yind + xind
			        # start = binstarts[ind]
			        # end = binends[ind]
			        # bins = np.linspace(start, end, 101)
			        if varstoplot[ind]:
			            pass
			        else:
			            continue
			        simDist_dvcs, bins = np.histogram(dvcsSimOutbCDFT_opt[varstoplot[ind]], 100, density = True)
			        simDist_dvpi0, _ = np.histogram(bkgSimOutbCDFT_opt[varstoplot[ind]], bins, density = True)
			        simDist = (1-contOutb)*simDist_dvcs + contOutb*simDist_dvpi0
			        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
			        axs[yind, xind].hist(epgExpOutbCDFT_selected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
			        axs[yind, xind].set_title(title[ind])
			        # axs[yind, xind].set_xlim([start, end])
			        if (unit[ind]):
			            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
			        else:
			            axs[yind, xind].set_xlabel(title[ind])
			plt.tight_layout()
			# plt.savefig(outDir+"OutbCDFT{}_{:.3f}.pdf".format(i, sigma2_opt))
			plt.savefig(outDir+"OutbCDFT{}_{:.3f}_{:.3f}.pdf".format(i, sigma1_opt, sigma2_opt))
			plt.clf()

		print(sigmas_opt)#, corrections)


	def SmearingV1(self, df, sigma1, sigma2, sigma3):
		df_epg = copy(df)
		regulator = np.abs(2*(1/(1+np.exp(-(df_epg.loc[df_epg["Psector"]>7, "Pp"]-0.3)/0.01))-0.5))
		df_epg.loc[df_epg["Psector"]>7, "Pp"] = df_epg.loc[df_epg["Psector"]>7, "Pp"]*np.random.normal(1, regulator*sigma1, len(df_epg.loc[df_epg.Psector>7]))
		df_epg.loc[df_epg["Psector"]>7, "Ptheta"] = df_epg.loc[df_epg["Psector"]>7, "Ptheta"] + np.random.normal(0, sigma2, len(df_epg.loc[df_epg.Psector>7]))
		df_epg.loc[df_epg["Psector"]>7, "Pphi"] = df_epg.loc[df_epg["Psector"]>7, "Pphi"] + np.random.normal(0, sigma3, len(df_epg.loc[df_epg.Psector>7])) 
		df_epg.loc[:, 'Pe'] = np.sqrt( df_epg.Pp**2 + M**2)
		df_epg.loc[:, "Ppx"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.cos(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppy"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.sin(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppz"] = df_epg.loc[:, "Pp"]*np.cos(np.radians(df_epg.loc[:, "Ptheta"]))

		self.df_epg = df_epg

	def CorrectionV1(self, df, correction):
		df_epg = copy(df)
		# df_epg.loc[df_epg["Psector"]>7, "Pp"] = df_epg.loc[df_epg["Psector"]>7, "Pp"] + correction
		df_epg.loc[df_epg["Psector"]>7, "Ptheta"] = df_epg.loc[df_epg["Psector"]>7, "Ptheta"] + correction
		# df_epg.loc[df_epg["Psector"]>7, "Pphi"] = df_epg.loc[df_epg["Psector"]>7, "Pphi"] + correction
		df_epg.loc[:, 'Pe'] = np.sqrt( df_epg.Pp**2 + M**2)
		df_epg.loc[:, "Ppx"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.cos(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppy"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.sin(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppz"] = df_epg.loc[:, "Pp"]*np.cos(np.radians(df_epg.loc[:, "Ptheta"]))

		self.df_epg = df_epg

	def MakeV3(self, inDir = "SimtoDat/v2/", outDir = "SimtoDat/v3/"):

		# binsMEepgInb = np.linspace(-0.439, 0.484, 101)
		# binsMM2egInb = np.linspace(0.246, 1.569, 101)
		# binsMEepgOutb = np.linspace(-0.436, 0.481, 101) 
		# binsMM2egOutb = np.linspace(0.223, 1.602, 101) 

		def distance(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1, _ = np.histogram(df1.loc[:, var], bins = bins)
			unchist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			unchist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1 = np.sqrt(unchist1)
			unchist2 = np.sqrt(unchist2)
			unchist_exp = np.sqrt(unchist_exp)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			uncdist1 = unchist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			uncdist2 = unchist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist_exp = unchist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist = np.sqrt((1-cont)**2 * uncdist1**2 + cont**2 * uncdist2 **2 + unchist_exp**2)
			uncdist = np.where(uncdist>0, uncdist, np.inf)
			distance = np.sum(((1-cont)*dist1 + cont*dist2 -dist_exp)**2)
			return distance


		def corr(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			simDist = (1-cont)*dist1 + cont*dist2
			expDist = dist_exp
			simPeak = bincenters[np.argmax(simDist)]
			expPeak = bincenters[np.argmax(expDist)]
			return expPeak - simPeak

		# if isinstance(sigma, str):
		# 	sigma = float(sigma)

		GeEdges = [2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5]
		# GthetaEdges = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35]
		sigmas = [0, 0.01, 0.02, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06]
		corrections = []
		sigmas_opt = []

		epgExpInbCD = pd.read_pickle(inDir + "/epgExpInbCD")
		epgExpOutbCD = pd.read_pickle(inDir + "/epgExpOutbCD")

		pi0ExpInbCD = pd.read_pickle(inDir+"/pi0ExpInbCD")
		pi0ExpOutbCD = pd.read_pickle(inDir+"/pi0ExpOutbCD")

		dvcsSimInbCD = pd.read_pickle(inDir+"/dvcsSimInbCD")
		dvcsSimOutbCD = pd.read_pickle(inDir+"/dvcsSimOutbCD")

		bkgSimInbCD = pd.read_pickle(inDir+"/bkgSimInbCD")
		bkgSimOutbCD = pd.read_pickle(inDir+"/bkgSimOutbCD")

		pi0SimInbCD = pd.read_pickle(inDir+"/pi0SimInbCD")
		pi0SimOutbCD = pd.read_pickle(inDir+"/pi0SimOutbCD")


		# #performing correcting
		# correction = 0.186/(1+np.exp(1.254*(epgExpInbCD.Gp-4.5)))
		# self.CorrectionV2(epgExpInbCD, correction, mode = "epg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "inbending")
		# epgExpInbCD = self.df_epg
		# correction = 0.186/(1+np.exp(1.254*(pi0ExpInbCD.Gp-4.5)))
		# self.CorrectionV2(pi0ExpInbCD, correction, mode = "epgg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "inbending")
		# pi0ExpInbCD = self.df_epg
		# #performing correcting
		# correction = 0.186/(1+np.exp(1.254*(epgExpOutbCD.Gp-4.5)))
		# self.CorrectionV2(epgExpOutbCD, correction, mode = "epg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "outbending")
		# epgExpOutbCD = self.df_epg
		# correction = 0.186/(1+np.exp(1.254*(pi0ExpOutbCD.Gp-4.5)))
		# self.CorrectionV2(pi0ExpOutbCD, correction, mode = "epgg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "outbending")
		# pi0ExpOutbCD = self.df_epg

		# epgExpInbCD.to_pickle(outDir + "/epgExpInbCD")
		# epgExpOutbCD.to_pickle(outDir + "/epgExpOutbCD")

		# pi0ExpInbCD.to_pickle(outDir+"/pi0ExpInbCD")
		# pi0ExpOutbCD.to_pickle(outDir+"/pi0ExpOutbCD")


		for sector in range(1, 7):
			correction_s = []
			sigmas_opt_s = []
			for i in range(len(GeEdges)-1):
			# for i in range(len(GthetaEdges)-1):

				distances = []
				GeMin = GeEdges[i]
				GeMax = GeEdges[i+1]

				epgExpInbCD_selected = epgExpInbCD.loc[(epgExpInbCD.Ge>GeMin) & (epgExpInbCD.Ge<GeMax) & (epgExpInbCD.Gsector==sector)]
				pi0ExpInbCD_selected = pi0ExpInbCD.loc[(pi0ExpInbCD.Ge>GeMin) & (pi0ExpInbCD.Ge<GeMax) & (pi0ExpInbCD.Gsector==sector)]
				dvcsSimInbCD_selected = dvcsSimInbCD.loc[(dvcsSimInbCD.Ge>GeMin) & (dvcsSimInbCD.Ge<GeMax) & (dvcsSimInbCD.Gsector==sector)]
				pi0SimInbCD_selected = pi0SimInbCD.loc[(pi0SimInbCD.Ge>GeMin) & (pi0SimInbCD.Ge<GeMax) & (pi0SimInbCD.Gsector==sector)]
				bkgSimInbCD_selected = bkgSimInbCD.loc[(bkgSimInbCD.Ge>GeMin) & (bkgSimInbCD.Ge<GeMax) & (bkgSimInbCD.Gsector==sector)]

				epgExpOutbCD_selected = epgExpOutbCD.loc[(epgExpOutbCD.Ge>GeMin) & (epgExpOutbCD.Ge<GeMax) & (epgExpOutbCD.Gsector==sector)]
				pi0ExpOutbCD_selected = pi0ExpOutbCD.loc[(pi0ExpOutbCD.Ge>GeMin) & (pi0ExpOutbCD.Ge<GeMax) & (pi0ExpOutbCD.Gsector==sector)]
				dvcsSimOutbCD_selected = dvcsSimOutbCD.loc[(dvcsSimOutbCD.Ge>GeMin) & (dvcsSimOutbCD.Ge<GeMax) & (dvcsSimOutbCD.Gsector==sector)]
				pi0SimOutbCD_selected = pi0SimOutbCD.loc[(pi0SimOutbCD.Ge>GeMin) & (pi0SimOutbCD.Ge<GeMax) & (pi0SimOutbCD.Gsector==sector)]
				bkgSimOutbCD_selected = bkgSimOutbCD.loc[(bkgSimOutbCD.Ge>GeMin) & (bkgSimOutbCD.Ge<GeMax) & (bkgSimOutbCD.Gsector==sector)]

				# GthetaMin = GthetaEdges[i]
				# GthetaMax = GthetaEdges[i+1]

				# epgExpInbCD_selected = epgExpInbCD.loc[(epgExpInbCD.Gtheta>GthetaMin) & (epgExpInbCD.Gtheta<GthetaMax) & (epgExpInbCD.Gsector == sector)]
				# pi0ExpInbCD_selected = pi0ExpInbCD.loc[(pi0ExpInbCD.Gtheta>GthetaMin) & (pi0ExpInbCD.Gtheta<GthetaMax) & (pi0ExpInbCD.Gsector == sector)]
				# dvcsSimInbCD_selected = dvcsSimInbCD.loc[(dvcsSimInbCD.Gtheta>GthetaMin) & (dvcsSimInbCD.Gtheta<GthetaMax) & (dvcsSimInbCD.Gsector == sector)]
				# pi0SimInbCD_selected = pi0SimInbCD.loc[(pi0SimInbCD.Gtheta>GthetaMin) & (pi0SimInbCD.Gtheta<GthetaMax) & (pi0SimInbCD.Gsector == sector)]
				# bkgSimInbCD_selected = bkgSimInbCD.loc[(bkgSimInbCD.Gtheta>GthetaMin) & (bkgSimInbCD.Gtheta<GthetaMax) & (bkgSimInbCD.Gsector == sector)]

				# epgExpOutbCD_selected = epgExpOutbCD.loc[(epgExpOutbCD.Gtheta>GthetaMin) & (epgExpOutbCD.Gtheta<GthetaMax) & (epgExpOutbCD.Gsector == sector)]
				# pi0ExpOutbCD_selected = pi0ExpOutbCD.loc[(pi0ExpOutbCD.Gtheta>GthetaMin) & (pi0ExpOutbCD.Gtheta<GthetaMax) & (pi0ExpOutbCD.Gsector == sector)]
				# dvcsSimOutbCD_selected = dvcsSimOutbCD.loc[(dvcsSimOutbCD.Gtheta>GthetaMin) & (dvcsSimOutbCD.Gtheta<GthetaMax) & (dvcsSimOutbCD.Gsector == sector)]
				# pi0SimOutbCD_selected = pi0SimOutbCD.loc[(pi0SimOutbCD.Gtheta>GthetaMin) & (pi0SimOutbCD.Gtheta<GthetaMax) & (pi0SimOutbCD.Gsector == sector)]
				# bkgSimOutbCD_selected = bkgSimOutbCD.loc[(bkgSimOutbCD.Gtheta>GthetaMin) & (bkgSimOutbCD.Gtheta<GthetaMax) & (bkgSimOutbCD.Gsector == sector)]

				contInb = 0
				contOutb = 0
				if len(epgExpInbCD_selected)*len(pi0SimInbCD_selected) > 0:
					contInb = len(bkgSimInbCD_selected)/len(pi0SimInbCD_selected)*len(pi0ExpInbCD_selected)/len(epgExpInbCD_selected)
				if len(epgExpOutbCD_selected)*len(pi0SimOutbCD_selected) > 0:
					contOutb = len(bkgSimOutbCD_selected)/len(pi0SimOutbCD_selected)*len(pi0ExpOutbCD_selected)/len(epgExpOutbCD_selected)

				correction1 = epgExpInbCD_selected.ME_epg.mean() - (1-contInb)*dvcsSimInbCD_selected.ME_epg.mean() - contInb*bkgSimInbCD_selected.ME_epg.mean()
				correction2 = epgExpOutbCD_selected.ME_epg.mean() - (1-contOutb)*dvcsSimOutbCD_selected.ME_epg.mean() - contOutb*bkgSimOutbCD_selected.ME_epg.mean()
				# correction1 = corr(dvcsSimInbCD_selected, bkgSimInbCD_selected, epgExpInbCD_selected, cont = contInb, var = "ME_epg")
				# correction2 = corr(dvcsSimOutbCD_selected, bkgSimOutbCD_selected, epgExpOutbCD_selected, cont = contOutb, var = "ME_epg")
				correction = (correction1+correction2)/2
				print(correction1, contInb, correction2, contOutb)
				corrections.append(correction)
				correction_s.append(correction)

				#performing correcting
				self.CorrectionV2(epgExpInbCD, correction, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS()
				epgExpInbCD_corrected = self.df_epg
				self.CorrectionV2(pi0ExpInbCD, correction, mode = "epgg")
				self.saveDVpi0Pvars()
				self.makeDVpi0P()
				pi0ExpInbCD_corrected = self.df_epg

				#performing correcting
				self.CorrectionV2(epgExpOutbCD, correction, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "outbending")
				epgExpOutbCD_corrected = self.df_epg
				self.CorrectionV2(pi0ExpOutbCD, correction, mode = "epgg")
				self.saveDVpi0Pvars()
				self.makeDVpi0P(pol = "outbending")
				pi0ExpOutbCD_corrected = self.df_epg

				epgExpInbCD_corrected = epgExpInbCD_corrected.loc[(epgExpInbCD_corrected.Ge>GeMin) & (epgExpInbCD_corrected.Ge<GeMax) & (epgExpInbCD_corrected.Gsector == sector)]
				epgExpOutbCD_corrected = epgExpOutbCD_corrected.loc[(epgExpOutbCD_corrected.Ge>GeMin) & (epgExpOutbCD_corrected.Ge<GeMax) & (epgExpOutbCD_corrected.Gsector == sector)]

				pi0ExpInbCD_corrected = pi0ExpInbCD_corrected.loc[(pi0ExpInbCD_corrected.Ge>GeMin) & (pi0ExpInbCD_corrected.Ge<GeMax) & (pi0ExpInbCD_corrected.Gsector == sector)]
				pi0ExpOutbCD_corrected = pi0ExpOutbCD_corrected.loc[(pi0ExpOutbCD_corrected.Ge>GeMin) & (pi0ExpOutbCD_corrected.Ge<GeMax) & (pi0ExpOutbCD_corrected.Gsector == sector)]

				# epgExpInbCD_corrected = epgExpInbCD_corrected.loc[(epgExpInbCD_corrected.Gtheta>GthetaMin) & (epgExpInbCD_corrected.Gtheta<GthetaMax) & (epgExpInbCD_corrected.Gsector == sector)]
				# epgExpOutbCD_corrected = epgExpOutbCD_corrected.loc[(epgExpOutbCD_corrected.Gtheta>GthetaMin) & (epgExpOutbCD_corrected.Gtheta<GthetaMax) & (epgExpOutbCD_corrected.Gsector == sector)]

				# pi0ExpInbCD_corrected = pi0ExpInbCD_corrected.loc[(pi0ExpInbCD_corrected.Gtheta>GthetaMin) & (pi0ExpInbCD_corrected.Gtheta<GthetaMax) & (pi0ExpInbCD_corrected.Gsector == sector)]
				# pi0ExpOutbCD_corrected = pi0ExpOutbCD_corrected.loc[(pi0ExpOutbCD_corrected.Gtheta>GthetaMin) & (pi0ExpOutbCD_corrected.Gtheta<GthetaMax) & (pi0ExpOutbCD_corrected.Gsector == sector)]


				# epgExpInbCD_corrected = epgExpInbCD.loc[(epgExpInbCD.Ge>GeMin) & (epgExpInbCD.Ge<GeMax)]
				# epgExpOutbCD_corrected = epgExpOutbCD.loc[(epgExpOutbCD.Ge>GeMin) & (epgExpOutbCD.Ge<GeMax)]


				for sigma in sigmas:

					print("smearing with {:.3f}".format(sigma))

					#performing smearing
					self.SmearingV2(dvcsSimInbCD, sigma, mode = "epg")
					self.saveDVCSvars()
					self.makeDVCS(pol = "inbending")
					dvcsSimInbCD_smeared = self.df_epg

					self.SmearingV2(bkgSimInbCD, sigma, mode = "epg")
					self.saveDVCSvars()
					self.makeDVCS(pol = "inbending")
					bkgSimInbCD_smeared = self.df_epg

					self.SmearingV2(dvcsSimOutbCD, sigma, mode = "epg")
					self.saveDVCSvars()
					self.makeDVCS(pol = "outbending")
					dvcsSimOutbCD_smeared = self.df_epg

					self.SmearingV2(bkgSimOutbCD, sigma, mode = "epg")
					self.saveDVCSvars()
					self.makeDVCS(pol = "outbending")
					bkgSimOutbCD_smeared = self.df_epg

					dvcsSimInbCD_smeared = dvcsSimInbCD_smeared.loc[(dvcsSimInbCD_smeared.Ge>GeMin) & (dvcsSimInbCD_smeared.Ge<GeMax) & ( dvcsSimInbCD_smeared.Gsector == sector)]
					dvcsSimOutbCD_smeared = dvcsSimOutbCD_smeared.loc[(dvcsSimOutbCD_smeared.Ge>GeMin) & (dvcsSimOutbCD_smeared.Ge<GeMax) & ( dvcsSimOutbCD_smeared.Gsector == sector)]

					bkgSimInbCD_smeared = bkgSimInbCD_smeared.loc[(bkgSimInbCD_smeared.Ge>GeMin) & (bkgSimInbCD_smeared.Ge<GeMax) & ( bkgSimInbCD_smeared.Gsector == sector)]
					bkgSimOutbCD_smeared = bkgSimOutbCD_smeared.loc[(bkgSimOutbCD_smeared.Ge>GeMin) & (bkgSimOutbCD_smeared.Ge<GeMax) & ( bkgSimOutbCD_smeared.Gsector == sector)]

					# dvcsSimInbCD_smeared = dvcsSimInbCD_smeared.loc[(dvcsSimInbCD_smeared.Gtheta>GthetaMin) & (dvcsSimInbCD_smeared.Gtheta<GthetaMax) & ( dvcsSimInbCD_smeared.Gsector == sector)]
					# dvcsSimOutbCD_smeared = dvcsSimOutbCD_smeared.loc[(dvcsSimOutbCD_smeared.Gtheta>GthetaMin) & (dvcsSimOutbCD_smeared.Gtheta<GthetaMax) & ( dvcsSimOutbCD_smeared.Gsector == sector)]

					# bkgSimInbCD_smeared = bkgSimInbCD_smeared.loc[(bkgSimInbCD_smeared.Gtheta>GthetaMin) & (bkgSimInbCD_smeared.Gtheta<GthetaMax) & ( bkgSimInbCD_smeared.Gsector == sector)]
					# bkgSimOutbCD_smeared = bkgSimOutbCD_smeared.loc[(bkgSimOutbCD_smeared.Gtheta>GthetaMin) & (bkgSimOutbCD_smeared.Gtheta<GthetaMax) & ( bkgSimOutbCD_smeared.Gsector == sector)]

					distance1 = distance(dvcsSimInbCD_smeared, bkgSimInbCD_smeared, epgExpInbCD_corrected, cont = contInb, var = "ME_epg")
					distance2 = distance(dvcsSimOutbCD_smeared, bkgSimOutbCD_smeared, epgExpOutbCD_corrected, cont = contOutb, var = "ME_epg")

					distances.append(( distance1+ distance2 )/2) 

				sigma_opt = sigmas[np.argmin(distances)]
				sigmas_opt.append(sigma_opt)
				sigmas_opt_s.append(sigma_opt)

				self.SmearingV2(dvcsSimInbCD, sigma_opt, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "inbending")
				dvcsSimInbCD_opt = self.df_epg

				self.SmearingV2(dvcsSimOutbCD, sigma_opt, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "outbending")
				dvcsSimOutbCD_opt = self.df_epg

				self.SmearingV2(bkgSimInbCD, sigma_opt, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "inbending")
				bkgSimInbCD_opt = self.df_epg

				self.SmearingV2(bkgSimOutbCD, sigma_opt, mode = "epg")
				self.saveDVCSvars()
				self.makeDVCS(pol = "outbending")
				bkgSimOutbCD_opt = self.df_epg

				self.SmearingV2(pi0SimInbCD, sigma_opt, mode = "epgg")
				self.saveDVpi0Pvars()
				self.makeDVpi0P(pol = "inbending")
				pi0SimInbCD_opt = self.df_epg

				self.SmearingV2(pi0SimOutbCD, sigma_opt, mode = "epgg")
				self.saveDVpi0Pvars()
				self.makeDVpi0P(pol = "outbending")
				pi0SimOutbCD_opt = self.df_epg

				dvcsSimInbCD_opt = dvcsSimInbCD_opt.loc[(dvcsSimInbCD_opt.Ge>GeMin) & (dvcsSimInbCD_opt.Ge<GeMax) & (dvcsSimInbCD_opt.Gsector == sector)]
				dvcsSimOutbCD_opt = dvcsSimOutbCD_opt.loc[(dvcsSimOutbCD_opt.Ge>GeMin) & (dvcsSimOutbCD_opt.Ge<GeMax) & (dvcsSimOutbCD_opt.Gsector == sector)]

				bkgSimInbCD_opt = bkgSimInbCD_opt.loc[(bkgSimInbCD_opt.Ge>GeMin) & (bkgSimInbCD_opt.Ge<GeMax) & (bkgSimInbCD_opt.Gsector == sector)]
				bkgSimOutbCD_opt = bkgSimOutbCD_opt.loc[(bkgSimOutbCD_opt.Ge>GeMin) & (bkgSimOutbCD_opt.Ge<GeMax) & (bkgSimOutbCD_opt.Gsector == sector)]

				pi0SimInbCD_opt = pi0SimInbCD_opt.loc[(pi0SimInbCD_opt.Ge>GeMin) & (pi0SimInbCD_opt.Ge<GeMax) & (pi0SimInbCD_opt.Gsector == sector)]
				pi0SimOutbCD_opt = pi0SimOutbCD_opt.loc[(pi0SimOutbCD_opt.Ge>GeMin) & (pi0SimOutbCD_opt.Ge<GeMax) & (pi0SimOutbCD_opt.Gsector == sector)]

				# dvcsSimInbCD_opt = dvcsSimInbCD_opt.loc[(dvcsSimInbCD_opt.Gtheta>GthetaMin) & (dvcsSimInbCD_opt.Gtheta<GthetaMax) & (dvcsSimInbCD_opt.Gsector == sector)]
				# dvcsSimOutbCD_opt = dvcsSimOutbCD_opt.loc[(dvcsSimOutbCD_opt.Gtheta>GthetaMin) & (dvcsSimOutbCD_opt.Gtheta<GthetaMax) & (dvcsSimOutbCD_opt.Gsector == sector)]

				# bkgSimInbCD_opt = bkgSimInbCD_opt.loc[(bkgSimInbCD_opt.Gtheta>GthetaMin) & (bkgSimInbCD_opt.Gtheta<GthetaMax) & (bkgSimInbCD_opt.Gsector == sector)]
				# bkgSimOutbCD_opt = bkgSimOutbCD_opt.loc[(bkgSimOutbCD_opt.Gtheta>GthetaMin) & (bkgSimOutbCD_opt.Gtheta<GthetaMax) & (bkgSimOutbCD_opt.Gsector == sector)]

				# pi0SimInbCD_opt = pi0SimInbCD_opt.loc[(pi0SimInbCD_opt.Gtheta>GthetaMin) & (pi0SimInbCD_opt.Gtheta<GthetaMax) & (pi0SimInbCD_opt.Gsector == sector)]
				# pi0SimOutbCD_opt = pi0SimOutbCD_opt.loc[(pi0SimOutbCD_opt.Gtheta>GthetaMin) & (pi0SimOutbCD_opt.Gtheta<GthetaMax) & (pi0SimOutbCD_opt.Gsector == sector)]

				contInb = 0
				contOutb = 0
				if len(epgExpInbCD_corrected)*len(pi0SimInbCD_opt) > 0:
					contInb = len(bkgSimInbCD_opt)/len(pi0SimInbCD_opt)*len(pi0ExpInbCD_corrected)/len(epgExpInbCD_corrected)
				if len(epgExpOutbCD_corrected)*len(pi0SimOutbCD_opt) > 0:
					contOutb = len(bkgSimOutbCD_opt)/len(pi0SimOutbCD_opt)*len(pi0ExpOutbCD_corrected)/len(epgExpOutbCD_corrected)

				print(contInb, contOutb)

				print(len(dvcsSimInbCD_opt), len(dvcsSimOutbCD_opt))
				print(len(epgExpInbCD_corrected), len(epgExpOutbCD_corrected))
				print(GeMin, GeMax, distances, sigma_opt, correction) 
				# print(GthetaMin, GthetaMax, distances, sigma_opt, correction) 

				varstoplot = ["Pp", "Ptheta", "Pphi","Gp", "Gtheta", "Gphi",  "coneAngle", "MM2_eg", "", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
				title = [r"$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", r"$p_{\gamma}$", r"$\theta_{\gamma}$", r"$\phi_{\gamma}$", r"$\theta_{e'\gamma}$", r"$MM^2_{e'\gamma}$", "", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
				unit = [GeV, degree, degree, GeV, degree, degree, degree, GeV2, "", GeV2, degree, degree, degree, degree, GeV, GeV2, GeV2, GeV2, GeVc]
				# binstarts = [GeMin, 0, -180, 0, 0, 0, -0.5, -0.01, 0.1, 0]
				# binends = [GeMax, 7, 180, 30, 0.8, .1, 1.2, 0.01, 1.7, 10]

				fig, axs = plt.subplots(5, 3, figsize = (15,25))
				for yind in range(0, 5):
				    for xind in range(0, 3):
				        ind = 3*yind + xind
				        if varstoplot[ind]:
				            pass
				        else:
				            continue
				        simDist_dvcs, bins = np.histogram(dvcsSimInbCD_opt[varstoplot[ind]], 100, density = True)
				        simDist_dvpi0, _ = np.histogram(bkgSimInbCD_opt[varstoplot[ind]], bins, density = True)
				        simDist = (1-contInb)*simDist_dvcs + contInb*simDist_dvpi0
				        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
				        axs[yind, xind].hist(epgExpInbCD_corrected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
				        axs[yind, xind].set_title(title[ind])
				        # axs[yind, xind].set_xlim([start, end])
				        if (unit[ind]):
				            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				        else:
				            axs[yind, xind].set_xlabel(title[ind])
				plt.tight_layout()
				plt.savefig(outDir+"InbCD{}_{:.3f}_{:.3f}_{:d}.pdf".format(i, sigma_opt, correction, sector))
				plt.clf()

				fig, axs = plt.subplots(5, 3, figsize = (15,25))
				for yind in range(0, 5):
				    for xind in range(0, 3):
				        ind = 3*yind + xind
				        # start = binstarts[ind]
				        # end = binends[ind]
				        # bins = np.linspace(start, end, 101)
				        if varstoplot[ind]:
				            pass
				        else:
				            continue
				        simDist_dvcs, bins = np.histogram(dvcsSimOutbCD_opt[varstoplot[ind]], 100, density = True)
				        simDist_dvpi0, _ = np.histogram(bkgSimOutbCD_opt[varstoplot[ind]], bins, density = True)
				        simDist = (1-contOutb)*simDist_dvcs + contOutb*simDist_dvpi0
				        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
				        axs[yind, xind].hist(epgExpOutbCD_corrected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
				        axs[yind, xind].set_title(title[ind])
				        # axs[yind, xind].set_xlim([start, end])
				        if (unit[ind]):
				            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				        else:
				            axs[yind, xind].set_xlabel(title[ind])
				plt.tight_layout()
				plt.savefig(outDir+"OutbCD{}_{:.3f}_{:.3f}_{:d}.pdf".format(i, sigma_opt, correction, sector))
				plt.clf()
			print("correction in sector{}".format(sector), correction_s)
			print("smearings in sector{}".format(sector), sigmas_opt_s)

		print(sigmas_opt, corrections)


	def CorrectionV2(self, df, correction, mode = "epg"):
		df_epg = copy(df)
		if mode == "epg":
			df_epg.loc[df_epg.Gsector<7, 'Gp'] = df_epg.Gp + correction
			df_epg.loc[df_epg.Gsector<7, 'Ge'] = df_epg.loc[df_epg.Gsector<7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))
		elif mode == "epgg":
			df_epg.loc[df_epg.Gsector<7, 'Gp'] = df_epg.Gp + correction
			df_epg.loc[df_epg.Gsector<7, 'Ge'] = df_epg.loc[df_epg.Gsector<7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))			
		self.df_epg = df_epg

	def SmearingV2(self, df, sigma, mode = "epg"):
		df_epg = copy(df)
		if mode == "epg":
			df_epg.loc[df_epg.Gsector<7, 'Gp'] = np.random.normal(1, sigma, len(df_epg.loc[df_epg.Gsector<7]))*df_epg.loc[df_epg.Gsector<7, 'Gp']
			df_epg.loc[df_epg.Gsector<7, 'Ge'] = df_epg.loc[df_epg.Gsector<7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))
		if mode == "epgg":
			df_epg.loc[df_epg.Gsector<7, 'Gp'] = np.random.normal(1, sigma, len(df_epg.loc[df_epg.Gsector<7]))*df_epg.loc[df_epg.Gsector<7, 'Gp']
			df_epg.loc[df_epg.Gsector<7, 'Ge'] = df_epg.loc[df_epg.Gsector<7, 'Gp']
			df_epg.loc[:, "Gpx"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.cos(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpy"] = df_epg.loc[:, "Gp"]*np.sin(np.radians(df_epg.loc[:, "Gtheta"]))*np.sin(np.radians(df_epg.loc[:, "Gphi"]))
			df_epg.loc[:, "Gpz"] = df_epg.loc[:, "Gp"]*np.cos(np.radians(df_epg.loc[:, "Gtheta"]))

			df_epg.loc[df_epg.Gsector<7, 'Gp2'] = np.random.normal(1, 0.035, len(df_epg.loc[df_epg.Gsector<7]))*df_epg.loc[df_epg.Gsector<7, 'Gp2']
			df_epg.loc[df_epg.Gsector<7, 'Ge2'] = df_epg.loc[df_epg.Gsector<7, 'Gp2']
			df_epg.loc[:, "Gpx2"] = df_epg.loc[:, "Gp2"]*np.sin(np.radians(df_epg.loc[:, "Gtheta2"]))*np.cos(np.radians(df_epg.loc[:, "Gphi2"]))
			df_epg.loc[:, "Gpy2"] = df_epg.loc[:, "Gp2"]*np.sin(np.radians(df_epg.loc[:, "Gtheta2"]))*np.sin(np.radians(df_epg.loc[:, "Gphi2"]))
			df_epg.loc[:, "Gpz2"] = df_epg.loc[:, "Gp2"]*np.cos(np.radians(df_epg.loc[:, "Gtheta2"]))

		self.df_epg = df_epg

	def MakeV4(self, inDir = "SimtoDat/v3/", outDir = "SimtoDat/v4/", pol = "inbending"):

		# binsMEepgInb = np.linspace(-0.439, 0.484, 101)
		# binsMM2egInb = np.linspace(0.246, 1.569, 101)
		# binsMEepgOutb = np.linspace(-0.436, 0.481, 101) 
		# binsMM2egOutb = np.linspace(0.223, 1.602, 101) 

		def distance(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1, _ = np.histogram(df1.loc[:, var], bins = bins)
			unchist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			unchist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			unchist1 = np.sqrt(unchist1)
			unchist2 = np.sqrt(unchist2)
			unchist_exp = np.sqrt(unchist_exp)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			uncdist1 = unchist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			uncdist2 = unchist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist_exp = unchist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			uncdist = np.sqrt((1-cont)**2 * uncdist1**2 + cont**2 * uncdist2 **2 + unchist_exp**2)
			uncdist = np.where(uncdist>0, uncdist, np.inf)
			distance = np.sum(((1-cont)*dist1 + cont*dist2 -dist_exp)**2)
			return distance


		def corr(df1, df2, df_exp, cont = 0, var = "ME_epg"):
			hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
			hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
			hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
			bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
			dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
			dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
			dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
			simDist = (1-cont)*dist1 + cont*dist2
			expDist = dist_exp
			simPeak = bincenters[np.argmax(simDist)]
			expPeak = bincenters[np.argmax(expDist)]
			return expPeak - simPeak

		# if isinstance(sigma, str):
		# 	sigma = float(sigma)

		PpEdges = [0.42, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
		sigmas = [0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08]
		# PthetaEdges = [40, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65]
		sigmas_temp = []
		corrections = np.linspace(-0.04, 0.04, 21)
		sigmas_opt = []
		corrections_opt = []

		polarity = "Inb"
		if pol == "outbending":
			polarity = "Outb"
			PpEdges = PpEdges[1:]

		print("FD {} proton correction/ smearing...".format(pol))
		pi0ExpFD = pd.read_pickle(inDir+"/pi0Exp{}FD".format(polarity))
		pi0SimFD = pd.read_pickle(inDir+"/pi0Sim{}FD".format(polarity))

		for sector in range(1, 7):
			correction_s = []
			sigmas_opt_s = []
			for i in range(len(PpEdges)-1):
			# for i in range(len(PthetaEdges)-1):

				distances = []
				PpMin = PpEdges[i]
				PpMax = PpEdges[i+1]
				
				scores = []
				pi0SimFD_selected = pi0SimFD.loc[(pi0SimFD.Pp>PpMin) & (pi0SimFD.Pp<PpMax) & ( pi0SimFD.Psector == sector)]
				for correction in corrections:
					print("adding momentum {}".format(correction))

					self.CorrectionV3(pi0ExpFD, correction)
					self.saveDVpi0Pvars()
					self.makeDVpi0P(pol = pol)
					pi0ExpFD_corrected = self.df_epg
					pi0ExpFD_corrected = pi0ExpFD_corrected.loc[(pi0ExpFD_corrected.Pp>PpMin) & (pi0ExpFD_corrected.Pp<PpMax) & ( pi0ExpFD_corrected.Psector == sector)]
					scores.append(np.abs(pi0SimFD_selected.MM2_ep.mean() - pi0ExpFD_corrected.MM2_ep.mean()))
				print(scores)
				correction_opt = corrections[np.argmin(scores)]
				corrections_opt.append(correction_opt)
				correction_s.append(correction_opt)
				print("the additional, optimal momentum correction is {}".format(correction_opt))
				self.CorrectionV3(pi0ExpFD, correction_opt)
				self.saveDVpi0Pvars()
				self.makeDVpi0P(pol = pol)
				pi0ExpFD_corrected = self.df_epg
				pi0ExpFD_corrected = pi0ExpFD_corrected.loc[(pi0ExpFD_corrected.Pp>PpMin) & (pi0ExpFD_corrected.Pp<PpMax) & ( pi0ExpFD_corrected.Psector == sector)]

				for sigma in sigmas:

					print("smearing with {:.3f}".format(sigma))

					#performing smearing
					self.SmearingV3(pi0SimFD, sigma, pol = pol)
					self.saveDVpi0Pvars()
					self.makeDVpi0P(pol = pol)
					pi0SimFD_smeared = self.df_epg
					pi0SimFD_smeared = pi0SimFD_smeared.loc[(pi0SimFD_smeared.Pp>PpMin) & (pi0SimFD_smeared.Pp<PpMax) & ( pi0SimFD_smeared.Psector == sector)]

					distance1 = distance(pi0SimFD_smeared, pi0SimFD_smeared, pi0ExpFD_corrected, cont = 0, var = "MM2_ep")
					distances.append(distance1) 

				sigma_opt = sigmas[np.argmin(distances)]
				sigmas_opt.append(sigma_opt)
				sigmas_opt_s.append(sigma_opt)

				self.SmearingV3(pi0SimFD, sigma_opt, pol = pol)
				self.saveDVpi0Pvars()
				self.makeDVpi0P(pol = "inbending")
				pi0SimFD_opt = self.df_epg
				pi0SimFD_opt = pi0SimFD_opt.loc[(pi0SimFD_opt.Pp>PpMin) & (pi0SimFD_opt.Pp<PpMax) & (pi0SimFD_opt.Psector == sector)]

				print(len(pi0SimFD_opt))
				print(len(pi0ExpFD_corrected))
				print(PpMin, PpMax, distances, sigma_opt, correction_opt) 

				varstoplot = ["Pp", "Ptheta", "Pphi","Gp", "Gtheta", "Gphi",  "coneAngle1", "MM2_egg", "Mpi0", "reconPi", "coplanarity", "ME_epgg", "MM2_epgg", "MM2_ep", "MPt"]
				title = [r"$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", r"$p_{\gamma 1}$", r"$\theta_{\gamma 1}$", r"$\phi_{\gamma 1}$", r"$\theta_{e'\gamma 1}$", r"$MM^2_{e'\gamma 1}$", "", r"$\theta_{\pi^0_{det.}\pi^0_{rec.}}$", r"$\Delta\phi$" , "ME"+r"${}_{ep\pi^0}$", "MM"+r"${}^{2}_{ep\pi^0}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{ep\pi^0}$"]
				unit = [GeV, degree, degree, GeV, degree, degree, degree, GeV2, GeV, GeV2, degree, degree, degree, degree, GeV, GeV2, GeV2, GeV2, GeVc]
				# binstarts = [PpMin, 0, -180, 0, 0, 0, -0.5, -0.01, 0.1, 0]
				# binends = [PpMax, 7, 180, 30, 0.8, .1, 1.2, 0.01, 1.7, 10]

				fig, axs = plt.subplots(5, 3, figsize = (15,25))
				for yind in range(0, 5):
				    for xind in range(0, 3):
				        ind = 3*yind + xind
				        if varstoplot[ind]:
				            pass
				        else:
				            continue
				        simDist, bins = np.histogram(pi0SimFD_opt[varstoplot[ind]], 100, density = True)
				        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
				        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
				        axs[yind, xind].hist(pi0ExpFD_corrected[varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
				        axs[yind, xind].set_title(title[ind])
				        # axs[yind, xind].set_xlim([start, end])
				        if (unit[ind]):
				            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
				        else:
				            axs[yind, xind].set_xlabel(title[ind])
				plt.tight_layout()
				plt.savefig(outDir+"{}FD{}_{:.3f}_{:.3f}_{:d}.pdf".format(polarity, i, sigma_opt, correction, sector))
				plt.clf()

			print("corrections in sector{}".format(sector), correction_s)
			print("smearings in sector{}".format(sector), sigmas_opt_s)

		print("done")


	def SmearingV3(self, df, sigma1, pol):
		df_epg = copy(df)
		if pol == "inbending":
			regulator = (1/(1+np.exp(-(df_epg.loc[df_epg["Psector"]<7, "Pp"]-0.5)/0.05)))
		elif pol == "outbending":
			regulator = (1/(1+np.exp(-(df_epg.loc[df_epg["Psector"]<7, "Pp"]-0.6)/0.05)))
		df_epg.loc[df_epg["Psector"]<7, "Pp"] = df_epg.loc[df_epg["Psector"]<7, "Pp"]*np.random.normal(1, regulator*sigma1, len(df_epg.loc[df_epg.Psector<7]))
		df_epg.loc[:, 'Pe'] = np.sqrt( df_epg.Pp**2 + M**2)
		df_epg.loc[:, "Ppx"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.cos(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppy"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.sin(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppz"] = df_epg.loc[:, "Pp"]*np.cos(np.radians(df_epg.loc[:, "Ptheta"]))

		self.df_epg = df_epg

	def CorrectionV3(self, df, correction):
		df_epg = copy(df)
		df_epg.loc[df_epg["Psector"]<7, "Pp"] = df_epg.loc[df_epg["Psector"]<7, "Pp"] + correction
		df_epg.loc[:, 'Pe'] = np.sqrt( df_epg.Pp**2 + M**2)
		df_epg.loc[:, "Ppx"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.cos(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppy"] = df_epg.loc[:, "Pp"]*np.sin(np.radians(df_epg.loc[:, "Ptheta"]))*np.sin(np.radians(df_epg.loc[:, "Pphi"]))
		df_epg.loc[:, "Ppz"] = df_epg.loc[:, "Pp"]*np.cos(np.radians(df_epg.loc[:, "Ptheta"]))

		self.df_epg = df_epg


	def saveDVCSvars(self):
	    #set up dvcs variables
	    df_epg = self.df_epg

	    ele = [df_epg['Epx'], df_epg['Epy'], df_epg['Epz']]
	    df_epg.loc[:, 'Ep'] = mag(ele)
	    df_epg.loc[:, 'Ee'] = getEnergy(ele, me)
	    df_epg.loc[:, 'Etheta'] = getTheta(ele)
	    df_epg.loc[:, 'Ephi'] = getPhi(ele)

	    pro = [df_epg['Ppx'], df_epg['Ppy'], df_epg['Ppz']]

	    gam = [df_epg['Gpx'], df_epg['Gpy'], df_epg['Gpz']]
	    df_epg.loc[:, 'Gp'] = mag(gam)
	    df_epg.loc[:, 'Ge'] = getEnergy(gam, 0)
	    df_epg.loc[:, 'Gtheta'] = getTheta(gam)
	    df_epg.loc[:, 'Gphi'] = getPhi(gam)

	    Ppt = mag([df_epg['Ppx'], df_epg['Ppy'], 0])

	    VGS = [-df_epg['Epx'], -df_epg['Epy'], pbeam - df_epg['Epz']]
	    v3l = cross(beam, ele)
	    v3h = cross(pro, VGS)
	    v3g = cross(VGS, gam)
	    VmissG = [-df_epg["Epx"] - df_epg["Ppx"], -df_epg["Epy"] - df_epg["Ppy"],
	              pbeam - df_epg["Epz"] - df_epg["Ppz"]]
	    VmissP = [-(df_epg["Epx"] + df_epg["Gpx"]), -(df_epg["Epy"] + df_epg["Gpy"]),
	              -(-pbeam + df_epg["Epz"] + df_epg["Gpz"])]
	    Vmiss = [-(df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"]), -(df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"]),
	             -(-pbeam + df_epg["Epz"] + df_epg["Ppz"] + df_epg["Gpz"])]
	    costheta = cosTheta(VGS, gam)

	    df_epg.loc[:, 'Mpx'], df_epg.loc[:, 'Mpy'], df_epg.loc[:, 'Mpz'] = Vmiss

	    # binning kinematics
	    df_epg.loc[:,'Q2'] = -((ebeam - df_epg['Ee'])**2 - mag2(VGS))
	    df_epg.loc[:,'nu'] = (ebeam - df_epg['Ee'])
	    df_epg.loc[:,'y'] = df_epg['nu']/ebeam
	    df_epg.loc[:,'xB'] = df_epg['Q2'] / 2.0 / M / df_epg['nu']
	    df_epg.loc[:,'t1'] = 2 * M * (df_epg['Pe'] - M)
	    df_epg.loc[:,'t2'] = (M * df_epg['Q2'] + 2 * M * df_epg['nu'] * (df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta))\
	    / (M + df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta)
	    df_epg.loc[:,'W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epg['Ee'])**2 - mag2(VGS)))

	    # trento angles
	    df_epg.loc[:,'phi1'] = angle(v3l, v3h)
	    df_epg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
	                              df_epg['phi1'], df_epg['phi1'])
	    df_epg.loc[:,'phi2'] = angle(v3l, v3g)
	    df_epg.loc[:,'phi2'] = np.where(dot(v3l, gam) <
	                              0, 360.0 - df_epg['phi2'], df_epg['phi2'])

	    # exclusivity variables
	    df_epg.loc[:,'MM2_epg'] = (-M - ebeam + df_epg["Ee"] +
	                         df_epg["Pe"] + df_epg["Ge"])**2 - mag2(Vmiss)
	    df_epg.loc[:,'ME_epg'] = (M + ebeam - df_epg["Ee"] - df_epg["Pe"] - df_epg["Ge"])
	    df_epg.loc[:,'MM2_ep'] = (-M - ebeam + df_epg["Ee"] + df_epg["Pe"])**2 - mag2(VmissG)
	    df_epg.loc[:,'MM2_eg'] = (-M - ebeam + df_epg["Ee"] + df_epg["Ge"])**2 - mag2(VmissP)
	    df_epg.loc[:,'MPt'] = np.sqrt((df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"])**2 +
	                            (df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"])**2)
	    df_epg.loc[:,'coneAngle'] = angle(ele, gam)
	    df_epg.loc[:,'reconGam'] = angle(gam, VmissG)
	    df_epg.loc[:,'coplanarity'] = angle(v3h, v3g)

	    eps = 2*M*df_epg.xB / np.sqrt(df_epg.Q2)
	    df_epg.loc[:,'ycol1'] = (df_epg.Q2-df_epg.t2)/(df_epg.Q2-df_epg.xB*df_epg.t2)
	    df_epg.loc[:,'ycol2'] = 1 - (1-df_epg.xB)*df_epg.t2/df_epg.Q2
	    df_epg.loc[:,'ymax1'] = 2*(np.sqrt(1+eps**2)-1)/(eps**2)
	    df_epg.loc[:,'ymax2'] = 1 - (M**2)*(df_epg.xB**2)/df_epg.Q2
	    df_epg.loc[:,'tmin1'] = df_epg.Q2*(2*(1-df_epg.xB)*(1-np.sqrt(1+eps**2))+eps**2)/(4*df_epg.xB*(1-df_epg.xB) + eps**2)
	    df_epg.loc[:,'tmin2'] = M*M*(df_epg.xB**2)/(1-df_epg.xB+df_epg.xB*M*M/df_epg.Q2)
	    df_epg.loc[:,'tcol'] = df_epg.Q2*(df_epg.Q2-2*df_epg.xB*M*ebeam)/df_epg.xB/(df_epg.Q2-2*M*ebeam)

	    self.df_epg = df_epg

	def saveDVpi0Pvars(self):
	    #set up pi0 variables
	    df_epgg = copy(self.df_epg)

	    # useful objects
	    ele = [df_epgg['Epx'], df_epgg['Epy'], df_epgg['Epz']]
	    df_epgg.loc[:, 'Ep'] = mag(ele)
	    df_epgg.loc[:, 'Ee'] = getEnergy(ele, me)
	    df_epgg.loc[:, 'Etheta'] = getTheta(ele)
	    df_epgg.loc[:, 'Ephi'] = getPhi(ele)

	    pro = [df_epgg['Ppx'], df_epgg['Ppy'], df_epgg['Ppz']]

	    gam = [df_epgg['Gpx'], df_epgg['Gpy'], df_epgg['Gpz']]
	    df_epgg.loc[:, 'Gp'] = mag(gam)
	    df_epgg.loc[:, 'Ge'] = getEnergy(gam, 0)
	    df_epgg.loc[:, 'Gtheta'] = getTheta(gam)
	    df_epgg.loc[:, 'Gphi'] = getPhi(gam)

	    gam2 = [df_epgg['Gpx2'], df_epgg['Gpy2'], df_epgg['Gpz2']]
	    df_epgg.loc[:, 'Gp2'] = mag(gam2)
	    df_epgg.loc[:,'Ge2'] = getEnergy(gam2, 0)
	    df_epgg.loc[:, 'Gtheta2'] = getTheta(gam2)
	    df_epgg.loc[:, 'Gphi2'] = getPhi(gam2)

	    pi0 = vecAdd(gam, gam2)
	    VGS = [-df_epgg['Epx'], -df_epgg['Epy'], pbeam - df_epgg['Epz']]
	    v3l = cross(beam, ele)
	    v3h = cross(pro, VGS)
	    v3g = cross(VGS, gam)
	    v3pi0 = cross(VGS, pi0)

	    VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
	                df_epgg["Ppy"], pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]
	    VmissP = [-df_epgg["Epx"] - df_epgg["Gpx"] - df_epgg["Gpx2"], -df_epgg["Epy"] -
	                df_epgg["Gpy"] - df_epgg["Gpy2"], pbeam - df_epgg["Epz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
	    Vmiss = [-df_epgg["Epx"] - df_epgg["Ppx"] - df_epgg["Gpx"] - df_epgg["Gpx2"],
	                -df_epgg["Epy"] - df_epgg["Ppy"] - df_epgg["Gpy"] - df_epgg["Gpy2"],
	                pbeam - df_epgg["Epz"] - df_epgg["Ppz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
	    costheta = cosTheta(VGS, gam)

	    df_epgg.loc[:, 'Mpx'], df_epgg.loc[:, 'Mpy'], df_epgg.loc[:, 'Mpz'] = Vmiss

	    # binning kinematics
	    df_epgg.loc[:,'Q2'] = -((ebeam - df_epgg['Ee'])**2 - mag2(VGS))
	    df_epgg.loc[:,'nu'] = (ebeam - df_epgg['Ee'])
	    df_epgg.loc[:,'y'] = df_epgg['nu']/ebeam
	    df_epgg.loc[:,'xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
	    df_epgg.loc[:,'t1'] = 2 * M * (df_epgg['Pe'] - M)
	    df_epgg.loc[:,'t2'] = (M * df_epgg['Q2'] + 2 * M * df_epgg['nu'] * (df_epgg['nu'] - np.sqrt(df_epgg['nu'] * df_epgg['nu'] + df_epgg['Q2']) * costheta))\
	    / (M + df_epgg['nu'] - np.sqrt(df_epgg['nu'] * df_epgg['nu'] + df_epgg['Q2']) * costheta)
	    df_epgg.loc[:,'W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epgg['Ee'])**2 - mag2(VGS)))
	    df_epgg.loc[:,'MPt'] = np.sqrt((df_epgg["Epx"] + df_epgg["Ppx"] + df_epgg["Gpx"] + df_epgg["Gpx2"])**2 +
	                             (df_epgg["Epy"] + df_epgg["Ppy"] + df_epgg["Gpy"] + df_epgg["Gpy2"])**2)
	    # trento angles
	    df_epgg.loc[:,'phi1'] = angle(v3l, v3h)
	    df_epgg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
	                              df_epgg['phi1'], df_epgg['phi1'])
	    df_epgg.loc[:,'phi2'] = angle(v3l, v3g)
	    df_epgg.loc[:,'phi2'] = np.where(dot(v3l, gam) <
	                              0, 360.0 - df_epgg['phi2'], df_epgg['phi2'])

	    # exclusivity variables
	    df_epgg.loc[:,'MM2_ep'] = (-M - ebeam + df_epgg["Ee"] +
	                         df_epgg["Pe"])**2 - mag2(VmissPi0)
	    df_epgg.loc[:,'MM2_egg'] = (-M - ebeam + df_epgg["Ee"] +
	                         df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(VmissP)
	    df_epgg.loc[:,'MM2_epgg'] = (-M - ebeam + df_epgg["Ee"] + df_epgg["Pe"] +
	                         df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(Vmiss)
	    df_epgg.loc[:,'ME_epgg'] = (M + ebeam - df_epgg["Ee"] - df_epgg["Pe"] - df_epgg["Ge"] - df_epgg["Ge2"])
	    df_epgg.loc[:,'Mpi0'] = pi0InvMass(gam, gam2)
	    df_epgg.loc[:,'reconPi'] = angle(VmissPi0, pi0)
	    df_epgg.loc[:,"Pie"] = df_epgg['Ge'] + df_epgg['Ge2']
	    df_epgg.loc[:,'coplanarity'] = angle(v3h, v3pi0)
	    df_epgg.loc[:,'coneAngle1'] = angle(ele, gam)
	    df_epgg.loc[:,'coneAngle2'] = angle(ele, gam2)
	    
	    df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)

	    self.df_epg = df_epgg

	def makeDVCS(self, pol = "inbending"):
	    #make dvcs pairs
	    df_dvcs = self.df_epg

	    #common cuts
	    cut_xBupper = df_dvcs["xB"] < 1  # xB
	    cut_xBlower = df_dvcs["xB"] > 0  # xB
	    cut_Q2 = df_dvcs["Q2"] > 1  # Q2
	    cut_W = df_dvcs["W"] > 2  # W
	    cut_Ee = df_dvcs["Ee"] > 2  # Ee
	    cut_Ge = df_dvcs["Ge"] > 1.8  # Ge
	    cut_Esector = (df_dvcs["Esector"]!=df_dvcs["Gsector"])
	    cut_Psector = ~( ((df_dvcs["Pstat"]//10)%10>0) & (df_dvcs["Psector"]==df_dvcs["Gsector"]))
	    cut_Ppmax = df_dvcs.Pp < 2  # Pp
	    # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
	    cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Esector & cut_Psector & cut_Ppmax

	    df_dvcs = df_dvcs[cut_common]

	    # proton reconstruction quality
	    # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) & (df_dvcs.loc[:, "Ptheta"]<35)
	    # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) & (df_dvcs.loc[:, "Ptheta"]>45) & (df_dvcs.loc[:, "Ptheta"]<65)
	    # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) #& (df_dvcs.loc[:, "Ptheta"]<37)
	    # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) #& (df_dvcs.loc[:, "Ptheta"]<66) #& (df_dvcs.loc[:, "Ptheta"]>40) 
	    # cut_proton = (cut_FD_proton)|(cut_CD_proton)
	    #(cut_FD_proton)|(cut_CD_proton)

	    df_dvcs.loc[:, "config"] = 0

	    if pol == "inbending":
	        #CDFT
	        cut_Pp1_CDFT = df_dvcs.Pp > 0.3  # Pp
	        cut_Psector_CDFT = df_dvcs.Psector>7
	        cut_Ptheta1_CDFT = df_dvcs.Ptheta<65
	        cut_Ptheta2_CDFT = df_dvcs.Ptheta>30
	        cut_Gsector_CDFT = df_dvcs.Gsector>7
	        cut_GFid_CDFT = df_dvcs.GFid==1
	        cut_mmep1_CDFT = df_dvcs["MM2_ep"] < 0.601  # mmep
	        cut_mmep2_CDFT = df_dvcs["MM2_ep"] > -0.528  # mmep
	        cut_mmeg1_CDFT = df_dvcs["MM2_eg"] < 1.569  # mmeg
	        cut_mmeg2_CDFT = df_dvcs["MM2_eg"] > 0.246  # mmeg
	        cut_meepg1_CDFT = df_dvcs["ME_epg"] < 0.484 # meepg
	        cut_meepg2_CDFT = df_dvcs["ME_epg"] > -0.439  # meepg
	        cut_cone1_CDFT = df_dvcs["coneAngle"] < 0.133 * df_dvcs.Gp**2 + 1.140*df_dvcs.Gp + 10.129  # coneangle
	        cut_cone2_CDFT = df_dvcs["coneAngle"] > 0.115 * df_dvcs.Gp**2 + 0.702*df_dvcs.Gp + 5.359 #12.159  # coneangle
	        cut_mpt_CDFT = df_dvcs["MPt"] < 0.0958  # mpt
	        cut_recon_CDFT = df_dvcs["reconGam"] < 0.659 + 0.4*(df_dvcs.Pp-0.3)  # recon gam angle
	        cut_coplanarity_CDFT = df_dvcs["coplanarity"] < 8.251  # coplanarity angle
	        cut_mmepg1_CDFT = np.abs(df_dvcs["MM2_epg"]) < 0.0122  # mmepg
	        cut_mmepg2_CDFT = np.abs(df_dvcs["MM2_epg"]) > -0.0151  # mmepg

	        cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT &
	                    cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mmeg1_CDFT & cut_mmeg2_CDFT &
	                    cut_meepg1_CDFT & cut_meepg2_CDFT & cut_cone1_CDFT & cut_cone2_CDFT &
	                    cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepg1_CDFT & cut_mmepg2_CDFT)


	        #CD
	        cut_Pp1_CD = df_dvcs.Pp > 0.3  # Pp
	        cut_Psector_CD = df_dvcs.Psector>7
	        cut_Ptheta_CD = df_dvcs.Ptheta<60
	        cut_Gsector_CD = df_dvcs.Gsector<7
	        cut_GFid_CD = df_dvcs.GFid==1
	        cut_mmep1_CD = df_dvcs["MM2_ep"] < 0.428  # mmep
	        cut_mmep2_CD = df_dvcs["MM2_ep"] > -0.418  # mmep
	        cut_mmeg1_CD = df_dvcs["MM2_eg"] < 2.321  # mmeg
	        cut_mmeg2_CD = df_dvcs["MM2_eg"] > -0.428 # mmeg
	        cut_meepg1_CD = df_dvcs["ME_epg"] < 1.045  # meepg
	        cut_meepg2_CD = df_dvcs["ME_epg"] > -0.922  # meepg
	        cut_cone1_CD = df_dvcs["coneAngle"] < 0.416 * df_dvcs.Gp**2 - 5.401*df_dvcs.Gp + 50.980  # coneangle
	        cut_cone2_CD = df_dvcs["coneAngle"] > 0.278 * df_dvcs.Gp**2 - 2.711*df_dvcs.Gp + 23.021  # coneangle
	        cut_cone1_CD = df_dvcs["coneAngle"] < 32.817  # coneangle
	        cut_cone2_CD = df_dvcs["coneAngle"] > 5#12.791  # coneangle
	        cut_mpt_CD = df_dvcs["MPt"] < 0.182  # mpt
	        cut_recon_CD = df_dvcs["reconGam"] < 1.5#0.876  # recon gam angle
	        cut_coplanarity_CD = df_dvcs["coplanarity"] < 8.352  # coplanarity angle
	        cut_mmepg1_CD = np.abs(df_dvcs["MM2_epg"]) < 0.0241  # mmepg
	        cut_mmepg2_CD = np.abs(df_dvcs["MM2_epg"]) > -0.0275  # mmepg

	        cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta_CD & cut_Gsector_CD & cut_GFid_CD &
	                    cut_mmep1_CD & cut_mmep2_CD & cut_mmeg1_CD & cut_mmeg2_CD &
	                    cut_meepg1_CD & cut_meepg2_CD & cut_cone1_CD & cut_cone2_CD &
	                    cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepg1_CD & cut_mmepg2_CD)

	        #FD
	        cut_Pp1_FD = df_dvcs.Pp > 0.42  # Pp
	        cut_Psector_FD = df_dvcs.Psector<7
	        cut_Ptheta_FD = df_dvcs.Ptheta>2.477
	        cut_Gsector_FD = df_dvcs.Gsector<7
	        cut_GFid_FD = df_dvcs.GFid==1
	        cut_mmep1_FD = df_dvcs["MM2_ep"] < 0.406  # mmep
	        cut_mmep2_FD = df_dvcs["MM2_ep"] > -0.396  # mmep
	        cut_mmeg1_FD = df_dvcs["MM2_eg"] < 2.045  # mmeg
	        cut_mmeg2_FD = df_dvcs["MM2_eg"] > -0.166  # mmeg
	        cut_meepg1_FD = df_dvcs["ME_epg"] < 0.996  # meepg
	        cut_meepg2_FD = df_dvcs["ME_epg"] > -0.864 # meepg
	        cut_cone1_FD = df_dvcs["coneAngle"] < 0.0673 * df_dvcs.Gp**2 - 0.752*df_dvcs.Gp + 41.678  # coneangle
	        cut_cone2_FD = df_dvcs["coneAngle"] > 0.712 * df_dvcs.Gp**2 - 7.676*df_dvcs.Gp + 46.720  # coneangle
	        cut_mpt_FD = df_dvcs["MPt"] < 0.184  # mpt
	        cut_recon_FD = df_dvcs["reconGam"] < 2.5#1.181  # recon gam angle
	        cut_coplanarity_FD = df_dvcs["coplanarity"] < 13.178  # coplanarity angle - no cut
	        cut_mmepg1_FD = np.abs(df_dvcs["MM2_epg"]) < 0.0267  # mmepg
	        cut_mmepg2_FD = np.abs(df_dvcs["MM2_epg"]) > -0.0302  # mmepg

	        cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta_FD & cut_Gsector_FD & cut_GFid_FD &
	                    cut_mmep1_FD & cut_mmep2_FD & cut_mmeg1_FD & cut_mmeg2_FD &
	                    cut_meepg1_FD & cut_meepg2_FD & cut_cone1_FD & cut_cone2_FD &
	                    cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepg1_FD & cut_mmepg2_FD)

	    elif pol == "outbending":
	        #CDFT
	        cut_Pp1_CDFT = df_dvcs.Pp > 0.3  # Pp
	        cut_Psector_CDFT = df_dvcs.Psector>7
	        cut_Ptheta1_CDFT = df_dvcs.Ptheta<65
	        cut_Ptheta2_CDFT = df_dvcs.Ptheta>40
	        cut_Gsector_CDFT = df_dvcs.Gsector>7
	        cut_GFid_CDFT = df_dvcs.GFid==1
	        cut_mmep1_CDFT = df_dvcs["MM2_ep"] < 0.567  # mmep
	        cut_mmep2_CDFT = df_dvcs["MM2_ep"] > -0.479  # mmep
	        cut_mmeg1_CDFT = df_dvcs["MM2_eg"] < 1.602  # mmeg
	        cut_mmeg2_CDFT = df_dvcs["MM2_eg"] > 0.223  # mmeg
	        cut_meepg1_CDFT = df_dvcs["ME_epg"] < 0.481 # meepg
	        cut_meepg2_CDFT = df_dvcs["ME_epg"] > -0.436  # meepg
	        cut_cone1_CDFT = df_dvcs["coneAngle"] < 0.133 * df_dvcs.Gp**2 + 1.140*df_dvcs.Gp + 10.129  # coneangle
	        cut_cone2_CDFT = df_dvcs["coneAngle"] > 0.115 * df_dvcs.Gp**2 + 0.702*df_dvcs.Gp + 5.359 #12.159  # coneangle
	        cut_mpt_CDFT = df_dvcs["MPt"] < 0.100  # mpt
	        cut_recon_CDFT = df_dvcs["reconGam"] < 0.697 + 0.4*(df_dvcs.Pp-0.3) # recon gam angle
	        cut_coplanarity_CDFT = df_dvcs["coplanarity"] < 8.479  # coplanarity angle
	        cut_mmepg1_CDFT = np.abs(df_dvcs["MM2_epg"]) < 0.0128  # mmepg
	        cut_mmepg2_CDFT = np.abs(df_dvcs["MM2_epg"]) > -0.0155  # mmepg

	        cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT &
	                    cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mmeg1_CDFT & cut_mmeg2_CDFT &
	                    cut_meepg1_CDFT & cut_meepg2_CDFT & cut_cone1_CDFT & cut_cone2_CDFT &
	                    cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepg1_CDFT & cut_mmepg2_CDFT)


	        #CD
	        cut_Pp1_CD = df_dvcs.Pp > 0.3  # Pp
	        cut_Psector_CD = df_dvcs.Psector>7
	        cut_Ptheta_CD = df_dvcs.Ptheta<60
	        cut_Gsector_CD = df_dvcs.Gsector<7
	        cut_GFid_CD = df_dvcs.GFid==1
	        cut_mmep1_CD = df_dvcs["MM2_ep"] < 0.358  # mmep
	        cut_mmep2_CD = df_dvcs["MM2_ep"] > -0.336  # mmep
	        cut_mmeg1_CD = df_dvcs["MM2_eg"] < 2.207  # mmeg
	        cut_mmeg2_CD = df_dvcs["MM2_eg"] > -0.387 # mmeg
	        cut_meepg1_CD = df_dvcs["ME_epg"] < 0.888  # meepg
	        cut_meepg2_CD = df_dvcs["ME_epg"] > -0.825  # meepg
	        cut_cone1_CD = df_dvcs["coneAngle"] < 0.416 * df_dvcs.Gp**2 - 5.401*df_dvcs.Gp + 50.980  # coneangle
	        cut_cone2_CD = df_dvcs["coneAngle"] > 0.278 * df_dvcs.Gp**2 - 2.711*df_dvcs.Gp + 23.021  # coneangle
	        cut_mpt_CD = df_dvcs["MPt"] < 0.183  # mpt
	        cut_recon_CD = df_dvcs["reconGam"] < 1.5#0.954  # recon gam angle
	        cut_coplanarity_CD = df_dvcs["coplanarity"] < 7.347  # coplanarity angle
	        cut_mmepg1_CD = np.abs(df_dvcs["MM2_epg"]) < 0.0211  # mmepg
	        cut_mmepg2_CD = np.abs(df_dvcs["MM2_epg"]) > -0.0241  # mmepg

	        cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta_CD & cut_Gsector_CD & cut_GFid_CD &
	                    cut_mmep1_CD & cut_mmep2_CD & cut_mmeg1_CD & cut_mmeg2_CD &
	                    cut_meepg1_CD & cut_meepg2_CD & cut_cone1_CD & cut_cone2_CD &
	                    cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepg1_CD & cut_mmepg2_CD)

	        #FD
	        cut_Pp1_FD = df_dvcs.Pp > 0.5  # Pp
	        cut_Psector_FD = df_dvcs.Psector<7
	        cut_Ptheta_FD = df_dvcs.Ptheta>2.477
	        cut_Gsector_FD = df_dvcs.Gsector<7
	        cut_GFid_FD = df_dvcs.GFid==1
	        cut_mmep1_FD = df_dvcs["MM2_ep"] < 0.439  # mmep
	        cut_mmep2_FD = df_dvcs["MM2_ep"] > -0.418  # mmep
	        cut_mmeg1_FD = df_dvcs["MM2_eg"] < 2.049  # mmeg
	        cut_mmeg2_FD = df_dvcs["MM2_eg"] > -0.205  # mmeg
	        cut_meepg1_FD = df_dvcs["ME_epg"] < 0.947  # meepg
	        cut_meepg2_FD = df_dvcs["ME_epg"] > -0.796  # meepg
	        cut_cone1_FD = df_dvcs["coneAngle"] < 0.0673 * df_dvcs.Gp**2 - 0.752*df_dvcs.Gp + 41.678  # coneangle
	        cut_cone2_FD = df_dvcs["coneAngle"] > 0.712 * df_dvcs.Gp**2 - 7.676*df_dvcs.Gp + 46.720  # coneangle
	        cut_mpt_FD = df_dvcs["MPt"] < 0.186  # mpt
	        cut_recon_FD = df_dvcs["reconGam"] < 2.5#1.662  # recon gam angle
	        cut_coplanarity_FD = df_dvcs["coplanarity"] < 11.685  # coplanarity angle - no cut
	        cut_mmepg1_FD = np.abs(df_dvcs["MM2_epg"]) < 0.0294  # mmepg
	        cut_mmepg2_FD = np.abs(df_dvcs["MM2_epg"]) > -0.0346  # mmepg

	        cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta_FD & cut_Gsector_FD & cut_GFid_FD &
	                    cut_mmep1_FD & cut_mmep2_FD & cut_mmeg1_FD & cut_mmeg2_FD &
	                    cut_meepg1_FD & cut_meepg2_FD & cut_cone1_FD & cut_cone2_FD &
	                    cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepg1_FD & cut_mmepg2_FD)            

	    df_dvcs.loc[cut_CDFT, "config"] = 3
	    df_dvcs.loc[cut_CD, "config"] = 2
	    df_dvcs.loc[cut_FD, "config"] = 1

	    df_dvcs = df_dvcs[df_dvcs.config>0]

	    # dealing with duplicates
	    df_dvcs = df_dvcs.sort_values(by=['reconGam', 'Ge', 'Pe'], ascending = [True, False, False])
	    df_dvcs = df_dvcs.loc[~df_dvcs.event.duplicated(), :]
	    df_dvcs = df_dvcs.sort_values(by='event')

	    self.df_epg = df_dvcs                        

	def makeDVpi0P(self, pol = "inbending"):
	    #make dvpi0 pairs
	    df_dvpi0p = copy(self.df_epg)

	    #common cuts
	    cut_xBupper = df_dvpi0p.loc[:, "xB"] < 1  # xB
	    cut_xBlower = df_dvpi0p.loc[:, "xB"] > 0  # xB
	    cut_Q2 = df_dvpi0p.loc[:, "Q2"] > 1  # Q2
	    cut_W = df_dvpi0p.loc[:, "W"] > 2  # W
	    cut_Ee = df_dvpi0p["Ee"] > 2  # Ee
	    cut_Ge2 = df_dvpi0p["Ge2"] > 0.6  # Ge cut. Ge>3 for DVCS module.
	    cut_Esector = (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector"]) & (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector2"])
	    cut_Psector = ~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector"]) ) & ~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector2"]) )
	    cut_Ppmax = df_dvpi0p.Pp < 2  # Pp
	    # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
	    cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge2 & cut_Esector & cut_Psector & cut_Ppmax

	    df_dvpi0p = df_dvpi0p[cut_common]

	    # proton reconstruction quality
	    # cut_FD_proton = (df_epgg.loc[:, "Psector"]<7) & (df_epgg.loc[:, "Ptheta"]<35)
	    # cut_CD_proton = (df_epgg.loc[:, "Psector"]>7) & (df_epgg.loc[:, "Ptheta"]>45) & (df_epgg.loc[:, "Ptheta"]<65)
	    # cut_proton = (cut_FD_proton)|(cut_CD_proton)
	    cut_proton = 1

	    df_dvpi0p.loc[:, "config"] = 0

	    if pol == "inbending":
	        #CDFT
	        cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
	        cut_Psector_CDFT = df_dvpi0p.Psector>7
	        cut_Ptheta_CDFT = df_dvpi0p.Ptheta<60
	        cut_Gsector_CDFT = df_dvpi0p.Gsector>7
	        cut_GFid2_CDFT = df_dvpi0p.GFid2==1
	        cut_mpi01_CDFT = df_dvpi0p["Mpi0"] < 0.157  # mpi0
	        cut_mpi02_CDFT = df_dvpi0p["Mpi0"] > 0.118  # mpi0
	        cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] < 0.914  # mmep
	        cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] > -0.715  # mmep
	        cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] < 2.155  # mmegg
	        cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] > -0.417  # mmegg
	        cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] < 0.799  # meepgg
	        cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] > -0.792  # meepgg
	        cut_mpt_CDFT = df_dvpi0p["MPt"] < 0.189  # mpt
	        cut_recon_CDFT = df_dvpi0p["reconPi"] < 1.468  # recon gam angle
	        cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] < 15.431  # coplanarity angle
	        cut_mmepgg1_CDFT = np.abs(df_dvpi0p["MM2_epgg"]) < 0.0440  # mmepgg
	        cut_mmepgg2_CDFT = np.abs(df_dvpi0p["MM2_epgg"]) > -0.0478  # mmepgg

	        cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta_CDFT & cut_Gsector_CDFT & cut_GFid2_CDFT &
	                    cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
	                    cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
	                    cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


	        #CD
	        cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
	        cut_Psector_CD = df_dvpi0p.Psector>7
	        cut_Ptheta_CD = df_dvpi0p.Ptheta<60
	        cut_Gsector_CD = df_dvpi0p.Gsector<7
	        cut_Gsector2_CD = df_dvpi0p.Gsector2<7
	        cut_GFid_CD = df_dvpi0p.GFid==1
	        cut_GFid2_CD = df_dvpi0p.GFid2==1
	        cut_mpi01_CD = df_dvpi0p["Mpi0"] < 0.162  # mpi0
	        cut_mpi02_CD = df_dvpi0p["Mpi0"] > 0.107  # mpi0
	        cut_mmep1_CD = df_dvpi0p["MM2_ep"] < 0.354  # mmep
	        cut_mmep2_CD = df_dvpi0p["MM2_ep"] > -0.283  # mmep
	        cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < 1.922  # mmegg
	        cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > 0.007  # mmegg
	        cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < 0.822  # meepgg
	        cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > -0.677  # meepgg
	        cut_mpt_CD = df_dvpi0p["MPt"] < 0.176  # mpt
	        cut_recon_CD = df_dvpi0p["reconPi"] < 1.476  # recon gam angle
	        cut_coplanarity_CD = df_dvpi0p["coplanarity"] < 10.203  # coplanarity angle
	        cut_mmepgg1_CD = np.abs(df_dvpi0p["MM2_epgg"]) < 0.0208  # mmepgg
	        cut_mmepgg2_CD = np.abs(df_dvpi0p["MM2_epgg"]) > -0.0250  # mmepgg

	        cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta_CD & cut_Gsector_CD & cut_Gsector2_CD & 
	                    cut_GFid_CD & cut_GFid2_CD &
	                    cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
	                    cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
	                    cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

	        #FD
	        cut_Pp1_FD = df_dvpi0p.Pp > 0.42  # Pp
	        cut_Psector_FD = df_dvpi0p.Psector<7
	        cut_Ptheta_FD = df_dvpi0p.Ptheta>2.477
	        cut_Gsector_FD = df_dvpi0p.Gsector<7
	        cut_Gsector2_FD = df_dvpi0p.Gsector2<7
	        cut_GFid_FD = df_dvpi0p.GFid==1
	        cut_GFid2_FD = df_dvpi0p.GFid2==1
	        cut_mpi01_FD = df_dvpi0p["Mpi0"] < 0.178  # mpi0
	        cut_mpi02_FD = df_dvpi0p["Mpi0"] > 0.0910  # mpi0
	        cut_mmep1_FD = df_dvpi0p["MM2_ep"] < 0.335  # mmep
	        cut_mmep2_FD = df_dvpi0p["MM2_ep"] > -0.271  # mmep
	        cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < 1.762  # mmegg
	        cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > 0.117  # mmegg
	        cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < 0.816 # meepgg
	        cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > -0.685  # meepgg
	        cut_mpt_FD = df_dvpi0p["MPt"] < 0.180  # mpt
	        cut_recon_FD = df_dvpi0p["reconPi"] < 1.363  # recon gam angle
	        cut_coplanarity_FD = df_dvpi0p["coplanarity"] < 9.190  # coplanarity angle
	        cut_mmepgg1_FD = np.abs(df_dvpi0p["MM2_epgg"]) < 0.0189  # mmepgg
	        cut_mmepgg2_FD = np.abs(df_dvpi0p["MM2_epgg"]) > -0.0224  # mmepgg

	        cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta_FD & cut_Gsector_FD & cut_Gsector2_FD &
	                    cut_GFid_FD & cut_GFid2_FD &
	                    cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
	                    cut_mmegg1_FD & cut_mmegg2_FD & cut_meepgg1_FD & cut_meepgg2_FD &
	                    cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepgg1_FD & cut_mmepgg2_FD)

	    elif pol == "outbending":
	        #CDFT
	        cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
	        cut_Psector_CDFT = df_dvpi0p.Psector>7
	        cut_Ptheta_CDFT = df_dvpi0p.Ptheta<60
	        cut_Gsector_CDFT = df_dvpi0p.Gsector>7
	        cut_GFid2_CDFT = df_dvpi0p.GFid2==1
	        cut_mpi01_CDFT = df_dvpi0p["Mpi0"] < 0.160  # mpi0
	        cut_mpi02_CDFT = df_dvpi0p["Mpi0"] > 0.115  # mpi0
	        cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] < 0.892  # mmep
	        cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] > -0.694  # mmep
	        cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] < 2.184  # mmegg
	        cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] > -0.412  # mmegg
	        cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] < 0.844  # meepgg
	        cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] > -0.806  # meepgg
	        cut_mpt_CDFT = df_dvpi0p["MPt"] < 0.210  # mpt
	        cut_recon_CDFT = df_dvpi0p["reconPi"] < 1.630  # recon gam angle
	        cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] < 17.817  # coplanarity angle
	        cut_mmepgg1_CDFT = np.abs(df_dvpi0p["MM2_epgg"]) < 0.0549  # mmepgg
	        cut_mmepgg2_CDFT = np.abs(df_dvpi0p["MM2_epgg"]) > -0.0575  # mmepgg

	        cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta_CDFT & cut_Gsector_CDFT & cut_GFid2_CDFT &
	                    cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
	                    cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
	                    cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


	        #CD
	        cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
	        cut_Psector_CD = df_dvpi0p.Psector>7
	        cut_Ptheta_CD = df_dvpi0p.Ptheta<60
	        cut_Gsector_CD = df_dvpi0p.Gsector<7
	        cut_Gsector2_CD = df_dvpi0p.Gsector2<7
	        cut_GFid_CD = df_dvpi0p.GFid==1
	        cut_GFid2_CD = df_dvpi0p.GFid2==1
	        cut_mpi01_CD = df_dvpi0p["Mpi0"] < 0.163  # mpi0
	        cut_mpi02_CD = df_dvpi0p["Mpi0"] > 0.106  # mpi0
	        cut_mmep1_CD = df_dvpi0p["MM2_ep"] < 0.294  # mmep
	        cut_mmep2_CD = df_dvpi0p["MM2_ep"] > -0.218  # mmep
	        cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < 1.876  # mmegg
	        cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > -0.0142  # mmegg
	        cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < 0.700  # meepgg
	        cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > -0.597  # meepgg
	        cut_mpt_CD = df_dvpi0p["MPt"] < 0.194  # mpt
	        cut_recon_CD = df_dvpi0p["reconPi"] < 1.761  # recon gam angle
	        cut_coplanarity_CD = df_dvpi0p["coplanarity"] < 9.530  # coplanarity angle
	        cut_mmepgg1_CD = np.abs(df_dvpi0p["MM2_epgg"]) < 0.0182  # mmepgg
	        cut_mmepgg2_CD = np.abs(df_dvpi0p["MM2_epgg"]) > -0.0219  # mmepgg

	        cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta_CD & cut_Gsector_CD & cut_Gsector2_CD & 
	                    cut_GFid_CD & cut_GFid2_CD &
	                    cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
	                    cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
	                    cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

	        #FD
	        cut_Pp1_FD = df_dvpi0p.Pp > 0.5  # Pp
	        cut_Psector_FD = df_dvpi0p.Psector<7
	        cut_Ptheta_FD = df_dvpi0p.Ptheta>2.477
	        cut_Gsector_FD = df_dvpi0p.Gsector<7
	        cut_Gsector2_FD = df_dvpi0p.Gsector2<7
	        cut_GFid_FD = df_dvpi0p.GFid==1
	        cut_GFid2_FD = df_dvpi0p.GFid2==1
	        cut_mpi01_FD = df_dvpi0p["Mpi0"] < 0.164  # mpi0
	        cut_mpi02_FD = df_dvpi0p["Mpi0"] > 0.105  # mpi0
	        cut_mmep1_FD = df_dvpi0p["MM2_ep"] < 0.323  # mmep
	        cut_mmep2_FD = df_dvpi0p["MM2_ep"] > -0.256  # mmep
	        cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < 1.828  # mmegg
	        cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > 0.0491  # mmegg
	        cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < 0.754  # meepgg
	        cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > -0.583  # meepgg
	        cut_mpt_FD = df_dvpi0p["MPt"] < 0.177  # mpt
	        cut_recon_FD = df_dvpi0p["reconPi"] < 1.940  # recon gam angle
	        cut_coplanarity_FD = df_dvpi0p["coplanarity"] < 7.498  # coplanarity angle
	        cut_mmepgg1_FD = np.abs(df_dvpi0p["MM2_epgg"]) < 0.0195  # mmepgg
	        cut_mmepgg2_FD = np.abs(df_dvpi0p["MM2_epgg"]) > -0.0240  # mmepgg

	        cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta_FD & cut_Gsector_FD & cut_Gsector2_FD &
	                    cut_GFid_FD & cut_GFid2_FD &
	                    cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
	                    cut_mmegg1_FD & cut_mmegg2_FD & cut_meepgg1_FD & cut_meepgg2_FD &
	                    cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepgg1_FD & cut_mmepgg2_FD)

	    df_dvpi0p.loc[cut_CDFT, "config"] = 3
	    df_dvpi0p.loc[cut_CD, "config"] = 2
	    df_dvpi0p.loc[cut_FD, "config"] = 1

	    df_dvpi0p = df_dvpi0p[df_dvpi0p.config>0]

	    #For an event, there can be two gg's passed conditions above.
	    #Take only one gg's that makes pi0 invariant mass
	    #This case is very rare.
	    #For now, duplicated proton is not considered.
	    df_dvpi0p = df_dvpi0p.sort_values(by=['closeness', 'Psector', 'Gsector'], ascending = [True, True, True])
	    df_dvpi0p = df_dvpi0p.loc[~df_dvpi0p.event.duplicated(), :]
	    df_dvpi0p = df_dvpi0p.sort_values(by='event')        
	    self.df_epg = df_dvpi0p #done with saving x



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--version", help="version", default="v2")
    parser.add_argument("-o","--outDir", help="output directory", default="/SimtoDat/v2")
    parser.add_argument("-p","--polarity", help="polarity", default="inbending")
    
    args = parser.parse_args()

    smearingDist = smearingDist(version = args.version, outDir = args.outDir, pol = args.polarity)
