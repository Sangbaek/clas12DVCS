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

	def __init__(self, version, outDir = "SimtoDat/v0/"):
		self.GetVersion(version)
		if version == "v0":
			self.MakeV0(outDir = outDir)
		if version == "v1":
			self.MakeV1()
		if version == "v2":
			self.MakeV2()

	def GetVersion(self, version):
		self.version = version
		self.outDir = "SimtoDat/{}".format(version)

	def MakeV0(self, outDir = "SimtoDat/v0/"):

		if self.version == "v0":
			inDirInb = "../nov2021/convPkl/"
			inDirOutb = "../nov2021/convPkl_outb/"
		else:
			inDirInb = "SimtoDat/v0/inb/"
			inDirOutb = "SimtoDat/v0/outb/"

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

		if isinstance(sigma, str):
			sigma = float(sigma)

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
		# self.CorrectingV0(epgExpInbCDFT, correction, mode = "epg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "inbending")
		# epgExpInbCDFT = self.df_epg
		# correction = 0.186/(1+np.exp(1.254*(pi0ExpInbCDFT.Gp-4.5)))
		# self.CorrectingV0(pi0ExpInbCDFT, correction, mode = "epgg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "inbending")
		# pi0ExpInbCDFT = self.df_epg
		# #performing correcting
		# correction = 0.186/(1+np.exp(1.254*(epgExpOutbCDFT.Gp-4.5)))
		# self.CorrectingV0(epgExpOutbCDFT, correction, mode = "epg")
		# self.saveDVCSvars()
		# self.makeDVCS(pol = "outbending")
		# epgExpOutbCDFT = self.df_epg
		# correction = 0.186/(1+np.exp(1.254*(pi0ExpOutbCDFT.Gp-4.5)))
		# self.CorrectingV0(pi0ExpOutbCDFT, correction, mode = "epgg")
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
			self.CorrectingV0(epgExpInbCDFT, correction, mode = "epg")
			self.saveDVCSvars()
			self.makeDVCS()
			epgExpInbCDFT_corrected = self.df_epg
			self.CorrectingV0(pi0ExpInbCDFT, correction, mode = "epgg")
			self.saveDVCSvars()
			self.makeDVCS()
			pi0ExpInbCDFT_corrected = self.df_epg

			#performing correcting
			self.CorrectingV0(epgExpOutbCDFT, correction, mode = "epg")
			self.saveDVCSvars()
			self.makeDVCS()
			epgExpOutbCDFT_corrected = self.df_epg
			self.CorrectingV0(pi0ExpOutbCDFT, correction, mode = "epgg")
			self.saveDVCSvars()
			self.makeDVCS()
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

	def CorrectingV0(self, df, correction, mode = "epg"):
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


		# def corr(df1, df2, df_exp, cont = 0, var = "ME_epg"):
		# 	hist1, bins = np.histogram(df1.loc[:, var], bins = 101)
		# 	hist2, _ = np.histogram(df2.loc[:, var], bins = bins)
		# 	hist_exp, _ = np.histogram(df_exp.loc[:, var], bins = bins)
		# 	bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
		# 	dist1 = hist1/np.sum(hist1)/(np.diff(bincenters)[0])
		# 	dist2 = hist2/np.sum(hist2)/(np.diff(bincenters)[0])
		# 	dist_exp = hist_exp/np.sum(hist_exp)/(np.diff(bincenters)[0])
		# 	simDist = (1-cont)*dist1 + cont*dist2
		# 	expDist = dist_exp
		# 	simPeak = bincenters[np.argmax(simDist)]
		# 	expPeak = bincenters[np.argmax(expDist)]
		# 	return expPeak - simPeak

		PpEdges = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
		sigma1s = np.linspace(0.05, 0.1, 1)
		sigma2s = np.linspace(0.5, 1, 1)
		sigma3s = np.linspace(1.8, 2.5, 1)
		# corrections = []
		sigma1s_temp = []
		sigma2s_temp = []
		sigma3s_temp = []

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

		for i in range(len(PpEdges)-1):

			distances = []
			PpMin = PpEdges[i]
			PpMax = PpEdges[i+1]

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

			contInb = 0
			contOutb = 0
			if len(epgExpInbCDFT_selected)*len(pi0SimInbCDFT_selected) > 0:
				contInb = len(bkgSimInbCDFT_selected)/len(pi0SimInbCDFT_selected)*len(pi0ExpInbCDFT_selected)/len(epgExpInbCDFT_selected)
			if len(epgExpInbCDFT_selected)*len(pi0SimOutbCDFT_selected) > 0:
				contOutb = len(bkgSimOutbCDFT_selected)/len(pi0SimOutbCDFT_selected)*len(pi0ExpOutbCDFT_selected)/len(epgExpOutbCDFT_selected)

			# correction1 = epgExpInbCDFT_selected.ME_epg.mean() - (1-contInb)*dvcsSimInbCDFT_selected.ME_epg.mean() - contInb*bkgSimInbCDFT_selected.ME_epg.mean()
			# correction2 = epgExpOutbCDFT_selected.ME_epg.mean() - (1-contOutb)*dvcsSimOutbCDFT_selected.ME_epg.mean() - contOutb*bkgSimOutbCDFT_selected.ME_epg.mean()
			## correction1 = corr(dvcsSimInbCDFT_selected, bkgSimInbCDFT_selected, epgExpInbCDFT_selected, cont = contInb, var = "ME_epg")
			## correction2 = corr(dvcsSimOutbCDFT_selected, bkgSimOutbCDFT_selected, epgExpOutbCDFT_selected, cont = contOutb, var = "ME_epg")
			# correction = (correction1+correction2)/2
			# print(correction1, contInb, correction2, contOutb)
			# corrections.append(correction)

			# #performing correcting
			# self.CorrectingV0(epgExpInbCDFT, correction, mode = "epg")
			# self.saveDVCSvars()
			# self.makeDVCS()
			# epgExpInbCDFT_corrected = self.df_epg
			# self.CorrectingV0(pi0ExpInbCDFT, correction, mode = "epgg")
			# self.saveDVCSvars()
			# self.makeDVCS()
			# pi0ExpInbCDFT_corrected = self.df_epg

			# #performing correcting
			# self.CorrectingV0(epgExpOutbCDFT, correction, mode = "epg")
			# self.saveDVCSvars()
			# self.makeDVCS()
			# epgExpOutbCDFT_corrected = self.df_epg
			# self.CorrectingV0(pi0ExpOutbCDFT, correction, mode = "epgg")
			# self.saveDVCSvars()
			# self.makeDVCS()
			# pi0ExpOutbCDFT_corrected = self.df_epg

			# epgExpInbCDFT_corrected = epgExpInbCDFT_corrected.loc[(epgExpInbCDFT_corrected.Ge>GeMin) & (epgExpInbCDFT_corrected.Ge<GeMax)]
			# epgExpOutbCDFT_corrected = epgExpOutbCDFT_corrected.loc[(epgExpOutbCDFT_corrected.Ge>GeMin) & (epgExpOutbCDFT_corrected.Ge<GeMax)]

			# epgExpInbCDFT_corrected = epgExpInbCDFT.loc[(epgExpInbCDFT.Ge>GeMin) & (epgExpInbCDFT.Ge<GeMax)]
			# epgExpOutbCDFT_corrected = epgExpOutbCDFT.loc[(epgExpOutbCDFT.Ge>GeMin) & (epgExpOutbCDFT.Ge<GeMax)]

			for sigma1 in sigma1s:
				for sigma2 in sigma2s:
					for sigma3 in sigma3s:

						print("smearing with {:.3f}{:.3f}{:.3f}".format(sigma1, sigma2, sigma3))

						#performing smearing
						self.SmearingV1(dvcsSimInbCDFT, sigma1, sigma2, sigma3)
						self.saveDVCSvars()
						self.makeDVCS(pol = "inbending")
						dvcsSimInbCDFT_smeared = self.df_epg

						self.SmearingV1(bkgSimInbCDFT, sigma1, sigma2, sigma3)
						self.saveDVCSvars()
						self.makeDVCS(pol = "inbending")
						bkgSimInbCDFT_smeared = self.df_epg

						self.SmearingV1(dvcsSimOutbCDFT, sigma1, sigma2, sigma3)
						self.saveDVCSvars()
						self.makeDVCS(pol = "outbending")
						dvcsSimOutbCDFT_smeared = self.df_epg

						self.SmearingV1(bkgSimOutbCDFT, sigma1, sigma2, sigma3)
						self.saveDVCSvars()
						self.makeDVCS(pol = "outbending")
						bkgSimOutbCDFT_smeared = self.df_epg

						dvcsSimInbCDFT_smeared = dvcsSimInbCDFT_smeared.loc[(dvcsSimInbCDFT_smeared.Pp>PpMin) & (dvcsSimInbCDFT_smeared.Pp<PpMax)]
						dvcsSimOutbCDFT_smeared = dvcsSimOutbCDFT_smeared.loc[(dvcsSimOutbCDFT_smeared.Pp>PpMin) & (dvcsSimOutbCDFT_smeared.Pp<PpMax)]

						distance1 = distance(dvcsSimInbCDFT_smeared, bkgSimInbCDFT_smeared, epgExpInbCDFT_selected, cont = contInb, var = "MM2_ep")
						distance2 = distance(dvcsSimOutbCDFT_smeared, bkgSimOutbCDFT_smeared, epgExpOutbCDFT_selected, cont = contOutb, var = "MM2_ep")

						distances.append(( distance1+ distance2 )/2) 
						sigma1s_temp.append(sigma1) 
						sigma2s_temp.append(sigma2) 
						sigma3s_temp.append(sigma3) 

			sigma1_opt = sigma1s_temp[np.argmin(distances)]
			sigma2_opt = sigma2s_temp[np.argmin(distances)]
			sigma3_opt = sigma3s_temp[np.argmin(distances)]
			sigmas_opt.append([sigma1_opt, sigma2_opt, sigma3_opt])
			self.SmearingV1(dvcsSimInbCDFT, sigma1_opt, sigma2_opt, sigma3_opt)
			self.saveDVCSvars()
			self.makeDVCS(pol = "inbending")
			dvcsSimInbCDFT_opt = self.df_epg

			self.SmearingV1(dvcsSimOutbCDFT, sigma1_opt, sigma2_opt, sigma3_opt)
			self.saveDVCSvars()
			self.makeDVCS(pol = "outbending")
			dvcsSimOutbCDFT_opt = self.df_epg
			dvcsSimInbCDFT_opt = dvcsSimInbCDFT_opt.loc[(dvcsSimInbCDFT_opt.Pp>PpMin) & (dvcsSimInbCDFT_opt.Pp<PpMax)]
			dvcsSimOutbCDFT_opt = dvcsSimOutbCDFT_opt.loc[(dvcsSimOutbCDFT_opt.Pp>PpMin) & (dvcsSimOutbCDFT_opt.Pp<PpMax)]
			print(len(dvcsSimInbCDFT_opt), len(dvcsSimOutbCDFT_opt))
			print(len(epgExpInbCDFT_selected), len(epgExpOutbCDFT_selected))
			print(PpMin, PpMax, distances, sigma_opt, correction) 

			varstoplot = ["t1", "Ptheta", "Pphi",  "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
			title = [r"$|t|$", r"$\theta_{p}$", r"$\phi_{p}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
			unit = [GeV2, degree, degree, degree, degree, GeV, GeV, GeV2, GeV2, degree]
			# binstarts = [PpMin, 0, -180, 0, 0, 0, -0.5, -0.01, 0.1, 0]
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
			plt.savefig(outDir+"InbCDFT{}_{:.3f}_{:.3f}_{:.3f}pdf".format(i, sigma1_opt, sigma2_opt, sigma3_opt))
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
			plt.savefig(outDir+"OutbCDFT{}_{:.3f}_{:.3f}_{:.3f}pdf".format(i, sigma1_opt, sigma2_opt, sigma3_opt))

		print(sigmas_opt)

	def SmearingV1(self, df, sigma1, sigma2, sigma3):
		df_epg = copy(df)
		df_epg.loc[df_epg["Psector"]>7, "Pp"] = df_epg.loc[df_epg["Psector"]>7, "Pp"]*np.random.normal(1, sigma1, len(df_epg.loc[df_epg.Psector>7]))
		df_epg.loc[df_epg["Psector"]>7, "Ptheta"] = df_epg.loc[df_epg["Psector"]>7, "Ptheta"] + np.random.normal(0, sigma2, len(df_epg.loc[df_epg.Psector>7]))
		df_epg.loc[df_epg["Psector"]>7, "Pphi"] = df_epg.loc[df_epg["Psector"]>7, "Pphi"] + np.random.normal(0, sigma3, len(df_epg.loc[df_epg.Psector>7])) 
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

	def saveDVpi0vars(self):
	    #set up pi0 variables
	    df_epgg = self.df_epgg

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

	    self.df_epgg = df_epgg

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
	        cut_Ptheta_CDFT = df_dvcs.Ptheta<60
	        cut_Gsector_CDFT = df_dvcs.Gsector>7
	        cut_mmep1_CDFT = df_dvcs["MM2_ep"] < 0.601  # mmep
	        cut_mmep2_CDFT = df_dvcs["MM2_ep"] > -0.528  # mmep
	        cut_mmeg1_CDFT = df_dvcs["MM2_eg"] < 1.569  # mmeg
	        cut_mmeg2_CDFT = df_dvcs["MM2_eg"] > 0.246  # mmeg
	        cut_meepg1_CDFT = df_dvcs["ME_epg"] < 0.484 # meepg
	        cut_meepg2_CDFT = df_dvcs["ME_epg"] > -0.439  # meepg
	        cut_cone1_CDFT = df_dvcs["coneAngle"] < 29.652  # coneangle
	        cut_cone2_CDFT = df_dvcs["coneAngle"] > 12.159  # coneangle
	        cut_mpt_CDFT = df_dvcs["MPt"] < 0.0958  # mpt
	        cut_recon_CDFT = df_dvcs["reconGam"] < 0.659  # recon gam angle
	        cut_coplanarity_CDFT = df_dvcs["coplanarity"] < 8.251  # coplanarity angle
	        cut_mmepg1_CDFT = np.abs(df_dvcs["MM2_epg"]) < 0.0122  # mmepg
	        cut_mmepg2_CDFT = np.abs(df_dvcs["MM2_epg"]) > -0.0151  # mmepg

	        cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta_CDFT & cut_Gsector_CDFT &
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
	        cut_cone1_CD = df_dvcs["coneAngle"] < 32.817  # coneangle
	        cut_cone2_CD = df_dvcs["coneAngle"] > 12.791  # coneangle
	        cut_mpt_CD = df_dvcs["MPt"] < 0.182  # mpt
	        cut_recon_CD = df_dvcs["reconGam"] < 0.876  # recon gam angle
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
	        cut_cone1_FD = df_dvcs["coneAngle"] < 35.794  # coneangle
	        cut_cone2_FD = df_dvcs["coneAngle"] > 23.175  # coneangle
	        cut_mpt_FD = df_dvcs["MPt"] < 0.184  # mpt
	        cut_recon_FD = df_dvcs["reconGam"] < 1.181  # recon gam angle
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
	        cut_Ptheta_CDFT = df_dvcs.Ptheta<60
	        cut_Gsector_CDFT = df_dvcs.Gsector>7
	        cut_mmep1_CDFT = df_dvcs["MM2_ep"] < 0.567  # mmep
	        cut_mmep2_CDFT = df_dvcs["MM2_ep"] > -0.479  # mmep
	        cut_mmeg1_CDFT = df_dvcs["MM2_eg"] < 1.602  # mmeg
	        cut_mmeg2_CDFT = df_dvcs["MM2_eg"] > 0.223  # mmeg
	        cut_meepg1_CDFT = df_dvcs["ME_epg"] < 0.481 # meepg
	        cut_meepg2_CDFT = df_dvcs["ME_epg"] > -0.436  # meepg
	        cut_cone1_CDFT = df_dvcs["coneAngle"] < 28.800  # coneangle
	        cut_cone2_CDFT = df_dvcs["coneAngle"] > 12.573  # coneangle
	        cut_mpt_CDFT = df_dvcs["MPt"] < 0.100  # mpt
	        cut_recon_CDFT = df_dvcs["reconGam"] < 0.697  # recon gam angle
	        cut_coplanarity_CDFT = df_dvcs["coplanarity"] < 8.479  # coplanarity angle
	        cut_mmepg1_CDFT = np.abs(df_dvcs["MM2_epg"]) < 0.0128  # mmepg
	        cut_mmepg2_CDFT = np.abs(df_dvcs["MM2_epg"]) > -0.0155  # mmepg

	        cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta_CDFT & cut_Gsector_CDFT &
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
	        cut_cone1_CD = df_dvcs["coneAngle"] < 29.471  # coneangle
	        cut_cone2_CD = df_dvcs["coneAngle"] > 7.591  # coneangle
	        cut_mpt_CD = df_dvcs["MPt"] < 0.183  # mpt
	        cut_recon_CD = df_dvcs["reconGam"] < 0.954  # recon gam angle
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
	        cut_cone1_FD = df_dvcs["coneAngle"] < 34.592  # coneangle
	        cut_cone2_FD = df_dvcs["coneAngle"] > 24.390  # coneangle
	        cut_mpt_FD = df_dvcs["MPt"] < 0.186  # mpt
	        cut_recon_FD = df_dvcs["reconGam"] < 1.662  # recon gam angle
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
	    df_dvpi0p = self.df_epgg

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
	    self.df_epgg = df_dvpi0p #done with saving x



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--version", help="version", default="v0")
    parser.add_argument("-o","--outDir", help="output directory", default="/SimtoDat/v0")
    
    args = parser.parse_args()

    smearingDist = smearingDist(version = args.version, outDir = args.outDir)