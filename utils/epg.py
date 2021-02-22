#!/usr/bin/env python3
"""
Modules set up pandas DataFrame for epg business, DVCS and DVpi0P.
"""

import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from utils.const import *
from utils.physics import *

class epg:
	"""
	Parent class setting up dvcs and dvpi0p
	"""
	def __init(self):
		pass

	def setDVCSvars(self):
		#set up dvcs variables
		df_epg = self.df_epg
		if ('Evz' in df_epg.index):
			df_epg = df_epg[np.abs(df_epg["Evz"] - df_epg["Pvz"]) < 2.5 +
			                2.5 / mag([df_epg["Ppx"], df_epg["Ppy"], df_epg["Ppz"]])]
		df_epg['Ee'] = np.sqrt(me**2 + df_epg["Epx"]**2 +
		                       df_epg["Epy"]**2 + df_epg["Epz"]**2)
		df_epg['Pe'] = np.sqrt(M**2 + df_epg["Ppx"]**2 +
		                       df_epg["Ppy"]**2 + df_epg["Ppz"]**2)
		df_epg['Ge'] = np.sqrt(df_epg["Gpx"]**2 + df_epg["Gpy"]**2 + df_epg["Gpz"]**2)
		ele = [df_epg['Epx'], df_epg['Epy'], df_epg['Epz']]
		pro = [df_epg['Ppx'], df_epg['Ppy'], df_epg['Ppz']]
		gam = [df_epg['Gpx'], df_epg['Gpy'], df_epg['Gpz']]
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

		# binning kinematics
		df_epg['Q2'] = -((ebeam - df_epg['Ee'])**2 - mag2(VGS))
		df_epg['nu'] = (ebeam - df_epg['Ee'])
		df_epg['y'] = df_epg['nu']/ebeam
		df_epg['xB'] = df_epg['Q2'] / 2.0 / M / df_epg['nu']
		df_epg['t1'] = 2 * M * (df_epg['Pe'] - M)
		df_epg['t2'] = (M * df_epg['Q2'] + 2 * M * df_epg['nu'] * (df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta))\
		/ (M + df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta)
		df_epg['W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epg['Ee'])**2 - mag2(VGS)))

		# trento angles
		df_epg['phi1'] = angle(v3l, v3h)
		df_epg['phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
		                          df_epg['phi1'], df_epg['phi1'])
		df_epg['phi2'] = angle(v3l, v3g)
		df_epg['phi2'] = np.where(dot(VGS, cross(v3l, v3g)) <
		                          0, 360.0 - df_epg['phi2'], df_epg['phi2'])

		# exclusivity variables
		df_epg['MM2_epg'] = (-M - ebeam + df_epg["Ee"] +
		                     df_epg["Pe"] + df_epg["Ge"])**2 - mag2(Vmiss)
		df_epg['ME_epg'] = (M + ebeam - df_epg["Ee"] - df_epg["Pe"] - df_epg["Ge"])
		df_epg['MM2_ep'] = (-M - ebeam + df_epg["Ee"] + df_epg["Pe"])**2 - mag2(VmissG)
		df_epg['MM2_eg'] = (-M - ebeam + df_epg["Ee"] + df_epg["Ge"])**2 - mag2(VmissP)
		df_epg['MPt'] = np.sqrt((df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"])**2 +
		                        (df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"])**2)
		df_epg['coneAngle'] = angle(ele, gam)
		df_epg['reconGam'] = angle(gam, VmissG)
		df_epg['coplanarity'] = angle(v3h, v3g)
		self.df_epg = df_epg

	def makeDVCS(self):
		#make dvcs pairs
		df_dvcs = self.df_epg
		df_dvcs = df_dvcs[df_dvcs["MM2_eg"] > 0]  # mmeg

		cut_xBupper = df_dvcs["xB"] < 1  # xB
		cut_xBlower = df_dvcs["xB"] > 0  # xB
		cut_Q2 = df_dvcs["Q2"] > 1  # Q2
		cut_W = df_dvcs["W"] > 2  # W
		cut_Ee = df_dvcs["Ee"] > 2  # Ee
		cut_Ge = df_dvcs["Ge"] > 2  # Ge
		cut_Pp = mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]]) > 0.12  # Pp

		#	Exclusivity cuts
		cut_mmepg = np.abs(df_dvcs["MM2_epg"]) < 0.04  # mmepg
		cut_mmegupper = np.sqrt(df_dvcs["MM2_eg"]) < 1.7  # mmeg
		cut_mmeglower = np.sqrt(df_dvcs["MM2_eg"]) > 0.1  # mmeg
		cut_meepgupper = df_dvcs["ME_epg"] < 1.2  # meepg
		cut_meepglower = df_dvcs["ME_epg"] > -0.5  # meepg
		cut_mpt = df_dvcs["MPt"] < 0.12  # mpt
		cut_cone = df_dvcs["coneAngle"] > 10  # coneangle
		cut_recon = df_dvcs["reconGam"] < 1.1  # recon gam angle
		if ("Esector" in df_dvcs.index):
			cut_sector = df_dvcs["Esector"]!=df_dvcs["Gsector"]
		else:
			cut_sector = 1

		df_dvcs = df_dvcs[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Pp & cut_mmepg &
		                 cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon & cut_sector]
		self.df_dvcs = df_dvcs               
		# #cut by detector status
		# cut_pFD = df_dvcs["Pstat"] < 4000  # FD
		# cut_pCD = df_dvcs["Pstat"] > 4000  # CD
		# cut_gFD = df_dvcs["Gstat"] > 2000  # FD
		# cut_gFT = df_dvcs["Gstat"] < 2000  # FT

	def setDVpi0vars(self):
		#set up pi0 variables
		df_epgg = self.df_epgg
		if ('Evz' in df_epgg.index):
			df_epgg = df_epgg[np.abs(df_epgg["Evz"] - df_epgg["Pvz"]) < 2.5 +
			                  2.5 / mag([df_epgg["Ppx"], df_epgg["Ppy"], df_epgg["Ppz"]])]
		df_epgg['Ee'] = np.sqrt(me**2 + df_epgg["Epx"]**2 + df_epgg["Epy"]**2 + df_epgg["Epz"]**2)
		df_epgg['Pe'] = np.sqrt(M**2 + df_epgg["Ppx"]**2 + df_epgg["Ppy"]**2 + df_epgg["Ppz"]**2)
		df_epgg['Ge'] = np.sqrt(df_epgg["Gpx"]**2 + df_epgg["Gpy"]**2 + df_epgg["Gpz"]**2)
		df_epgg['Ge2'] = np.sqrt(df_epgg["Gpx2"]**2 + df_epgg["Gpy2"]**2 + df_epgg["Gpz2"]**2)
		
		# useful objects
		ele = [df_epgg['Epx'], df_epgg['Epy'], df_epgg['Epz']]
		pro = [df_epgg['Ppx'], df_epgg['Ppy'], df_epgg['Ppz']]
		gam = [df_epgg['Gpx'], df_epgg['Gpy'], df_epgg['Gpz']]
		gam2 = [df_epgg['Gpx2'], df_epgg['Gpy2'], df_epgg['Gpz2']]
		pi0 = vecAdd(gam, gam2)
		VGS = [-df_epgg['Epx'], -df_epgg['Epy'], pbeam - df_epgg['Epz']]
		v3l = cross(beam, ele)
		v3h = cross(pro, VGS)
		v3g = cross(VGS, gam)
		VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
		            df_epgg["Ppy"], pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]

		# binning kinematics
		df_epgg['Q2'] = -((ebeam - df_epgg['Ee'])**2 - mag2(VGS))
		df_epgg['nu'] = (ebeam - df_epgg['Ee'])
		df_epgg['xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
		df_epgg['W'] = np.sqrt((ebeam + M - df_epgg['Ee'])**2 - mag2(VGS))
		df_epgg['MPt'] = np.sqrt((df_epgg["Epx"] + df_epgg["Ppx"] + df_epgg["Gpx"] + df_epgg["Gpx2"])**2 +
		                         (df_epgg["Epy"] + df_epgg["Ppy"] + df_epgg["Gpy"] + df_epgg["Gpy2"])**2)

		# exclusivity variables
		df_epgg['MM2_ep'] = (-M - ebeam + df_epgg["Ee"] +
		                     df_epgg["Pe"])**2 - mag2(VmissPi0)
		df_epgg['ME_epgg'] = (M + ebeam - df_epgg["Ee"] - df_epgg["Pe"] - df_epgg["Ge"] - df_epgg["Ge2"])
		df_epgg['Mpi0'] = pi0InvMass(gam, gam2)
		df_epgg['reconPi'] = angle(VmissPi0, pi0)
		df_epgg["Pie"] = df_epgg['Ge'] + df_epgg['Ge2']
		self.df_epgg = df_epgg

	def makeDVpi0(self):
		#make dvpi0 pairs
		df_epgg = self.df_epgg
		cut_xBupper = df_epgg["xB"] < 1  # xB
		cut_xBlower = df_epgg["xB"] > 0  # xB
		cut_Q2 = df_epgg["Q2"] > 1  # Q2
		cut_W = df_epgg["W"] > 2  # W

		# Exclusivity cuts
		cut_mmep = df_epgg["MM2_ep"] < 0.7  # mmep
		cut_meepgg = df_epgg["ME_epgg"] < 0.7  # meepgg
		cut_mpt = df_epgg["MPt"] < 0.2  # mpt
		cut_recon = df_epgg["reconPi"] < 2  # recon gam angle
		cut_pi0upper = df_epgg["Mpi0"] < 0.2
		cut_pi0lower = df_epgg["Mpi0"] > 0.07
		if ("Esector" in df_epgg.index):
			cut_sector = (df_epgg["Esector"]!=df_epgg["Gsector"]) & (df_epgg["Esector"]!=df_epgg["Gsector2"])
		else:
			cut_sector = 1

		df_dvpi0 = df_epgg[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_mmep & cut_meepgg &
		                   cut_mpt & cut_recon & cut_pi0upper & cut_pi0lower & cut_sector]
		self.df_dvpi0 = df_dvpi0

		# #cut by detector status
		# cut_pFD = df_dvpi0["Pstat"] < 4000  # FD
		# cut_pCD = df_dvpi0["Pstat"] > 4000  # CD
		# cut_gFD = df_dvpi0["Gstat"] > 2000  # FD
		# cut_gFT = df_dvpi0["Gstat"] < 2000  # FT
		# cut_g2FD = df_dvpi0["Gstat2"] > 2000  # FD
		# cut_g2FT = df_dvpi0["Gstat2"] < 2000  # FT
		# cut_FTFT = cut_gFT & cut_g2FT
		# cut_FDFT = (cut_gFD & cut_g2FD) | (cut_gFT & cut_g2FT)
		# cut_FDFD = cut_gFD & cut_g2FD

	def getEPGG(self):
		#returning pd.DataFrame of epgg
		self.setDVpi0vars()
		return self.df_epgg

	def getDVCS(self, sub2g = None):
		#returning pd.DataFrame of epg
		self.setDVCSvars()
		return self.df_epg

	def getDVpi0(self):
		#returning pd.DataFrame of dvpi0
		self.setDVpi0vars()
		self.makeDVpi0()
		return self.df_dvpi0

	def getDVCS(self, sub2g = None):
		#returning pd.DataFrame of dvcs
		self.setDVCSvars()
		self.makeDVCS()
		if(sub2g):
			self.pi02gSubtraction()
		return self.df_dvcs

	def pi02gSubtraction(self):
		#exclude dvpi0 from dvcs. use only when both set up.
		df_dvcs = self.df_dvcs
		pi0to2gammas = df_dvcs["event"].isin(self.df_dvpi0["event"])
		df_dvcs = df_dvcs[~pi0to2gammas]
		self.df_dvcs = df_dvcs

class epgFromROOT(epg):
	#class to read root to make epg pairs, inherited from epg
	def __init__(self, fname, entry_stop = None, mc = False):
		self.fname = fname
		self.mc = mc
		self.readEPG(entry_stop)

	def readFile(self):
		#read root using uproot
		self.file = uproot.open(self.fname)
		self.tree = self.file["T"]

	def closeFile(self):
		#close file for saving memory
		self.file = None
		self.tree = None

	def readEPG(self, entry_stop = None):
		#save data into df_epg, df_epgg for parent class epg
		self.readFile()

		df_electron = pd.DataFrame()
		df_proton = pd.DataFrame()
		df_gamma = pd.DataFrame()

		eleKeys = ["Epx", "Epy", "Epz", "Evz", "Esector"]
		proKeys = ["Ppx", "Ppy", "Ppz", "Pvz", "PorigIndx", "Pstat", "Psector"]
		gamKeys = ["Gpx", "Gpy", "Gpz", "GorigIndx", "Gstat", "Gsector"]
		if (self.mc):
			eleKeys = eleKeys[:-1]
			proKeys = proKeys[:-2]
			gamKeys = gamKeys[:-2]

		for key in eleKeys:
			df_electron[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

		for key in proKeys:
			df_proton[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

		for key in gamKeys:
			df_gamma[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

		self.closeFile()

		df_electron = df_electron.astype({"Epx": float, "Epy": float, "Epz": float})
		df_proton = df_proton.astype({"Ppx": float, "Ppy": float, "Ppz": float})
		df_gamma = df_gamma.astype({"Gpx": float, "Gpy": float, "Gpz": float})

		df_electron['event'] = df_electron.index.get_level_values('entry')
		df_proton['event'] = df_proton.index.get_level_values('entry')
		df_gamma['event'] = df_gamma.index.get_level_values('entry')

		df_gg = pd.merge(df_gamma, df_gamma,
		                 how='outer', on='event', suffixes=("", "2"))
		df_gg = df_gg[df_gg["GorigIndx"] < df_gg["GorigIndx2"]]
		df_ep = pd.merge(df_electron, df_proton, how='outer', on='event')
		df_epg = pd.merge(df_ep, df_gamma, how='outer', on='event')
		df_epgg = pd.merge(df_ep, df_gg, how='outer', on='event')
		df_epgg = df_epgg[~np.isnan(df_epgg["Gpx"])]

		self.df_epg = df_epg
		self.df_epgg = df_epgg

class epgFromLund(epg):
	#class to read lund format to make epg pairs, inherited from epg
	def __init__(self, fname, entry_stop = None):
		self.fname = fname
		self.readEPG(entry_stop)

	def readFile(self):
		#read tsv file using python built-in functions
		with open(self.fname,"r") as file:
		    self.data = file.read()
	
	def closeFile(self):
		#close file for saving memory
		self.data = None

	def readEPG(self, entry_stop = None):
		#save data into df_epg, df_epgg for parent class epg
		self.readFile()
		partArray = []

		txtlst = self.data.split("\n")
		for ind, line in enumerate(txtlst):
			if entry_stop:
				if ind==entry_stop:
					break
			if ind %400000 == 0:
				print("On event {}".format(ind/4))
			if ind % 4 == 0:
				header = line
				eleLine = txtlst[ind+1]
				eleQuantities = eleLine.split()
				Epx = eleQuantities[6]
				Epy = eleQuantities[7]
				Epz = eleQuantities[8]
				proLine = txtlst[ind+2]
				proQuantities = proLine.split()
				Ppx = proQuantities[6]
				Ppy = proQuantities[7]
				Ppz = proQuantities[8]
				gamLine = txtlst[ind+3]
				gamQuantities = gamLine.split()
				Gpx = gamQuantities[6]
				Gpy = gamQuantities[7]
				Gpz = gamQuantities[8]
				partArray.append([float(Epx), float(Epy), float(Epz), float(Ppx), float(Ppy), float(Ppz), float(Gpx), float(Gpy), float(Gpz)])

		self.df_epg = pd.DataFrame(partArray, columns = ["Epx", "Epy", "Epz", "Ppx", "Ppy", "Ppz", "Gpx", "Gpy", "Gpz"])
		self.df_epgg = None
		self.closeFile()