import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from lmfit import Model, Parameter, report_fit
from scipy.optimize import curve_fit
from copy import copy
from utils.physics import *
import icecream as ic

M = 0.938272081
me = 0.5109989461 * 0.001
ebeam = 10.604
pbeam = np.sqrt(ebeam * ebeam - me * me)
beam = [0, 0, pbeam]
target = [0, 0, 0]

class epgFromData:
	def __init__(self, fname, entry_stop = None):
		self.fname = fname
		self.readEPG(entry_stop)

	def readFile(self):
		self.file = uproot.open(self.fname)
		self.tree = self.file["T"]

	def closeFile(self):
		self.file = None
		self.tree = None

	def readEPG(self, entry_stop = None):

		self.readFile()

		df_electron = pd.DataFrame()
		df_proton = pd.DataFrame()
		df_gamma = pd.DataFrame()

		for key in ["Epx", "Epy", "Epz", "Evz", "Esector"]:
			df_electron[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

		for key in ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "PorigIndx", "Psector"]:
			df_proton[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

		for key in ["Gpx", "Gpy", "Gpz", "Gstat", "GorigIndx", "Gsector"]:
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

		df_epg = df_epg[np.abs(df_epg["Evz"] - df_epg["Pvz"]) < 2.5 +
		                2.5 / mag([df_epg["Epx"], df_epg["Epy"], df_epg["Epz"]])]
		df_epgg = df_epgg[np.abs(df_epgg["Evz"] - df_epgg["Pvz"]) < 2.5 +
		                  2.5 / mag([df_epgg["Epx"], df_epgg["Epy"], df_epgg["Epz"]])]

		self.df_epg = df_epg
		self.df_epgg = df_epgg


	def setDVCSvars(self):
		#useful objects
		df_epg = self.df_epg
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
		df_epg['W'] = np.sqrt((ebeam + M - df_epg['Ee'])**2 - mag2(VGS))

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
		df_epg = df_epg[df_epg["MM2_eg"] > 0]  # mmeg
		df_epg['MPt'] = np.sqrt((df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"])**2 +
		                        (df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"])**2)
		df_epg['coneAngle'] = angle(ele, gam)
		df_epg['reconGam'] = angle(gam, VmissG)
		df_epg['coplanarity'] = angle(v3h, v3g)
		self.df_epg = df_epg

	def makeDVCS(self):
		df_dvcs = self.df_epg
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
		cut_sector = df_dvcs["Esector"]!=df_dvcs["Gsector"]

		df_dvcs = df_dvcs[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Pp & cut_mmepg &
		                 cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon & cut_sector]
		self.df_dvcs = df_dvcs               
		# #cut by detector status
		# cut_pFD = df_dvcs["Pstat"] < 4000  # FD
		# cut_pCD = df_dvcs["Pstat"] > 4000  # CD
		# cut_gFD = df_dvcs["Gstat"] > 2000  # FD
		# cut_gFT = df_dvcs["Gstat"] < 2000  # FT

	def setDVpi0vars(self):
		df_epgg = self.df_epgg
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
		cut_sector = (df_epgg["Esector"]!=df_epgg["Gsector"]) & (df_epgg["Esector"]!=df_epgg["Gsector2"])
		
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

	def getDVpi0(self):
		self.setDVpi0vars()
		self.makeDVpi0()
		return self.df_dvpi0

	def getDVCS(self):
		self.setDVCSvars()
		self.makeDVCS()
		return self.df_dvcs

	def pi02gSubtraction(self):
		df_dvcs = self.df_dvcs
		pi0to2gammas = df_dvcs["event"].isin(self.df_dvpi0["event"])
		df_dvcs = df_dvcs[~pi0to2gammas]
		self.df_dvcs = df_dvcs
		return self.df_dvcs