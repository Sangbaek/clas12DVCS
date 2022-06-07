#!/usr/bin/env python3
"""
A simple script to save data in pickle.
"""

import pandas as pd
import numpy as np
import argparse
from copy import copy
from utils.const import *
from utils.physics import *

class lund2pickle():
	#class to read lund format to make epg pairs, inherited from epg
	def __init__(self, fname, entry_stop = None, Q2xBtbin = None):
		self.fname = fname
		self.readEPG(entry_stop, Q2xBtbin)

	def readFile(self):
		#read tsv file using python built-in functions
		with open(self.fname,"r") as file:
		    self.data = file.read()

	def closeFile(self):
		#close file for saving memory
		self.data = None

	def readEPG(self, entry_stop = None, Q2xBtbin = None):
		#save data into df_epg, df_epgg for parent class epg
		self.readFile()
		kinArray = []

		txtlst = self.data.split("\n")
		txts = iter(txtlst)
		skip = 0
		for ind, line in enumerate(txtlst[:-1]):
			num_particles = line[11]
			if skip > 0:
				skip = skip - 1
				continue
			if num_particles == "3":
				logLine = txtlst[ind]
				eleLine = txtlst[ind+1]
				proLine = txtlst[ind+2]
				gamLine = txtlst[ind+3]
				skip = 3 
			elif num_particles == "4":
				logLine = txtlst[ind]
				eleLine = txtlst[ind+1]
				proLine = txtlst[ind+2]
				gamLine = txtlst[ind+3]
				radLine = txtlst[ind+4]
				skip = 4
			logQuantities = logLine.split()
			eleQuantities = eleLine.split()
			proQuantities = proLine.split()
			gamQuantities = gamLine.split()
			dsigma = logQuantities[-1]
			helicity = logQuantities[4]
			xB = eleQuantities[1]
			radMode = eleQuantities[5]
			Epx = eleQuantities[6]
			Epy = eleQuantities[7]
			Epz = eleQuantities[8]
			Q2 = eleQuantities[9]
			t1 = eleQuantities[10]
			phi1 = proQuantities[1]
			Ppx = proQuantities[6]
			Ppy = proQuantities[7]
			Ppz = proQuantities[8]
			Gpx = gamQuantities[6]
			Gpy = gamQuantities[7]
			Gpz = gamQuantities[8]
			if num_particles == "3":
				Gpx2 = 0
				Gpy2 = 0
				Gpz2 = 0
			elif num_particles == "4":
				radQuantities = radLine.split()
				Gpx2 = radQuantities[6]
				Gpy2 = radQuantities[7]
				Gpz2 = radQuantities[8]

			Ee = np.sqrt(Epx**2 + Epy**2 + Epz**2 + me**2)
			Pe = np.sqrt(Ppx**2 + Ppy**2 + Ppz**2 + M**2)
			Ge = np.sqrt(Gpx**2 + Gpy**2 + Gpz**2)
			Ge2 = np.sqrt(Gpx2**2 + Gpy2**2 + Gpz2**2)

			kinArray.append([float(Epx), float(Epy), float(Epz), float(Ee), float(Ppx), float(Ppy), float(Ppz), float(Pe),
			 float(Gpx), float(Gpy), float(Gpz), float(Ge), float(Gpx2), float(Gpy2), float(Gpz2),float(Ge2), 
			 float(xB), float(Q2), float(t1), np.degrees(float(phi1)), float(dsigma), int(radMode), int(float(helicity))])

		df_epgg = pd.DataFrame(kinArray, columns = ["Epx", "Epy", "Epz", "Ee", "Ppx", "Ppy", "Ppz", "Pe", "Gpx", "Gpy", "Gpz", "Ge", "Gpx2", "Gpy2", "Gpz2", "Ge2",
			"xB", "Q2", "t1", "phi1", "dsigma", "radMode", "helicity"])

		if Q2xBtbin:
			# encode unassigned bin as -1
			df_epgg.loc[:, "Q2bin"] = -1
			df_epgg.loc[:, "xBbin"] = -1
			df_epgg.loc[:, "tbin"] = -1
			# df_epgg.loc[:, "tbin2"] = -1
			df_epgg.loc[:, "phibin"] = -1
			# df_epgg.loc[:, "phibin2"] = -1
			df_epgg.loc[:, "Q2xBtbin"] = -1
			df_epgg.loc[:, "Q2xBtbin"] = -1
			# df_epgg.loc[:, "Q2xBtbin2"] = -1
			df_epgg.loc[:, "Q2xBtphibin"] = -1
			Q2xBbin = 0

			# encode all binning
			for Q2bin in range(len(Q2bin_i)):
				#square Q2 binning
				df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]), "Q2bin"] = Q2bin
				#adaptive xB binning
				for xBbin in range(len(xBbin_i[Q2bin])):
					if Q2bin < len(Q2bin_i) -1:
						if xBbin == 0:
							df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin #0
							df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin #0
						elif xBbin < len(xBbin_i[Q2bin])-1:
							df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin
							df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin
						else:
							df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
							df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin
					else:
						df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
						df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin #0

			Q2xBbin = Q2xBbin + 1
			for tbin in range(len(tbin_i)):
				#square t binning
				df_epgg.loc[(df_epgg.t1>=tbin_i[tbin]) & (df_epgg.t1<tbin_f[tbin]), "tbin"] = tbin
				# df_epgg.loc[(df_epgg.t2>=tbin_i[tbin]) & (df_epgg.t2<tbin_f[tbin]), "tbin2"] = tbin
			for phibin in range(len(phibin_i)):
				#square phi binning
				df_epgg.loc[(df_epgg.phi1>=phibin_i[phibin]) & (df_epgg.phi1<phibin_f[phibin]), "phibin"] = phibin
				# df_epgg.loc[(df_epgg.phi2>=phibin_i[phibin]) & (df_epgg.phi2<phibin_f[phibin]), "phibin2"] = phibin

			df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtbin"] = len(tbin_i) * df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBbin"] + df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "tbin"]
			# df_epgg.loc[(df_epgg.Q2bin>0)&(df_epgg.xBbin>0)&(df_epgg.tbin2>0), "Q2xBtbin2"] = df_epgg.Q2bin.astype(str) + df_epgg.xBbin.astype(str) + df_epgg.tbin2.astype(str)
			df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtphibin"] = len(phibin_i) * df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtbin"] + df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "phibin"]

			df_epgg = df_epgg.astype({"Q2bin": int, "xBbin": int, "tbin": int, "phibin": int, "Q2xBbin": int, "Q2xBtbin": int, "Q2xBtphibin": int})

			df_epgg = df_epgg.loc[df_epgg.Q2xBtbin == Q2xBtbin, :]

		self.df = df_epgg
		self.closeFile()

if __name__ == "__main__":
	print("run lund2pickle for radiative DVCS generator")

	parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
	parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
	parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
	parser.add_argument("-Q2xBtbin","--Q2xBtbin", help="Q2xBtbin, None for full", default = None)

	args = parser.parse_args()

	if args.Q2xBtbin:
		converter = lund2pickle(args.fname, entry_stop = args.entry_stop, Q2xBtbin = int(args.Q2xBtbin))
	else:
		converter = lund2pickle(args.fname, entry_stop = args.entry_stop)
	df = converter.df

	df.to_pickle(args.out)