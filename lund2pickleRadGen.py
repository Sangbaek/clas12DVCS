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
		txts = iter(txtlst)
		skip = 0
		for ind, line in enumerate(txtlst[:-1]):
			num_particles = line[11]
			if skip > 0:
				skip = skip - 1
				continue
			if num_particles == "3":
				eleLine = txtlst[ind+1]
				proLine = txtlst[ind+2]
				gamLine = txtlst[ind+3]
				skip = 3 
			elif num_particles == "4":
				eleLine = txtlst[ind+1]
				proLine = txtlst[ind+2]
				gamLine = txtlst[ind+3]
				radLine = txtlst[ind+4]
				skip = 4
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
			if num_particles == "3":
				Gpx2 = 0
				Gpy2 = 0
				Gpz2 = 0
			elif num_particles == "4":
				radQuantities = radLine.split()
				Gpx2 = radQuantities[6]
				Gpy2 = radQuantities[7]
				Gpz2 = radQuantities[8]

			partArray.append([float(Epx), float(Epy), float(Epz), float(Ppx), float(Ppy), float(Ppz), float(Gpx), float(Gpy), float(Gpz), float(Gpx2), float(Gpy2), float(Gpz2)])

		self.df = pd.DataFrame(partArray, columns = ["Epx", "Epy", "Epz", "Ppx", "Ppy", "Ppz", "Gpx", "Gpy", "Gpz", "Gpx2", "Gpy2", "Gpz2"])
		self.closeFile()

if __name__ == "__main__":
	print("run lund2pickle for radiative DVCS generator")

	parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
	parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
	parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)

	args = parser.parse_args()

	converter = lund2pickle(args.fname, entry_stop = args.entry_stop)
	df = converter.df

	df.to_pickle(args.out)