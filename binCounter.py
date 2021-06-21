#!/usr/bin/env python3
"""
An independent script that takes pandas dataframe
(e.g.) df_epg) to create a df_xBQ2t, a binned data that has 8 other columns

lower_xB, upper_xB, lower_Q2, upper_Q2, -lower_t, -upper_t, lower_phi, upper_phi

should useful for cross section study, and especially for counting N_gen.
Often, the entire lund file doesn't need to be investigated.
Each count in each bin is only in interest for acceptance study, or in tighter binning for bin volume effect study.
"""
import pandas as pd
import numpy as np
import argparse
from copy import copy

class binCounter():
	def __init__(self, df):
		self.df = df
		binedges_xB = [0.29, 0.32, 0.35, 0.38] # start with small range of xB to ignore finite bin volume effects.
		binedges_Q2 = [2.2, 2.7, 3.2] # start with small range of xB to ignore finite bin volume effects.
		binedges_t = [0.23, 0.3, 0.39, 0.52]
		# binedges_phi = np.linspace(0, 360, 21) # 21 = 360/30 + 1

		for ind_xbq2 in range(4):
			for ind_t in range(3):
				ind_bin = ind_xbq2*3 + ind_t
				if ind_xbq2 == 0:
					lower_xB = 0.29
					upper_xB = 0.32
					lower_Q2 = 2.2
					upper_Q2 = 2.7
				if ind_xbq2 == 1:
					lower_xB = 0.32
					upper_xB = 0.35
					lower_Q2 = 2.2
					upper_Q2 = 2.7
				if ind_xbq2 == 2:
					lower_xB = 0.35
					upper_xB = 0.38
					lower_Q2 = 2.2
					upper_Q2 = 2.7
				if ind_xbq2 == 3:
					lower_xB = 0.35
					upper_xB = 0.38
					lower_Q2 = 2.7
					upper_Q2 = 3.2
				lower_t = binedges_t[ind_t]
				upper_t = binedges_t[ind_t+1]

				print(ind_bin, self.dfinOneBin(lower_xB, upper_xB, lower_Q2, upper_Q2, lower_t, upper_t)[0])

		# dict_xB = {'lower_xB': np.array(binedges_xB[:-1]), 'upper_xB': np.array(binedges_xB[1:])}
		# dict_Q2 = {'lower_Q2': np.array(binedges_Q2[:-1]), 'upper_Q2': np.array(binedges_Q2[1:])}
		# dict_t = {'lower_t': np.array(binedges_t[:-1]), 'upper_t': np.array(binedges_t[1:])}
		# dict_phi = {'lower_phi': np.array(binedges_phi[:-1]), 'upper_phi': np.array(binedges_phi[1:])}

		# df_xB = pd.DataFrame(data = dict_xB, index = np.zeros(len(binedges_xB)-1))
		# df_Q2 = pd.DataFrame(data = dict_Q2, index = np.zeros(len(binedges_Q2)-1))
		# df_t  = pd.DataFrame(data = dict_t,  index = np.zeros(len(binedges_t)-1))
		# df_phi = pd.DataFrame(data = dict_phi, index = np.zeros(len(binedges_phi)-1))

		# df_xBQ2 = pd.merge(df_xB, df_Q2, left_index = True, right_index = True)
		# df_xBQ2t = pd.merge(df_xBQ2, df_t, left_index = True, right_index = True)
		# df_xBQ2tphi = pd.merge(df_xBQ2t, df_phi, left_index = True, right_index = True)

		# df_xBQ2tphi.index = np.linspace(0, len(df_xBQ2tphi)-1, len(df_xBQ2tphi), dtype = int)
		# for i in range(len(df_xBQ2tphi)):
		# 	df_xBQ2tphi.loc[i, "Count"] = self.countinOneBin(df_xBQ2tphi.loc[i, "lower_xB"], df_xBQ2tphi.loc[i, "upper_xB"], df_xBQ2tphi.loc[i, "lower_Q2"], df_xBQ2tphi.loc[i, "upper_Q2"], 
		# 										df_xBQ2tphi.loc[i, "lower_t"], df_xBQ2tphi.loc[i, "upper_t"], df_xBQ2tphi.loc[i, "lower_phi"], df_xBQ2tphi.loc[i, "upper_phi"])
		
		# print(np.sum(df_xBQ2tphi.Count))
		# self.df_xBQ2tphi = df_xBQ2tphi

	def dfinOneBin(self, lower_xB, upper_xB, lower_Q2, upper_Q2, lower_t, upper_t):

		df = copy(self.df)

		cond_lower_Q2 = df["Q2"] > lower_Q2
		cond_upper_Q2 = df["Q2"] < upper_Q2
		cond_lower_xB = df["xB"] > lower_xB
		cond_upper_xB = df["xB"] < upper_xB
		if "t2" in df.columns:
			cond_lower_t = df["t2"] > lower_t
			cond_upper_t = df["t2"] < upper_t
		else:	
			cond_lower_t = df["t"] > lower_t
			cond_upper_t = df["t"] < upper_t
		if "phi2" in df.columns:
			phi = df[cond_lower_xB & cond_upper_xB & cond_lower_Q2 & cond_upper_Q2 & cond_lower_t & cond_upper_t]
			# cond_lower_phi = df["phi2"] > lower_phi
			# cond_upper_phi = df["phi2"] < upper_phi
		else:
			phi = df[cond_lower_xB & cond_upper_xB & cond_lower_Q2 & cond_upper_Q2 & cond_lower_t & cond_upper_t]
		# 	cond_lower_phi = df["phi"] > lower_phi
		# 	cond_upper_phi = df["phi"] < upper_phi


		return np.histogram(phi, bins = np.linspace(0, 360, 21))



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a pickle file contains xB, Q2, t, phi", default="pi0.pkl")
    parser.add_argument("-o","--out", help="a single pickle file that only contains xBQ2tphi binned data", default="xBQ2tphi.pkl")
    parser.add_argument("-d","--detector", help="detector components/FDFD/CDFD/CDFT", default= "all")
    
    args = parser.parse_args()

    df = pd.read_pickle(args.fname)
    if args.detector == "all":
    	pass
    elif args.detector == "FD":
    	df = df[(df.Psector<7)]

    elif args.detector == "CD":
    	df = df[(df.Psector>6)]

    # elif args.detector == "CDFT":
    # 	df = df[(df.Psector<7) & (df.Gsector>7)]

    # elif args.detector == "FDFT":
    # 	df = df[(df.Psector<7) & (df.Gsector>7)]

    # else:
    # 	print("check the desired detector component.")
    # 	exit()
    binCounter(df)
    # df_xBQ2tphi = binCounter(df).df_xBQ2tphi
    # df_xBQ2tphi.to_pickle(args.out)