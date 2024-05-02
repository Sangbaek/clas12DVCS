#!/usr/bin/env python3
"""
A simple script to save data in pickle.
"""

import uproot
import pandas as pd
import numpy as np
import argparse
from copy import copy
from utils.const import *
from utils.physics import *
import awkward as ak


class root2pickle():
    #class to read root to make epg pairs, inherited from epg
    def __init__(self, fname, entry_stop = None):
        self.fname = fname

        self.readBinScheme()
        self.readEPGG(entry_stop)
        # self.saveHistogram()
        self.saveRaw()

    def readBinScheme(self):
        self.bin_scheme = np.loadtxt('/work/clas12/sangbaek/km15gen/bin_scheme.csv', delimiter = ',')
        self.fringe_bin_scheme = np.loadtxt('/work/clas12/sangbaek/km15gen/fringe_bin_scheme.csv', delimiter = ',')
    def readFile(self):
        #read root using uproot
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        #close file for saving memory
        self.file = None
        self.tree = None

    def readEPGG(self, entry_stop = None):
        #save data into df_epg, df_epgg for parent class epg
        self.readFile()

        # data frames and their keys to read Z part
        df_epg = pd.DataFrame()

        eleKeysGen = ["GenxB", "GenQ2", "Gent", "Genphi", "GenWeight", "BornWeight", "crossRef", "helicity", "radMode", "config", "beamEnergy"]
        # read keys
        for key in eleKeysGen:
            df_epg[key] = ak.to_dataframe(self.tree[key].array(library="ak"))

        # df_epg = df_epg.rename(columns ={"GenxB": "xB", "GenQ2": "Q2", "Gent": "t1", "Genphi": "phi1", "beamEnergy": "beamE"})
        df_epg = df_epg.rename(columns ={"beamEnergy": "beamE"})
        df_epg.loc[:, "event"] = df_epg.index

        df_epg.loc[:, "integrated_binnum_gen"] = 0

        for binnum, bin in enumerate(self.bin_scheme):
            xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin
            try:
                assert np.sum(df_epg.loc[ (df_epg.GenxB>=xBmin) & (df_epg.GenxB<xBmax) & (df_epg.GenQ2>=Q2min) & (df_epg.GenQ2<Q2max)  & (df_epg.Gent>=tmin) & (df_epg.Gent<tmax), "integrated_binnum_gen"] != 0) == 0        
            except:
                print("This bin overlaps with others. Check the geometry. {}".format(bin))
            df_epg.loc[ (df_epg.GenxB>=xBmin) & (df_epg.GenxB<xBmax) & (df_epg.GenQ2>=Q2min) & (df_epg.GenQ2<Q2max)  & (df_epg.Gent>=tmin) & (df_epg.Gent<tmax), "integrated_binnum_gen"] = binnum + 1

        for binnum, bin in enumerate(self.fringe_bin_scheme):
            xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin
            try:
                assert np.sum(df_epg.loc[ (df_epg.GenxB>=xBmin) & (df_epg.GenxB<xBmax) & (df_epg.GenQ2>=Q2min) & (df_epg.GenQ2<Q2max)  & (df_epg.Gent>=tmin) & (df_epg.Gent<tmax), "integrated_binnum_gen"] != 0) == 0        
            except:
                print("This bin overlaps with others. Check the geometry. {}".format(bin))
            df_epg.loc[ (df_epg.GenxB>=xBmin) & (df_epg.GenxB<xBmax) & (df_epg.GenQ2>=Q2min) & (df_epg.GenQ2<Q2max)  & (df_epg.Gent>=tmin) & (df_epg.Gent<tmax), "integrated_binnum_gen"] = binnum + 1 + len(bin_scheme)

        phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
        for phi_binnum in range(24):
            phimin = phibins[phi_binnum]
            phimax = phibins[phi_binnum+1]
            df_epg.loc[ (df_epg.Genphi>=phimin) & (df_epg.Genphi<phimax), "phi_binnum_gen"] = phi_binnum

        df_epg = df_epg.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})


        self.df_epg = df_epg

    def saveRaw(self):
        self.df = self.df_epg


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    
    args = parser.parse_args()

    converter = root2pickle(args.fname, entry_stop = args.entry_stop)
    df = converter.df
    df.to_pickle(args.out)