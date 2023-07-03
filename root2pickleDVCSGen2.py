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

        self.readEPGG(entry_stop)
        # self.saveHistogram()
        self.saveRaw()

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

        eleKeysGen = ["GenxB", "GenQ2", "Gent", "Genphi", "GenWeight", "BornWeight", "crossRef", "helicity", "radMode", "config"]
        # read keys
        for key in eleKeysGen:
            df_epg[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))

        df_epg = df_epg.rename(columns ={"GenxB": "xB", "GenQ2": "Q2", "Gent": "t1", "Genphi": "phi1"})
        df_epg.loc[:, "event"] = df_epg.index

        df_epg.loc[df_epg.phi1<0, "phi1"] = 10**(-4)
        df_epg.loc[df_epg.phi1>=360, "phi1"] = 360-10**(-4)

        df_epg.loc[:, 'xBbin'] = np.zeros(len(df_epg.xB), dtype = 'int') - 1
        df_epg.loc[:, 'Q2bin'] = np.zeros(len(df_epg.Q2), dtype = 'int') - 1
        df_epg.loc[:, 'tbin'] = np.zeros(len(df_epg.t1), dtype = 'int') - 1
        df_epg.loc[:, 'phibin'] = np.zeros(len(df_epg.phi1), dtype = 'int') - 1
        for xB in newxBbins2:
            df_epg.xBbin = df_epg.xBbin + (df_epg.xB>xB).astype("int").to_numpy(dtype = 'int')

        for Q2 in newQ2bins2:
            df_epg.Q2bin = df_epg.Q2bin + (df_epg.Q2>Q2).astype("int").to_numpy(dtype = 'int')

        for t1 in newtbins:
            df_epg.tbin = df_epg.tbin + (df_epg.t1>t1).astype("int").to_numpy(dtype = 'int')

        for phi1 in phibins:
            df_epg.phibin = df_epg.phibin + (df_epg.phi1>phi1).astype("int").to_numpy(dtype = 'int')

        # # encode unassigned bin as -1
        # df_epg.loc[:, "Q2bin"] = -1
        # df_epg.loc[:, "xBbin"] = -1
        # df_epg.loc[:, "tbin"] = -1
        # # df_epg.loc[:, "tbin2"] = -1
        # df_epg.loc[:, "phibin"] = -1
        # # df_epg.loc[:, "phibin2"] = -1
        # df_epg.loc[:, "Q2xBbin"] = -1
        # df_epg.loc[:, "Q2xBtbin"] = -1
        # # df_epg.loc[:, "Q2xBtbin2"] = -1
        # df_epg.loc[:, "Q2xBtphibin"] = -1
        # Q2xBbin = 0

        # # encode all binning
        # for Q2bin in range(len(Q2bin_i)):
        #     #square Q2 binning
        #     df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]), "Q2bin"] = Q2bin
        #     #adaptive xB binning
        #     for xBbin in range(len(xBbin_i[Q2bin])):
        #         if Q2bin < len(Q2bin_i) -1:
        #             if xBbin == 0:
        #                 df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin #0
        #                 df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin #0
        #             elif xBbin < len(xBbin_i[Q2bin])-1:
        #                 df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin
        #                 df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin
        #             else:
        #                 df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "xBbin"] = xBbin
        #                 df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "Q2xBbin"] = Q2xBbin
        #         else:
        #             df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB)& (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "xBbin"] = xBbin
        #             df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB)& (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "Q2xBbin"] = Q2xBbin #0

        #         Q2xBbin = Q2xBbin + 1
        # for tbin in range(len(tbin_i)):
        #     #square t binning
        #     df_epg.loc[(df_epg.t1>=tbin_i[tbin]) & (df_epg.t1<tbin_f[tbin]), "tbin"] = tbin
        #     # df_epg.loc[(df_epg.t2>=tbin_i[tbin]) & (df_epg.t2<tbin_f[tbin]), "tbin2"] = tbin
        # for phibin in range(len(phibin_i)):
        #     #square phi binning
        #     df_epg.loc[(df_epg.phi1>=phibin_i[phibin]) & (df_epg.phi1<phibin_f[phibin]), "phibin"] = phibin
        #     # df_epg.loc[(df_epg.phi2>=phibin_i[phibin]) & (df_epg.phi2<phibin_f[phibin]), "phibin2"] = phibin

        # df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBtbin"] = len(tbin_i) * df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBbin"] + df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "tbin"]
        # # df_epg.loc[(df_epg.Q2bin>0)&(df_epg.xBbin>0)&(df_epg.tbin2>0), "Q2xBtbin2"] = df_epg.Q2bin.astype(str) + df_epg.xBbin.astype(str) + df_epg.tbin2.astype(str)
        # df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBtphibin"] = len(phibin_i) * df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBtbin"] + df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "phibin"]

        # df_epg = df_epg.astype({"Q2bin": int, "xBbin": int, "tbin": int, "phibin": int, "Q2xBbin": int, "Q2xBtbin": int, "Q2xBtphibin": int, "helicity": int})


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