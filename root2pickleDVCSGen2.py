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

        eleKeysGen = ["xB", "Q2", "t2", "phi2"]
        # read keys
        for key in eleKeysGen:
            df_epg[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

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