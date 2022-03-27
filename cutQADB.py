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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single exp pickle file before clasqa", default="/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/exp/dvcs.pkl")
    parser.add_argument("-q","--qadb", help="a qadb query", default="/volatile/clas12/sangbaek/clasqaDB/src/outs/qadb.pkl")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="qadb_applied.pkl")

    args = parser.parse_args()

    qadb = pd.read_pickle(args.qadb)
    qadb = qadb.loc[qadb.QAresult == 1, ["EventNum", "RunNum", "beamCurrent", "runConfig"]]
    exp = pd.read_pickle(args.fname)

    exp = pd.merge([exp, qadb], how = 'inner', on = ['RunNum', 'EventNum'])

    exp.to_pickle(args.out)