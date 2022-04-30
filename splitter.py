import pandas as pd
import numpy as np
import argparse
from utils.const import *
from utils.physics import *
import subprocess, os



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="file to be split", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-b","--binning", help="Q2, Q2xB, Q2xBt, Q2xBtphi", default="Q2xBt")
    
    args = parser.parse_args()

    fname = args.fname
    df = pd.read_pickle(args.fname)


    dirname = fname.rsplit('/', 1)[0] + "/{}/".format(args.binning)
    jobn = fname.rsplit('/', 1)[1][:-4]

    os.mkdir(dirname)

    if args.binning == "Q2xBt":
        for Q2xBtbin in np.sort((df.Q2xBtbin).unique()):
            if Q2xBtbin == -1:
                continue
            df_sub = df.loc[df.Q2xBtbin == Q2xBtbin, :]
            ofname = dirname + jobn+"_{}.pkl".format(Q2xBtbin)
            df_sub.to_pickle(ofname)

    # df.loc[:, "BHrad"] = XsecObs
    # df.loc[:, "BHborn"] = XsecBorn
    # df.to_pickle(args.ofname)