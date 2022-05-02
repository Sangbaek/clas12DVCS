import pandas as pd
import numpy as np
import argparse
from utils.const import *
from utils.physics import *


def countDF(total, df_global):
    # numbers1 = []
    # numbers2 = []
    # numbers3 = []
    numbers = []
    for i in df_global.Q2xBtphibin:
        if i%50 == 0:
            print(i)
        for j in df_global.Q2xBtphibin:
            number = sum((total.Q2xBtphibin == i) & (total.GenQ2xBtphibin == j))
            df_global.loc[df_gloab.Q2xBtphibin == i, "Gen{}".format(j)] = number

    return df_global


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-if","--ifname", help="a single pickle to be numbered", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-of","--ofname", help="output with numbered", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-g","--gen", help="gen or rec", action = "store_true")
    parser.add_argument("-i","--inglobal", help="a single pickle file name as an input", default="df_global_Feb.pkl")
    parser.add_argument("-o","--outglobal", help="a single pickle file name as an output", default="df_global_Feb_out.pkl")
    
    args = parser.parse_args()

    print("reading..")
    df = pd.read_pickle(args.ifname)
    print("reading seed..")
    df_global = pd.read_pickle(args.inglobal)
    print("count Rec..")
    df_global = countDF(df, df_global)
    print("done with counting..")
    df_global.to_pickle(args.outglobal)