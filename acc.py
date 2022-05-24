import pandas as pd
import numpy as np
import argparse
from utils.const import *
from utils.physics import *


def countDF(total, df_global, colName = "new"):

    hist, _ = np.histogram(total.Q2xBtphibin, bins = np.linspace(-0.5, 4199.5, 4201))
    hist = hist[activebins]
    hist1, _ = np.histogram(total.loc[total.config==1, "Q2xBtphibin"], bins = np.linspace(-0.5, 4199.5, 4201))
    hist1 = hist1[activebins]
    hist2, _ = np.histogram(total.loc[total.config==2, "Q2xBtphibin"], bins = np.linspace(-0.5, 4199.5, 4201))
    hist2 = hist2[activebins]
    hist3, _ = np.histogram(total.loc[total.config==3, "Q2xBtphibin"], bins = np.linspace(-0.5, 4199.5, 4201))
    hist3 = hist3[activebins]

    df_global.loc[:, colName+"1"] = hist1
    df_global.loc[:, colName+"2"] = hist2
    df_global.loc[:, colName+"3"] = hist3
    df_global.loc[:, colName] = hist
    return df_global

# def countGenDF(total, df_global, colName = "new"):
#     numbers = []
#     for i in df_global.Q2xBtphibin:
#         if i%10==0:
#             print(i)
#         number = sum((total.Q2xBtphibin == i))
#         numbers.append(number)
#     df_global.loc[:, colName] = numbers
#     return df_global

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-if","--ifname", help="a single pickle to be numbered", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-of","--ofname", help="output with numbered", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-g","--gen", help="gen or rec", action = "store_true")
    parser.add_argument("-i","--inglobal", help="a single pickle file name as an input", default="df_global_Feb.pkl")
    parser.add_argument("-o","--outglobal", help="a single pickle file name as an output", default="df_global_Feb_out.pkl")
    parser.add_argument("-N","--nonrad", help="to make fractional", action="store_true")
    parser.add_argument("-S","--speak", help="to make fractional", action="store_true")
    parser.add_argument("-P","--ppeak", help="to make fractional", action="store_true")
    parser.add_argument("-c","--colName", help="columnName", default="dvcsSimInb50nA")
    
    args = parser.parse_args()

    print("reading..")
    df = pd.read_pickle(args.ifname)
    print("reading seed..")
    df_global = pd.read_pickle(args.inglobal)
    # if args.gen:
    #     print("count Gen..")
    #     if args.nonrad:
    #         print("select non rad events only")
    #         df = df.loc[df.radMode == 1, :]
    #     if args.speak:
    #         print("select s-peak events only")
    #         df = df.loc[df.radMode == 2, :]
    #     if args.ppeak:
    #         print("select p-peak events only")
    #         df = df.loc[df.radMode == 3, :]
    #     df_global = countGenDF(df, df_global, args.colName)
    #     print("done with counting..")
    #     df_global.to_pickle(args.outglobal)
    # else:
    #     if args.nonrad:
    #         print("select non rad events only")
    #         df = df.loc[df.radMode == 1, :]
    #     if args.speak:
    #         print("select s-peak events only")
    #         df = df.loc[df.radMode == 2, :]
    #     if args.ppeak:
    #         print("select p-peak events only")
    #         df = df.loc[df.radMode == 3, :]
    #     print("count Rec..")
    df_global = countDF(df, df_global, args.colName)
    print("done with counting..")
    df_global.to_pickle(args.outglobal)