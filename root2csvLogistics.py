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

df_logistics1 = pd.DataFrame()
file1 = uproot.open("/volatile/clas12/sangbaek/nov2021/convROOT/exp/inb/dvcs.root")
df_logistics2 = pd.DataFrame()
file2 = uproot.open("/volatile/clas12/sangbaek/nov2021/convROOT/exp/inb/pi0.root")
df_logistics3 = pd.DataFrame()
file3 = uproot.open("/volatile/clas12/sangbaek/nov2021/convROOT/exp/outb/dvcs.root")
df_logistics4 = pd.DataFrame()
file4 = uproot.open("/volatile/clas12/sangbaek/nov2021/convROOT/exp/outb/pi0.root")
tree1 = file1["T"]
tree2 = file2["T"]
tree3 = file3["T"]
tree4 = file4["T"]

for key in ["RunNum", "EventNum", "beamQ"]:
    df_logistics1[key] = tree1[key].array(library="pd")
    df_logistics2[key] = tree2[key].array(library="pd")
    df_logistics3[key] = tree3[key].array(library="pd")
    df_logistics4[key] = tree4[key].array(library="pd")

df = pd.concat([df_logistics1, df_logistics2, df_logistics3, df_logistics4]).drop_duplicates().reset_index(drop=True)
df = df.sort_values(by=['RunNum', 'EventNum'], ascending = [True, True])

for run in df.RunNum.unique():
    df.loc[df.RunNum == run, :].to_csv("/volatile/clas12/sangbaek/clasqaDB/src/csvs/exp{}.csv".format(run))
# df.index = np.linspace(0, len(df)-1, len(df), dtype = int)
# for i in range(len(df)//100000):
#     df[i*100000: (i+1)*100000].to_csv("/volatile/clas12/sangbaek/clasqaDB/src/csvs/exp{}.csv".format(i))
# df[100000*(len(df)//100000)].to_csv("/volatile/clas12/sangbaek/clasqaDB/src/csvs/exp{}.csv".format(i))