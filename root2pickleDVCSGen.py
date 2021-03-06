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
        self.saveDVCSvars()
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
        df_electronGen = pd.DataFrame()
        df_protonGen = pd.DataFrame()
        df_gammaGen = pd.DataFrame()
        eleKeysGen = ["GenEpx", "GenEpy", "GenEpz"]
        proKeysGen = ["GenPpx", "GenPpy", "GenPpz"]
        gamKeysGen = ["GenGpx", "GenGpy", "GenGpz"]
        # read keys
        for key in eleKeysGen:
            df_electronGen[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        for key in proKeysGen:
            df_protonGen[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        for key in gamKeysGen:
            df_gammaGen[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

        #convert data type to standard double
        df_electronGen = df_electronGen.astype({"GenEpx": float, "GenEpy": float, "GenEpz": float})
        df_protonGen = df_protonGen.astype({"GenPpx": float, "GenPpy": float, "GenPpz": float})
        df_gammaGen = df_gammaGen.astype({"GenGpx": float, "GenGpy": float, "GenGpz": float})

        #set up a dummy index for merging
        df_electronGen.loc[:,'event'] = df_electronGen.index
        df_protonGen.loc[:,'event'] = df_protonGen.index
        df_gammaGen.loc[:,'event'] = df_gammaGen.index.get_level_values('entry')

        #sort columns for readability
        df_electronGen = df_electronGen.loc[:, ["event", "GenEpx", "GenEpy", "GenEpz"]]

        #spherical coordinates
        eleGen = [df_electronGen["GenEpx"], df_electronGen["GenEpy"], df_electronGen["GenEpz"]]
        df_electronGen.loc[:, 'GenEp'] = mag(eleGen)
        df_electronGen.loc[:, 'GenEtheta'] = getTheta(eleGen)
        df_electronGen.loc[:, 'GenEphi'] = getPhi(eleGen)

        proGen = [df_protonGen["GenPpx"], df_protonGen["GenPpy"], df_protonGen["GenPpz"]]
        df_protonGen.loc[:, 'GenPp'] = mag(proGen)
        df_protonGen.loc[:, 'GenPtheta'] = getTheta(proGen)
        df_protonGen.loc[:, 'GenPphi'] = getPhi(proGen)

        df_z = pd.merge(df_electronGen, df_protonGen, how='inner', on='event')

        df_gammaGen = df_gammaGen[df_gammaGen.index.get_level_values('subentry')==0]
        gamGen = [df_gammaGen["GenGpx"], df_gammaGen["GenGpy"], df_gammaGen["GenGpz"]]
        df_gammaGen.loc[:, 'GenGp'] = mag(gamGen)
        df_gammaGen.loc[:, 'GenGtheta'] = getTheta(gamGen)
        df_gammaGen.loc[:, 'GenGphi'] = getPhi(gamGen)

        df_z = pd.merge(df_z, df_gammaGen, how='inner', on='event')

        keys = {key:key[3:] for key in df_z.columns}
        keys["event"] = "event"
        self.df_epg = df_z.rename(columns=keys)

    def saveDVCSvars(self, correction=None):
        #set up dvcs variables
        df_epg = self.df_epg

        ele = [df_epg['Epx'], df_epg['Epy'], df_epg['Epz']]
        df_epg.loc[:, 'Ep'] = mag(ele)
        df_epg.loc[:, 'Ee'] = getEnergy(ele, me)
        df_epg.loc[:, 'Etheta'] = getTheta(ele)
        df_epg.loc[:, 'Ephi'] = getPhi(ele)

        pro = [df_epg['Ppx'], df_epg['Ppy'], df_epg['Ppz']]
        df_epg.loc[:, 'Pp'] = mag(pro)
        df_epg.loc[:, 'Pe'] = getEnergy(pro, M)
        df_epg.loc[:, 'Ptheta'] = getTheta(pro)
        df_epg.loc[:, 'Pphi'] = getPhi(pro)

        gam = [df_epg['Gpx'], df_epg['Gpy'], df_epg['Gpz']]
        df_epg.loc[:, 'Gp'] = mag(gam)
        df_epg.loc[:, 'Ge'] = getEnergy(gam, 0)
        df_epg.loc[:, 'Gtheta'] = getTheta(gam)
        df_epg.loc[:, 'Gphi'] = getPhi(gam)

        Ppt = mag([df_epg['Ppx'], df_epg['Ppy'], 0])

        VGS = [-df_epg['Epx'], -df_epg['Epy'], pbeam - df_epg['Epz']]
        v3l = cross(beam, ele)
        v3h = cross(pro, VGS)
        v3g = cross(VGS, gam)
        VmissG = [-df_epg["Epx"] - df_epg["Ppx"], -df_epg["Epy"] - df_epg["Ppy"],
                  pbeam - df_epg["Epz"] - df_epg["Ppz"]]
        VmissP = [-(df_epg["Epx"] + df_epg["Gpx"]), -(df_epg["Epy"] + df_epg["Gpy"]),
                  -(-pbeam + df_epg["Epz"] + df_epg["Gpz"])]
        Vmiss = [-(df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"]), -(df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"]),
                 -(-pbeam + df_epg["Epz"] + df_epg["Ppz"] + df_epg["Gpz"])]
        costheta = cosTheta(VGS, gam)

        df_epg.loc[:, 'Mpx'], df_epg.loc[:, 'Mpy'], df_epg.loc[:, 'Mpz'] = Vmiss

        # binning kinematics
        df_epg.loc[:,'Q2'] = -((ebeam - df_epg['Ee'])**2 - mag2(VGS))
        df_epg.loc[:,'nu'] = (ebeam - df_epg['Ee'])
        df_epg.loc[:,'y'] = df_epg['nu']/ebeam
        df_epg.loc[:,'xB'] = df_epg['Q2'] / 2.0 / M / df_epg['nu']
        df_epg.loc[:,'t1'] = 2 * M * (df_epg['Pe'] - M)
        df_epg.loc[:,'t2'] = (M * df_epg['Q2'] + 2 * M * df_epg['nu'] * (df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta))\
        / (M + df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta)
        df_epg.loc[:,'W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epg['Ee'])**2 - mag2(VGS)))

        # trento angles
        df_epg.loc[:,'phi1'] = angle(v3l, v3h)
        df_epg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_epg['phi1'], df_epg['phi1'])
        df_epg.loc[:,'phi2'] = angle(v3l, v3g)
        df_epg.loc[:,'phi2'] = np.where(dot(v3l, gam) <
                                  0, 360.0 - df_epg['phi2'], df_epg['phi2'])

        # exclusivity variables
        df_epg.loc[:,'MM2_epg'] = (-M - ebeam + df_epg["Ee"] +
                             df_epg["Pe"] + df_epg["Ge"])**2 - mag2(Vmiss)
        df_epg.loc[:,'ME_epg'] = (M + ebeam - df_epg["Ee"] - df_epg["Pe"] - df_epg["Ge"])
        df_epg.loc[:,'MM2_ep'] = (-M - ebeam + df_epg["Ee"] + df_epg["Pe"])**2 - mag2(VmissG)
        df_epg.loc[:,'MM2_eg'] = (-M - ebeam + df_epg["Ee"] + df_epg["Ge"])**2 - mag2(VmissP)
        df_epg.loc[:,'MPt'] = np.sqrt((df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"])**2 +
                                (df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"])**2)
        df_epg.loc[:,'coneAngle'] = angle(ele, gam)
        df_epg.loc[:,'reconGam'] = angle(gam, VmissG)
        df_epg.loc[:,'coplanarity'] = angle(v3h, v3g)
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