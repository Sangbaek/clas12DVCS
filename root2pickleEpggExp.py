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
    def __init__(self, fname, entry_stop = None, pol = "inbending"):
        self.fname = fname
        self.readEPGG(entry_stop, pol = pol)
        self.saveDVpi0vars()
        self.makeDVpi0()
        self.saveRaw()

    def readFile(self):
        #read root using uproot
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        #close file for saving memory
        self.file = None
        self.tree = None

    def readEPGG(self, entry_stop = None, pol = "inbending"):
        #save data into df_epg, df_epgg for parent class epg
        self.readFile()
        print(entry_stop)
        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Evx", "Evy", "Evz", "Esector", "RunNum", "beamQ", "liveTime", "helicity"]
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Psector"]
        gamKeysRec = ["Gpx", "Gpy", "Gpz", "Gsector"]
        # read them
        for key in eleKeysRec:
            df_electronRec[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        for key in proKeysRec:
            df_protonRec[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        for key in gamKeysRec:
            df_gammaRec[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float, "Evx": float, "Evy": float, "Evz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float, "Pvz": float})
        df_gammaRec = df_gammaRec.astype({"Gpx": float, "Gpy": float, "Gpz": float})

        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'event'] = df_gammaRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'GIndex'] = df_gammaRec.index.get_level_values('subentry')

        #save only FD protons and photons
        # df_protonRec = df_protonRec[df_protonRec["Psector"]<7]
        #proton momentum correction
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)

        if pol == "inbending":
            const = np.select([df_protonRec.Ptheta<27, (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [-0.0123049 + 0.00028887*df_protonRec.Ptheta, -0.138227479 + 8.07557430*0.001*df_protonRec.Ptheta -1.34807927*0.0001*df_protonRec.Ptheta*df_protonRec.Ptheta, -0.0275235])
            coeff = np.select([df_protonRec.Ptheta<27, (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [0.01528006 - 0.00024079*df_protonRec.Ptheta, 5.65817597*0.01 -2.36903348*0.001*df_protonRec.Ptheta + 4.93780046*0.00001*df_protonRec.Ptheta*df_protonRec.Ptheta, 0.03998975])    

            CorrectedPp = const + coeff/df_protonRec.loc[:, "Pp"] + df_protonRec.loc[:, "Pp"]

            const = np.select([df_protonRec.Ptheta<19.5, (df_protonRec.Ptheta>=19.5) & (df_protonRec.Ptheta<27), (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<39), (df_protonRec.Ptheta>=39) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [2.63643690*0.01, 0.50047232 -0.03834672 *df_protonRec.Ptheta + 0.00071967*df_protonRec.Ptheta*df_protonRec.Ptheta, 6.91308654 - 0.439839300*df_protonRec.Ptheta +6.83075548*0.001*df_protonRec.Ptheta*df_protonRec.Ptheta, 1.59424606, 1.47198581*10])
            coeff = np.select([df_protonRec.Ptheta<19.5, (df_protonRec.Ptheta>=19.5) & (df_protonRec.Ptheta<27), (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<39), (df_protonRec.Ptheta>=39) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [-1.46440415, 74.99891704  -6.1576777*df_protonRec.Ptheta + 0.11469137*df_protonRec.Ptheta*df_protonRec.Ptheta, 682.909471 - 43.9551177 * df_protonRec.Ptheta + 0.682383790 * df_protonRec.Ptheta * df_protonRec.Ptheta, -8.19627119, -23.55701865])    
            coeff2 = np.select([df_protonRec.Ptheta<19.5, (df_protonRec.Ptheta>=19.5) & (df_protonRec.Ptheta<27), (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<39), (df_protonRec.Ptheta>=39) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [-3.47690993, 47.71351973 -4.34918241*df_protonRec.Ptheta + 0.08841191*df_protonRec.Ptheta*df_protonRec.Ptheta, 100.33995753 - 6.96600416*df_protonRec.Ptheta + 0.11223046*df_protonRec.Ptheta*df_protonRec.Ptheta, -1.25261927, -0.40113733])    

            CorrectedPtheta = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Ptheta"]

            const = np.select([df_protonRec.Ptheta<16.5, (df_protonRec.Ptheta>=16.5) & (df_protonRec.Ptheta<27), (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [-0.190662844, -0.20725736 -0.00675627 *df_protonRec.Ptheta + 0.0007863*df_protonRec.Ptheta*df_protonRec.Ptheta, 12.1881698 - 0.78906294*df_protonRec.Ptheta +0.01297898*df_protonRec.Ptheta*df_protonRec.Ptheta, -4.59743066*10])
            coeff = np.select([df_protonRec.Ptheta<16.5, (df_protonRec.Ptheta>=16.5) & (df_protonRec.Ptheta<27), (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [6.48745941, 142.96379788  -16.66339055*df_protonRec.Ptheta + 0.51311212*df_protonRec.Ptheta*df_protonRec.Ptheta, 2.1853046 + 5.78521226 * df_protonRec.Ptheta - 0.09727796 * df_protonRec.Ptheta * df_protonRec.Ptheta, 7.46969457*10])    
            coeff2 = np.select([df_protonRec.Ptheta<16.5, (df_protonRec.Ptheta>=16.5) & (df_protonRec.Ptheta<27), (df_protonRec.Ptheta>=27) & (df_protonRec.Ptheta<42), df_protonRec.Ptheta>=42],
                              [-3.14646608, 17.39529095 -1.78403359*df_protonRec.Ptheta + 0.0335692*df_protonRec.Ptheta*df_protonRec.Ptheta, -1.03655317*10 + 0.161333213*df_protonRec.Ptheta -1.29625675*0.001*df_protonRec.Ptheta*df_protonRec.Ptheta, -4.41246899*0.1])    

            CorrectedPphi = const + coeff*np.exp(coeff2*df_protonRec.loc[:, "Pp"]) + df_protonRec.loc[:, "Pphi"]

            df_protonRec.loc[df_protonRec["Psector"]<7, "Pp"] = CorrectedPp
            df_protonRec.loc[df_protonRec["Psector"]<7, "Ptheta"] = CorrectedPtheta
            df_protonRec.loc[df_protonRec["Psector"]<7, "Pphi"] = CorrectedPphi

            df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
            df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
            df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
            pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

        # df_gammaRec = df_gammaRec[df_gammaRec["Gsector"]<7]
        #photon momentum correction
        newGpz = df_gammaRec.Gpz*0 + np.select([df_gammaRec.Gpz>=2, (df_gammaRec.Gpz<2) & (df_gammaRec.Gpz>1), df_gammaRec.Gpz<=1],[df_gammaRec.Gpz+0.13, df_gammaRec.Gpz+0.13*(df_gammaRec.Gpz-1), df_gammaRec.Gpz])
        newGpx = df_gammaRec.Gpx*0 + np.select([df_gammaRec.Gpz>=2, (df_gammaRec.Gpz<2) & (df_gammaRec.Gpz>1), df_gammaRec.Gpz<=1],[df_gammaRec.Gpx+0.13*df_gammaRec.Gpx/df_gammaRec.Gpz, df_gammaRec.Gpx+0.13*(df_gammaRec.Gpz-1)*df_gammaRec.Gpx/df_gammaRec.Gpz, df_gammaRec.Gpx])
        newGpy = df_gammaRec.Gpy*0 + np.select([df_gammaRec.Gpz>=2, (df_gammaRec.Gpz<2) & (df_gammaRec.Gpz>1), df_gammaRec.Gpz<=1],[df_gammaRec.Gpy+0.13*df_gammaRec.Gpy/df_gammaRec.Gpz, df_gammaRec.Gpy+0.13*(df_gammaRec.Gpz-1)*df_gammaRec.Gpy/df_gammaRec.Gpz, df_gammaRec.Gpy])
        
        df_gammaRec.loc[df_gammaRec["Gsector"]<7, "Gpx"] = newGpx
        df_gammaRec.loc[df_gammaRec["Gsector"]<7, "Gpy"] = newGpy
        df_gammaRec.loc[df_gammaRec["Gsector"]<7, "Gpz"] = newGpz

        df_gg = pd.merge(df_gammaRec, df_gammaRec,
                         how='outer', on='event', suffixes=("", "2"))
        df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
        df_ep = pd.merge(df_electronRec, df_protonRec, how='outer', on='event')

        df_epgg = pd.merge(df_ep, df_gg, how='outer', on='event')
        df_epgg = df_epgg[~np.isnan(df_epgg["Ppx"])]
        df_epgg = df_epgg[~np.isnan(df_epgg["Gpx"])]
        df_epgg = df_epgg[~np.isnan(df_epgg["Gpx2"])]

        print(len(df_gg))
        print(len(df_ep))
        print(len(df_epgg))

        self.df_epgg = df_epgg #temporarily save df_epgg

    def saveDVpi0vars(self):
        #set up pi0 variables
        df_epgg = self.df_epgg

        # useful objects
        ele = [df_epgg['Epx'], df_epgg['Epy'], df_epgg['Epz']]
        df_epgg.loc[:, 'Ep'] = mag(ele)
        df_epgg.loc[:, 'Ee'] = getEnergy(ele, me)
        df_epgg.loc[:, 'Etheta'] = getTheta(ele)
        df_epgg.loc[:, 'Ephi'] = getPhi(ele)

        pro = [df_epgg['Ppx'], df_epgg['Ppy'], df_epgg['Ppz']]

        gam = [df_epgg['Gpx'], df_epgg['Gpy'], df_epgg['Gpz']]
        df_epgg.loc[:, 'Gp'] = mag(gam)
        df_epgg.loc[:, 'Ge'] = getEnergy(gam, 0)
        df_epgg.loc[:, 'Gtheta'] = getTheta(gam)
        df_epgg.loc[:, 'Gphi'] = getPhi(gam)

        gam2 = [df_epgg['Gpx2'], df_epgg['Gpy2'], df_epgg['Gpz2']]
        df_epgg.loc[:, 'Gp2'] = mag(gam2)
        df_epgg.loc[:,'Ge2'] = getEnergy(gam2, 0)
        df_epgg.loc[:, 'Gtheta2'] = getTheta(gam2)
        df_epgg.loc[:, 'Gphi2'] = getPhi(gam2)

        pi0 = vecAdd(gam, gam2)
        VGS = [-df_epgg['Epx'], -df_epgg['Epy'], pbeam - df_epgg['Epz']]
        v3l = cross(beam, ele)
        v3h = cross(pro, VGS)
        v3g = cross(VGS, gam)
        VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
                    df_epgg["Ppy"], pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]
        VmissP = [-df_epgg["Epx"] - df_epgg["Gpx"] - df_epgg["Gpx2"], -df_epgg["Epy"] -
                    df_epgg["Gpy"] - df_epgg["Gpy2"], pbeam - df_epgg["Epz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        Vmiss = [-df_epgg["Epx"] - df_epgg["Ppx"] - df_epgg["Gpx"] - df_epgg["Gpx2"],
                    -df_epgg["Epy"] - df_epgg["Ppy"] - df_epgg["Gpy"] - df_epgg["Gpy2"],
                    pbeam - df_epgg["Epz"] - df_epgg["Ppz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]

        df_epgg.loc[:, 'Mpx'], df_epgg.loc[:, 'Mpy'], df_epgg.loc[:, 'Mpz'] = Vmiss

        # binning kinematics
        df_epgg.loc[:,'Q2'] = -((ebeam - df_epgg['Ee'])**2 - mag2(VGS))
        df_epgg.loc[:,'nu'] = (ebeam - df_epgg['Ee'])
        df_epgg.loc[:,'xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
        df_epgg.loc[:,'t'] = 2 * M * (df_epgg['Pe'] - M)
        df_epgg.loc[:,'W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epgg['Ee'])**2 - mag2(VGS)))
        df_epgg.loc[:,'MPt'] = np.sqrt((df_epgg["Epx"] + df_epgg["Ppx"] + df_epgg["Gpx"] + df_epgg["Gpx2"])**2 +
                                 (df_epgg["Epy"] + df_epgg["Ppy"] + df_epgg["Gpy"] + df_epgg["Gpy2"])**2)
        # trento angles
        df_epgg.loc[:,'phi1'] = angle(v3l, v3h)
        df_epgg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_epgg['phi1'], df_epgg['phi1'])

        # exclusivity variables
        df_epgg.loc[:,'MM2_ep'] = (-M - ebeam + df_epgg["Ee"] +
                             df_epgg["Pe"])**2 - mag2(VmissPi0)
        df_epgg.loc[:,'MM2_egg'] = (-M - ebeam + df_epgg["Ee"] +
                             df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(VmissP)
        df_epgg.loc[:,'MM2_epgg'] = (-M - ebeam + df_epgg["Ee"] + df_epgg["Pe"] +
                             df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(Vmiss)
        df_epgg.loc[:,'ME_epgg'] = (M + ebeam - df_epgg["Ee"] - df_epgg["Pe"] - df_epgg["Ge"] - df_epgg["Ge2"])
        df_epgg.loc[:,'Mpi0'] = pi0InvMass(gam, gam2)
        df_epgg.loc[:,'reconPi'] = angle(VmissPi0, pi0)
        df_epgg.loc[:,"Pie"] = df_epgg['Ge'] + df_epgg['Ge2']
        self.df_epgg = df_epgg

    def makeDVpi0(self):
        #make dvpi0 pairs
        df_epgg = self.df_epgg

        df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)

        cut_xBupper = df_epgg.loc[:, "xB"] < 1  # xB
        cut_xBlower = df_epgg.loc[:, "xB"] > 0  # xB
        cut_Q2 = df_epgg.loc[:, "Q2"] > 1  # Q2
        cut_W = df_epgg.loc[:, "W"] > 2  # W

        # Exclusivity cuts
        cut_mmep = df_epgg.loc[:, "MM2_ep"] < 0.7  # mmep
        cut_meepgg = df_epgg.loc[:, "ME_epgg"] < 0.7  # meepgg
        cut_mpt = df_epgg.loc[:, "MPt"] < 0.2  # mpt
        cut_recon = df_epgg.loc[:, "reconPi"] < 2  # recon gam angle
        cut_pi0upper = df_epgg.loc[:, "Mpi0"] < 0.2
        cut_pi0lower = df_epgg.loc[:, "Mpi0"] > 0.07
        cut_sector = (df_epgg.loc[:, "Esector"]!=df_epgg.loc[:, "Gsector"]) & (df_epgg.loc[:, "Esector"]!=df_epgg.loc[:, "Gsector2"])
        cut_Vz = np.abs(df_epgg["Evz"] - df_epgg["Pvz"]) < 2.5 + 2.5 / mag([df_epgg["Ppx"], df_epgg["Ppy"], df_epgg["Ppz"]])

        df_dvpi0 = df_epgg.loc[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_mmep & cut_meepgg & cut_Vz &
                           cut_mpt & cut_recon & cut_pi0upper & cut_pi0lower & cut_sector, :]

        print(len(df_dvpi0))

        #For an event, there can be two gg's passed conditions above.
        #Take only one gg's that makes pi0 invariant mass
        #This case is very rare.
        #For now, duplicated proton is not considered.
        df_dvpi0 = df_dvpi0.sort_values(by=['closeness', 'Psector', 'Gsector'], ascending = [True, True, True])
        df_dvpi0 = df_dvpi0.loc[~df_dvpi0.event.duplicated(), :]
        df_dvpi0 = df_dvpi0.sort_values(by='event')        
        print(len(df_dvpi0))
        self.df_x = df_dvpi0 #done with saving x

    def saveRaw(self):
        df_x = self.df_x
        self.df = df_x


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-p","--polarity", help="polarity", default = "inbending")

    args = parser.parse_args()

    converter = root2pickle(args.fname, entry_stop = args.entry_stop, pol = args.polarity)
    df = converter.df

    df.to_pickle(args.out)