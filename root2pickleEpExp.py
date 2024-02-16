#!/usr/bin/env python3
"""
A simple script to save data in pickle.
"""

import uproot
import argparse
from copy import copy
from utils.const import *
from utils.physics import *
from utils.fiducial import *
from utils.kinCorrection import *
import awkward as ak
pd.options.mode.chained_assignment = None 

class root2pickle():
    #class to read root to make epg pairs, inherited from epg
    def __init__(self, fname, entry_start = None, entry_stop = None, pol = "inbending",
     detRes = False, raw = False, logistics = False, width = "mid", nofid = False, nocorr = False, noeloss = False, nopcorr = False,
     fidlevel = 'mid', allowsamesector = False, allowduplicates = False, ebeam = 10.604, ebar = False):
        '''
            clas init.
            Args
            --------------------
            fname: root file name to be read
            entry_start: the lower bound of root entry
            entry_stop: the upper bound of root entry
            pol: polarity
            gen: generator
            raw: no exclusivity cuts
            detRes: include detector responses in the output (full skim)
            width: data selection window width
            nofid: do not apply fid cuts.
            nocorr: do not apply the momentum correction
            allowsamesector: allow same sectors.

            Attributes
            --------------------
            fname: root file name to be read
            Methods
            --------------------
            determineWidth: to determine event selection window
            readFile: read root file to pandas
            closeFile: nullify the files to save memory
            readEP: read and post-process data. fiducial cut/ proton energy loss correction/ reconstruction bias correction.
            saveDVCSvars: 4 momentum algebra to save DVCS vars
            saveDVpi0vars: 4 momentum algebra to save DVpi0P vars
            makeDVpi0P_DVCS: select pi0->2g events that are overlapped w/ DVCS.
            pi02gSubtraction: exclude pi0->2g events first
            makeDVCS: select BH-DVCS candidates
            save: save output
        '''
        self.fname = fname
        self.ebeam = ebeam # beam energy
        self.pbeam = np.sqrt(ebeam * ebeam - me * me) # beam electron momentum
        self.beam = [0, 0, self.pbeam] # beam vector

        self.readEP(entry_start = entry_start, entry_stop = entry_stop, pol = pol, 
            detRes = detRes, logistics = logistics, nofid = nofid, nocorr = nocorr, noeloss = noeloss, nopcorr = nopcorr,
            fidlevel = fidlevel, ebar = ebar)
        self.saveEPvars(ebar = ebar)

    def readFile(self):
        '''read root using uproot'''
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        '''close file for saving memory'''
        self.file = None
        self.tree = None

    def readEP(self, entry_start = None, entry_stop = None, pol = "inbending", 
        detRes = False, logistics = False, nofid = False, 
        nocorr = False, noeloss = False, nopcorr = False, fidlevel = 'mid', ebar = False):
        '''save data into df_ep for parent class epg'''
        self.readFile()

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Epa", "Estat"]
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pstat"]

        if ebar:
            df_positronRec = pd.DataFrame()
            posKeysRec = ["Ebarpx", "Ebarpy", "Ebarpz"]
            for key in posKeysRec:
                df_positronRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
            df_positronRec = df_positronRec.astype({"Ebarpx": float, "Ebarpy": float, "Ebarpz": float})
            df_positronRec.loc[:, "event"] = df_positronRec.index.get_level_values('entry')
            pos = [df_positronRec.Ebarpx, df_positronRec.Ebarpy, df_positronRec.Ebarpz]
            df_positronRec.loc[:, "Ebarp"] = mag(pos)
            df_positronRec.loc[:, "Ebare"] = getEnergy(pos, me)
            df_positronRec.loc[:, "Ebartheta"] = getTheta(pos)
            df_positronRec.loc[:, "Ebarphi"] = getPhi(pos)

        # read them
        for key in eleKeysRec:
            df_electronRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in proKeysRec:
            df_protonRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        if logistics:
            df_logisticsRec = pd.DataFrame()
            logKeysRec = ["nmlbar", "nma", "nmc", "TriggerPid", "TriggerBit", "EventNum", "RunNum", "beamQ", "liveTime", "helicity"]
            for key in logKeysRec:
                df_logisticsRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
            df_logisticsRec.loc[:,'event'] = df_logisticsRec.index

        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        ele = [df_electronRec['Epx'], df_electronRec['Epy'], df_electronRec['Epz']]
        df_electronRec.loc[:, 'Ep'] = mag(ele)
        df_electronRec.loc[:, 'Ee'] = getEnergy(ele, me)
        df_electronRec.loc[:, 'Etheta'] = getTheta(ele)
        df_electronRec.loc[:, 'Ephi'] = getPhi(ele)
        VGS = [-df_electronRec['Epx'], -df_electronRec['Epy'], self.pbeam - df_electronRec['Epz']]
        df_electronRec.loc[:,'Q2'] = -((self.ebeam - df_electronRec['Ee'])**2 - mag2(VGS))
        df_electronRec.loc[:,'nu'] = (self.ebeam - df_electronRec['Ee'])
        df_electronRec.loc[:,'y'] = df_electronRec['nu']/self.ebeam
        df_electronRec.loc[:,'xB'] = df_electronRec['Q2'] / 2.0 / M / df_electronRec['nu']
        df_electronRec.loc[:,'W'] = np.sqrt(np.maximum(0, (self.ebeam + M - df_electronRec['Ee'])**2 - mag2(VGS)))

        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index.get_level_values('entry')
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')
        #apply fiducial cuts
        print(len(df_electronRec), len(df_protonRec))

        #prepare for proton energy loss corrections correction
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)
        
        df_ep = pd.merge(df_electronRec, df_protonRec, how='outer', on='event')
        df_ep = df_ep.loc[df_ep.Estat < 2000, :] #Save FT only

        if ebar:
            df_eebar = pd.merge(df_electronRec, df_positronRec, how = 'inner', on = 'event')
            eebarInvmass2 = (df_eebar.Ee+df_eebar.Ebare)**2 - (df_eebar.Epx+df_eebar.Ebarpx)**2 - (df_eebar.Epy+df_eebar.Ebarpy)**2 - (df_eebar.Epz+df_eebar.Ebarpz)**2
            eebarInvmass2 = np.where(eebarInvmass2 >= 0, eebarInvmass2, 10**6)
            eebarInvmass  = np.sqrt(eebarInvmass2)
            eebarInvmass = np.where(eebarInvmass > 100, -1000, eebarInvmass)
            df_eebar.loc[:, "IM_eebar"] = eebarInvmass
            df_eebar = pd.merge(df_eebar, df_protonRec, how = 'inner', on = 'event')
            df_eebar.loc[:, "Ge"] = df_eebar.Ee + df_eebar.Ebare + df_eebar.Pe - M
            df_eebar.loc[:, "Mpx"] = - (df_eebar.Epx + df_eebar.Ebarpx + df_eebar.Ppx)
            df_eebar.loc[:, "Mpy"] = - (df_eebar.Epy + df_eebar.Ebarpy + df_eebar.Ppy)
            df_eebar.loc[:, "Mpz"] = self.pbeam - (df_eebar.Epz + df_eebar.Ebarpz + df_eebar.Ppz)
            df_eebar.loc[:, "Mp"] = mag([df_eebar.Mpx, df_eebar.Mpy, df_eebar.Mpz])
            df_eebar.loc[:, "ME"] = self.ebeam + M - (df_eebar.Ee + df_eebar.Ebare + df_eebar.Pe)
            df_eebar.loc[:, "MM2_X"] = df_eebar.ME**2 - df_eebar.Mp**2
            df_eebar.loc[:, "costheta_miss"] = df_eebar.Mpz/ df_eebar.Mp
            df_eebar.loc[:, "Q2"] = 2*self.ebeam*df_eebar.Mp*(1-df_eebar.costheta_miss)
            scat_ele = [df_eebar.Mpx, df_eebar.Mpy, df_eebar. Mpz]
            df_eebar.loc[:, "SEp"] = mag(scat_ele)
            df_eebar.loc[:, "SEe"] = getEnergy(scat_ele, me)
            df_eebar.loc[:, "SEtheta"] = getTheta(scat_ele)
            df_eebar.loc[:, "SEphi"] = getPhi(scat_ele)
            Vmiss = [-df_eebar["Mpx"] - df_eebar["Ppx"], -df_eebar["Mpy"] - df_eebar["Ppy"],
                      10.604 - df_eebar["Mpz"] - df_eebar["Ppz"]]
            df_eebar.loc[:,'MM2_ep'] = (-M - 10.604 + df_eebar["SEe"] + df_eebar["Pe"])**2 - mag2(Vmiss)
            VGS = [-df_eebar['Mpx'], -df_eebar['Mpy'], 10.604 - df_eebar['Mpz']]
            df_eebar.loc[:,'Q2_new'] = -((10.604 - df_eebar['SEe'])**2 - mag2(VGS))
            df_eebar.loc[:,'nu'] = (10.604 - df_eebar['SEe'])
            df_eebar.loc[:,'y'] = df_eebar['nu']/10.604
            df_eebar.loc[:,'xB'] = df_eebar['Q2_new'] / 2.0 / M / df_eebar['nu']
            df_eebar.loc[:,'W'] = np.sqrt(np.maximum(0, (10.604 + M - df_eebar['SEe'])**2 - mag2(VGS)))

            df_eeebar = pd.merge(df_electronRec, df_eebar, how = 'inner', on = 'event', suffixes=("", "jpsi"))
            df_eeebar = df_eeebar.loc[(df_eeebar.Epa != df_eeebar.Epajpsi) & (df_eeebar.Estat < 2000), :]

            df_peebar = pd.merge(df_protonRec, df_eebar, how = 'inner', on = 'event')
            df_epeebar = pd.merge(df_protonRec, df_eeebar, how = 'inner', on = 'event')

        if logistics:
            df_ep = pd.merge(df_ep, df_logisticsRec, how='outer', on='event')
            df_eebar = pd.merge(df_eebar, df_logisticsRec, how='outer', on='event')
            df_eeebar = pd.merge(df_eeebar, df_logisticsRec, how='outer', on='event')
            df_peebar = pd.merge(df_peebar, df_logisticsRec, how='outer', on='event')
            df_epeebar = pd.merge(df_epeebar, df_logisticsRec, how='outer', on='event')

        self.df_ep = df_ep
        self.df_eebar = df_eebar
        self.df_eeebar = df_eeebar
        self.df_peebar = df_peebar
        self.df_epeebar = df_epeebar


    def saveEPvars(self, correction=None, ebar = False):
        #set up dvcs variables
        df_ep = self.df_ep

        ele = [df_ep['Epx'], df_ep['Epy'], df_ep['Epz']]
        pro = [df_ep['Ppx'], df_ep['Ppy'], df_ep['Ppz']]

        # Ppt = mag([df_ep['Ppx'], df_ep['Ppy'], 0])

        # v3l = cross(self.beam, ele)
        # v3h = cross(pro, VGS)
        Vmiss = [-df_ep["Epx"] - df_ep["Ppx"], -df_ep["Epy"] - df_ep["Ppy"],
                  self.pbeam - df_ep["Epz"] - df_ep["Ppz"]]

        # binning kinematics

        # exclusivity variables
        df_ep.loc[:,'t1'] = 2 * M * (df_ep['Pe'] - M)
        df_ep.loc[:,'MM2_ep'] = (-M - self.ebeam + df_ep["Ee"] + df_ep["Pe"])**2 - mag2(Vmiss)

        self.df_ep = df_ep


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-S","--entry_start", help="entry_start to start reading the root file", default = None)
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-p","--polarity", help="polarity", default = "inbending")
    parser.add_argument("-d","--detRes", help="include detector response", action = "store_true")
    parser.add_argument("-l","--logistics", help="include logistics", action = "store_true")
    parser.add_argument("-w","--width", help="width of selection cuts", default = "default")
    parser.add_argument("-r","--raw", help="save raw only", default = False, action = "store_true")
    parser.add_argument("-nf","--nofid", help="no additional fiducial cuts", action = "store_true")
    parser.add_argument("-nc","--nocorr", help="no momentum correction", action = "store_true")
    parser.add_argument("-ne","--noeloss", help="no energy loss correction", action = "store_true")
    parser.add_argument("-np","--nopcorr", help="no proton correction at all", action = "store_true")
    parser.add_argument("-fl","--fidlevel", help="fiducial cut level", default = 'mid')
    parser.add_argument("-as","--allowsamesector", help="allow same sector conditions", action = "store_true")
    parser.add_argument("-ad","--allowduplicates", help="allow duplicates", action = "store_true")
    parser.add_argument("-be","--beam", help="beam energy", default = "10.604")
    parser.add_argument("-ebar","--ebar", help="use positron or not", action = "store_true")

    args = parser.parse_args()

    if args.entry_start:
        args.entry_start = int(args.entry_start)
    if args.entry_stop:
        args.entry_stop = int(args.entry_stop)

    be = float(args.beam)
    converter = root2pickle(args.fname, entry_start = args.entry_start,
     entry_stop = args.entry_stop, pol = args.polarity, detRes = args.detRes, raw = args.raw,
     logistics = args.logistics, width = args.width, nofid = args.nofid, nocorr = args.nocorr, noeloss = args.noeloss,nopcorr = args.nopcorr,
     fidlevel = args.fidlevel, allowsamesector = args.allowsamesector, allowduplicates = args.allowduplicates, ebeam = be, ebar = args.ebar)

    filename_ep      = args.out.replace(".pkl", "_ep.pkl")
    filename_eebar   = args.out.replace(".pkl", "_eebar.pkl")
    filename_eeebar  = args.out.replace(".pkl", "_eeebar.pkl")
    filename_peebar  = args.out.replace(".pkl", "_peebar.pkl")
    filename_epeebar = args.out.replace(".pkl", "_epeebar.pkl")

    df = converter.df_ep
    df.to_pickle(filename_ep)

    if args.ebar:
        converter.df_eebar.to_pickle("filename_eebar") 
        converter.df_eeebar.to_pickle("filename_eeebar")
        converter.df_peebar.to_pickle("filename_peebar")
        converter.df_epeebar.to_pickle("filename_epeebar")