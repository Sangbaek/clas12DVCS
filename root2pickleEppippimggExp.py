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
     fidlevel = 'mid', allowsamesector = False, allowduplicates = False, ebeam = 10.604):
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
            fidlevel = fidlevel)

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
        nocorr = False, noeloss = False, nopcorr = False, fidlevel = 'mid'):
        '''save data into df_ep for parent class epg'''
        self.readFile()

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_pipRec = pd.DataFrame()
        df_pimRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()

        eleKeysRec = ["Epx", "Epy", "Epz", "Estat"]
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pstat"]
        pipKeysRec = ["Pippx", "Pippy", "Pippz"]
        pimKeysRec = ["Pimpx", "Pimpy", "Pimpz"]
        gamKeysRec = ["Gpx", "Gpy", "Gpz"]

        # read them
        for key in eleKeysRec:
            df_electronRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in proKeysRec:
            df_protonRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in pipKeysRec:
            df_pipRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in pimKeysRec:
            df_pimRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in gamKeysRec:
            df_gammaRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        if logistics:
            df_logisticsRec = pd.DataFrame()
            logKeysRec = ["nma", "nmc", "TriggerBit", "EventNum", "RunNum", "beamQ", "liveTime", "helicity"]
            for key in logKeysRec:
                df_logisticsRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
            df_logisticsRec.loc[:,'event'] = df_logisticsRec.index

        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index.get_level_values('entry')
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')

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

        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)
        df_protonRec.loc[:, 't1'] = 2 * M * (df_protonRec['Pe'] - M)

        df_ep = pd.merge(df_electronRec, df_protonRec, how='outer', on='event')
        ep_inner_cond = (~np.isnan(df_ep.Epx)) & (~np.isnan(df_ep.Ppx))
        Vmiss = [-df_ep.loc[ep_inner_cond, "Epx"] - df_ep.loc[ep_inner_cond, "Ppx"], -df_ep.loc[ep_inner_cond, "Epy"] - df_ep.loc[ep_inner_cond, "Ppy"],
                  self.pbeam - df_ep.loc[ep_inner_cond, "Epz"] - df_ep.loc[ep_inner_cond, "Ppz"]]
        df_ep.loc[ep_inner_cond, 'MM2_ep'] = (-M - self.ebeam + df_ep.loc[ep_inner_cond, "Ee"] + df_ep.loc[ep_inner_cond, "Pe"])**2 - mag2(Vmiss)
        if logistics:
            df_ep = pd.merge(df_ep, df_logisticsRec, how='outer', on='event')

        print("Electron, proton, ep cases:")
        print(len(df_electronRec), len(df_protonRec), len(self.df_ep))

        #convert data type to standard double
        df_pipRec = df_pipRec.astype({"Pippx": float, "Pippy": float, "Pippz": float})
        df_pimRec = df_pimRec.astype({"Pimpx": float, "Pimpy": float, "Pimpz": float})
        df_gammaRec = df_gammaRec.astype({"Gpx": float, "Gpy": float, "Gpz": float})
        #set up a dummy index for merging
        df_pipRec.loc[:,'event'] = df_pipRec.index.get_level_values('entry')
        df_pimRec.loc[:,'event'] = df_pimRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'event'] = df_gammaRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'GIndex'] = df_gammaRec.index.get_level_values('subentry')

        df_gg = pd.merge(df_gammaRec, df_gammaRec,
                         how='outer', on='event', suffixes=("", "2"))
        df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
        df_gg = df_gg.drop(['GIndex', 'GIndex2'], axis = 1)

        gam = [df_gg['Gpx'], df_gg['Gpy'], df_gg['Gpz']]
        df_gg.loc[:, 'Gp'] = mag(gam)
        df_gg.loc[:, 'Ge'] = getEnergy(gam, 0)
        df_gg.loc[:, 'Gtheta'] = getTheta(gam)
        df_gg.loc[:, 'Gphi'] = getPhi(gam)

        gam2 = [df_gg['Gpx2'], df_gg['Gpy2'], df_gg['Gpz2']]
        df_gg.loc[:, 'Gp2'] = mag(gam2)
        df_gg.loc[:,'Ge2'] = getEnergy(gam2, 0)
        df_gg.loc[:, 'Gtheta2'] = getTheta(gam2)
        df_gg.loc[:, 'Gphi2'] = getPhi(gam2)

        pip = [df_pipRec.Pipx, df_pipRec.Pipy, df_pipRec.Pipz]
        df_pipRec.loc[:, "Pipp"] = mag(pip)
        df_pipRec.loc[:, "Pipe"] = getEnergy(pip, 0.13957039)
        df_pipRec.loc[:, "Piptheta"] = getTheta(pip)
        df_pipRec.loc[:, "Pipphi"] = getPhi(pip)

        pim = [df_pimRec.Pipx, df_pimRec.Pipy, df_pimRec.Pipz]
        df_pimRec.loc[:, "Pimp"] = mag(pim)
        df_pimRec.loc[:, "Pime"] = getEnergy(pim, 0.13957039)
        df_pimRec.loc[:, "Pimtheta"] = getTheta(pim)
        df_pimRec.loc[:, "Pimphi"] = getPhi(pim)

        df_gg.loc[:,'Mgg'] = pi0InvMass(gam, gam2)

        df_pippimeta = pd.merge(df_pipRec, df_pimRec, how='inner', on='event')
        df_pippimeta = pd.merge(df_pippimeta, df_gg, how='inner', on='event')
        df_pippimeta.loc[:, "IM2_pippimgg"] = (df_pippimeta.Ge + df_pippimeta.Ge2 + df_pippimeta.Pipe + df_pippimeta.Pime)**2 - \
        (df_pippimeta.Gpx + df_pippimeta.Gpx2 + df_pippimeta.Pippx + df_pippimeta.Pimpx)**2 -\
        (df_pippimeta.Gpy + df_pippimeta.Gpy2 + df_pippimeta.Pippy + df_pippimeta.Pimpy)**2 -\
        (df_pippimeta.Gpz + df_pippimeta.Gpz2 + df_pippimeta.Pippz + df_pippimeta.Pimpz)**2 
        df_pippimeta.loc[:, "IM2_pimgg"] = (df_pippimeta.Ge + df_pippimeta.Ge2 + df_pippimeta.Pime)**2 - \
        (df_pippimeta.Gpx + df_pippimeta.Gpx2 + df_pippimeta.Pimpx)**2 -\
        (df_pippimeta.Gpy + df_pippimeta.Gpy2 + df_pippimeta.Pimpy)**2 -\
        (df_pippimeta.Gpz + df_pippimeta.Gpz2 + df_pippimeta.Pimpz)**2 
        df_pippimeta.loc[:, "IM2_pipgg"] = (df_pippimeta.Ge + df_pippimeta.Ge2 + df_pippimeta.Pipe)**2 - \
        (df_pippimeta.Gpx + df_pippimeta.Gpx2 + df_pippimeta.Pippx)**2 -\
        (df_pippimeta.Gpy + df_pippimeta.Gpy2 + df_pippimeta.Pippy)**2 -\
        (df_pippimeta.Gpz + df_pippimeta.Gpz2 + df_pippimeta.Pippz)**2 
        df_pippimeta.loc[:, "IM2_pippim"] = (df_pippimeta.Pipe + df_pippimeta.Pime)**2 - \
        (df_pippimeta.Pippx + df_pippimeta.Pimpx)**2 -\
        (df_pippimeta.Pippy + df_pippimeta.Pimpy)**2 -\
        (df_pippimeta.Pippz + df_pippimeta.Pimpz)**2 

        self.df_eppippimeta = pd.merge(df_ep, df_pippimeta, how = 'right', on = 'event')



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

    args = parser.parse_args()

    if args.entry_start:
        args.entry_start = int(args.entry_start)
    if args.entry_stop:
        args.entry_stop = int(args.entry_stop)

    be = float(args.beam)
    converter = root2pickle(args.fname, entry_start = args.entry_start,
     entry_stop = args.entry_stop, pol = args.polarity, detRes = args.detRes, raw = args.raw,
     logistics = args.logistics, width = args.width, nofid = args.nofid, nocorr = args.nocorr, noeloss = args.noeloss,nopcorr = args.nopcorr,
     fidlevel = args.fidlevel, allowsamesector = args.allowsamesector, allowduplicates = args.allowduplicates, ebeam = be)

    filename = args.out

    df = converter.df_eppippimeta
    df.to_pickle(filename)