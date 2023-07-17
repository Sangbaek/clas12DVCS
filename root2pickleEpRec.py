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
        self.saveEPvars()
        # self.save(raw = raw, pol = pol)


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

        # data frames and their keys to read Z part
        df_electronGen = pd.DataFrame()
        df_protonGen = pd.DataFrame()

        eleKeysGen = ["GenEpx", "GenEpy", "GenEpz", "nmEBAR"]
        proKeysGen = ["GenPpx", "GenPpy", "GenPpz"]

        # read keys
        for key in eleKeysGen:
            df_electronGen[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in proKeysGen:
            df_protonGen[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))

        df_electronGen = df_electronGen.astype({"GenEpx": float, "GenEpy": float, "GenEpz": float})
        df_protonGen = df_protonGen.astype({"GenPpx": float, "GenPpy": float, "GenPpz": float})

        #set up a dummy index for merging
        df_electronGen.loc[:,'event'] = df_electronGen.index
        df_protonGen.loc[:,'event'] = df_protonGen.index

        #spherical coordinates
        eleGen = [df_electronGen["GenEpx"], df_electronGen["GenEpy"], df_electronGen["GenEpz"]]
        df_electronGen.loc[:, 'GenEp'] = mag(eleGen)
        df_electronGen.loc[:, 'GenEtheta'] = getTheta(eleGen)
        df_electronGen.loc[:, 'GenEphi'] = getPhi(eleGen)

        proGen = [df_protonGen["GenPpx"], df_protonGen["GenPpy"], df_protonGen["GenPpz"]]
        df_protonGen.loc[:, 'GenPp'] = mag(proGen)
        df_protonGen.loc[:, 'GenPtheta'] = getTheta(proGen)
        df_protonGen.loc[:, 'GenPphi'] = getPhi(proGen)

        df_MC = pd.merge(df_electronGen, df_protonGen, how='inner', on='event')
        self.df_MC = df_MC

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Estat", "nmebar"]
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pstat"]
        if logistics:
            eleKeysRec.extend(["EventNum", "RunNum", "beamQ", "liveTime", "helicity"])

        # read them
        for key in eleKeysRec:
            df_electronRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in proKeysRec:
            df_protonRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        ele = [df_electronRec['Epx'], df_electronRec['Epy'], df_electronRec['Epz']]
        df_electronRec.loc[:, 'Ep'] = mag(ele)

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

        self.df_ep = df_ep # saves df_ep

    def saveEPvars(self, correction=None):
        #set up dvcs variables
        df_ep = self.df_ep

        ele = [df_ep['Epx'], df_ep['Epy'], df_ep['Epz']]
        df_ep.loc[:, 'Ep'] = mag(ele)
        df_ep.loc[:, 'Ee'] = getEnergy(ele, me)
        df_ep.loc[:, 'Etheta'] = getTheta(ele)
        df_ep.loc[:, 'Ephi'] = getPhi(ele)

        pro = [df_ep['Ppx'], df_ep['Ppy'], df_ep['Ppz']]

        Ppt = mag([df_ep['Ppx'], df_ep['Ppy'], 0])

        VGS = [-df_ep['Epx'], -df_ep['Epy'], self.pbeam - df_ep['Epz']]
        v3l = cross(self.beam, ele)
        v3h = cross(pro, VGS)
        Vmiss = [-df_ep["Epx"] - df_ep["Ppx"], -df_ep["Epy"] - df_ep["Ppy"],
                  self.pbeam - df_ep["Epz"] - df_ep["Ppz"]]

        # binning kinematics
        df_ep.loc[:,'Q2'] = -((self.ebeam - df_ep['Ee'])**2 - mag2(VGS))
        df_ep.loc[:,'nu'] = (self.ebeam - df_ep['Ee'])
        df_ep.loc[:,'y'] = df_ep['nu']/self.ebeam
        df_ep.loc[:,'xB'] = df_ep['Q2'] / 2.0 / M / df_ep['nu']
        df_ep.loc[:,'t1'] = 2 * M * (df_ep['Pe'] - M)
        df_ep.loc[:,'W'] = np.sqrt(np.maximum(0, (self.ebeam + M - df_ep['Ee'])**2 - mag2(VGS)))

        # trento angles
        df_ep.loc[:,'phi1'] = angle(v3l, v3h)
        df_ep.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_ep['phi1'], df_ep['phi1'])

        # exclusivity variables
        df_ep.loc[:,'MM2_ep'] = (-M - self.ebeam + df_ep["Ee"] + df_ep["Pe"])**2 - mag2(Vmiss)

        df_ep = pd.merge(df_ep, self.df_MC, how = 'inner', on='event')
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
    df = converter.df_ep

    df.to_pickle(args.out)