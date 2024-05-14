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
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

class root2pickle():
    #class to read root to make epg pairs, inherited from epg
    def __init__(self, fname, entry_start = None, entry_stop = None, pol = "inbending",
     detRes = False, raw = False, logistics = False, width = "mid", nofid = False, nocorr = False, noeloss = False, nopcorr = False,
     fidlevel = 'mid', allowsamesector = False, allowduplicates = False, ebeam = 10.604, efficiency = False):
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
            readEPGG: read and post-process data. fiducial cut/ proton energy loss correction/ reconstruction bias correction.
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

        self.determineWidth(width = width)
        self.readBinScheme()
        self.readEPGG(entry_start = entry_start, entry_stop = entry_stop, pol = pol, 
            detRes = detRes, logistics = logistics, nofid = nofid, nocorr = nocorr, noeloss = noeloss, nopcorr = nopcorr,
            fidlevel = fidlevel, efficiency = efficiency)
        self.saveDVCSvars()
        self.saveDVpi0vars()
        if not raw:
            self.makeDVpi0P_DVCS(pol = pol, nofid = nofid)
            self.pi02gSubtraction()
            self.makeDVCS(pol = pol, nofid = nofid, allowsamesector = allowsamesector, allowduplicates = allowduplicates)
        self.save(raw = raw, pol = pol, efficiency = efficiency)


    def readBinScheme(self):
        self.bin_scheme = np.loadtxt('/work/clas12/sangbaek/km15gen/bin_scheme.csv', delimiter = ',')
        self.fringe_bin_scheme = np.loadtxt('/work/clas12/sangbaek/km15gen/fringe_bin_scheme.csv', delimiter = ',')

    def readFile(self):
        '''read root using uproot'''
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        '''close file for saving memory'''
        self.file = None
        self.tree = None

    def determineWidth(self, width = "mid"):
        '''determine event selection window'''
        print("determine width level: {}".format(width))
        if width == "default":
            self.cuts_dvcs_CDFT_Inb  = cuts_dvcs_default
            self.cuts_dvcs_CD_Inb    = cuts_dvcs_default
            self.cuts_dvcs_FD_Inb    = cuts_dvcs_default
            self.cuts_dvcs_CDFT_Outb = cuts_dvcs_default
            self.cuts_dvcs_CD_Outb   = cuts_dvcs_default
            self.cuts_dvcs_FD_Outb   = cuts_dvcs_default
        if width == "mid":
            self.cuts_dvcs_CDFT_Inb = cuts_dvcs_CDFT_Inb_3sigma
            self.cuts_dvcs_CD_Inb = cuts_dvcs_CD_Inb_3sigma
            self.cuts_dvcs_FD_Inb = cuts_dvcs_FD_Inb_3sigma
            self.cuts_dvcs_CDFT_Outb = cuts_dvcs_CDFT_Outb_3sigma
            self.cuts_dvcs_CD_Outb = cuts_dvcs_CD_Outb_3sigma
            self.cuts_dvcs_FD_Outb = cuts_dvcs_FD_Outb_3sigma
        if width == "tight":
            self.cuts_dvcs_CDFT_Inb = cuts_dvcs_CDFT_Inb_2sigma
            self.cuts_dvcs_CD_Inb = cuts_dvcs_CD_Inb_2sigma
            self.cuts_dvcs_FD_Inb = cuts_dvcs_FD_Inb_2sigma
            self.cuts_dvcs_CDFT_Outb = cuts_dvcs_CDFT_Outb_2sigma
            self.cuts_dvcs_CD_Outb = cuts_dvcs_CD_Outb_2sigma
            self.cuts_dvcs_FD_Outb = cuts_dvcs_FD_Outb_2sigma
        if width == "loose":
            self.cuts_dvcs_CDFT_Inb = cuts_dvcs_CDFT_Inb_4sigma
            self.cuts_dvcs_CD_Inb = cuts_dvcs_CD_Inb_4sigma
            self.cuts_dvcs_FD_Inb = cuts_dvcs_FD_Inb_4sigma
            self.cuts_dvcs_CDFT_Outb = cuts_dvcs_CDFT_Outb_4sigma
            self.cuts_dvcs_CD_Outb = cuts_dvcs_CD_Outb_4sigma
            self.cuts_dvcs_FD_Outb = cuts_dvcs_FD_Outb_4sigma

    def readEPGG(self, entry_start = None, entry_stop = None, pol = "inbending", 
        detRes = False, logistics = False, nofid = False, 
        nocorr = False, noeloss = False, nopcorr = False, fidlevel = 'mid', efficiency = False):
        '''save data into df_epg, df_epgg for parent class epg'''
        self.readFile()

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Eedep", "Evz", "Esector", "Estat", "Epa"]
        eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
        eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
        eleKeysRec.extend(["EcalHx1", "EcalHy1"])
        eleKeysRec.extend(["EcalHx3", "EcalHy3"])
        eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc2Hitx", "EDc2Hity", "EDc2Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
        eleKeysRec.extend(["EFtof1bSector", "EFtof1bComponent"])
        eleKeysRec.extend(["Enphe", "EhtccX", "EhtccY"])
        if efficiency:
            eleKeysRec.extend(["EhtcctrajX", "EhtcctrajY"])
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "Psector", "Pchi2pid"]
        proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"])
        proKeysRec.extend(["PDc2Hitx", "PDc2Hity", "PDc2Hitz", "PDc3Hitx", "PDc3Hity", "PDc3Hitz"])
        proKeysRec.extend(["PFtof1aTime", "PFtof1bTime", "PFtof2Time", "PCtofTime"])
        proKeysRec.extend(["PFtof1bSector", "PFtof1bComponent"])
        # proKeysRec.extend(["Pchi2pid", "Pchi2track", "PNDFtrack"])
        gamKeysRec = ["Gpx", "Gpy", "Gpz", "Gedep", "GcX", "GcY", "Gsector"]
        gamKeysRec.extend(["GcalU1", "GcalV1", "GcalW1", "Gbeta"])
        gamKeysRec.extend(["GcalX1", "GcalY1"])
        gamKeysRec.extend(["GcalX3", "GcalY3"])

        if detRes:
            eleKeysRec.extend(["Evx", "Evy"])
            # eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
            # eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
            # eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
            eleKeysRec.extend(["EcalU2", "EcalV2", "EcalW2"])
            eleKeysRec.extend(["EcalU3", "EcalV3", "EcalW3"])
            eleKeysRec.extend(["EcalHx2", "EcalHy2"])
            # eleKeysRec.extend(["Enphe"])
            eleKeysRec.extend(["EhtccZ"])
            gamKeysRec.extend(["Gedep1", "Gedep2", "Gedep3"])
            # gamKeysRec.extend(["GcalU1", "GcalV1", "GcalW1"])
            gamKeysRec.extend(["GcalU2", "GcalV2", "GcalW2"])
            gamKeysRec.extend(["GcalU3", "GcalV3", "GcalW3"])
            gamKeysRec.extend(["GcalX2", "GcalY2"])
            # gamKeysRec.extend(["Gbeta"])
            # proKeysRec.extend(["Pvz"])
            proKeysRec.extend(["PCvt1Hitx", "PCvt1Hity", "PCvt1Hitz", "PCvt3Hitx", "PCvt3Hity", "PCvt3Hitz", "PCvt5Hitx", "PCvt5Hity", "PCvt5Hitz", "PCvt7Hitx", "PCvt7Hity", "PCvt7Hitz"])
            # proKeysRec.extend(["PDc2Hitx", "PDc2Hity", "PDc2Hitz", "PDc3Hitx", "PDc3Hity", "PDc3Hitz"])
            eleKeysRec.extend(["startTime"])
            proKeysRec.extend(["PFtof1aTime", "PFtof1bTime", "PFtof2Time", "PCtofTime"])
            proKeysRec.extend(["PFtof1aHitx", "PFtof1bHitx", "PFtof2Hitx", "PCtofHitx"])
            proKeysRec.extend(["PFtof1aHity", "PFtof1bHity", "PFtof2Hity", "PCtofHity"])
            proKeysRec.extend(["PFtof1aHitz", "PFtof1bHitz", "PFtof2Hitz", "PCtofHitz"])
            proKeysRec.extend(["Pchi2track", "PNDFtrack"])

        # read them
        for key in eleKeysRec:
            df_electronRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in proKeysRec:
            df_protonRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in gamKeysRec:
            df_gammaRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        if logistics:
            df_logisticsRec = pd.DataFrame()
            logKeysRec = ["TriggerBit", "EventNum", "RunNum", "beamQ", "liveTime", "helicity"]
            for key in logKeysRec:
                df_logisticsRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
            df_logisticsRec.loc[:,'event'] = df_logisticsRec.index
        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        df_gammaRec = df_gammaRec.astype({"Gpx": float, "Gpy": float, "Gpz": float, "Gedep": float, "GcX": float, "GcY": float})
        ele = [df_electronRec['Epx'], df_electronRec['Epy'], df_electronRec['Epz']]
        df_electronRec.loc[:, 'Ep'] = mag(ele)
        df_electronRec.loc[:, 'Ee'] = getEnergy(ele, me)
        df_electronRec.loc[:, 'Etheta'] = getTheta(ele)
        df_electronRec.loc[:, 'Ephi'] = getPhi(ele)
        df_electronRec.loc[:,'ESamplFrac'] = df_electronRec.Eedep/ df_electronRec.Ep

        #proton momentum preparation for the fiducial cut
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)

        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index.get_level_values('entry')
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'event'] = df_gammaRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'GIndex'] = df_gammaRec.index.get_level_values('subentry')

        #create df_gg for pi0 exclusion
        gam = [df_gammaRec['Gpx'], df_gammaRec['Gpy'], df_gammaRec['Gpz']]
        df_gammaRec.loc[:, 'Gp'] = mag(gam)
        df_gammaRec.loc[:, 'Ge'] = getEnergy(gam, 0)
        df_gammaRec.loc[:, 'Gtheta'] = getTheta(gam)
        df_gammaRec.loc[:, 'Gphi'] = getPhi(gam)
        df_gammaRec.loc[:,'GSamplFrac'] = df_gammaRec.Gedep/ df_gammaRec.Gp

        #apply fiducial cuts
        print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))
        if nofid:
            # skip the fiducial cuts
            df_electronRec.loc[:, "EFid"] = 1
            df_protonRec.loc[:, "PFid"] = 1
            df_gammaRec.loc[:, "GFid"] = 1
            df_protonRec.loc[:, "PCvt12theta"] = -100000
            df_protonRec.loc[:, "PCvt12phi"] = -100000
            df_protonRec.loc[df_protonRec.Psector > 7, "PCvt12theta"] = getTheta([df_protonRec.loc[df_protonRec.Psector > 7].PCvt12Hitx, df_protonRec.loc[df_protonRec.Psector > 7].PCvt12Hity, df_protonRec.loc[df_protonRec.Psector > 7].PCvt12Hitz])
            df_protonRec.loc[df_protonRec.Psector > 7, "PCvt12phi"] = getPhi([df_protonRec.loc[df_protonRec.Psector > 7].PCvt12Hitx, df_protonRec.loc[df_protonRec.Psector > 7].PCvt12Hity, df_protonRec.loc[df_protonRec.Psector > 7].PCvt12Hitz])
        else:
            # perform the fiducial cuts
            df_electronRec = electronFiducial(df_electronRec, mc = False, fidlevel = fidlevel)
            df_protonRec = protonFiducial(df_protonRec, fidlevel = fidlevel)
            df_gammaRec = gammaFiducial(df_gammaRec, fidlevel = fidlevel)
            print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))
            coincidence = reduce(np.intersect1d, (df_electronRec.event, df_protonRec.event, df_gammaRec.event))
            df_electronRec = df_electronRec.loc[df_electronRec.event.isin(coincidence), :]
            df_protonRec = df_protonRec.loc[df_protonRec.event.isin(coincidence), :]
            df_gammaRec = df_gammaRec.loc[df_gammaRec.event.isin(coincidence), :]
            print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))

        # Done with the fiducial cuts.

        # Post-processing
        # e1: Electron correction  - exp only
        # e2: Electron smearing    - mc only 
        # p1: Proton energy loss   - exp, mc both
        # p2: Proton correction    - exp only
        # p3: Proton smearing      - mc only
        # g1: Gamma correction     - exp only
        # g2: Gamma smearing       - mc only
        if nocorr:
            print("no correction applied")
            pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
            df_protonRec.loc[:, 'Pp'] = mag(pro)
            df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
            df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
            df_protonRec.loc[:, 'Pphi'] = getPhi(pro)
            df_gg = pd.merge(df_gammaRec, df_gammaRec,
                             how='inner', on='event', suffixes=("", "2"))
            df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
            df_gg = df_gg.drop(['GIndex', 'GIndex2'], axis = 1)
        else:
            #e1
            df_electronRec = electronMomentumCorrection(pol, df_electronRec)
            # #e2
            # df_electronRec = electronMomentumSmearing(df_electronRec)
            #p1
            if not args.noeloss:
                df_protonRec = protonEnergyLossCorr(pol, df_protonRec)
            # #p2
            df_protonRec = protonMomentumCorrection(pol, df_protonRec)
            # #p3
            # df_protonRec = protonMomentumSmearing(pol, df_protonRec, smearing = smearing)
            # #g2
            # df_gammaRec  = gammaMomentumSmearing(df_gammaRec, smearing = smearing)
            df_gg = pd.merge(df_gammaRec, df_gammaRec,
                             how='inner', on='event', suffixes=("", "2"))
            df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
            df_gg = df_gg.drop(['GIndex', 'GIndex2'], axis = 1)
            #g1
            df_gg, df_gammaRec  = gammaMomentumCorrection(pol, df_gg, df_gammaRec)

        if detRes:
            df_protonRec.loc[:, "PDc3theta"] = -100000

            df_electronRec.loc[:, "EDc1theta"] = getTheta([df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity, df_electronRec.EDc1Hitz])
            df_electronRec.loc[:, "EDc3theta"] = getTheta([df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, df_electronRec.EDc3Hitz])
            df_electronRec.loc[:, "EAngleDiff"] = df_electronRec.loc[:, "EDc3theta"] - df_electronRec.loc[:, "EDc1theta"]

            df_protonRec.loc[:, "PCvt1r"] = -100000
            df_protonRec.loc[:, "PCvt1theta"] = -100000
            df_protonRec.loc[:, "PCvt1phi"] = -100000
            df_protonRec.loc[:, "PCvt3r"] = -100000
            df_protonRec.loc[:, "PCvt3theta"] = -100000
            df_protonRec.loc[:, "PCvt3phi"] = -100000
            df_protonRec.loc[:, "PCvt5r"] = -100000
            df_protonRec.loc[:, "PCvt5theta"] = -100000
            df_protonRec.loc[:, "PCvt5phi"] = -100000
            df_protonRec.loc[:, "PCvt7r"] = -100000
            df_protonRec.loc[:, "PCvt7theta"] = -100000
            df_protonRec.loc[:, "PCvt7phi"] = -100000
            df_protonRec.loc[:, "PCvt12r"] = -100000
            df_protonRec.loc[:, "PCvt12theta"] = -100000
            df_protonRec.loc[:, "PCvt12phi"] = -100000
            cutCD = df_protonRec.Psector>7
            df_protonRec.loc[cutCD, "PCvt1r"] = mag([df_protonRec.loc[cutCD].PCvt1Hitx, df_protonRec.loc[cutCD].PCvt1Hity, df_protonRec.loc[cutCD].PCvt1Hitz])
            df_protonRec.loc[cutCD, "PCvt1theta"] = getTheta([df_protonRec.loc[cutCD].PCvt1Hitx, df_protonRec.loc[cutCD].PCvt1Hity, df_protonRec.loc[cutCD].PCvt1Hitz])
            df_protonRec.loc[cutCD, "PCvt1phi"] = getPhi([df_protonRec.loc[cutCD].PCvt1Hitx, df_protonRec.loc[cutCD].PCvt1Hity, df_protonRec.loc[cutCD].PCvt1Hitz])
            df_protonRec.loc[cutCD, "PCvt3r"] = mag([df_protonRec.loc[cutCD].PCvt3Hitx, df_protonRec.loc[cutCD].PCvt3Hity, df_protonRec.loc[cutCD].PCvt3Hitz])
            df_protonRec.loc[cutCD, "PCvt3theta"] = getTheta([df_protonRec.loc[cutCD].PCvt3Hitx, df_protonRec.loc[cutCD].PCvt3Hity, df_protonRec.loc[cutCD].PCvt3Hitz])
            df_protonRec.loc[cutCD, "PCvt3phi"] = getPhi([df_protonRec.loc[cutCD].PCvt3Hitx, df_protonRec.loc[cutCD].PCvt3Hity, df_protonRec.loc[cutCD].PCvt3Hitz])
            df_protonRec.loc[cutCD, "PCvt5r"] = mag([df_protonRec.loc[cutCD].PCvt5Hitx, df_protonRec.loc[cutCD].PCvt5Hity, df_protonRec.loc[cutCD].PCvt5Hitz])
            df_protonRec.loc[cutCD, "PCvt5theta"] = getTheta([df_protonRec.loc[cutCD].PCvt5Hitx, df_protonRec.loc[cutCD].PCvt5Hity, df_protonRec.loc[cutCD].PCvt5Hitz])
            df_protonRec.loc[cutCD, "PCvt5phi"] = getPhi([df_protonRec.loc[cutCD].PCvt5Hitx, df_protonRec.loc[cutCD].PCvt5Hity, df_protonRec.loc[cutCD].PCvt5Hitz])
            df_protonRec.loc[cutCD, "PCvt7r"] = mag([df_protonRec.loc[cutCD].PCvt7Hitx, df_protonRec.loc[cutCD].PCvt7Hity, df_protonRec.loc[cutCD].PCvt7Hitz])
            df_protonRec.loc[cutCD, "PCvt7theta"] = getTheta([df_protonRec.loc[cutCD].PCvt7Hitx, df_protonRec.loc[cutCD].PCvt7Hity, df_protonRec.loc[cutCD].PCvt7Hitz])
            df_protonRec.loc[cutCD, "PCvt7phi"] = getPhi([df_protonRec.loc[cutCD].PCvt7Hitx, df_protonRec.loc[cutCD].PCvt7Hity, df_protonRec.loc[cutCD].PCvt7Hitz])
            df_protonRec.loc[cutCD, "PCvt12r"] = mag([df_protonRec.loc[cutCD].PCvt12Hitx, df_protonRec.loc[cutCD].PCvt12Hity, df_protonRec.loc[cutCD].PCvt12Hitz])
            df_protonRec.loc[cutCD, "PCvt12theta"] = getTheta([df_protonRec.loc[cutCD].PCvt12Hitx, df_protonRec.loc[cutCD].PCvt12Hity, df_protonRec.loc[cutCD].PCvt12Hitz])
            df_protonRec.loc[cutCD, "PCvt12phi"] = getPhi([df_protonRec.loc[cutCD].PCvt12Hitx, df_protonRec.loc[cutCD].PCvt12Hity, df_protonRec.loc[cutCD].PCvt12Hitz])


        #moduli proton phi
        df_protonRec.loc[:, "Pphi"] = np.where(df_protonRec.loc[:, "Pphi"]%360<180, df_protonRec.loc[:, "Pphi"]%360, df_protonRec.loc[:, "Pphi"]%360-360)

        df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
        df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
        df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

        df_gg.loc[:, "Gpx"] = df_gg.loc[:, "Gp"]*np.sin(np.radians(df_gg.loc[:, "Gtheta"]))*np.cos(np.radians(df_gg.loc[:, "Gphi"]))
        df_gg.loc[:, "Gpy"] = df_gg.loc[:, "Gp"]*np.sin(np.radians(df_gg.loc[:, "Gtheta"]))*np.sin(np.radians(df_gg.loc[:, "Gphi"]))
        df_gg.loc[:, "Gpz"] = df_gg.loc[:, "Gp"]*np.cos(np.radians(df_gg.loc[:, "Gtheta"]))
        df_gg.loc[:,'GSamplFrac'] = df_gg.Gedep/ df_gg.Gp

        df_gammaRec.loc[:, "Gpx"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.cos(np.radians(df_gammaRec.loc[:, "Gphi"]))
        df_gammaRec.loc[:, "Gpy"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.sin(np.radians(df_gammaRec.loc[:, "Gphi"]))
        df_gammaRec.loc[:, "Gpz"] = df_gammaRec.loc[:, "Gp"]*np.cos(np.radians(df_gammaRec.loc[:, "Gtheta"]))
        df_gammaRec.loc[:,'GSamplFrac'] = df_gammaRec.Gedep/ df_gammaRec.Gp

        if detRes:
            df_gg = df_gg.loc[:, ~df_gg.columns.duplicated()]
            df_gg.loc[:, "Gedep2_tot"] = df_gg.Gedep12 + df_gg.Gedep22 + df_gg.Gedep32
        else:
            # df_protonRec = df_protonRec.drop(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PDc1theta", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"], axis = 1)
            df_gammaRec = df_gammaRec.drop(["GcX", "GcY"], axis = 1)
            df_gg = df_gg.drop(["GcX", "GcY", "GcX2", "GcY2"], axis = 1)
        
        df_ep = pd.merge(df_electronRec, df_protonRec, how='inner', on='event')
        if logistics:
            df_ep = pd.merge(df_ep, df_logisticsRec, how='inner', on='event')

        df_epgg = pd.merge(df_ep, df_gg, how='inner', on='event')
        df_epgg = df_epgg.loc[~np.isnan(df_epgg["Ppx"]), :]
        df_epgg = df_epgg.loc[~np.isnan(df_epgg["Gpx"]), :]
        df_epgg = df_epgg.loc[~np.isnan(df_epgg["Gpx2"]), :]

        self.df_epgg = df_epgg # saves df_epgg

        df_epg = pd.merge(df_ep, df_gammaRec, how='inner', on='event')
        df_epg = df_epg.loc[~np.isnan(df_epg["Ppx"]), :]
        df_epg = df_epg.loc[~np.isnan(df_epg["Gpx"]), :]

        self.df_epg = df_epg # saves df_epgg

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
        df_epgg.loc[:, 'Pp'] = mag(pro)
        df_epgg.loc[:, 'Pe'] = getEnergy(pro, M)
        df_epgg.loc[:, 'Ptheta'] = getTheta(pro)
        df_epgg.loc[:, 'Pphi'] = getPhi(pro)

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
        VGS = [-df_epgg['Epx'], -df_epgg['Epy'], self.pbeam - df_epgg['Epz']]
        v3l = cross(self.beam, ele)
        v3h = cross(pro, VGS)
        v3g = cross(VGS, gam)
        v3pi0 = cross(VGS, pi0)

        VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
                    df_epgg["Ppy"], self.pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]
        VmissP = [-df_epgg["Epx"] - df_epgg["Gpx"] - df_epgg["Gpx2"], -df_epgg["Epy"] -
                    df_epgg["Gpy"] - df_epgg["Gpy2"], self.pbeam - df_epgg["Epz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        Vmiss = [-df_epgg["Epx"] - df_epgg["Ppx"] - df_epgg["Gpx"] - df_epgg["Gpx2"],
                    -df_epgg["Epy"] - df_epgg["Ppy"] - df_epgg["Gpy"] - df_epgg["Gpy2"],
                    self.pbeam - df_epgg["Epz"] - df_epgg["Ppz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        costheta = cosTheta(VGS, gam)

        df_epgg.loc[:, 'Mpx'], df_epgg.loc[:, 'Mpy'], df_epgg.loc[:, 'Mpz'] = Vmiss

        # binning kinematics
        df_epgg.loc[:,'Q2'] = -((self.ebeam - df_epgg['Ee'])**2 - mag2(VGS))
        df_epgg.loc[:,'nu'] = (self.ebeam - df_epgg['Ee'])
        df_epgg.loc[:,'xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
        df_epgg.loc[:,'t1'] = 2 * M * (df_epgg['Pe'] - M)
        # df_epgg.loc[:,'t1Orig'] = 2 * M * (df_epgg['PeOrig'] - M)
        df_epgg.loc[:,'t2'] = (M * df_epgg['Q2'] + 2 * M * df_epgg['nu'] * (df_epgg['nu'] - np.sqrt(df_epgg['nu'] * df_epgg['nu'] + df_epgg['Q2']) * costheta))\
        / (M + df_epgg['nu'] - np.sqrt(df_epgg['nu'] * df_epgg['nu'] + df_epgg['Q2']) * costheta)
        df_epgg.loc[:,'W'] = np.sqrt(np.maximum(0, (self.ebeam + M - df_epgg['Ee'])**2 - mag2(VGS)))
        df_epgg.loc[:,'MPt'] = np.sqrt((df_epgg["Epx"] + df_epgg["Ppx"] + df_epgg["Gpx"] + df_epgg["Gpx2"])**2 +
                                 (df_epgg["Epy"] + df_epgg["Ppy"] + df_epgg["Gpy"] + df_epgg["Gpy2"])**2)
        # trento angles
        df_epgg.loc[:,'phi1'] = angle(v3l, v3h)
        df_epgg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_epgg['phi1'], df_epgg['phi1'])
        df_epgg.loc[:,'phi2'] = angle(v3l, v3g)
        df_epgg.loc[:,'phi2'] = np.where(dot(v3l, gam) <
                                  0, 360.0 - df_epgg['phi2'], df_epgg['phi2'])

        # exclusivity variables
        df_epgg.loc[:,'MM2_ep'] = (-M - self.ebeam + df_epgg["Ee"] +
                             df_epgg["Pe"])**2 - mag2(VmissPi0)
        df_epgg.loc[:,'MM2_egg'] = (-M - self.ebeam + df_epgg["Ee"] +
                             df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(VmissP)
        df_epgg.loc[:,'MM2_epgg'] = (-M - self.ebeam + df_epgg["Ee"] + df_epgg["Pe"] +
                             df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(Vmiss)
        df_epgg.loc[:,'ME_epgg'] = (M + self.ebeam - df_epgg["Ee"] - df_epgg["Pe"] - df_epgg["Ge"] - df_epgg["Ge2"])
        df_epgg.loc[:,'Mpi0'] = pi0InvMass(gam, gam2)
        df_epgg.loc[:,'reconPi'] = angle(VmissPi0, pi0)
        df_epgg.loc[:,"Pie"] = df_epgg['Ge'] + df_epgg['Ge2']
        df_epgg.loc[:,'coplanarity'] = angle(v3h, v3pi0)
        df_epgg.loc[:,'coneAngle1'] = angle(ele, gam)
        df_epgg.loc[:,'coneAngle2'] = angle(ele, gam2)
        
        df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)

        self.df_epgg = df_epgg

    def makeDVpi0P_DVCS(self, pol = "inbending", nofid = False):
        #make dvpi0 pairs
        df_dvpi0p = self.df_epgg

        #common cuts
        cut_xBupper = df_dvpi0p["xB"] < 1  # xB
        cut_xBlower = df_dvpi0p["xB"] > 0  # xB
        cut_Q2 = df_dvpi0p["Q2"] > 1  # Q2
        cut_W = df_dvpi0p["W"] > 2  # W
        cut_Ee = df_dvpi0p["Ee"] > 2  # Ee
        cut_Ge = df_dvpi0p["Ge"] > 2  # Ge
        cut_Esector = 1#(df_dvpi0p["Esector"]!=df_dvpi0p["Gsector"]) & (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector2"])
        cut_Psector = 1#~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector"])) & ~( ((df_dvpi0p["Pstat"]//10)%10>0) & df_dvpi0p["Psector"]!=df_dvpi0p["Gsector2"])
        cut_Ppmax = df_dvpi0p.Pp < 1.6  # Pp
        cut_Pthetamin = df_dvpi0p.Ptheta > 0  # Ptheta
        cut_Trigger = ((df_dvpi0p.TriggerBit & 1 << 1) > 0) | ((df_dvpi0p.TriggerBit & 1 << 2) > 0) | ((df_dvpi0p.TriggerBit & 1 << 3) > 0) | ((df_dvpi0p.TriggerBit & 1 << 4) > 0) | ((df_dvpi0p.TriggerBit & 1 << 5) > 0) | ((df_dvpi0p.TriggerBit & 1 << 6) > 0)
        # cut_Vz = np.abs(df_dvpi0p["Evz"] - df_dvpi0p["Pvz"]) < 2.5 + 2.5 / mag([df_dvpi0p["Ppx"], pi0SimInb_forDVCS["Ppy"], pi0SimInb_forDVCS["Ppz"]])
        cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin & cut_Trigger

        df_dvpi0p = df_dvpi0p[cut_common]

        if len(df_dvpi0p):
            # proton reconstruction quality
            # cut_FD_proton = (df_epgg.loc[:, "Psector"]<7) & (df_epgg.loc[:, "Ptheta"]<35)
            # cut_CD_proton = (df_epgg.loc[:, "Psector"]>7) & (df_epgg.loc[:, "Ptheta"]>45) & (df_epgg.loc[:, "Ptheta"]<65)
            # cut_proton = (cut_FD_proton)|(cut_CD_proton)
            cut_proton = 1

            CD_Ptheta_ub = CD_Ptheta_ub_nominal
            CD_Ptheta_lb = CD_Ptheta_lb_nominal
            FD_Ptheta_inb_ub = FD_Ptheta_inb_ub_nominal
            FD_Ptheta_outb_ub = FD_Ptheta_outb_ub_nominal
            FD_Ptheta_lb = FD_Ptheta_lb_nominal
            if nofid:
                CD_Ptheta_lb = 0
                FD_Ptheta_inb_ub = 90
            
            if pol == 'inbending':
                #CDFT
                cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
                cut_Psector_CDFT = df_dvpi0p.Psector>7
                cut_Ptheta1_CDFT = df_dvpi0p.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CDFT = df_dvpi0p.Ptheta>CD_Ptheta_lb
                cut_Gsector_CDFT = df_dvpi0p.Gsector>7
                cut_GFid_CDFT = df_dvpi0p.GFid==1
                cut_PFid_CDFT = df_dvpi0p.PFid==1
                cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] < 0.6  # mmep
                cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] > -0.6  # mmep
                cut_mpi01_CDFT = df_dvpi0p["Mpi0"] < 0.161  # mpi0
                cut_mpi02_CDFT = df_dvpi0p["Mpi0"] > 0.114  # mpi0
                cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] < 2.181  # mmegg
                cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] > -0.525  # mmegg
                cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] < 0.902  # meepgg
                cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] > -0.914  # meepgg
                cut_mpt_CDFT = df_dvpi0p["MPt"] < 0.205  # mpt
                cut_recon_CDFT = df_dvpi0p["reconPi"] < 1.648  # recon gam angle
                cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] < 14.444  # coplanarity angle
                cut_mmepgg1_CDFT = df_dvpi0p["MM2_epgg"] < 0.0401  # mmepgg
                cut_mmepgg2_CDFT = df_dvpi0p["MM2_epgg"] > -0.0438  # mmepgg

                cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT &
                            cut_PFid_CDFT & cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
                            cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
                            cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


                #CD
                cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
                cut_Psector_CD = df_dvpi0p.Psector>7
                cut_Ptheta1_CD = df_dvpi0p.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CD = df_dvpi0p.Ptheta>CD_Ptheta_lb
                cut_Gsector_CD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
                cut_GFid_CD = df_dvpi0p.GFid==1
                cut_PFid_CD = df_dvpi0p.PFid==1
                cut_mmep1_CD = df_dvpi0p["MM2_ep"] < 0.6  # mmep
                cut_mmep2_CD = df_dvpi0p["MM2_ep"] > -0.6  # mmep
                cut_mpi01_CD = df_dvpi0p["Mpi0"] < 0.166  # mpi0
                cut_mpi02_CD = df_dvpi0p["Mpi0"] > 0.104  # mpi0
                cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < 2.112  # mmegg
                cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > -0.335  # mmegg
                cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < 0.882  # meepgg
                cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > -0.853  # meepgg
                cut_mpt_CD = df_dvpi0p["MPt"] < 0.177  # mpt
                cut_recon_CD = df_dvpi0p["reconPi"] < 1.111  # recon gam angle
                cut_coplanarity_CD = df_dvpi0p["coplanarity"] < 9.719  # coplanarity angle
                cut_mmepgg1_CD = df_dvpi0p["MM2_epgg"] < 0.0266  # mmepgg
                cut_mmepgg2_CD = df_dvpi0p["MM2_epgg"] > -0.0311
                  # mmepgg

                cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_GFid_CD &
                            cut_PFid_CD & cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
                            cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
                            cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

                #FD
                cut_Pp1_FD = df_dvpi0p.Pp > 0.42  # Pp
                cut_Psector_FD = df_dvpi0p.Psector<7
                cut_Ptheta1_FD = df_dvpi0p.Ptheta<FD_Ptheta_inb_ub
                cut_Ptheta2_FD = df_dvpi0p.Ptheta>FD_Ptheta_lb
                cut_Gsector_FD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
                cut_GFid_FD = df_dvpi0p.GFid==1
                cut_PFid_FD = df_dvpi0p.PFid==1
                cut_mmep1_FD = df_dvpi0p["MM2_ep"] < 0.6  # mmep
                cut_mmep2_FD = df_dvpi0p["MM2_ep"] > -0.6  # mmep
                cut_mpi01_FD = df_dvpi0p["Mpi0"] < 0.166  # mpi0
                cut_mpi02_FD = df_dvpi0p["Mpi0"] > 0.105  # mpi0
                cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < 1.806  # mmegg
                cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > -0.0905  # mmegg
                cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < 0.812  # meepgg
                cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > -0.804  # meepgg
                cut_mpt_FD = df_dvpi0p["MPt"] < 0.171  # mpt
                cut_recon_FD = df_dvpi0p["reconPi"] < 1.024  # recon gam angle
                cut_coplanarity_FD = df_dvpi0p["coplanarity"] < 9.519  # coplanarity angle
                cut_mmepgg1_FD = df_dvpi0p["MM2_epgg"] < 0.0242  # mmepgg
                cut_mmepgg2_FD = df_dvpi0p["MM2_epgg"] > -0.0279  # mmepgg

                cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_GFid_FD &
                            cut_PFid_FD & cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
                            cut_mmegg1_FD & cut_mmegg2_FD & cut_meepgg1_FD & cut_meepgg2_FD &
                            cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepgg1_FD & cut_mmepgg2_FD)

            elif pol == "outbending":
                #CDFT
                cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
                cut_Psector_CDFT = df_dvpi0p.Psector>7
                cut_Ptheta1_CDFT = df_dvpi0p.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CDFT = df_dvpi0p.Ptheta>CD_Ptheta_lb
                cut_Gsector_CDFT = df_dvpi0p.Gsector>7
                cut_GFid_CDFT = df_dvpi0p.GFid==1
                cut_PFid_CDFT = df_dvpi0p.PFid==1
                cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] < 0.6  # mmep
                cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] > -0.6  # mmep
                cut_mpi01_CDFT = df_dvpi0p["Mpi0"] < 0.161  # mpi0
                cut_mpi02_CDFT = df_dvpi0p["Mpi0"] > 0.114  # mpi0
                cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] < 2.151  # mmegg
                cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] > -0.390  # mmegg
                cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] < 0.882  # meepgg
                cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] > -0.861  # meepgg
                cut_mpt_CDFT = df_dvpi0p["MPt"] < 0.220  # mpt
                cut_recon_CDFT = df_dvpi0p["reconPi"] < 1.150  # recon gam angle
                cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] < 13.571  # coplanarity angle
                cut_mmepgg1_CDFT = df_dvpi0p["MM2_epgg"] < 0.0413  # mmepgg
                cut_mmepgg2_CDFT = df_dvpi0p["MM2_epgg"] > -0.0440  # mmepgg

                cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT &
                            cut_PFid_CDFT & cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
                            cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
                            cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


                #CD
                cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
                cut_Psector_CD = df_dvpi0p.Psector>7
                cut_Ptheta1_CD = df_dvpi0p.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CD = df_dvpi0p.Ptheta>CD_Ptheta_lb
                cut_Gsector_CD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
                cut_GFid_CD = df_dvpi0p.GFid==1
                cut_PFid_CD = df_dvpi0p.PFid==1
                cut_mmep1_CD = df_dvpi0p["MM2_ep"] < 0.6  # mmep
                cut_mmep2_CD = df_dvpi0p["MM2_ep"] > -0.6  # mmep
                cut_mpi01_CD = df_dvpi0p["Mpi0"] < 0.168  # mpi0
                cut_mpi02_CD = df_dvpi0p["Mpi0"] > 0.103  # mpi0
                cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < 2.079  # mmegg
                cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > -0.348  # mmegg
                cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < 0.814  # meepgg
                cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > -0.810  # meepgg
                cut_mpt_CD = df_dvpi0p["MPt"] < 0.202  # mpt
                cut_recon_CD = df_dvpi0p["reconPi"] < 1.287  # recon gam angle
                cut_coplanarity_CD = df_dvpi0p["coplanarity"] < 8.643  # coplanarity angle
                cut_mmepgg1_CD = df_dvpi0p["MM2_epgg"] < 0.0257  # mmepgg
                cut_mmepgg2_CD = df_dvpi0p["MM2_epgg"] > -0.0302  # mmepgg

                cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_GFid_CD &
                            cut_PFid_CD & cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
                            cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
                            cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

                #FD
                cut_Pp1_FD = df_dvpi0p.Pp > 0.5  # Pp
                cut_Psector_FD = df_dvpi0p.Psector<7
                cut_Ptheta1_FD = df_dvpi0p.Ptheta<FD_Ptheta_outb_ub
                cut_Ptheta2_FD = df_dvpi0p.Ptheta>FD_Ptheta_lb
                cut_Gsector_FD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
                cut_GFid_FD = df_dvpi0p.GFid==1
                cut_PFid_FD = df_dvpi0p.PFid==1
                cut_mmep1_FD = df_dvpi0p["MM2_ep"] < 0.6  # mmep
                cut_mmep2_FD = df_dvpi0p["MM2_ep"] > -0.6  # mmep
                cut_mpi01_FD = df_dvpi0p["Mpi0"] < 0.166  # mpi0
                cut_mpi02_FD = df_dvpi0p["Mpi0"] > 0.106  # mpi0
                cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < 1.888  # mmegg
                cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > -0.191  # mmegg
                cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < 0.864  # meepgg
                cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > -0.801  # meepgg
                cut_mpt_FD = df_dvpi0p["MPt"] < 0.195  # mpt
                cut_recon_FD = df_dvpi0p["reconPi"] < 1.299  # recon gam angle
                cut_coplanarity_FD = df_dvpi0p["coplanarity"] < 11.316  # coplanarity angle
                cut_mmepgg1_FD = df_dvpi0p["MM2_epgg"] < 0.0354  # mmepgg
                cut_mmepgg2_FD = df_dvpi0p["MM2_epgg"] > -0.0480  # mmepgg

                cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_GFid_FD &
                            cut_PFid_FD & cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
                            cut_mmegg1_FD & cut_mmegg2_FD & cut_meepgg1_FD & cut_meepgg2_FD &
                            cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepgg1_FD & cut_mmepgg2_FD)

            df_dvpi0p.loc[cut_CDFT, "config"] = 3
            df_dvpi0p.loc[cut_CD, "config"] = 2
            df_dvpi0p.loc[cut_FD, "config"] = 1

            df_dvpi0p = df_dvpi0p[df_dvpi0p.config>0]
    
        self.df_dvpi0p = df_dvpi0p #no need to reduce duplicates of pi0. remove the event if any.

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

        VGS = [-df_epg['Epx'], -df_epg['Epy'], self.pbeam - df_epg['Epz']]
        v3l = cross(self.beam, ele)
        v3h = cross(pro, VGS)
        v3g = cross(VGS, gam)
        VmissG = [-df_epg["Epx"] - df_epg["Ppx"], -df_epg["Epy"] - df_epg["Ppy"],
                  self.pbeam - df_epg["Epz"] - df_epg["Ppz"]]
        VmissP = [-(df_epg["Epx"] + df_epg["Gpx"]), -(df_epg["Epy"] + df_epg["Gpy"]),
                  -(-self.pbeam + df_epg["Epz"] + df_epg["Gpz"])]
        Vmiss = [-(df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"]), -(df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"]),
                 -(-self.pbeam + df_epg["Epz"] + df_epg["Ppz"] + df_epg["Gpz"])]
        costheta = cosTheta(VGS, gam)

        df_epg.loc[:, 'Mpx'], df_epg.loc[:, 'Mpy'], df_epg.loc[:, 'Mpz'] = Vmiss

        # binning kinematics
        df_epg.loc[:,'Q2'] = -((self.ebeam - df_epg['Ee'])**2 - mag2(VGS))
        df_epg.loc[:,'nu'] = (self.ebeam - df_epg['Ee'])
        df_epg.loc[:,'y'] = df_epg['nu']/self.ebeam
        df_epg.loc[:,'xB'] = df_epg['Q2'] / 2.0 / M / df_epg['nu']
        df_epg.loc[:,'t1'] = 2 * M * (df_epg['Pe'] - M)
        # df_epg.loc[:,'t1Orig'] = 2 * M * (df_epg['PeOrig'] - M)
        df_epg.loc[:,'t2'] = (M * df_epg['Q2'] + 2 * M * df_epg['nu'] * (df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta))\
        / (M + df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta)
        df_epg.loc[:,'W'] = np.sqrt(np.maximum(0, (self.ebeam + M - df_epg['Ee'])**2 - mag2(VGS)))

        # trento angles
        df_epg.loc[:,'phi1'] = angle(v3l, v3h)
        df_epg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_epg['phi1'], df_epg['phi1'])
        df_epg.loc[:,'phi2'] = angle(v3l, v3g)
        df_epg.loc[:,'phi2'] = np.where(dot(v3l, gam) <
                                  0, 360.0 - df_epg['phi2'], df_epg['phi2'])

        # exclusivity variables
        df_epg.loc[:,'MM2_epg'] = (-M - self.ebeam + df_epg["Ee"] +
                             df_epg["Pe"] + df_epg["Ge"])**2 - mag2(Vmiss)
        df_epg.loc[:,'ME_epg'] = (M + self.ebeam - df_epg["Ee"] - df_epg["Pe"] - df_epg["Ge"])
        df_epg.loc[:,'MM2_ep'] = (-M - self.ebeam + df_epg["Ee"] + df_epg["Pe"])**2 - mag2(VmissG)
        df_epg.loc[:,'MM2_eg'] = (-M - self.ebeam + df_epg["Ee"] + df_epg["Ge"])**2 - mag2(VmissP)
        df_epg.loc[:,'MPt'] = np.sqrt((df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"])**2 +
                                (df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"])**2)
        df_epg.loc[:,'coneAngle'] = angle(ele, gam)
        df_epg.loc[:,'reconGam'] = angle(gam, VmissG)
        df_epg.loc[:,'coplanarity'] = angle(v3h, v3g)

        df_epg.loc[:,'closeness2'] = np.abs(df_epg.MM2_ep)

        eps = 2*M*df_epg.xB / np.sqrt(df_epg.Q2)
        df_epg.loc[:,'ycol1'] = (df_epg.Q2-df_epg.t2)/(df_epg.Q2-df_epg.xB*df_epg.t2)
        df_epg.loc[:,'ycol2'] = 1 - (1-df_epg.xB)*df_epg.t2/df_epg.Q2
        df_epg.loc[:,'ymax1'] = 2*(np.sqrt(1+eps**2)-1)/(eps**2)
        df_epg.loc[:,'ymax2'] = 1 - (M**2)*(df_epg.xB**2)/df_epg.Q2
        df_epg.loc[:,'tmin1'] = df_epg.Q2*(2*(1-df_epg.xB)*(1-np.sqrt(1+eps**2))+eps**2)/(4*df_epg.xB*(1-df_epg.xB) + eps**2)
        df_epg.loc[:,'tmin2'] = M*M*(df_epg.xB**2)/(1-df_epg.xB+df_epg.xB*M*M/df_epg.Q2)
        df_epg.loc[:,'tcol'] = df_epg.Q2*(df_epg.Q2-2*df_epg.xB*M*self.ebeam)/df_epg.xB/(df_epg.Q2-2*M*self.ebeam)

        df_epg.loc[:, 'vzdiff'] = df_epg.Evz - df_epg.Pvz

        self.df_epg = df_epg

    def makeDVCS(self, pol = "inbending", nofid = False, allowsamesector = False, allowduplicates = False):
        #make dvcs pairs
        df_dvcs = self.df_epg

        #common cuts
        cut_xBupper = df_dvcs["xB"] < 1  # xB
        cut_xBlower = df_dvcs["xB"] > 0  # xB
        cut_Q2 = df_dvcs["Q2"] > 1  # Q2
        cut_W = df_dvcs["W"] > 2  # W
        cut_Ee = df_dvcs["Ee"] > 2  # Ee
        cut_Ge = df_dvcs["Ge"] > 2  # Ge
        if allowsamesector:
            cut_Esector = 1
            cut_Psector = 1
        else:
            cut_Esector = (df_dvcs["Esector"]!=df_dvcs["Gsector"])
            cut_Psector = ~( ((df_dvcs["Pstat"]//10)%10>0) & (df_dvcs["Psector"]==df_dvcs["Gsector"]))
        cut_Ppmax = df_dvcs.Pp < 1.6  # Pp
        cut_Pthetamin = df_dvcs.Ptheta > 0 #Ptheta
        cut_Trigger = ((df_dvcs.TriggerBit & 1 << 1) > 0) | ((df_dvcs.TriggerBit & 1 << 2) > 0) | ((df_dvcs.TriggerBit & 1 << 3) > 0) | ((df_dvcs.TriggerBit & 1 << 4) > 0) | ((df_dvcs.TriggerBit & 1 << 5) > 0) | ((df_dvcs.TriggerBit & 1 << 6) > 0)
        cut_EFid = df_dvcs.EFid == 1
        # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
        cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin & cut_Trigger & cut_EFid

        df_dvcs = df_dvcs[cut_common]

        if len(df_dvcs):

            # proton reconstruction quality
            # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) & (df_dvcs.loc[:, "Ptheta"]<35)
            # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) & (df_dvcs.loc[:, "Ptheta"]>45) & (df_dvcs.loc[:, "Ptheta"]<65)
            # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) #& (df_dvcs.loc[:, "Ptheta"]<37)
            # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) #& (df_dvcs.loc[:, "Ptheta"]<66) #& (df_dvcs.loc[:, "Ptheta"]>40) 
            # cut_proton = (cut_FD_proton)|(cut_CD_proton)
            #(cut_FD_proton)|(cut_CD_proton)

            cuts_dvcs_CDFT_Inb = self.cuts_dvcs_CDFT_Inb 
            cuts_dvcs_CD_Inb = self.cuts_dvcs_CD_Inb 
            cuts_dvcs_FD_Inb = self.cuts_dvcs_FD_Inb
            cuts_dvcs_CDFT_Outb = self.cuts_dvcs_CDFT_Outb 
            cuts_dvcs_CD_Outb = self.cuts_dvcs_CD_Outb 
            cuts_dvcs_FD_Outb = self.cuts_dvcs_FD_Outb

            CD_Ptheta_ub = CD_Ptheta_ub_nominal
            CD_Ptheta_lb = CD_Ptheta_lb_nominal
            FD_Ptheta_inb_ub = FD_Ptheta_inb_ub_nominal
            FD_Ptheta_outb_ub = FD_Ptheta_outb_ub_nominal
            FD_Ptheta_lb = FD_Ptheta_lb_nominal
            if nofid:
                CD_Ptheta_lb = 0
                FD_Ptheta_inb_ub = 90

            if pol == "inbending":
                #CDFT
                cut_Pp1_CDFT = df_dvcs.Pp > 0.3  # Pp
                cut_Psector_CDFT = df_dvcs.Psector>7
                cut_Ptheta1_CDFT = df_dvcs.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CDFT = df_dvcs.Ptheta>CD_Ptheta_lb
                cut_Gsector_CDFT = df_dvcs.Gsector>7
                cut_GFid_CDFT = df_dvcs.GFid==1
                cut_PFid_CDFT = df_dvcs.PFid==1
                cut_mmep1_CDFT = df_dvcs["MM2_ep"] < cuts_dvcs_CDFT_Inb["MM2_ep_ub"]  # mmep
                cut_mmep2_CDFT = df_dvcs["MM2_ep"] > cuts_dvcs_CDFT_Inb["MM2_ep_lb"]  # mmep
                cut_mmeg1_CDFT = df_dvcs["MM2_eg"] < cuts_dvcs_CDFT_Inb["MM2_eg_ub"]  # mmeg
                cut_mmeg2_CDFT = df_dvcs["MM2_eg"] > cuts_dvcs_CDFT_Inb["MM2_eg_lb"]  # mmeg
                cut_meepg1_CDFT = df_dvcs["ME_epg"] < cuts_dvcs_CDFT_Inb["ME_epg_ub"] # meepg
                cut_meepg2_CDFT = df_dvcs["ME_epg"] > cuts_dvcs_CDFT_Inb["ME_epg_lb"]  # meepg
                cut_cone1_CDFT = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_CDFT_Inb["coneAngle_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_CDFT = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_CDFT_Inb["coneAngle_lb"])(df_dvcs.Etheta) # coneangle
                cut_mpt_CDFT = df_dvcs["MPt"] < cuts_dvcs_CDFT_Inb["MPt_ub"]  # mpt
                cut_recon_CDFT = df_dvcs["reconGam"] < cuts_dvcs_CDFT_Inb["reconGam_ub"]  # recon gam angle
                cut_coplanarity_CDFT = df_dvcs["coplanarity"] < cuts_dvcs_CDFT_Inb["coplanarity_ub"]  # coplanarity angle
                cut_mmepg1_CDFT = df_dvcs["MM2_epg"] < cuts_dvcs_CDFT_Inb["MM2_epg_ub"]  # mmepg
                cut_mmepg2_CDFT = df_dvcs["MM2_epg"] > cuts_dvcs_CDFT_Inb["MM2_epg_lb"]  # mmepg

                cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT & 
                            cut_PFid_CDFT & cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mmeg1_CDFT & cut_mmeg2_CDFT &
                            cut_meepg1_CDFT & cut_meepg2_CDFT & cut_cone1_CDFT & cut_cone2_CDFT &
                            cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepg1_CDFT & cut_mmepg2_CDFT)


                cut_cone1_CR = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_CDFT_Inb["coneAngleCR_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_CR = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_CDFT_Inb["coneAngleCR_lb"])(df_dvcs.Etheta) # coneangle

                cut_CR = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT & 
                            cut_PFid_CDFT & cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mmeg1_CDFT & cut_mmeg2_CDFT &
                            cut_meepg1_CDFT & cut_meepg2_CDFT & cut_cone1_CR & cut_cone2_CR &
                            cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepg1_CDFT & cut_mmepg2_CDFT)

                #CD
                cut_Pp1_CD = df_dvcs.Pp > 0.3  # Pp
                cut_Psector_CD = df_dvcs.Psector>7
                cut_Ptheta1_CD = df_dvcs.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CD = df_dvcs.Ptheta>CD_Ptheta_lb
                cut_Gsector_CD = (df_dvcs.Gsector<7) & (df_dvcs.Gsector>0)
                cut_GFid_CD = df_dvcs.GFid==1
                cut_PFid_CD = df_dvcs.PFid==1
                cut_mmep1_CD = df_dvcs["MM2_ep"] < cuts_dvcs_CD_Inb["MM2_ep_ub"]  # mmep
                cut_mmep2_CD = df_dvcs["MM2_ep"] > cuts_dvcs_CD_Inb["MM2_ep_lb"]  # mmep
                cut_mmeg1_CD = df_dvcs["MM2_eg"] < cuts_dvcs_CD_Inb["MM2_eg_ub"]  # mmeg
                cut_mmeg2_CD = df_dvcs["MM2_eg"] > cuts_dvcs_CD_Inb["MM2_eg_lb"]  # mmeg
                cut_meepg1_CD = df_dvcs["ME_epg"] < cuts_dvcs_CD_Inb["ME_epg_ub"] # meepg
                cut_meepg2_CD = df_dvcs["ME_epg"] > cuts_dvcs_CD_Inb["ME_epg_lb"]  # meepg
                cut_cone1_CD = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_CD_Inb["coneAngle_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_CD = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_CD_Inb["coneAngle_lb"])(df_dvcs.Etheta) # coneangle
                cut_mpt_CD = df_dvcs["MPt"] < cuts_dvcs_CD_Inb["MPt_ub"]  # mpt
                cut_recon_CD = df_dvcs["reconGam"] < cuts_dvcs_CD_Inb["reconGam_ub"]  # recon gam angle
                cut_coplanarity_CD = df_dvcs["coplanarity"] < cuts_dvcs_CD_Inb["coplanarity_ub"]  # coplanarity angle
                cut_mmepg1_CD = df_dvcs["MM2_epg"] < cuts_dvcs_CD_Inb["MM2_epg_ub"]  # mmepg
                cut_mmepg2_CD = df_dvcs["MM2_epg"] > cuts_dvcs_CD_Inb["MM2_epg_lb"]  # mmepg

                cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_GFid_CD &
                            cut_PFid_CD & cut_mmep1_CD & cut_mmep2_CD & cut_mmeg1_CD & cut_mmeg2_CD &
                            cut_meepg1_CD & cut_meepg2_CD & cut_cone1_CD & cut_cone2_CD &
                            cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepg1_CD & cut_mmepg2_CD)

                #FD
                cut_Pp1_FD = df_dvcs.Pp > 0.42  # Pp
                cut_Psector_FD = df_dvcs.Psector<7
                cut_Ptheta1_FD = df_dvcs.Ptheta<FD_Ptheta_inb_ub
                cut_Ptheta2_FD = df_dvcs.Ptheta>FD_Ptheta_lb
                cut_Gsector_FD = (df_dvcs.Gsector<7) & (df_dvcs.Gsector>0)
                cut_GFid_FD = df_dvcs.GFid==1
                cut_PFid_FD = df_dvcs.PFid==1
                cut_mmep1_FD = df_dvcs["MM2_ep"] < cuts_dvcs_FD_Inb["MM2_ep_ub"]  # mmep
                cut_mmep2_FD = df_dvcs["MM2_ep"] > cuts_dvcs_FD_Inb["MM2_ep_lb"]  # mmep
                cut_mmeg1_FD = df_dvcs["MM2_eg"] < cuts_dvcs_FD_Inb["MM2_eg_ub"]  # mmeg
                cut_mmeg2_FD = df_dvcs["MM2_eg"] > cuts_dvcs_FD_Inb["MM2_eg_lb"]  # mmeg
                cut_meepg1_FD = df_dvcs["ME_epg"] < cuts_dvcs_FD_Inb["ME_epg_ub"] # meepg
                cut_meepg2_FD = df_dvcs["ME_epg"] > cuts_dvcs_FD_Inb["ME_epg_lb"]  # meepg
                cut_cone1_FD = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_FD_Inb["coneAngle_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_FD = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_FD_Inb["coneAngle_lb"])(df_dvcs.Etheta) # coneangle
                cut_mpt_FD = df_dvcs["MPt"] < cuts_dvcs_FD_Inb["MPt_ub"]  # mpt
                cut_recon_FD = df_dvcs["reconGam"] < cuts_dvcs_FD_Inb["reconGam_ub"]  # recon gam angle
                cut_coplanarity_FD = df_dvcs["coplanarity"] < cuts_dvcs_FD_Inb["coplanarity_ub"]  # coplanarity angle
                cut_mmepg1_FD = df_dvcs["MM2_epg"] < cuts_dvcs_FD_Inb["MM2_epg_ub"]  # mmepg
                cut_mmepg2_FD = df_dvcs["MM2_epg"] > cuts_dvcs_FD_Inb["MM2_epg_lb"]  # mmepg

                cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_GFid_FD &
                            cut_PFid_FD & cut_mmep1_FD & cut_mmep2_FD & cut_mmeg1_FD & cut_mmeg2_FD &
                            cut_meepg1_FD & cut_meepg2_FD & cut_cone1_FD & cut_cone2_FD &
                            cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepg1_FD & cut_mmepg2_FD)

            elif pol == "outbending":
                #CDFT
                cut_Pp1_CDFT = df_dvcs.Pp > 0.3  # Pp
                cut_Psector_CDFT = df_dvcs.Psector>7
                cut_Ptheta1_CDFT = df_dvcs.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CDFT = df_dvcs.Ptheta>CD_Ptheta_lb
                cut_Gsector_CDFT = df_dvcs.Gsector>7
                cut_GFid_CDFT = df_dvcs.GFid==1
                cut_PFid_CDFT = df_dvcs.PFid==1
                cut_mmep1_CDFT = df_dvcs["MM2_ep"] < cuts_dvcs_CDFT_Outb["MM2_ep_ub"]  # mmep
                cut_mmep2_CDFT = df_dvcs["MM2_ep"] > cuts_dvcs_CDFT_Outb["MM2_ep_lb"]  # mmep
                cut_mmeg1_CDFT = df_dvcs["MM2_eg"] < cuts_dvcs_CDFT_Outb["MM2_eg_ub"]  # mmeg
                cut_mmeg2_CDFT = df_dvcs["MM2_eg"] > cuts_dvcs_CDFT_Outb["MM2_eg_lb"]  # mmeg
                cut_meepg1_CDFT = df_dvcs["ME_epg"] < cuts_dvcs_CDFT_Outb["ME_epg_ub"] # meepg
                cut_meepg2_CDFT = df_dvcs["ME_epg"] > cuts_dvcs_CDFT_Outb["ME_epg_lb"]  # meepg
                cut_cone1_CDFT = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_CDFT_Outb["coneAngle_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_CDFT = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_CDFT_Outb["coneAngle_lb"])(df_dvcs.Etheta) # coneangle
                cut_mpt_CDFT = df_dvcs["MPt"] < cuts_dvcs_CDFT_Outb["MPt_ub"]  # mpt
                cut_recon_CDFT = df_dvcs["reconGam"] < cuts_dvcs_CDFT_Outb["reconGam_ub"]  # recon gam angle
                cut_coplanarity_CDFT = df_dvcs["coplanarity"] < cuts_dvcs_CDFT_Outb["coplanarity_ub"]  # coplanarity angle
                cut_mmepg1_CDFT = df_dvcs["MM2_epg"] < cuts_dvcs_CDFT_Outb["MM2_epg_ub"]  # mmepg
                cut_mmepg2_CDFT = df_dvcs["MM2_epg"] > cuts_dvcs_CDFT_Outb["MM2_epg_lb"]  # mmepg

                cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT & 
                            cut_PFid_CDFT & cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mmeg1_CDFT & cut_mmeg2_CDFT &
                            cut_meepg1_CDFT & cut_meepg2_CDFT & cut_cone1_CDFT & cut_cone2_CDFT &
                            cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepg1_CDFT & cut_mmepg2_CDFT)


                cut_cone1_CR = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_CDFT_Outb["coneAngleCR_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_CR = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_CDFT_Outb["coneAngleCR_lb"])(df_dvcs.Etheta) # coneangle

                cut_CR = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT & 
                            cut_PFid_CDFT & cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mmeg1_CDFT & cut_mmeg2_CDFT &
                            cut_meepg1_CDFT & cut_meepg2_CDFT & cut_cone1_CR & cut_cone2_CR &
                            cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepg1_CDFT & cut_mmepg2_CDFT)

                #CD
                cut_Pp1_CD = df_dvcs.Pp > 0.3  # Pp
                cut_Psector_CD = df_dvcs.Psector>7
                cut_Ptheta1_CD = df_dvcs.Ptheta<CD_Ptheta_ub
                cut_Ptheta2_CD = df_dvcs.Ptheta>CD_Ptheta_lb
                cut_Gsector_CD = (df_dvcs.Gsector<7) & (df_dvcs.Gsector>0)
                cut_GFid_CD = df_dvcs.GFid==1
                cut_PFid_CD = df_dvcs.PFid==1
                cut_mmep1_CD = df_dvcs["MM2_ep"] < cuts_dvcs_CD_Outb["MM2_ep_ub"]  # mmep
                cut_mmep2_CD = df_dvcs["MM2_ep"] > cuts_dvcs_CD_Outb["MM2_ep_lb"]  # mmep
                cut_mmeg1_CD = df_dvcs["MM2_eg"] < cuts_dvcs_CD_Outb["MM2_eg_ub"]  # mmeg
                cut_mmeg2_CD = df_dvcs["MM2_eg"] > cuts_dvcs_CD_Outb["MM2_eg_lb"]  # mmeg
                cut_meepg1_CD = df_dvcs["ME_epg"] < cuts_dvcs_CD_Outb["ME_epg_ub"] # meepg
                cut_meepg2_CD = df_dvcs["ME_epg"] > cuts_dvcs_CD_Outb["ME_epg_lb"]  # meepg
                cut_cone1_CD = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_CD_Outb["coneAngle_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_CD = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_CD_Outb["coneAngle_lb"])(df_dvcs.Etheta) # coneangle
                cut_mpt_CD = df_dvcs["MPt"] < cuts_dvcs_CD_Outb["MPt_ub"]  # mpt
                cut_recon_CD = df_dvcs["reconGam"] < cuts_dvcs_CD_Outb["reconGam_ub"]  # recon gam angle
                cut_coplanarity_CD = df_dvcs["coplanarity"] < cuts_dvcs_CD_Outb["coplanarity_ub"]  # coplanarity angle
                cut_mmepg1_CD = df_dvcs["MM2_epg"] < cuts_dvcs_CD_Outb["MM2_epg_ub"]  # mmepg
                cut_mmepg2_CD = df_dvcs["MM2_epg"] > cuts_dvcs_CD_Outb["MM2_epg_lb"]  # mmepg

                cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_GFid_CD &
                            cut_PFid_CD & cut_mmep1_CD & cut_mmep2_CD & cut_mmeg1_CD & cut_mmeg2_CD &
                            cut_meepg1_CD & cut_meepg2_CD & cut_cone1_CD & cut_cone2_CD &
                            cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepg1_CD & cut_mmepg2_CD)

                #FD
                cut_Pp1_FD = df_dvcs.Pp > 0.5  # Pp
                cut_Psector_FD = df_dvcs.Psector<7
                cut_Ptheta1_FD = df_dvcs.Ptheta<FD_Ptheta_outb_ub
                cut_Ptheta2_FD = df_dvcs.Ptheta>FD_Ptheta_lb
                cut_Gsector_FD = (df_dvcs.Gsector<7) & (df_dvcs.Gsector>0)
                cut_GFid_FD = df_dvcs.GFid==1
                cut_PFid_FD = df_dvcs.PFid==1
                cut_mmep1_FD = df_dvcs["MM2_ep"] < cuts_dvcs_FD_Outb["MM2_ep_ub"]  # mmep
                cut_mmep2_FD = df_dvcs["MM2_ep"] > cuts_dvcs_FD_Outb["MM2_ep_lb"]  # mmep
                cut_mmeg1_FD = df_dvcs["MM2_eg"] < cuts_dvcs_FD_Outb["MM2_eg_ub"]  # mmeg
                cut_mmeg2_FD = df_dvcs["MM2_eg"] > cuts_dvcs_FD_Outb["MM2_eg_lb"]  # mmeg
                cut_meepg1_FD = df_dvcs["ME_epg"] < cuts_dvcs_FD_Outb["ME_epg_ub"] # meepg
                cut_meepg2_FD = df_dvcs["ME_epg"] > cuts_dvcs_FD_Outb["ME_epg_lb"]  # meepg
                cut_cone1_FD = df_dvcs["coneAngle"] < np.poly1d(cuts_dvcs_FD_Outb["coneAngle_ub"])(df_dvcs.Etheta) # coneangle
                cut_cone2_FD = df_dvcs["coneAngle"] > np.poly1d(cuts_dvcs_FD_Outb["coneAngle_lb"])(df_dvcs.Etheta) # coneangle
                cut_mpt_FD = df_dvcs["MPt"] < cuts_dvcs_FD_Outb["MPt_ub"]  # mpt
                cut_recon_FD = df_dvcs["reconGam"] < cuts_dvcs_FD_Outb["reconGam_ub"]  # recon gam angle
                cut_coplanarity_FD = df_dvcs["coplanarity"] < cuts_dvcs_FD_Outb["coplanarity_ub"]  # coplanarity angle
                cut_mmepg1_FD = df_dvcs["MM2_epg"] < cuts_dvcs_FD_Outb["MM2_epg_ub"]  # mmepg
                cut_mmepg2_FD = df_dvcs["MM2_epg"] > cuts_dvcs_FD_Outb["MM2_epg_lb"]  # mmepg

                cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_GFid_FD &
                            cut_PFid_FD & cut_mmep1_FD & cut_mmep2_FD & cut_mmeg1_FD & cut_mmeg2_FD &
                            cut_meepg1_FD & cut_meepg2_FD & cut_cone1_FD & cut_cone2_FD &
                            cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepg1_FD & cut_mmepg2_FD)            

            df_dvcs.loc[cut_CR, "config"] = 4
            df_dvcs.loc[cut_CDFT, "config"] = 3
            df_dvcs.loc[cut_CD, "config"] = 2
            df_dvcs.loc[cut_FD, "config"] = 1

            df_dvcs = df_dvcs[df_dvcs.config>0]

            if allowduplicates:
                pass
            else:
                # dealing with duplicates
                df_dvcs = df_dvcs.sort_values(by=['reconGam', 'closeness2'], ascending = [True, True])
                df_dvcs = df_dvcs.loc[~df_dvcs.event.duplicated(), :]
                df_dvcs = df_dvcs.sort_values(by='event')

        self.df_dvcs = df_dvcs                        

    def pi02gSubtraction(self):
        #exclude dvpi0 from dvcs. use only when both set up.
        df_epg = self.df_epg
        pi0to2gammas = df_epg["event"].isin(self.df_dvpi0p["event"])
        df_epg = df_epg[~pi0to2gammas]
        self.df_epg = df_epg

    def save(self, raw = False, pol = 'inbending', efficiency = False):
        # df_x = self.df_dvcs
        # df_protonDet = self.df_protonDet
        # df = pd.merge(df_x, df_protonDet, how = 'inner', on ='event')
        if raw:
            CD_Ptheta_ub = CD_Ptheta_ub_nominal
            CD_Ptheta_lb = CD_Ptheta_lb_nominal
            FD_Ptheta_inb_ub = FD_Ptheta_inb_ub_nominal
            FD_Ptheta_outb_ub = FD_Ptheta_outb_ub_nominal
            FD_Ptheta_lb = FD_Ptheta_lb_nominal

            df_Rec = self.df_epg
            #common cuts
            cut_xBupper = df_Rec["xB"] < 1  # xB
            cut_xBlower = df_Rec["xB"] > 0  # xB
            cut_Q2 = df_Rec["Q2"] > 1  # Q2
            cut_W = df_Rec["W"] > 2  # W
            cut_Ee = df_Rec["Ee"] > 2  # Ee
            cut_Ge = df_Rec["Ge"] > 2  # Ge
            cut_Esector = (df_Rec["Esector"]!=df_Rec["Gsector"])
            cut_Psector = ~( ((df_Rec["Pstat"]//10)%10>0) & (df_Rec["Psector"]==df_Rec["Gsector"]))
            cut_Ppmax = df_Rec.Pp < 1.6  # Pp
            cut_Pthetamin = df_Rec.Ptheta > 0 # Ptheta
            # cut_Vz = np.abs(df_Rec["Evz"] - df_Rec["Pvz"]) < 2.5 + 2.5 / mag([df_Rec["Ppx"], df_Rec["Ppy"], df_Rec["Ppz"]])
            cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin

            df_Rec = df_Rec[cut_common]

            #CDFT
            cut_Pp1_CDFT = df_Rec.Pp > 0.3  # Pp
            cut_Psector_CDFT = df_Rec.Psector>7
            cut_Ptheta1_CDFT = df_Rec.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CDFT = df_Rec.Ptheta>CD_Ptheta_lb
            cut_Gsector_CDFT = df_Rec.Gsector>7
            cut_GFid_CDFT = df_Rec.GFid==1
            cut_PFid_CDFT = df_Rec.PFid==1
            #CD
            cut_Pp1_CD = df_Rec.Pp > 0.3  # Pp
            cut_Psector_CD = df_Rec.Psector>7
            cut_Ptheta1_CD = df_Rec.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CD = df_Rec.Ptheta>CD_Ptheta_lb
            cut_Gsector_CD = (df_Rec.Gsector<7)&(df_Rec.Gsector>0)
            cut_GFid_CD = df_Rec.GFid==1
            cut_PFid_CD = df_Rec.PFid==1
            #FD
            if pol == "inbending":
                cut_Pp1_FD = df_Rec.Pp > 0.42  # Pp
                cut_Ptheta1_FD = df_Rec.Ptheta<FD_Ptheta_inb_ub
            elif pol == "outbending":
                cut_Pp1_FD = df_Rec.Pp > 0.5  # Pp
                cut_Ptheta1_FD = df_Rec.Ptheta<FD_Ptheta_outb_ub
            cut_Psector_FD = df_Rec.Psector<7
            cut_Ptheta2_FD = df_Rec.Ptheta>FD_Ptheta_lb
            cut_Gsector_FD = (df_Rec.Gsector<7)&(df_Rec.Gsector>0)
            cut_GFid_FD = df_Rec.GFid==1
            cut_PFid_FD = df_Rec.PFid==1

            cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_GFid_CDFT & cut_PFid_CDFT)
            cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_GFid_CD & cut_PFid_CD)
            cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_GFid_FD & cut_PFid_FD)

            df_Rec.loc[cut_CDFT, "config"] = 3
            df_Rec.loc[cut_CD, "config"] = 2
            df_Rec.loc[cut_FD, "config"] = 1

            df_Rec = df_Rec[df_Rec.config>0]

            df_Rec = df_Rec.sort_values(by=['reconGam', 'closeness2'], ascending = [True, True])
            df_Rec = df_Rec.loc[~df_Rec.event.duplicated(), ]
            df_Rec = df_Rec.sort_values(by='event')
        else:
            df_Rec = self.df_dvcs

        df_Rec.loc[:, "integrated_binnum"] = 0

        for binnum, bin in enumerate(self.bin_scheme):
            xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin
            try:
                assert np.sum(df_Rec.loc[ (df_Rec.xB>=xBmin) & (df_Rec.xB<xBmax) & (df_Rec.Q2>=Q2min) & (df_Rec.Q2<Q2max)  & (df_Rec.t1>=tmin) & (df_Rec.t1<tmax), "integrated_binnum"] != 0) == 0        
            except:
                print("This bin overlaps with others. Check the geometry. {}".format(bin))
            df_Rec.loc[ (df_Rec.xB>=xBmin) & (df_Rec.xB<xBmax) & (df_Rec.Q2>=Q2min) & (df_Rec.Q2<Q2max)  & (df_Rec.t1>=tmin) & (df_Rec.t1<tmax), "integrated_binnum"] = binnum + 1

        for binnum, bin in enumerate(self.fringe_bin_scheme):
            xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin
            try:
                assert np.sum(df_Rec.loc[ (df_Rec.xB>=xBmin) & (df_Rec.xB<xBmax) & (df_Rec.Q2>=Q2min) & (df_Rec.Q2<Q2max)  & (df_Rec.t1>=tmin) & (df_Rec.t1<tmax), "integrated_binnum"] != 0) == 0        
            except:
                print("This bin overlaps with others. Check the geometry. {}".format(bin))
            df_Rec.loc[ (df_Rec.xB>=xBmin) & (df_Rec.xB<xBmax) & (df_Rec.Q2>=Q2min) & (df_Rec.Q2<Q2max)  & (df_Rec.t1>=tmin) & (df_Rec.t1<tmax), "integrated_binnum"] = binnum + 1 + len(self.bin_scheme)

        phibins = [-1] + list(np.linspace(0, 360, 24+1)[1:-1]) + [361]
        for phi_binnum in range(24):
            phimin = phibins[phi_binnum]
            phimax = phibins[phi_binnum+1]
            df_Rec.loc[ (df_Rec.phi1>=phimin) & (df_Rec.phi1<phimax), "phi_binnum"] = phi_binnum

        df_Rec = df_Rec.astype({"integrated_binnum": int, "config": int, "phi_binnum": int})

        if efficiency:
            df_Rec = assign_efficiency(df_Rec)

        self.df = df_Rec

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
    parser.add_argument("-binbybin","--binbybin", action = "store_true")
    parser.add_argument("-efficiency","--efficiency", action = "store_true")

    args = parser.parse_args()

    if args.entry_start:
        args.entry_start = int(args.entry_start)
    if args.entry_stop:
        args.entry_stop = int(args.entry_stop)

    be = float(args.beam)
    converter = root2pickle(args.fname, entry_start = args.entry_start,
     entry_stop = args.entry_stop, pol = args.polarity, detRes = args.detRes, raw = args.raw,
     logistics = args.logistics, width = args.width, nofid = args.nofid, nocorr = args.nocorr, noeloss = args.noeloss,nopcorr = args.nopcorr,
     fidlevel = args.fidlevel, allowsamesector = args.allowsamesector, allowduplicates = args.allowduplicates, ebeam = be, efficiency = args.efficiency)
    df = converter.df
    if args.binbybin:
        for integrated_binnum in range(len(converter.bin_scheme) + len(converter.fringe_bin_scheme) + 1):
            df.loc[ (df.integrated_binnum == integrated_binnum), :].to_pickle("{}.{}.pkl".format(args.out, integrated_binnum))
    else:
        df.to_pickle(args.out)