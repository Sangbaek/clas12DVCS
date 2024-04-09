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
        detRes = False, raw = False, width = "mid", logistics = False, nofid = False, nocorr = False, noeloss = False, nopcorr = False,
        fidlevel = 'mid', allowsamesector = False, ebeam = 10.604):
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
            saveDVpi0vars: 4 momentum algebra to save DVpi0P vars
            makeDVpi0P: select pi0->2g events.
            makeDVpi0P: select DVpi0P candidates
            save: save output
        '''
        self.fname = fname
        self.ebeam = ebeam # beam energy
        self.pbeam = np.sqrt(ebeam * ebeam - me * me) # beam electron momentum
        self.beam = [0, 0, self.pbeam] # beam vector
 
        self.determineWidth(width = width)
        self.readEPGG(entry_start = entry_start, entry_stop = entry_stop, pol = pol, 
            detRes = detRes, logistics = logistics, nofid = nofid, nocorr = nocorr, noeloss = noeloss, nopcorr =nopcorr,
            fidlevel = fidlevel)
        self.saveDVpi0vars()
        if not raw:
            self.makeDVpi0P(pol = pol, nofid = nofid, allowsamesector = allowsamesector)
        self.save(raw = raw, pol = pol)

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
            self.Ge2Threshold = Ge2Threshold_default
            self.cuts_dvpi0p_CDFT_Inb = cuts_dvpi0p_default
            self.cuts_dvpi0p_CD_Inb = cuts_dvpi0p_default
            self.cuts_dvpi0p_FD_Inb = cuts_dvpi0p_default
            self.cuts_dvpi0p_CDFT_Outb = cuts_dvpi0p_default
            self.cuts_dvpi0p_CD_Outb = cuts_dvpi0p_default
            self.cuts_dvpi0p_FD_Outb = cuts_dvpi0p_default
        if width == "mid":
            self.Ge2Threshold = Ge2Threshold_mid
            self.cuts_dvpi0p_CDFT_Inb = cuts_dvpi0p_CDFT_Inb_3sigma
            self.cuts_dvpi0p_CD_Inb = cuts_dvpi0p_CD_Inb_3sigma
            self.cuts_dvpi0p_FD_Inb = cuts_dvpi0p_FD_Inb_3sigma
            self.cuts_dvpi0p_CDFT_Outb = cuts_dvpi0p_CDFT_Outb_3sigma
            self.cuts_dvpi0p_CD_Outb = cuts_dvpi0p_CD_Outb_3sigma
            self.cuts_dvpi0p_FD_Outb = cuts_dvpi0p_FD_Outb_3sigma
        if width == "tight":
            self.Ge2Threshold = Ge2Threshold_tight
            self.cuts_dvpi0p_CDFT_Inb = cuts_dvpi0p_CDFT_Inb_3sigma
            self.cuts_dvpi0p_CD_Inb = cuts_dvpi0p_CD_Inb_3sigma
            self.cuts_dvpi0p_FD_Inb = cuts_dvpi0p_FD_Inb_3sigma
            self.cuts_dvpi0p_CDFT_Outb = cuts_dvpi0p_CDFT_Outb_3sigma
            self.cuts_dvpi0p_CD_Outb = cuts_dvpi0p_CD_Outb_3sigma
            self.cuts_dvpi0p_FD_Outb = cuts_dvpi0p_FD_Outb_3sigma
        if width == "loose":
            self.Ge2Threshold = Ge2Threshold_loose
            self.cuts_dvpi0p_CDFT_Inb = cuts_dvpi0p_CDFT_Inb_3sigma
            self.cuts_dvpi0p_CD_Inb = cuts_dvpi0p_CD_Inb_3sigma
            self.cuts_dvpi0p_FD_Inb = cuts_dvpi0p_FD_Inb_3sigma
            self.cuts_dvpi0p_CDFT_Outb = cuts_dvpi0p_CDFT_Outb_3sigma
            self.cuts_dvpi0p_CD_Outb = cuts_dvpi0p_CD_Outb_3sigma
            self.cuts_dvpi0p_FD_Outb = cuts_dvpi0p_FD_Outb_3sigma

    def readEPGG(self, entry_start = None, entry_stop = None, pol = "inbending", 
        detRes = False, logistics = False, nofid = False, 
        nocorr = False, noeloss = False, nopcorr = False, fidlevel = 'mid'):
        '''save data into df_epg, df_epgg for parent class epg'''
        self.readFile()

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Eedep", "Evz", "Esector", "Estat", "Epa"]
        eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
        eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
        eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc2Hitx", "EDc2Hity", "EDc2Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
        eleKeysRec.extend(["EFtof1bSector", "EFtof1bComponent"])
        eleKeysRec.extend(["Enphe", "EhtccX", "EhtccY"])
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "Psector", "Pchi2pid"]
        proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"])
        proKeysRec.extend(["PDc2Hitx", "PDc2Hity", "PDc2Hitz", "PDc3Hitx", "PDc3Hity", "PDc3Hitz"])
        proKeysRec.extend(["PFtof1aTime", "PFtof1bTime", "PFtof2Time", "PCtofTime"])
        proKeysRec.extend(["PFtof1bSector", "PFtof1bComponent"])
        # proKeysRec.extend(["Pchi2pid", "Pchi2track", "PNDFtrack"])
        gamKeysRec = ["Gpx", "Gpy", "Gpz", "Gedep", "GcX", "GcY", "Gsector"]
        gamKeysRec.extend(["GcalU1", "GcalV1", "GcalW1", "Gbeta"])

        if detRes:
            eleKeysRec.extend(["Evx", "Evy"])
            # eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
            # eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
            # eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
            eleKeysRec.extend(["EcalU2", "EcalV2", "EcalW2"])
            eleKeysRec.extend(["EcalU3", "EcalV3", "EcalW3"])
            # eleKeysRec.extend(["Enphe"])
            eleKeysRec.extend(["EhtccZ"])
            gamKeysRec.extend(["Gedep1", "Gedep2", "Gedep3"])
            # gamKeysRec.extend(["GcalU1", "GcalV1", "GcalW1"])
            gamKeysRec.extend(["GcalU2", "GcalV2", "GcalW2"])
            gamKeysRec.extend(["GcalU3", "GcalV3", "GcalW3"])
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
        else:
            # perform the fiducial cuts
            df_electronRec = electronFiducial(df_electronRec, mc = False, fidlevel = fidlevel)
            df_protonRec = protonFiducial(df_protonRec, pol = pol)
            df_gammaRec = gammaFiducial(df_gammaRec)
            print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))
            coincidence = reduce(np.intersect1d, (df_electronRec.event, df_protonRec.event, df_gammaRec.event))
            df_electronRec = df_electronRec.loc[df_electronRec.event.isin(coincidence), :]
            df_protonRec = df_protonRec.loc[df_protonRec.event.isin(coincidence), :]
            df_gammaRec = df_gammaRec.loc[df_gammaRec.event.isin(coincidence), :]
            print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))
        df_gg = pd.merge(df_gammaRec, df_gammaRec,
                         how='inner', on='event', suffixes=("", "2"))
        df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
        df_gg = df_gg.drop(['GIndex', 'GIndex2'], axis = 1)

        # Done with the fiducial cuts.

        # Post-processing
        # e1: Electron correction  - exp only
        # e2: Electron smearing    - mc only 
        # p1: Proton energy loss   - exp, mc both
        # p2: Proton correction    - exp only
        # p3: Proton smearing      - mc only
        # g1: Gamma correction     - exp only
        # g1: Gamma smearing       - mc only
        if nocorr:
            print("no correction applied")
            pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
            df_protonRec.loc[:, 'Pp'] = mag(pro)
            df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
            df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
            df_protonRec.loc[:, 'Pphi'] = getPhi(pro)
        else:
            #e1
            df_electronRec = electronMomentumCorrection(pol, df_electronRec)
            # #e2
            # df_electronRec = electronMomentumSmearing(df_electronRec)
            #p1
            df_protonRec = protonEnergyLossCorr(pol, df_protonRec)
            # #p2
            # df_protonRec = protonMomentumCorrection(pol, df_protonRec)
            # # #p3
            # # df_protonRec = protonMomentumSmearing(pol, df_protonRec, smearing = 1)
            # #g1
            # df_gg, df_gammaRec  = gammaMomentumCorrection(pol, df_gg, df_gammaRec)
            # # #g2
            # # df_gammaRec  = gammaMomentumSmearing(df_gammaRec, smearing = 1)

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
            df_protonRec = df_protonRec.drop(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PDc1theta", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"], axis = 1)
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
        df_epgg.loc[:,'y'] = df_epgg['nu']/self.ebeam
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
        df_epgg.loc[:,'openingAngle'] = angle(gam, gam2)
        df_epgg.loc[:, 'vzdiff'] = df_epgg.Evz - df_epgg.Pvz
        
        df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)
        df_epgg.loc[:, "closeness2"] = np.abs(df_epgg.loc[:, "MM2_ep"] - .1349766**2)

        df_epgg.loc[:, 'xBbin'] = np.zeros(len(df_epgg.xB), dtype = 'int') - 1
        df_epgg.loc[:, 'Q2bin'] = np.zeros(len(df_epgg.Q2), dtype = 'int') - 1
        df_epgg.loc[:, 'tbin'] = np.zeros(len(df_epgg.t1), dtype = 'int') - 1
        df_epgg.loc[:, 'phibin'] = np.zeros(len(df_epgg.phi1), dtype = 'int') - 1
        for xB in newxBbins2:
            df_epgg.xBbin = df_epgg.xBbin + (df_epgg.xB>xB).astype("int").to_numpy(dtype = 'int')

        for Q2 in newQ2bins2:
            df_epgg.Q2bin = df_epgg.Q2bin + (df_epgg.Q2>Q2).astype("int").to_numpy(dtype = 'int')

        for t1 in newtbins:
            df_epgg.tbin = df_epgg.tbin + (df_epgg.t1>t1).astype("int").to_numpy(dtype = 'int')

        for phi1 in phibins:
            df_epgg.phibin = df_epgg.phibin + (df_epgg.phi1>phi1).astype("int").to_numpy(dtype = 'int')

        # # encode unassigned bin as -1
        # df_epgg.loc[:, "Q2bin"] = -1
        # df_epgg.loc[:, "xBbin"] = -1
        # df_epgg.loc[:, "tbin"] = -1
        # # df_epgg.loc[:, "tbin2"] = -1
        # df_epgg.loc[:, "phibin"] = -1
        # # df_epgg.loc[:, "phibin2"] = -1
        # df_epgg.loc[:, "Q2xBbin"] = -1
        # df_epgg.loc[:, "Q2xBtbin"] = -1
        # # df_epgg.loc[:, "Q2xBtbin2"] = -1
        # df_epgg.loc[:, "Q2xBtphibin"] = -1
        # Q2xBbin = 0

        # # encode all binning
        # for Q2bin in range(len(Q2bin_i)):
        #     #square Q2 binning
        #     df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]), "Q2bin"] = Q2bin
        #     #adaptive xB binning
        #     for xBbin in range(len(xBbin_i[Q2bin])):
        #         if Q2bin < len(Q2bin_i) -1:
        #             if xBbin == 0:
        #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin #0
        #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin #0
        #             elif xBbin < len(xBbin_i[Q2bin])-1:
        #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin
        #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin
        #             else:
        #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
        #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin
        #         else:
        #             df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
        #             df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin #0

        #         Q2xBbin = Q2xBbin + 1
        # for tbin in range(len(tbin_i)):
        #     #square t binning
        #     df_epgg.loc[(df_epgg.t1>=tbin_i[tbin]) & (df_epgg.t1<tbin_f[tbin]), "tbin"] = tbin
        #     # df_epgg.loc[(df_epgg.t2>=tbin_i[tbin]) & (df_epgg.t2<tbin_f[tbin]), "tbin2"] = tbin
        # for phibin in range(len(phibin_i)):
        #     #square phi binning
        #     df_epgg.loc[(df_epgg.phi1>=phibin_i[phibin]) & (df_epgg.phi1<phibin_f[phibin]), "phibin"] = phibin
        #     # df_epgg.loc[(df_epgg.phi2>=phibin_i[phibin]) & (df_epgg.phi2<phibin_f[phibin]), "phibin2"] = phibin

        # df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtbin"] = len(tbin_i) * df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBbin"] + df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "tbin"]
        # # df_epgg.loc[(df_epgg.Q2bin>0)&(df_epgg.xBbin>0)&(df_epgg.tbin2>0), "Q2xBtbin2"] = df_epgg.Q2bin.astype(str) + df_epgg.xBbin.astype(str) + df_epgg.tbin2.astype(str)
        # df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtphibin"] = len(phibin_i) * df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtbin"] + df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "phibin"]

        # df_epgg = df_epgg.astype({"Q2bin": int, "xBbin": int, "tbin": int, "phibin": int, "Q2xBbin": int, "Q2xBtbin": int, "Q2xBtphibin": int})

        self.df_epgg = df_epgg

    def makeDVpi0P(self, pol = "inbending", nofid = False, allowsamesector = False):
        #make dvpi0 pairs
        df_dvpi0p = self.df_epgg

        #common cuts
        cut_xBupper = df_dvpi0p.loc[:, "xB"] < 1  # xB
        cut_xBlower = df_dvpi0p.loc[:, "xB"] > 0  # xB
        cut_Q2 = df_dvpi0p.loc[:, "Q2"] > 1  # Q2
        cut_W = df_dvpi0p.loc[:, "W"] > 2  # W
        cut_Ee = df_dvpi0p["Ee"] > 2  # Ee
        cut_Ge2 = df_dvpi0p["Ge2"] > self.Ge2Threshold  # Ge2 Threshold.
        if allowsamesector:
            cut_Esector = 1
            cut_Psector = 1
        else:
            cut_Esector = (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector"]) & (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector2"])
            cut_Psector = ~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector"]) ) & ~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector2"]) )
        cut_Ppmax = df_dvpi0p.Pp < 1.6  # Pp
        cut_Pthetamin = df_dvpi0p.Ptheta > 0  # Ptheta
        cut_Trigger = ((df_dvpi0p.TriggerBit & 1 << 1) > 0) | ((df_dvpi0p.TriggerBit & 1 << 2) > 0) | ((df_dvpi0p.TriggerBit & 1 << 3) > 0) | ((df_dvpi0p.TriggerBit & 1 << 4) > 0) | ((df_dvpi0p.TriggerBit & 1 << 5) > 0) | ((df_dvpi0p.TriggerBit & 1 << 6) > 0)
        cut_EFid = df_dvpi0p.EFid == 1
        # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
        cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge2 & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin & cut_Trigger & cut_EFid

        df_dvpi0p = df_dvpi0p[cut_common]

        # proton reconstruction quality
        # cut_FD_proton = (df_epgg.loc[:, "Psector"]<7) & (df_epgg.loc[:, "Ptheta"]<35)
        # cut_CD_proton = (df_epgg.loc[:, "Psector"]>7) & (df_epgg.loc[:, "Ptheta"]>45) & (df_epgg.loc[:, "Ptheta"]<65)
        # cut_proton = (cut_FD_proton)|(cut_CD_proton)
        cut_proton = 1

        df_dvpi0p.loc[:, "config"] = 0

        cuts_dvpi0p_CDFT_Inb = self.cuts_dvpi0p_CDFT_Inb 
        cuts_dvpi0p_CD_Inb = self.cuts_dvpi0p_CD_Inb 
        cuts_dvpi0p_FD_Inb = self.cuts_dvpi0p_FD_Inb
        cuts_dvpi0p_CDFT_Outb = self.cuts_dvpi0p_CDFT_Outb 
        cuts_dvpi0p_CD_Outb = self.cuts_dvpi0p_CD_Outb 
        cuts_dvpi0p_FD_Outb = self.cuts_dvpi0p_FD_Outb

        CD_Ptheta_ub = CD_Ptheta_ub_nominal
        CD_Ptheta_lb = CD_Ptheta_lb_nominal
        FD_Ptheta_inb_ub = FD_Ptheta_inb_ub_nominal
        FD_Ptheta_outb_ub = FD_Ptheta_outb_ub_nominal
        FD_Ptheta_lb = FD_Ptheta_lb_nominal
        if nofid:
            CD_Ptheta_lb = 0
            FD_Ptheta_inb_ub = 90

        if pol == "inbending":
            # if nofid:
            #     pass
            # else:
            #     vzdiffCD_lb,    vzdiffCD_ub    = -2.011, 2.314
            #     vzdiffFD_S1_lb, vzdiffFD_S1_ub = -3.209, 4.017
            #     vzdiffFD_S2_lb, vzdiffFD_S2_ub = -3.612, 4.139
            #     vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.328, 4.287
            #     vzdiffFD_S4_lb, vzdiffFD_S4_ub = -3.411, 4.108
            #     vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.607, 4.246
            #     vzdiffFD_S6_lb, vzdiffFD_S6_ub = -2.999, 3.927

            #     df_dvpi0p.loc[ (df_dvpi0p.Psector>4000) & ((df_dvpi0p.vzdiff<vzdiffCD_lb)   | (df_dvpi0p.vzdiff>vzdiffCD_ub)  ), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==1)   & ((df_dvpi0p.vzdiff<vzdiffFD_S1_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S1_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==2)   & ((df_dvpi0p.vzdiff<vzdiffFD_S2_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S2_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==3)   & ((df_dvpi0p.vzdiff<vzdiffFD_S3_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S3_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==4)   & ((df_dvpi0p.vzdiff<vzdiffFD_S4_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S4_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==5)   & ((df_dvpi0p.vzdiff<vzdiffFD_S5_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S5_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==6)   & ((df_dvpi0p.vzdiff<vzdiffFD_S6_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S6_ub)), "PFid"] = 0
            #CDFT
            cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CDFT = df_dvpi0p.Psector>7
            cut_Ptheta1_CDFT = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CDFT = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CDFT = df_dvpi0p.Gsector>7
            cut_GFid_CDFT = df_dvpi0p.GFid==1
            cut_GFid2_CDFT = df_dvpi0p.GFid2==1
            cut_PFid_CDFT = df_dvpi0p.PFid==1
            cut_mpi01_CDFT = df_dvpi0p["Mpi0"]                   < cuts_dvpi0p_CDFT_Inb["Mpi0_ub"]            #0.149#0.157  # mpi0
            cut_mpi02_CDFT = df_dvpi0p["Mpi0"]                   > cuts_dvpi0p_CDFT_Inb["Mpi0_lb"]            #0.126#0.118  # mpi0
            cut_mmep1_CDFT = df_dvpi0p["MM2_ep"]                 < cuts_dvpi0p_CDFT_Inb["MM2_ep_ub"]      #0.610#0.914  # mmep
            cut_mmep2_CDFT = df_dvpi0p["MM2_ep"]                 > cuts_dvpi0p_CDFT_Inb["MM2_ep_lb"]      #-0.384#-0.715  # mmep
            cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"]               < cuts_dvpi0p_CDFT_Inb["MM2_egg_ub"]         #1.641#2.155  # mmegg
            cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"]               > cuts_dvpi0p_CDFT_Inb["MM2_egg_lb"]         #0.0974#-0.417  # mmegg
            cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"]              < cuts_dvpi0p_CDFT_Inb["ME_epgg_ub"]         #0.481#0.799  # meepgg
            cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"]              > cuts_dvpi0p_CDFT_Inb["ME_epgg_lb"]         #-0.474#-0.792  # meepgg
            cut_mpt_CDFT = df_dvpi0p["MPt"]                      < cuts_dvpi0p_CDFT_Inb["MPt_ub"]             #0.1272#0.189  # mpt
            cut_recon_CDFT = df_dvpi0p["reconPi"]                < cuts_dvpi0p_CDFT_Inb["reconPi_ub"]         #0.955  # recon gam angle
            cut_coplanarity_CDFT = df_dvpi0p["coplanarity"]      < cuts_dvpi0p_CDFT_Inb["coplanarity_ub"]     #9.259#15.431  # coplanarity angle
            cut_mmepgg1_CDFT = df_dvpi0p["MM2_epgg"]             < cuts_dvpi0p_CDFT_Inb["MM2_epgg_ub"]        #0.02564#0.0440  # mmepgg
            cut_mmepgg2_CDFT = df_dvpi0p["MM2_epgg"]             > cuts_dvpi0p_CDFT_Inb["MM2_epgg_lb"]        #-0.02944#-0.0478  # mmepgg

            cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & 
                        cut_GFid_CDFT & cut_GFid2_CDFT & cut_PFid_CDFT &
                        cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
                        cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
                        cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


            #CD
            cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CD = df_dvpi0p.Psector>7
            cut_Ptheta1_CD = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CD = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
            cut_Gsector2_CD = (df_dvpi0p.Gsector2<7) & (df_dvpi0p.Gsector2>0)
            cut_GFid_CD = df_dvpi0p.GFid==1
            cut_GFid2_CD = df_dvpi0p.GFid2==1
            cut_PFid_CD = df_dvpi0p.PFid==1
            cut_mpi01_CD = df_dvpi0p["Mpi0"] < cuts_dvpi0p_CD_Inb["Mpi0_ub"]       #0.162  # mpi0
            cut_mpi02_CD = df_dvpi0p["Mpi0"] > cuts_dvpi0p_CD_Inb["Mpi0_lb"]       #0.107  # mpi0
            cut_mmep1_CD = df_dvpi0p["MM2_ep"] < cuts_dvpi0p_CD_Inb["MM2_ep_ub"]   #0.354  # mmep
            cut_mmep2_CD = df_dvpi0p["MM2_ep"] > cuts_dvpi0p_CD_Inb["MM2_ep_lb"]   #-0.283  # mmep
            cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < cuts_dvpi0p_CD_Inb["MM2_egg_ub"]    #1.922  # mmegg
            cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > cuts_dvpi0p_CD_Inb["MM2_egg_lb"]    #0.007  # mmegg
            cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < cuts_dvpi0p_CD_Inb["ME_epgg_ub"]   #0.822  # meepgg
            cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > cuts_dvpi0p_CD_Inb["ME_epgg_lb"]   #-0.677  # meepgg
            cut_mpt_CD = df_dvpi0p["MPt"] < cuts_dvpi0p_CD_Inb["MPt_ub"]       #0.176  # mpt
            cut_recon_CD = df_dvpi0p["reconPi"] < cuts_dvpi0p_CD_Inb["reconPi_ub"] #1.476  # recon gam angle
            cut_coplanarity_CD = df_dvpi0p["coplanarity"] < cuts_dvpi0p_CD_Inb["coplanarity_ub"]#10.203  # coplanarity angle
            cut_mmepgg1_CD = df_dvpi0p["MM2_epgg"] < cuts_dvpi0p_CD_Inb["MM2_epgg_ub"] #0.0208  # mmepgg
            cut_mmepgg2_CD = df_dvpi0p["MM2_epgg"] > cuts_dvpi0p_CD_Inb["MM2_epgg_lb"] #-0.0250  # mmepgg

            cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_Gsector2_CD & 
                        cut_GFid_CD & cut_GFid2_CD & cut_PFid_CD &
                        cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
                        cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
                        cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

            #FD
            cut_Pp1_FD = df_dvpi0p.Pp > 0.42  # Pp
            cut_Psector_FD = df_dvpi0p.Psector<7
            cut_Ptheta1_FD = df_dvpi0p.Ptheta<FD_Ptheta_inb_ub
            cut_Ptheta2_FD = df_dvpi0p.Ptheta>FD_Ptheta_lb
            cut_Gsector_FD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
            cut_Gsector2_FD = (df_dvpi0p.Gsector2<7) & (df_dvpi0p.Gsector2>0)
            cut_GFid_FD = df_dvpi0p.GFid==1
            cut_GFid2_FD = df_dvpi0p.GFid2==1
            cut_PFid_FD = df_dvpi0p.PFid==1
            cut_mpi01_FD = df_dvpi0p["Mpi0"] < cuts_dvpi0p_FD_Inb["Mpi0_ub"]       #                         0.178  # mpi0
            cut_mpi02_FD = df_dvpi0p["Mpi0"] > cuts_dvpi0p_FD_Inb["Mpi0_lb"]       #                         0.0910  # mpi0
            cut_mmep1_FD = df_dvpi0p["MM2_ep"] <   cuts_dvpi0p_FD_Inb["MM2_ep_ub"] #                         0.335  # mmep
            cut_mmep2_FD = df_dvpi0p["MM2_ep"] >   cuts_dvpi0p_FD_Inb["MM2_ep_lb"] #                         -0.271  # mmep
            cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < cuts_dvpi0p_FD_Inb["MM2_egg_ub"]    #                         1.762  # mmegg
            cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > cuts_dvpi0p_FD_Inb["MM2_egg_lb"]    #                         0.117  # mmegg
            cut_meepgg1_FD = df_dvpi0p["ME_epgg"] <    cuts_dvpi0p_FD_Inb["ME_epgg_ub"]    #                         0.816 # meepgg
            cut_meepgg2_FD = df_dvpi0p["ME_epgg"] >    cuts_dvpi0p_FD_Inb["ME_epgg_lb"]    #                         -0.685  # meepgg
            cut_mpt_FD = df_dvpi0p["MPt"] <    cuts_dvpi0p_FD_Inb["MPt_ub"]        #                         0.180  # mpt
            cut_recon_FD = df_dvpi0p["reconPi"] <  cuts_dvpi0p_FD_Inb["reconPi_ub"]    #                         1.363  # recon gam angle
            cut_coplanarity_FD = df_dvpi0p["coplanarity"] <    cuts_dvpi0p_FD_Inb["coplanarity_ub"]#                         9.190  # coplanarity angle
            cut_mmepgg1_FD = df_dvpi0p["MM2_epgg"] <   cuts_dvpi0p_FD_Inb["MM2_epgg_ub"]   #                         0.0189  # mmepgg
            cut_mmepgg2_FD = df_dvpi0p["MM2_epgg"] >   cuts_dvpi0p_FD_Inb["MM2_epgg_lb"]   #                         -0.0224  # mmepgg

            cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_Gsector2_FD &
                        cut_GFid_FD & cut_GFid2_FD & cut_PFid_FD &
                        cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
                        cut_mmegg1_FD & cut_mmegg2_FD & cut_meepgg1_FD & cut_meepgg2_FD &
                        cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepgg1_FD & cut_mmepgg2_FD)

        elif pol == "outbending":
            # if nofid:
            #     pass
            # else:
            #     vzdiffCD_lb,    vzdiffCD_ub    = -2.737, 2.096
            #     vzdiffFD_S1_lb, vzdiffFD_S1_ub = -4.435, 3.429
            #     vzdiffFD_S2_lb, vzdiffFD_S2_ub = -4.646, 2.978
            #     vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.922, 3.040
            #     vzdiffFD_S4_lb, vzdiffFD_S4_ub = -4.646, 3.493
            #     vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.901, 3.750
            #     vzdiffFD_S6_lb, vzdiffFD_S6_ub = -3.846, 3.623

            #     df_dvpi0p.loc[ (df_dvpi0p.Psector>4000) & ((df_dvpi0p.vzdiff<vzdiffCD_lb)   | (df_dvpi0p.vzdiff>vzdiffCD_ub)  ), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==1)   & ((df_dvpi0p.vzdiff<vzdiffFD_S1_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S1_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==2)   & ((df_dvpi0p.vzdiff<vzdiffFD_S2_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S2_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==3)   & ((df_dvpi0p.vzdiff<vzdiffFD_S3_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S3_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==4)   & ((df_dvpi0p.vzdiff<vzdiffFD_S4_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S4_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==5)   & ((df_dvpi0p.vzdiff<vzdiffFD_S5_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S5_ub)), "PFid"] = 0
            #     df_dvpi0p.loc[ (df_dvpi0p.Psector==6)   & ((df_dvpi0p.vzdiff<vzdiffFD_S6_lb) | (df_dvpi0p.vzdiff>vzdiffFD_S6_ub)), "PFid"] = 0
            #CDFT
            cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CDFT = df_dvpi0p.Psector>7
            cut_Ptheta1_CDFT = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CDFT = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CDFT = df_dvpi0p.Gsector>7
            cut_GFid_CDFT = df_dvpi0p.GFid==1
            cut_GFid2_CDFT = df_dvpi0p.GFid2==1
            cut_PFid_CDFT = df_dvpi0p.PFid==1
            cut_mpi01_CDFT = df_dvpi0p["Mpi0"] <  cuts_dvpi0p_CDFT_Outb["Mpi0_ub"]         #   0.151#0.160  # mpi0
            cut_mpi02_CDFT = df_dvpi0p["Mpi0"] >  cuts_dvpi0p_CDFT_Outb["Mpi0_lb"]         #   0.124#0.115  # mpi0
            cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] <  cuts_dvpi0p_CDFT_Outb["MM2_ep_ub"]     #   0.575#0.892  # mmep
            cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] >  cuts_dvpi0p_CDFT_Outb["MM2_ep_lb"]     #   -0.378#-0.694  # mmep
            cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] <  cuts_dvpi0p_CDFT_Outb["MM2_egg_ub"]  #   1.665#2.184  # mmegg
            cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] >  cuts_dvpi0p_CDFT_Outb["MM2_egg_lb"]  #   0.107#-0.412  # mmegg
            cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] <  cuts_dvpi0p_CDFT_Outb["ME_epgg_ub"]     #   0.514#0.844  # meepgg
            cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] >  cuts_dvpi0p_CDFT_Outb["ME_epgg_lb"]     #   -0.476#-0.806  # meepgg
            cut_mpt_CDFT = df_dvpi0p["MPt"] <  cuts_dvpi0p_CDFT_Outb["MPt_ub"]         #   0.146#0.210  # mpt
            cut_recon_CDFT = df_dvpi0p["reconPi"] <  cuts_dvpi0p_CDFT_Outb["reconPi_ub"]   #   1.114#1.630  # recon gam angle
            cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] <  cuts_dvpi0p_CDFT_Outb["coplanarity_ub"] #   10.69#17.817  # coplanarity angle
            cut_mmepgg1_CDFT = df_dvpi0p["MM2_epgg"] <  cuts_dvpi0p_CDFT_Outb["MM2_epgg_ub"]   #   0.0324#0.0549  # mmepgg
            cut_mmepgg2_CDFT = df_dvpi0p["MM2_epgg"] >  cuts_dvpi0p_CDFT_Outb["MM2_epgg_lb"]   #   -0.035#-0.0575  # mmepgg

            cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & 
                        cut_GFid_CDFT & cut_GFid2_CDFT & cut_PFid_CDFT &
                        cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
                        cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
                        cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


            #CD
            cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CD = df_dvpi0p.Psector>7
            cut_Ptheta1_CD = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CD = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
            cut_Gsector2_CD = (df_dvpi0p.Gsector2<7) & (df_dvpi0p.Gsector2>0)
            cut_GFid_CD = df_dvpi0p.GFid==1
            cut_GFid2_CD = df_dvpi0p.GFid2==1
            cut_PFid_CD = df_dvpi0p.PFid==1
            cut_mpi01_CD = df_dvpi0p["Mpi0"] < cuts_dvpi0p_CD_Outb["Mpi0_ub"]                 #0.163  # mpi0
            cut_mpi02_CD = df_dvpi0p["Mpi0"] > cuts_dvpi0p_CD_Outb["Mpi0_lb"]                 #0.106  # mpi0
            cut_mmep1_CD = df_dvpi0p["MM2_ep"] < cuts_dvpi0p_CD_Outb["MM2_ep_ub"]             #0.294  # mmep
            cut_mmep2_CD = df_dvpi0p["MM2_ep"] > cuts_dvpi0p_CD_Outb["MM2_ep_lb"]             #-0.218  # mmep
            cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < cuts_dvpi0p_CD_Outb["MM2_egg_ub"]              #1.876  # mmegg
            cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > cuts_dvpi0p_CD_Outb["MM2_egg_lb"]              #-0.0142  # mmegg
            cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < cuts_dvpi0p_CD_Outb["ME_epgg_ub"]             #0.700  # meepgg
            cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > cuts_dvpi0p_CD_Outb["ME_epgg_lb"]             #-0.597  # meepgg
            cut_mpt_CD = df_dvpi0p["MPt"] < cuts_dvpi0p_CD_Outb["MPt_ub"]                 #0.194  # mpt
            cut_recon_CD = df_dvpi0p["reconPi"] < cuts_dvpi0p_CD_Outb["reconPi_ub"]           #1.761  # recon gam angle
            cut_coplanarity_CD = df_dvpi0p["coplanarity"] < cuts_dvpi0p_CD_Outb["coplanarity_ub"]          #9.530  # coplanarity angle
            cut_mmepgg1_CD = df_dvpi0p["MM2_epgg"] < cuts_dvpi0p_CD_Outb["MM2_epgg_ub"]           #0.0182  # mmepgg
            cut_mmepgg2_CD = df_dvpi0p["MM2_epgg"] > cuts_dvpi0p_CD_Outb["MM2_epgg_lb"]           #-0.0219  # mmepgg

            cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_Gsector2_CD & 
                        cut_GFid_CD & cut_GFid2_CD & cut_PFid_CD &
                        cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
                        cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
                        cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

            #FD
            cut_Pp1_FD = df_dvpi0p.Pp > 0.5  # Pp
            cut_Psector_FD = df_dvpi0p.Psector<7
            cut_Ptheta1_FD = df_dvpi0p.Ptheta<FD_Ptheta_outb_ub
            cut_Ptheta2_FD = df_dvpi0p.Ptheta>FD_Ptheta_lb
            cut_Gsector_FD = (df_dvpi0p.Gsector<7) & (df_dvpi0p.Gsector>0)
            cut_Gsector2_FD = (df_dvpi0p.Gsector2<7) & (df_dvpi0p.Gsector2>0)
            cut_GFid_FD = df_dvpi0p.GFid==1
            cut_GFid2_FD = df_dvpi0p.GFid2==1
            cut_PFid_FD = df_dvpi0p.PFid==1
            cut_mpi01_FD = df_dvpi0p["Mpi0"] < cuts_dvpi0p_FD_Outb["Mpi0_ub"]                    # 0.164  # mpi0
            cut_mpi02_FD = df_dvpi0p["Mpi0"] > cuts_dvpi0p_FD_Outb["Mpi0_lb"]                    # 0.105  # mpi0
            cut_mmep1_FD = df_dvpi0p["MM2_ep"] < cuts_dvpi0p_FD_Outb["MM2_ep_ub"]                # 0.323  # mmep
            cut_mmep2_FD = df_dvpi0p["MM2_ep"] > cuts_dvpi0p_FD_Outb["MM2_ep_lb"]                # -0.256  # mmep
            cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < cuts_dvpi0p_FD_Outb["MM2_egg_ub"]                 # 1.828  # mmegg
            cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > cuts_dvpi0p_FD_Outb["MM2_egg_lb"]                 # 0.0491  # mmegg
            cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < cuts_dvpi0p_FD_Outb["ME_epgg_ub"]                # 0.754  # meepgg
            cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > cuts_dvpi0p_FD_Outb["ME_epgg_lb"]                # -0.583  # meepgg
            cut_mpt_FD = df_dvpi0p["MPt"] < cuts_dvpi0p_FD_Outb["MPt_ub"]                    # 0.177  # mpt
            cut_recon_FD = df_dvpi0p["reconPi"] < cuts_dvpi0p_FD_Outb["reconPi_ub"]              # 1.940  # recon gam angle
            cut_coplanarity_FD = df_dvpi0p["coplanarity"] < cuts_dvpi0p_FD_Outb["coplanarity_ub"]             # 7.498  # coplanarity angle
            cut_mmepgg1_FD = df_dvpi0p["MM2_epgg"] < cuts_dvpi0p_FD_Outb["MM2_epgg_ub"]              # 0.0195  # mmepgg
            cut_mmepgg2_FD = df_dvpi0p["MM2_epgg"] > cuts_dvpi0p_FD_Outb["MM2_epgg_lb"]              # -0.0240  # mmepgg

            cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_Gsector2_FD &
                        cut_GFid_FD & cut_GFid2_FD & cut_PFid_FD & 
                        cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
                        cut_mmegg1_FD & cut_mmegg2_FD & cut_meepgg1_FD & cut_meepgg2_FD &
                        cut_mpt_FD & cut_recon_FD & cut_coplanarity_FD & cut_mmepgg1_FD & cut_mmepgg2_FD)

        df_dvpi0p.loc[cut_CDFT, "config"] = 3
        df_dvpi0p.loc[cut_CD, "config"] = 2
        df_dvpi0p.loc[cut_FD, "config"] = 1

        df_dvpi0p = df_dvpi0p[df_dvpi0p.config>0]

        #For an event, there can be two gg's passed conditions above.
        #Take only one gg's that makes pi0 invariant mass
        #This case is very rare.
        #For now, duplicated proton is not considered.
        df_dvpi0p = df_dvpi0p.sort_values(by=['closeness', 'closeness2'], ascending = [True, True])
        df_dvpi0p = df_dvpi0p.loc[~df_dvpi0p.event.duplicated(), :]
        df_dvpi0p = df_dvpi0p.sort_values(by='event')        
        self.df_dvpi0p = df_dvpi0p #done with saving x

    def save(self, raw = False, pol = 'inbending'):
        if raw:
            CD_Ptheta_ub = CD_Ptheta_ub_nominal
            CD_Ptheta_lb = CD_Ptheta_lb_nominal
            FD_Ptheta_inb_ub = FD_Ptheta_inb_ub_nominal
            FD_Ptheta_outb_ub = FD_Ptheta_outb_ub_nominal
            FD_Ptheta_lb = FD_Ptheta_lb_nominal

            print("saving raw with common cuts")
            df_x = self.df_epgg
            #common cuts
            cut_xBupper = df_x.loc[:, "xB"] < 1  # xB
            cut_xBlower = df_x.loc[:, "xB"] > 0  # xB
            cut_Q2 = df_x.loc[:, "Q2"] > 1  # Q2
            cut_W = df_x.loc[:, "W"] > 2  # W
            cut_Ee = df_x["Ee"] > 2  # Ee
            cut_Ge2 = df_x["Ge2"] > self.Ge2Threshold  # Ge cut. Ge>3 for DVCS module.
            cut_Esector = (df_x["Esector"]!=df_x["Gsector"]) & (df_x["Esector"]!=df_x["Gsector2"]) 
            cut_Psector = ~( ((df_x["Pstat"]//10)%10>0) & (df_x["Psector"]==df_x["Gsector"]) ) & ~( ((df_x["Pstat"]//10)%10>0) & (df_x["Psector"]==df_x["Gsector2"]) )
            cut_Ppmax = df_x.Pp < 1.6  # Pp
            cut_Pthetamin = df_x.Ptheta > 0  # Ptheta
            # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
            cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge2 & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin
            df_x = df_x[cut_common]

            #CDFT
            cut_Pp1_CDFT = df_x.Pp > 0.3  # Pp
            cut_Psector_CDFT = df_x.Psector>7
            cut_Ptheta1_CDFT = df_x.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CDFT = df_x.Ptheta>CD_Ptheta_lb
            cut_Gsector_CDFT = df_x.Gsector>7
            cut_Gsector2_CDFT = df_x.Gsector2>7
            cut_GFid_CDFT = df_x.GFid==1
            cut_GFid2_CDFT = df_x.GFid2==1
            cut_PFid_CDFT = df_x.PFid==1
            #CD
            cut_Pp1_CD = df_x.Pp > 0.3  # Pp
            cut_Psector_CD = df_x.Psector>7
            cut_Ptheta1_CD = df_x.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CD = df_x.Ptheta>CD_Ptheta_lb
            cut_Gsector_CD = (df_x.Gsector<7) & (df_x.Gsector>0)
            cut_Gsector2_CD = (df_x.Gsector2<7) & (df_x.Gsector2>0)
            cut_GFid_CD = df_x.GFid==1
            cut_GFid2_CD = df_x.GFid2==1
            cut_PFid_CD = df_x.PFid==1
            #FD
            if pol == "inbending":
                cut_Pp1_FD = df_x.Pp > 0.42  # Pp
                cut_Ptheta1_FD = df_x.Ptheta<FD_Ptheta_inb_ub
            elif pol == "outbending":
                cut_Pp1_FD = df_x.Pp > 0.5  # Pp
                cut_Ptheta1_FD = df_x.Ptheta<FD_Ptheta_outb_ub
            cut_Psector_FD = df_x.Psector<7
            cut_Ptheta2_FD = df_x.Ptheta>FD_Ptheta_lb
            cut_Gsector_FD = (df_x.Gsector<7) & (df_x.Gsector>0)
            cut_Gsector2_FD = (df_x.Gsector2<7) & (df_x.Gsector2>0)
            cut_GFid_FD = df_x.GFid==1
            cut_GFid2_FD = df_x.GFid2==1
            cut_PFid_FD = df_x.PFid==1

            cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & cut_Gsector2_CDFT & cut_GFid_CDFT & cut_GFid2_CDFT & cut_PFid_CDFT)
            cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_Gsector2_CD & cut_GFid_CD & cut_GFid2_CD & cut_PFid_CD)
            cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_Gsector2_FD & cut_GFid_FD & cut_GFid2_FD & cut_PFid_FD)

            df_x.loc[cut_CDFT, "config"] = 3
            df_x.loc[cut_CD, "config"] = 2
            df_x.loc[cut_FD, "config"] = 1

            df_x = df_x[df_x.config>0]

            df_x = df_x.sort_values(by=['closeness', 'closeness2'], ascending = [True, True])
            df_x = df_x.loc[~df_x.event.duplicated(), :]
            df_x = df_x.sort_values(by='event')

        else:
            df_x = self.df_dvpi0p
        self.df = df_x


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-S","--entry_start", help="entry_start to start reading the root file", default = None)
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-p","--polarity", help="polarity", default = "inbending")
    parser.add_argument("-d","--detRes", help="include detector response", action = "store_true")
    parser.add_argument("-w","--width", help="width of selection cuts", default = "default")
    parser.add_argument("-r","--raw", help="save raw only", default = False, action = "store_true")
    parser.add_argument("-l","--logistics", help="include logistics", action = "store_true")
    parser.add_argument("-nf","--nofid", help="no additional fiducial cuts", action = "store_true")
    parser.add_argument("-nc","--nocorr", help="no momentum correction", action = "store_true")
    parser.add_argument("-ne","--noeloss", help="no energy loss correction", action = "store_true")
    parser.add_argument("-np","--nopcorr", help="no proton correction at all", action = "store_true")
    parser.add_argument("-fl","--fidlevel", help="fiducial cut level", default = 'mid')
    parser.add_argument("-as","--allowsamesector", help="allow same sector conditions", action = "store_true")
    parser.add_argument("-be","--beam", help="beam energy", default = "10.604")

    args = parser.parse_args()

    if args.entry_start:
        args.entry_start = int(args.entry_start)
    if args.entry_stop:
        args.entry_stop = int(args.entry_stop)

    be = float(args.beam)
    converter = root2pickle(args.fname, entry_start = args.entry_start, 
        entry_stop = args.entry_stop, pol = args.polarity, detRes = args.detRes, raw = args.raw,
        width = args.width, logistics = args.logistics, nofid = args.nofid, nocorr = args.nocorr, noeloss = args.noeloss, nopcorr = args.nopcorr,
        fidlevel = args.fidlevel, allowsamesector = args.allowsamesector, ebeam = be)
    df = converter.df

    df.to_pickle(args.out)