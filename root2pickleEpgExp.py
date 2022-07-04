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
from utils.fiducial import *
from scipy.stats import skewnorm

class root2pickle():
    #class to read root to make epg pairs, inherited from epg
    def __init__(self, fname, entry_start = None, entry_stop = None, pol = "inbending", detRes = False, logistics = False, width = "mid", nofid = False):
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

        self.determineWidth(width = width)
        self.readEPGG(entry_start = entry_start, entry_stop = entry_stop, pol = pol, detRes = detRes, logistics = logistics, nofid = nofid)
        self.saveDVCSvars()
        self.saveDVpi0vars()
        self.makeDVpi0P_DVCS(pol = pol)
        self.pi02gSubtraction()
        self.makeDVCS(pol = pol)
        self.save()

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
            self.cuts_dvcs_CDFT_Inb = cuts_dvcs_default_CDFT
            self.cuts_dvcs_CD_Inb = cuts_dvcs_default_CD
            self.cuts_dvcs_FD_Inb = cuts_dvcs_default_FD
            self.cuts_dvcs_CDFT_Outb = cuts_dvcs_default_CDFT
            self.cuts_dvcs_CD_Outb = cuts_dvcs_default_CD
            self.cuts_dvcs_FD_Outb = cuts_dvcs_default_FD
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

    def readEPGG(self, entry_start = None, entry_stop = None, pol = "inbending", detRes = False, logistics = False, nofid = False):
        '''save data into df_epg, df_epgg for parent class epg'''
        self.readFile()

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Eedep", "Evz", "Esector", "TriggerBit"]
        eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
        eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
        eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc2Hitx", "EDc2Hity", "EDc2Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
        eleKeysRec.extend(["Enphe"])
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "Psector", "Pchi2pid"]
        proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"])
        proKeysRec.extend(["PDc2Hitx", "PDc2Hity", "PDc2Hitz", "PDc3Hitx", "PDc3Hity", "PDc3Hitz"])
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
            eleKeysRec.extend(["EhtccX", "EhtccY", "EhtccZ"])
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
        if logistics:
            eleKeysRec.extend(["EventNum", "RunNum", "beamQ", "liveTime", "helicity"])

        # read them
        for key in eleKeysRec:
            df_electronRec[key] = self.tree[key].array(library="pd", entry_start = entry_start, entry_stop=entry_stop)
        for key in proKeysRec:
            df_protonRec[key] = self.tree[key].array(library="pd", entry_start = entry_start, entry_stop=entry_stop)
        for key in gamKeysRec:
            df_gammaRec[key] = self.tree[key].array(library="pd", entry_start = entry_start, entry_stop=entry_stop)
        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        df_gammaRec = df_gammaRec.astype({"Gpx": float, "Gpy": float, "Gpz": float, "Gedep": float, "GcX": float, "GcY": float})
        ele = [df_electronRec['Epx'], df_electronRec['Epy'], df_electronRec['Epz']]
        df_electronRec.loc[:, 'Ep'] = mag(ele)
        df_electronRec.loc[:,'ESamplFrac'] = df_electronRec.Eedep/ df_electronRec.Ep

        print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))
        df_electronRec = electronFiducial(df_electronRec, pol = pol, mc = False)
        df_protonRec = protonFiducial(df_protonRec, pol = pol)
        df_gammaRec = gammaFiducial(df_gammaRec)
        print(len(df_electronRec), len(df_protonRec), len(df_gammaRec))

        #apply photon fiducial cuts
        if nofid:
            df_gammaRec.loc[:, "GFid"] = 1
        else:
            #photon FD fiducial cuts by F.X. Girod
            df_gammaRec.loc[:, "GFid"] = 0

            sector_cond = [df_gammaRec.Gsector ==1, df_gammaRec.Gsector ==2, df_gammaRec.Gsector ==3, df_gammaRec.Gsector ==4, df_gammaRec.Gsector ==5, df_gammaRec.Gsector ==6]
            psplit = np.select(sector_cond, [87, 82, 85, 77, 78, 82])
            tleft = np.select(sector_cond, [58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822])
            tright = np.select(sector_cond, [58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914])
            sleft = np.select(sector_cond, [0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343])
            sright = np.select(sector_cond, [-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729])
            rleft = np.select(sector_cond, [64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997])
            rright = np.select(sector_cond, [65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641])
            qleft = np.select(sector_cond, [0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899])
            qright = np.select(sector_cond, [-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481])
            #first condition
            ang = np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
            GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
            GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)

            df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
            df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

            cond1_1 = df_gammaRec.GcX >= psplit
            cond1_2 = df_gammaRec.GcY < sleft * (df_gammaRec.GcX - tleft)
            cond1_3 = df_gammaRec.GcY > sright * (df_gammaRec.GcX - tright)
            cond1_4 = df_gammaRec.Gsector < 7
            cond1 = cond1_1 & cond1_2 & cond1_3 & cond1_4
            df_gammaRec.loc[cond1, "GFid"] = 1
            #second condition else if the first
            # cond2_0 = df_gammaRec.GFid == 0 # not necessary, because cond2_1 rules out the first (S. Lee)
            cond2_1 = df_gammaRec.GcX < psplit
            cond2_2 = df_gammaRec.GcY < qleft * (df_gammaRec.GcX - rleft)
            cond2_3 = df_gammaRec.GcY > qright * (df_gammaRec.GcX - rright)
            cond2_4 = df_gammaRec.Gsector < 7
            cond2 = cond2_1 & cond2_2 & cond2_3 & cond2_4
            df_gammaRec.loc[cond2, "GFid"] = 1

            df_gammaRec.loc[df_gammaRec.Gsector > 7, "GFid"] = 1

            #FT fiducial cuts
            circleCenterX1 = -8.419
            circleCenterY1 = 9.889
            circleRadius1 = 1.6

            circleCenterX2 = -9.89
            circleCenterY2 = -5.327
            circleRadius2 = 1.6

            circleCenterX3 = -6.15
            circleCenterY3 = -13
            circleRadius3 = 2.3

            circleCenterX4 = 3.7
            circleCenterY4 = -6.5
            circleRadius4 = 2
            
            circle1 = (df_gammaRec.GcX - circleCenterX1)**2 + (df_gammaRec.GcY - circleCenterY1)**2 < circleRadius1**2
            circle2 = (df_gammaRec.GcX - circleCenterX2)**2 + (df_gammaRec.GcY - circleCenterY2)**2 < circleRadius2**2
            circle3 = (df_gammaRec.GcX - circleCenterX3)**2 + (df_gammaRec.GcY - circleCenterY3)**2 < circleRadius3**2
            circle4 = (df_gammaRec.GcX - circleCenterX4)**2 + (df_gammaRec.GcY - circleCenterY4)**2 < circleRadius4**2

            df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle1, "GFid"] = 0
            df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle2, "GFid"] = 0
            df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle3, "GFid"] = 0
            df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle4, "GFid"] = 0

            exclusion1_1 = (df_gammaRec.GcalW1 > 74) & (df_gammaRec.GcalW1 < 79.8)
            exclusion1_2 = (df_gammaRec.GcalW1 > 83.6) & (df_gammaRec.GcalW1 < 92.2)
            exclusion1_3 = (df_gammaRec.GcalW1 > 212.5) & (df_gammaRec.GcalW1 < 230)
            exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
            df_gammaRec.loc[(df_gammaRec.Gsector == 1) & exclusion1, "GFid"] = 0
            exclusion2_1 = (df_gammaRec.GcalW1 < 14)
            exclusion2_2 = (df_gammaRec.GcalU1 > 111.2) & (df_gammaRec.GcalU1 < 119.3)
            exclusion2_3 = (df_gammaRec.GcalV1 > 113) & (df_gammaRec.GcalV1 < 118.7)
            exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
            df_gammaRec.loc[(df_gammaRec.Gsector == 2) & exclusion2, "GFid"] = 0
            exclusion3 = df_gammaRec.GcalW1 < 14
            df_gammaRec.loc[(df_gammaRec.Gsector == 3) & exclusion3, "GFid"] = 0
            exclusion4_1 = (df_gammaRec.GcalV1 < 14)
            exclusion4_2 = (df_gammaRec.GcalV1 > 229.4) & (df_gammaRec.GcalV1 < 240.7)
            exclusion4_3 = (df_gammaRec.GcalW1 > 135) & (df_gammaRec.GcalW1 < 150)
            exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
            df_gammaRec.loc[(df_gammaRec.Gsector == 4) & exclusion4, "GFid"] = 0
            exclusion6 = (df_gammaRec.GcalW1 > 170) & (df_gammaRec.GcalW1 < 192)
            df_gammaRec.loc[(df_gammaRec.Gsector == 6) & exclusion6, "GFid"] = 0

        #apply electron fiducial cuts
        if nofid:
            df_electronRec.loc[:, "EFid"] = 1
        else:
            df_electronRec.loc[:, "EFid"] = 1
            exclusion1_1 = (df_electronRec.EcalW1 > 74) & (df_electronRec.EcalW1 < 79.8)
            exclusion1_2 = (df_electronRec.EcalW1 > 83.6) & (df_electronRec.EcalW1 < 92.2)
            exclusion1_3 = (df_electronRec.EcalW1 > 212.5) & (df_electronRec.EcalW1 < 230)
            exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
            df_electronRec.loc[(df_electronRec.Esector == 1) & exclusion1, "EFid"] = 0
            exclusion2_1 = (df_electronRec.EcalW1 < 14)
            exclusion2_2 = (df_electronRec.EcalU1 > 111.2) & (df_electronRec.EcalU1 < 119.3)
            exclusion2_3 = (df_electronRec.EcalV1 > 113) & (df_electronRec.EcalV1 < 118.7)
            exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
            df_electronRec.loc[(df_electronRec.Esector == 2) & exclusion2, "EFid"] = 0
            exclusion3 = df_electronRec.EcalW1 < 14
            df_electronRec.loc[(df_electronRec.Esector == 3) & exclusion3, "EFid"] = 0
            exclusion4_1 = (df_electronRec.EcalV1 < 14)
            exclusion4_2 = (df_electronRec.EcalV1 > 229.4) & (df_electronRec.EcalV1 < 240.7)
            exclusion4_3 = (df_electronRec.EcalW1 > 135) & (df_electronRec.EcalW1 < 150)
            exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
            df_electronRec.loc[(df_electronRec.Esector == 4) & exclusion4, "EFid"] = 0
            exclusion6 = (df_electronRec.EcalW1 > 170) & (df_electronRec.EcalW1 < 192)
            df_electronRec.loc[(df_electronRec.Esector == 6) & exclusion6, "EFid"] = 0


        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'event'] = df_gammaRec.index.get_level_values('entry')
        df_gammaRec.loc[:,'GIndex'] = df_gammaRec.index.get_level_values('subentry')

        #prepare for proton energy loss corrections correction
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)

        df_protonRec.loc[:, "PDc1theta"] = -100000

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
        else:
            df_protonRec.loc[:, "PCvt12theta"] = -100000
            df_protonRec.loc[:, "PCvt12phi"] = -100000

        df_protonRecFD = df_protonRec.loc[df_protonRec.Psector<7, :]
        df_protonRecCD = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), :]
        df_protonRecOthers = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta>=75), :]

        correction = False

        #two band criterion
        def corr(x, t):
            x0, x1, x2, x3 = x
            return x0 + x1*np.power(t-np.ones(len(t))*0.3, x3)

        df_protonRecFD = df_protonRecFD.loc[df_protonRec.Pp > 0.3, :]
        df_protonRecFD.loc[:, "PDc1theta"] = getTheta([df_protonRecFD.PDc1Hitx, df_protonRecFD.PDc1Hity, df_protonRecFD.PDc1Hitz])
        if detRes:
            df_protonRecFD.loc[:, "PDc3theta"] = getTheta([df_protonRecFD.PDc3Hitx, df_protonRecFD.PDc3Hity, df_protonRecFD.PDc3Hitz])
        best_params = [-53.14680163254601, 79.61307254040804, 0.3, 0.05739232362022314]
        df_protonRecFD_1 = df_protonRecFD.loc[df_protonRecFD.PDc1theta < corr(best_params, df_protonRecFD.Pp), :]
        df_protonRecFD_2 = df_protonRecFD.loc[df_protonRecFD.PDc1theta >= corr(best_params, df_protonRecFD.Pp), :]

        if detRes:
            df_protonRecCD.loc[:, "PCvt1r"] = mag([df_protonRecCD.PCvt1Hitx, df_protonRecCD.PCvt1Hity, df_protonRecCD.PCvt1Hitz])
            df_protonRecCD.loc[:, "PCvt1theta"] = getTheta([df_protonRecCD.PCvt1Hitx, df_protonRecCD.PCvt1Hity, df_protonRecCD.PCvt1Hitz])
            df_protonRecCD.loc[:, "PCvt1phi"] = getPhi([df_protonRecCD.PCvt1Hitx, df_protonRecCD.PCvt1Hity, df_protonRecCD.PCvt1Hitz])
            df_protonRecCD.loc[:, "PCvt3r"] = mag([df_protonRecCD.PCvt3Hitx, df_protonRecCD.PCvt3Hity, df_protonRecCD.PCvt3Hitz])
            df_protonRecCD.loc[:, "PCvt3theta"] = getTheta([df_protonRecCD.PCvt3Hitx, df_protonRecCD.PCvt3Hity, df_protonRecCD.PCvt3Hitz])
            df_protonRecCD.loc[:, "PCvt3phi"] = getPhi([df_protonRecCD.PCvt3Hitx, df_protonRecCD.PCvt3Hity, df_protonRecCD.PCvt3Hitz])
            df_protonRecCD.loc[:, "PCvt5r"] = mag([df_protonRecCD.PCvt5Hitx, df_protonRecCD.PCvt5Hity, df_protonRecCD.PCvt5Hitz])
            df_protonRecCD.loc[:, "PCvt5theta"] = getTheta([df_protonRecCD.PCvt5Hitx, df_protonRecCD.PCvt5Hity, df_protonRecCD.PCvt5Hitz])
            df_protonRecCD.loc[:, "PCvt5phi"] = getPhi([df_protonRecCD.PCvt5Hitx, df_protonRecCD.PCvt5Hity, df_protonRecCD.PCvt5Hitz])
            df_protonRecCD.loc[:, "PCvt7r"] = mag([df_protonRecCD.PCvt7Hitx, df_protonRecCD.PCvt7Hity, df_protonRecCD.PCvt7Hitz])
            df_protonRecCD.loc[:, "PCvt7theta"] = getTheta([df_protonRecCD.PCvt7Hitx, df_protonRecCD.PCvt7Hity, df_protonRecCD.PCvt7Hitz])
            df_protonRecCD.loc[:, "PCvt7phi"] = getPhi([df_protonRecCD.PCvt7Hitx, df_protonRecCD.PCvt7Hity, df_protonRecCD.PCvt7Hitz])
            df_protonRecCD.loc[:, "PCvt12r"] = mag([df_protonRecCD.PCvt12Hitx, df_protonRecCD.PCvt12Hity, df_protonRecCD.PCvt12Hitz])
            df_protonRecCD.loc[:, "PCvt12theta"] = getTheta([df_protonRecCD.PCvt12Hitx, df_protonRecCD.PCvt12Hity, df_protonRecCD.PCvt12Hitz])
            df_protonRecCD.loc[:, "PCvt12phi"] = getPhi([df_protonRecCD.PCvt12Hitx, df_protonRecCD.PCvt12Hity, df_protonRecCD.PCvt12Hitz])
        else:
            df_protonRecCD.loc[:, "PCvt12theta"] = getTheta([df_protonRecCD.PCvt12Hitx, df_protonRecCD.PCvt12Hity, df_protonRecCD.PCvt12Hitz])
            df_protonRecCD.loc[:, "PCvt12phi"] = getPhi([df_protonRecCD.PCvt12Hitx, df_protonRecCD.PCvt12Hity, df_protonRecCD.PCvt12Hitz])

        #inbending proton energy loss correction
        if pol == "inbending":
            const_FD = -0.00051894 - 0.00018104 * df_protonRecFD_1.Ptheta
            coeff_FD = 3.29466917*10**(-3) +  5.73663160*10**(-4) * df_protonRecFD_1.Ptheta - 1.40807209 * 10**(-5) * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta
            CorrectedPp_FD_1 = np.select([df_protonRecFD_1.Pp<1, df_protonRecFD_1.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pp"], np.exp(-2.739 - 3.932*df_protonRecFD_1.Pp) + 0.002907+df_protonRecFD_1.Pp])

            const_FD = -0.16742969 + 0.00697925 * df_protonRecFD_1.Ptheta
            coeff_FD = 0.23352115 - 0.01338697 * df_protonRecFD_1.Ptheta
            CorrectedPtheta_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Ptheta"]

            const_FD = 0.21192125 -0.0115175 * df_protonRecFD_1.Ptheta
            coeff_FD = -8.94307411*0.1 + 1.66349766*0.1 * df_protonRecFD_1.Ptheta -8.90617559*0.001 * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta + 1.64803754*0.0001 * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta * df_protonRecFD_1.Ptheta
            CorrectedPphi_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pphi"]

            const_FD = -3.03346359*10**(-1) + 1.83368163*10**(-2)*df_protonRecFD_2.Ptheta - 2.86486404*10**(-4)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            coeff_FD =  2.01023276*10**(-1) - 1.13312215*10**(-2)*df_protonRecFD_2.Ptheta + 1.82487916*10**(-4)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            CorrectedPp_FD_2 = np.select([df_protonRecFD_2.Pp<1, df_protonRecFD_2.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pp"], np.exp(-1.2 - 4.228*df_protonRecFD_2.Pp) + 0.007502+df_protonRecFD_2.Pp])

            const_FD = 2.04334532 * 10 -1.81052405 * df_protonRecFD_2.Ptheta + 5.32556360*0.01 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta -5.23157558*0.0001 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta
            coeff_FD = 8.74233279 -7.63869344 * 0.1 * df_protonRecFD_2.Ptheta + 2.22376362*0.01 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta -2.16457260*0.0001 * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta * df_protonRecFD_2.Ptheta
            CorrectedPtheta_FD_2 = const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"]/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Ptheta"]

            const_FD = 0.54697831 -0.04896981*df_protonRecFD_2.Ptheta +  0.00111376*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            coeff_FD = -4.06733541*10**2 + 2.43696202*10*df_protonRecFD_2.Ptheta -3.36144736*10**(-1)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            coeff2_FD = 2.06378660*10 - 1.42866062*df_protonRecFD_2.Ptheta + 2.01085440*10**(-2)*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            CorrectedPphi_FD_2 = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD_2.loc[:, "Pp"]) + df_protonRecFD_2.loc[:, "Pphi"]

            #CD part
            const_CD = 1.93686914 - 0.116288824*df_protonRecCD.Ptheta + 0.00223685833*df_protonRecCD.Ptheta**2 - 1.40771969 * 10**(-5)*df_protonRecCD.Ptheta**3
            coeff_CD = -0.738047800 + 0.0443343685*df_protonRecCD.Ptheta - 8.50985972*10**(-4)*df_protonRecCD.Ptheta*df_protonRecCD.Ptheta + 5.36810280 * 10**(-6) * df_protonRecCD.Ptheta**3

            CorrectedPp_CD = const_CD + coeff_CD/df_protonRecCD.loc[:, "Pp"] + df_protonRecCD.loc[:, "Pp"]

            const_CD = -1.09849291*100 + 8.86664014 * df_protonRecCD.Ptheta - 0.26643881 * df_protonRecCD.Ptheta**2 + 3.53814210 * 10**(-3) * df_protonRecCD.Ptheta**3 - 1.75297107 * 10**(-5) * df_protonRecCD.Ptheta**4
            coeff_CD = 9.52034523*100 -5.74808292 * 10 * df_protonRecCD.Ptheta + 1.15386949 * df_protonRecCD.Ptheta**2 - 7.57970373 * 0.001 * df_protonRecCD.Ptheta**3
            coeff2_CD = -2.00387313*100 + 1.18979079 * 10 * df_protonRecCD.Ptheta - 2.37730217*0.1 * df_protonRecCD.Ptheta**2 + 1.55153003*0.001*df_protonRecCD.Ptheta**3

            CorrectedPtheta_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Ptheta"]

            const_CD = 4.94546178 -3.26662886*0.1 * df_protonRecCD.Ptheta +  7.39069603 * 0.001 * df_protonRecCD.Ptheta**2 -6.83599356*10**(-5) * df_protonRecCD.Ptheta**3 + 2.12303103*10**(-7) * df_protonRecCD.Ptheta**4
            coeff_CD = 1.72181613*10**(5) -1.36827111*10**(4) * df_protonRecCD.Ptheta + 4.00923146*10**(2) * df_protonRecCD.Ptheta**2 - 5.12792347 * df_protonRecCD.Ptheta**3 + 2.41793167*10**(-2) * df_protonRecCD.Ptheta**4
            coeff2_CD =  1.20477219*10**(2) -5.86630228 * df_protonRecCD.Ptheta + 7.44007875*10**(-2) * df_protonRecCD.Ptheta**2 -2.42652473*10**(-4) * df_protonRecCD.Ptheta**3
            CorrectedPphi_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Pphi"]

            correction = True
        #outbending proton energy loss correction
        elif pol == "outbending":
            #FD part
            const_FD = 0.05083242 -0.00469777*df_protonRecFD_1.Ptheta + 0.0001082*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
            coeff_FD = -1.47443264*0.01 + 1.58220893*0.001*df_protonRecFD_1.Ptheta -3.19490013*0.00001*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
            CorrectedPp_FD_1 = np.select([df_protonRecFD_1.Pp<1, df_protonRecFD_1.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pp"], np.exp(-2.739 - 3.932*df_protonRecFD_1.Pp) + 0.002907 + df_protonRecFD_1.Pp])

            const_FD = -2.56460305*10 + 3.29877542*df_protonRecFD_1.Ptheta -1.43106886*0.1*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta + 2.08341898*0.001*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
            coeff_FD =  9.12532740*10 -1.20100762*10*df_protonRecFD_1.Ptheta + 5.27654711*0.1*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta -7.72656759*0.001*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
            CorrectedPtheta_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Ptheta"]

            const_FD = -20.4780893 + 1.67020488*df_protonRecFD_1.Ptheta - 0.03419348*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
            coeff_FD = 35.02807194 - 2.9098043*df_protonRecFD_1.Ptheta +  0.06037906*df_protonRecFD_1.Ptheta*df_protonRecFD_1.Ptheta
            CorrectedPphi_FD_1 = const_FD + coeff_FD/df_protonRecFD_1.loc[:, "Pp"]/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pphi"]

            const_FD = 0.09832589 -0.0066463*df_protonRecFD_2.Ptheta + 0.00010312*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            coeff_FD = -9.61421691*0.01 + 6.85638807*0.001*df_protonRecFD_2.Ptheta -9.75766427*0.00001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            CorrectedPp_FD_2 = np.select([df_protonRecFD_2.Pp<1, df_protonRecFD_2.Pp>=1], [const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pp"], np.exp(-1.871 - 3.063*df_protonRecFD_2.Pp) + 0.007517 + df_protonRecFD_2.Pp])

            const_FD = -1.68873940 + 9.56867163*0.01*df_protonRecFD_2.Ptheta -1.43741464*0.001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            coeff_FD = 1.49978357*10 -1.40137094*df_protonRecFD_2.Ptheta + 4.38501543*0.01*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta -4.57982872*0.0001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            CorrectedPtheta_FD_2 = const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"]/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Ptheta"]

            const_FD = 6.75359137 - 0.43199851*df_protonRecFD_2.Ptheta + 0.0068995*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            coeff_FD = -1.68588219 + 1.05609627*0.1*df_protonRecFD_2.Ptheta -1.50452832*0.001*df_protonRecFD_2.Ptheta*df_protonRecFD_2.Ptheta
            CorrectedPphi_FD_2 = const_FD + coeff_FD/df_protonRecFD_2.loc[:, "Pp"]/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pphi"]
            #CD part
            const_CD = 1.92657376 - 0.113836734*df_protonRecCD.Ptheta + 0.00215038526*df_protonRecCD.Ptheta**2 - 1.32525053 * 10**(-5)*df_protonRecCD.Ptheta**3
            coeff_CD = -0.755650043 + 0.0445538936*df_protonRecCD.Ptheta - 8.38241864*10**(-4)*df_protonRecCD.Ptheta*df_protonRecCD.Ptheta + 5.16887255 * 10**(-6) * df_protonRecCD.Ptheta**3

            CorrectedPp_CD = const_CD + coeff_CD/df_protonRecCD.loc[:, "Pp"] + df_protonRecCD.loc[:, "Pp"]

            const_CD = -5.79024055*10 + 4.67197531 * df_protonRecCD.Ptheta - 0.140156897 * df_protonRecCD.Ptheta**2 + 1.85853057 * 10**(-3) * df_protonRecCD.Ptheta**3 - 9.19989908 * 10**(-6) * df_protonRecCD.Ptheta**4
            coeff_CD = 2.99700765*1000 - 2.18027982 * 10**2 * df_protonRecCD.Ptheta + 5.84757503 * df_protonRecCD.Ptheta**2 - 6.80409195 * 0.01 * df_protonRecCD.Ptheta**3 + 2.89244618 * 0.0001 * df_protonRecCD.Ptheta**4
            coeff2_CD = -1.82237904*100 + 1.10153549 * 10 * df_protonRecCD.Ptheta - 2.24699931*0.1 * df_protonRecCD.Ptheta**2 + 1.49390960*0.001*df_protonRecCD.Ptheta**3

            CorrectedPtheta_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Ptheta"]

            const_CD = 7.58761670 - 5.28224578*0.1 * df_protonRecCD.Ptheta +  1.31580117 * 0.01 * df_protonRecCD.Ptheta**2 -1.41738951*10**(-4) * df_protonRecCD.Ptheta**3 + 5.62884363*10**(-7) * df_protonRecCD.Ptheta**4
            coeff_CD = 1.07644097*10**(5) - 8.67994639*10**(3) * df_protonRecCD.Ptheta + 2.57187193*10**(2) * df_protonRecCD.Ptheta**2 - 3.31379317 * df_protonRecCD.Ptheta**3 + 1.56896621*10**(-2) * df_protonRecCD.Ptheta**4
            coeff2_CD =  1.92263184*10**(2) -1.00870704 * 10 * df_protonRecCD.Ptheta + 1.56575252*10**(-1) * df_protonRecCD.Ptheta**2 -7.71489734*10**(-4) * df_protonRecCD.Ptheta**3
            CorrectedPphi_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Pphi"]

            correction = True
        else:
            print("no correction applied")

        if correction:
            print("energy loss correction applied for " + pol)
            #proton correction
            df_protonRecCD.loc[:, "Pp"] = CorrectedPp_CD + 0.01
            df_protonRecCD.loc[:, "Ptheta"] = CorrectedPtheta_CD - 0.002129*CorrectedPtheta_CD**2 + 0.198*CorrectedPtheta_CD - 4.762 -0.2/(1+np.exp((CorrectedPp_CD-0.55)/(-0.05)))
            df_protonRecCD.loc[:, "Pphi"] = CorrectedPphi_CD

            df_protonRecFD_1.loc[:, "Pp"] = CorrectedPp_FD_1
            df_protonRecFD_1.loc[:, "Ptheta"] = CorrectedPtheta_FD_1
            df_protonRecFD_1.loc[:, "Pphi"] = CorrectedPphi_FD_1
            df_protonRecFD_1.loc[:, "Pband"] = "lower"

            df_protonRecFD_2.loc[:, "Pp"] = CorrectedPp_FD_2
            df_protonRecFD_2.loc[:, "Ptheta"] = CorrectedPtheta_FD_2
            df_protonRecFD_2.loc[:, "Pphi"] = CorrectedPphi_FD_2
            df_protonRecFD_2.loc[:, "Pband"] = "upper"

            df_protonRecFD = pd.concat([df_protonRecFD_1, df_protonRecFD_2])

            #correcting reconstruction bias after proton energy loss correction
            print("applying the proton kinematic corrections for " + pol)
            if pol == "inbending":
                corr = np.poly1d([1.671, -4.918, 5.151, -2.434])(df_protonRecFD.Pp)       
                corr = np.where(corr<0, corr, 0)
                df_protonRecFD.loc[:, "Ptheta"] = df_protonRecFD.Ptheta + corr #corr scale -1 to -0.3 degrees
            if pol == "outbending":
                df_protonRecFD.loc[df_protonRecFD.Psector<7, "Pp"] = df_protonRecFD.loc[df_protonRecFD.Psector<7, "Pp"] - 0.02
                df_protonRecFD.loc[:, "Ptheta"] = df_protonRecFD.Ptheta + 0.05*(np.abs(df_protonRecFD.Ptheta - 27) + (df_protonRecFD.Ptheta - 27))

            df_protonRec = pd.concat([df_protonRecFD, df_protonRecCD, df_protonRecOthers])

            #moduli proton phi
            df_protonRec.loc[:, "Pphi"] = np.where(df_protonRec.loc[:, "Pphi"]%360<180, df_protonRec.loc[:, "Pphi"]%360, df_protonRec.loc[:, "Pphi"]%360-360)

            df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
            df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
            df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
            pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

            df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

        gam = [df_gammaRec['Gpx'], df_gammaRec['Gpy'], df_gammaRec['Gpz']]
        df_gammaRec.loc[:, 'Gp'] = mag(gam)
        df_gammaRec.loc[:,'GSamplFrac'] = df_gammaRec.Gedep/ df_gammaRec.Gp

        df_gg = pd.merge(df_gammaRec, df_gammaRec,
                         how='outer', on='event', suffixes=("", "2"))
        df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
        df_gg = df_gg.drop(['GIndex', 'GIndex2'], axis = 1)

        if correction:
            print("applying the photon kinematic corrections for " + pol)

            #photon kinematic correction of df_gg
            gam = [df_gg['Gpx'], df_gg['Gpy'], df_gg['Gpz']]
            df_gg.loc[:, 'Gp'] = mag(gam)
            df_gg.loc[:, 'Gtheta'] = getTheta(gam)
            df_gg.loc[:, 'Gphi'] = getPhi(gam)
            #photon kinematic correction of df_gamma
            gam = [df_gammaRec['Gpx'], df_gammaRec['Gpy'], df_gammaRec['Gpz']]
            df_gammaRec.loc[:, 'Gp'] = mag(gam)
            df_gammaRec.loc[:, 'Gtheta'] = getTheta(gam)
            df_gammaRec.loc[:, 'Gphi'] = getPhi(gam)

            def cubic(args, x):
                a, b, c = args
                return a*x**3 + b*x**2 + c*x

            def quintic(args, x):
                a, b, c = args
                if b < 0:
                    return 0*x
                return a*x*(x-b)**3 * (x-c)

            def quartic(args, x):
                a, b, c, d = args
                x = np.array(x)
                return a*x**4 + b*x**3 + c*x**2 + d*x**1
            #FT - df_gg: perform correction for only one photon
            FT_phot_corr = (-0.00467*df_gg.loc[df_gg["Gsector"]>7, "Gp"]**2 + 0.0802 *df_gg.loc[df_gg["Gsector"]>7, "Gp"]  -0.352) + 0.25
            df_gg.loc[df_gg["Gsector"]>7, "Gp"] = df_gg.loc[df_gg["Gsector"]>7, "Gp"] + np.where(FT_phot_corr>0, FT_phot_corr, 0)
            #FT - df_gammaRec: perform every photon
            FT_phot_corr = (-0.00467*df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]**2 + 0.0802 *df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"]  -0.352) + 0.25
            df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"] = df_gammaRec.loc[df_gammaRec["Gsector"]>7, "Gp"] + np.where(FT_phot_corr>0, FT_phot_corr, 0)

            #FD
            if pol == "inbending":
                args = [[-0.0000732, 1.480, 9.344], [-0.000135, 3.070, 9.248], [-0.0000437, 0.719, 9.873], [-0.0000428, 0.00234, 0.0103], [0.000250, -0.00314, 0.0232], [-0.0000454, 0.517, 9.447]]
                funcs = [quintic, quintic, quintic, cubic, cubic, quintic]
                funcs_minor = [quintic, quintic, quintic, cubic, cubic, cubic]
                args_minor = [[-0.0000168, 0.821, 8.894], [-0.0000340, 2.720, 8.419], [-0.0000620, 2.793, 8.865], [ 0.000132, -0.00162,  0.00978], [-0.000135,  0.000282, 0.00650], [ 0.000263,  -0.00293,   0.0139]]
                for sector in range(1, 7):
                    cond = df_gg.Gsector == sector
                    FD_phot_corr_sector = funcs[sector-1](args[sector-1], df_gg.loc[cond, "Gp"])
                    df_gg.loc[cond, "Gp"] = df_gg.loc[cond, "Gp"] + FD_phot_corr_sector
                    FD_phot_corr_minor_sector = funcs_minor[sector-1](args_minor[sector-1], df_gg.loc[cond, "Gp"])
                    df_gg.loc[cond, "Gp"] = df_gg.loc[cond, "Gp"] + FD_phot_corr_minor_sector

                    cond = df_gammaRec.Gsector == sector
                    FD_phot_corr_sector = funcs[sector-1](args[sector-1], df_gammaRec.loc[cond, "Gp"])
                    df_gammaRec.loc[cond, "Gp"] = df_gammaRec.loc[cond, "Gp"] + FD_phot_corr_sector
                    FD_phot_corr_minor_sector = funcs_minor[sector-1](args_minor[sector-1], df_gammaRec.loc[cond, "Gp"])
                    df_gammaRec.loc[cond, "Gp"] = df_gammaRec.loc[cond, "Gp"] + FD_phot_corr_minor_sector

            if pol == "outbending":
                args = [[-0.000615,  0.0113, -0.0600,   0.115],[-0.000334,  0.00656, -0.0383,  0.0934],[-0.000911,  0.0157, -0.0806,  0.154],[ 0.000117, -0.000905,  0.00215,  0.0331],[-0.000119,  0.000979, -0.00400,  0.0499],[-0.000893,  0.0131, -0.0580,  0.111]]
                for sector in range(1, 7):
                    cond = df_gg.Gsector == sector
                    FD_phot_corr_sector = quartic(args[sector-1], df_gg.loc[cond, "Gp"])/(1+np.exp(-(df_gg.loc[cond, "Gp"]-2.2)/0.15))
                    df_gg.loc[cond, "Gp"] = df_gg.loc[cond, "Gp"] + FD_phot_corr_sector

                    cond = df_gammaRec.Gsector == sector
                    FD_phot_corr_sector = quartic(args[sector-1], df_gammaRec.loc[cond, "Gp"])/(1+np.exp(-(df_gammaRec.loc[cond, "Gp"]-2.2)/0.15))
                    df_gammaRec.loc[cond, "Gp"] = df_gammaRec.loc[cond, "Gp"] + FD_phot_corr_sector

            df_gg.loc[:, "Gpx"] = df_gg.loc[:, "Gp"]*np.sin(np.radians(df_gg.loc[:, "Gtheta"]))*np.cos(np.radians(df_gg.loc[:, "Gphi"]))
            df_gg.loc[:, "Gpy"] = df_gg.loc[:, "Gp"]*np.sin(np.radians(df_gg.loc[:, "Gtheta"]))*np.sin(np.radians(df_gg.loc[:, "Gphi"]))
            df_gg.loc[:, "Gpz"] = df_gg.loc[:, "Gp"]*np.cos(np.radians(df_gg.loc[:, "Gtheta"]))
            df_gg.loc[:,'GSamplFrac'] = df_gg.Gedep/ df_gg.Gp

            df_gammaRec.loc[:, "Gpx"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.cos(np.radians(df_gammaRec.loc[:, "Gphi"]))
            df_gammaRec.loc[:, "Gpy"] = df_gammaRec.loc[:, "Gp"]*np.sin(np.radians(df_gammaRec.loc[:, "Gtheta"]))*np.sin(np.radians(df_gammaRec.loc[:, "Gphi"]))
            df_gammaRec.loc[:, "Gpz"] = df_gammaRec.loc[:, "Gp"]*np.cos(np.radians(df_gammaRec.loc[:, "Gtheta"]))
            df_gammaRec.loc[:,'GSamplFrac'] = df_gammaRec.Gedep/ df_gammaRec.Gp

        # proton fiducial cuts
        if nofid:
            df_protonRec.loc[:, "PFid"] = 1
        else:
            df_protonRec.loc[:, "PFid"] = 0

            df_protonRec.loc[df_protonRec.Psector<7, "PFid"] = 1 #FD fid done by previous pipeline

            cut_CD = df_protonRec.Psector > 7
            cut_right = cut_CD & (df_protonRec.Ptheta<64.23)
            cut_bottom = cut_CD & (df_protonRec.PCvt12theta>44.5)
            cut_sidel = cut_CD & (df_protonRec.PCvt12theta<-2.942 + 1.274*df_protonRec.Ptheta)
            cut_sider = cut_CD & (df_protonRec.PCvt12theta>-3.523 + 1.046*df_protonRec.Ptheta)

            cut_trapezoid = cut_CD & cut_right & cut_bottom & cut_sidel & cut_sider

            cut_gaps1 = ~((df_protonRec.PCvt12phi>-95) & (df_protonRec.PCvt12phi<-80))
            cut_gaps2 = ~((df_protonRec.PCvt12phi>25) & (df_protonRec.PCvt12phi<40))
            cut_gaps3 = ~((df_protonRec.PCvt12phi>143) & (df_protonRec.PCvt12phi<158))
            cut_gaps = cut_CD & cut_gaps1 & cut_gaps2 & cut_gaps3
            cut_total = cut_gaps & cut_trapezoid

            df_protonRec.loc[cut_total, "PFid"] = 1 #CD fid
            
            if pol == "inbending":
                pchi2CD_lb,   pchi2CD_ub   = -5.000, 6.345
                pchi2FD_S1_lb, pchi2FD_S1_ub = -3.296, 3.508
                pchi2FD_S2_lb, pchi2FD_S2_ub = -3.552, 4.000
                pchi2FD_S3_lb, pchi2FD_S3_ub = -3.446, 3.937
                pchi2FD_S4_lb, pchi2FD_S4_ub = -2.747, 3.190
                pchi2FD_S5_lb, pchi2FD_S5_ub = -2.851, 3.418
                pchi2FD_S6_lb, pchi2FD_S6_ub = -3.174, 3.514
            elif pol == "outbending":
                pchi2CD_lb,   pchi2CD_ub   = -5.592,  6.785
                pchi2FD_S1_lb, pchi2FD_S1_ub = -3.905, 4.088
                pchi2FD_S2_lb, pchi2FD_S2_ub = -3.411, 3.939
                pchi2FD_S3_lb, pchi2FD_S3_ub = -4.042, 5.954
                pchi2FD_S4_lb, pchi2FD_S4_ub = -3.820, 5.065
                pchi2FD_S5_lb, pchi2FD_S5_ub = -3.384, 4.232
                pchi2FD_S6_lb, pchi2FD_S6_ub = -5.077, 5.100

            df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.Pchi2pid<pchi2CD_lb)   | (df_protonRec.Pchi2pid>pchi2CD_ub)  ), "PFid"] = 0
            df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.Pchi2pid<pchi2FD_S1_lb) | (df_protonRec.Pchi2pid>pchi2FD_S1_ub)), "PFid"] = 0
            df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.Pchi2pid<pchi2FD_S2_lb) | (df_protonRec.Pchi2pid>pchi2FD_S2_ub)), "PFid"] = 0
            df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.Pchi2pid<pchi2FD_S3_lb) | (df_protonRec.Pchi2pid>pchi2FD_S3_ub)), "PFid"] = 0
            df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.Pchi2pid<pchi2FD_S4_lb) | (df_protonRec.Pchi2pid>pchi2FD_S4_ub)), "PFid"] = 0
            df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.Pchi2pid<pchi2FD_S5_lb) | (df_protonRec.Pchi2pid>pchi2FD_S5_ub)), "PFid"] = 0
            df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.Pchi2pid<pchi2FD_S6_lb) | (df_protonRec.Pchi2pid>pchi2FD_S6_ub)), "PFid"] = 0

        if detRes:
            df_gg = df_gg.loc[:, ~df_gg.columns.duplicated()]
            df_gg.loc[:, "Gedep2_tot"] = df_gg.Gedep12 + df_gg.Gedep22 + df_gg.Gedep32
        # else:
        #     df_protonRec = df_protonRec.drop(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PDc1theta", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz", "PCvt12theta", "PCvt12phi"], axis = 1)
        #     df_gammaRec = df_gammaRec.drop(["GcX", "GcY"], axis = 1)
        #     df_gg = df_gg.drop(["GcX", "GcY", "GcX2", "GcY2"], axis = 1)
        
        df_ep = pd.merge(df_electronRec, df_protonRec, how='outer', on='event')

        df_epgg = pd.merge(df_ep, df_gg, how='outer', on='event')
        df_epgg = df_epgg.loc[~np.isnan(df_epgg["Ppx"]), :]
        df_epgg = df_epgg.loc[~np.isnan(df_epgg["Gpx"]), :]
        df_epgg = df_epgg.loc[~np.isnan(df_epgg["Gpx2"]), :]

        self.df_epgg = df_epgg # saves df_epgg

        df_epg = pd.merge(df_ep, df_gammaRec, how='outer', on='event')
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
        v3pi0 = cross(VGS, pi0)

        VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
                    df_epgg["Ppy"], pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]
        VmissP = [-df_epgg["Epx"] - df_epgg["Gpx"] - df_epgg["Gpx2"], -df_epgg["Epy"] -
                    df_epgg["Gpy"] - df_epgg["Gpy2"], pbeam - df_epgg["Epz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        Vmiss = [-df_epgg["Epx"] - df_epgg["Ppx"] - df_epgg["Gpx"] - df_epgg["Gpx2"],
                    -df_epgg["Epy"] - df_epgg["Ppy"] - df_epgg["Gpy"] - df_epgg["Gpy2"],
                    pbeam - df_epgg["Epz"] - df_epgg["Ppz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        costheta = cosTheta(VGS, gam)

        df_epgg.loc[:, 'Mpx'], df_epgg.loc[:, 'Mpy'], df_epgg.loc[:, 'Mpz'] = Vmiss

        # binning kinematics
        df_epgg.loc[:,'Q2'] = -((ebeam - df_epgg['Ee'])**2 - mag2(VGS))
        df_epgg.loc[:,'nu'] = (ebeam - df_epgg['Ee'])
        df_epgg.loc[:,'xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
        df_epgg.loc[:,'t1'] = 2 * M * (df_epgg['Pe'] - M)
        df_epgg.loc[:,'t2'] = (M * df_epgg['Q2'] + 2 * M * df_epgg['nu'] * (df_epgg['nu'] - np.sqrt(df_epgg['nu'] * df_epgg['nu'] + df_epgg['Q2']) * costheta))\
        / (M + df_epgg['nu'] - np.sqrt(df_epgg['nu'] * df_epgg['nu'] + df_epgg['Q2']) * costheta)
        df_epgg.loc[:,'W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epgg['Ee'])**2 - mag2(VGS)))
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
        df_epgg.loc[:,'coplanarity'] = angle(v3h, v3pi0)
        df_epgg.loc[:,'coneAngle1'] = angle(ele, gam)
        df_epgg.loc[:,'coneAngle2'] = angle(ele, gam2)
        
        df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)

        self.df_epgg = df_epgg

    def makeDVpi0P_DVCS(self, pol = "inbending"):
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

        # proton reconstruction quality
        # cut_FD_proton = (df_epgg.loc[:, "Psector"]<7) & (df_epgg.loc[:, "Ptheta"]<35)
        # cut_CD_proton = (df_epgg.loc[:, "Psector"]>7) & (df_epgg.loc[:, "Ptheta"]>45) & (df_epgg.loc[:, "Ptheta"]<65)
        # cut_proton = (cut_FD_proton)|(cut_CD_proton)
        cut_proton = 1

        df_dvpi0p.loc[:, "config"] = 0
        
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

        eps = 2*M*df_epg.xB / np.sqrt(df_epg.Q2)
        df_epg.loc[:,'ycol1'] = (df_epg.Q2-df_epg.t2)/(df_epg.Q2-df_epg.xB*df_epg.t2)
        df_epg.loc[:,'ycol2'] = 1 - (1-df_epg.xB)*df_epg.t2/df_epg.Q2
        df_epg.loc[:,'ymax1'] = 2*(np.sqrt(1+eps**2)-1)/(eps**2)
        df_epg.loc[:,'ymax2'] = 1 - (M**2)*(df_epg.xB**2)/df_epg.Q2
        df_epg.loc[:,'tmin1'] = df_epg.Q2*(2*(1-df_epg.xB)*(1-np.sqrt(1+eps**2))+eps**2)/(4*df_epg.xB*(1-df_epg.xB) + eps**2)
        df_epg.loc[:,'tmin2'] = M*M*(df_epg.xB**2)/(1-df_epg.xB+df_epg.xB*M*M/df_epg.Q2)
        df_epg.loc[:,'tcol'] = df_epg.Q2*(df_epg.Q2-2*df_epg.xB*M*ebeam)/df_epg.xB/(df_epg.Q2-2*M*ebeam)

        df_epg.loc[:, 'vzdiff'] = df_epg.Evz - df_epg.Pvz

        # encode unassigned bin as -1
        df_epg.loc[:, "Q2bin"] = -1
        df_epg.loc[:, "xBbin"] = -1
        df_epg.loc[:, "tbin"] = -1
        # df_epg.loc[:, "tbin2"] = -1
        df_epg.loc[:, "phibin"] = -1
        # df_epg.loc[:, "phibin2"] = -1
        df_epg.loc[:, "Q2xBbin"] = -1
        df_epg.loc[:, "Q2xBtbin"] = -1
        # df_epg.loc[:, "Q2xBtbin2"] = -1
        df_epg.loc[:, "Q2xBtphibin"] = -1
        Q2xBbin = 0

        # encode all binning
        for Q2bin in range(len(Q2bin_i)):
            #square Q2 binning
            df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]), "Q2bin"] = Q2bin
            #adaptive xB binning
            for xBbin in range(len(xBbin_i[Q2bin])):
                if Q2bin < len(Q2bin_i) -1:
                    if xBbin == 0:
                        df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin #0
                        df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin #0
                    elif xBbin < len(xBbin_i[Q2bin])-1:
                        df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin
                        df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin
                    else:
                        df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "xBbin"] = xBbin
                        df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "Q2xBbin"] = Q2xBbin
                else:
                    df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB)& (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "xBbin"] = xBbin
                    df_epg.loc[(df_epg.Q2>=Q2bin_i[Q2bin]) & (df_epg.Q2<Q2bin_f[Q2bin]) & (df_epg.Q2<=2*M*(10.604-2)*df_epg.xB)& (df_epg.Q2>=(4-M*M)*df_epg.xB/(1-df_epg.xB)), "Q2xBbin"] = Q2xBbin #0
                Q2xBbin = Q2xBbin + 1
        for tbin in range(len(tbin_i)):
            #square t binning
            df_epg.loc[(df_epg.t1>=tbin_i[tbin]) & (df_epg.t1<tbin_f[tbin]), "tbin"] = tbin
            # df_epg.loc[(df_epg.t2>=tbin_i[tbin]) & (df_epg.t2<tbin_f[tbin]), "tbin2"] = tbin
        for phibin in range(len(phibin_i)):
            #square phi binning
            df_epg.loc[(df_epg.phi1>=phibin_i[phibin]) & (df_epg.phi1<phibin_f[phibin]), "phibin"] = phibin
            # df_epg.loc[(df_epg.phi2>=phibin_i[phibin]) & (df_epg.phi2<phibin_f[phibin]), "phibin2"] = phibin

        df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBtbin"] = len(tbin_i) * df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBbin"] + df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "tbin"]
        # df_epg.loc[(df_epg.Q2bin>0)&(df_epg.xBbin>0)&(df_epg.tbin2>0), "Q2xBtbin2"] = df_epg.Q2bin.astype(str) + df_epg.xBbin.astype(str) + df_epg.tbin2.astype(str)
        df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBtphibin"] = len(phibin_i) * df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "Q2xBtbin"] + df_epg.loc[(df_epg.Q2xBbin>=0)&(df_epg.tbin>=0), "phibin"]

        df_epg = df_epg.astype({"Q2bin": int, "xBbin": int, "tbin": int, "phibin": int, "Q2xBbin": int, "Q2xBtbin": int, "Q2xBtphibin": int})

        self.df_epg = df_epg

    def makeDVCS(self, pol = "inbending"):
        #make dvcs pairs
        df_dvcs = self.df_epg

        #common cuts
        cut_xBupper = df_dvcs["xB"] < 1  # xB
        cut_xBlower = df_dvcs["xB"] > 0  # xB
        cut_Q2 = df_dvcs["Q2"] > 1  # Q2
        cut_W = df_dvcs["W"] > 2  # W
        cut_Ee = df_dvcs["Ee"] > 2  # Ee
        cut_Ge = df_dvcs["Ge"] > 2  # Ge
        cut_Esector = (df_dvcs["Esector"]!=df_dvcs["Gsector"])
        cut_Psector = ~( ((df_dvcs["Pstat"]//10)%10>0) & (df_dvcs["Psector"]==df_dvcs["Gsector"]))
        cut_Ppmax = df_dvcs.Pp < 1.6  # Pp
        cut_Pthetamin = df_dvcs.Ptheta > 0 #Ptheta
        cut_Trigger = ((df_dvcs.TriggerBit & 1 << 1) > 0) | ((df_dvcs.TriggerBit & 1 << 2) > 0) | ((df_dvcs.TriggerBit & 1 << 3) > 0) | ((df_dvcs.TriggerBit & 1 << 4) > 0) | ((df_dvcs.TriggerBit & 1 << 5) > 0) | ((df_dvcs.TriggerBit & 1 << 6) > 0)
        cut_EFid = df_dvcs.EFid == 1
        # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
        cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin & cut_Trigger & cut_EFid

        df_dvcs = df_dvcs[cut_common]

        # proton reconstruction quality
        # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) & (df_dvcs.loc[:, "Ptheta"]<35)
        # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) & (df_dvcs.loc[:, "Ptheta"]>45) & (df_dvcs.loc[:, "Ptheta"]<65)
        # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) #& (df_dvcs.loc[:, "Ptheta"]<37)
        # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) #& (df_dvcs.loc[:, "Ptheta"]<66) #& (df_dvcs.loc[:, "Ptheta"]>40) 
        # cut_proton = (cut_FD_proton)|(cut_CD_proton)
        #(cut_FD_proton)|(cut_CD_proton)

        df_dvcs.loc[:, "config"] = 0

        cuts_dvcs_CDFT_Inb = self.cuts_dvcs_CDFT_Inb 
        cuts_dvcs_CD_Inb = self.cuts_dvcs_CD_Inb 
        cuts_dvcs_FD_Inb = self.cuts_dvcs_FD_Inb
        cuts_dvcs_CDFT_Outb = self.cuts_dvcs_CDFT_Outb 
        cuts_dvcs_CD_Outb = self.cuts_dvcs_CD_Outb 
        cuts_dvcs_FD_Outb = self.cuts_dvcs_FD_Outb

        if pol == "inbending":
            vzdiffCD_lb,    vzdiffCD_ub    = -2.011, 2.314
            vzdiffFD_S1_lb, vzdiffFD_S1_ub = -3.209, 4.017
            vzdiffFD_S2_lb, vzdiffFD_S2_ub = -3.612, 4.139
            vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.328, 4.287
            vzdiffFD_S4_lb, vzdiffFD_S4_ub = -3.411, 4.108
            vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.607, 4.246
            vzdiffFD_S6_lb, vzdiffFD_S6_ub = -2.999, 3.927

            df_dvcs.loc[ (df_dvcs.Psector>4000) & ((df_dvcs.vzdiff<vzdiffCD_lb)   | (df_dvcs.vzdiff>vzdiffCD_ub)  ), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==1)   & ((df_dvcs.vzdiff<vzdiffFD_S1_lb) | (df_dvcs.vzdiff>vzdiffFD_S1_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==2)   & ((df_dvcs.vzdiff<vzdiffFD_S2_lb) | (df_dvcs.vzdiff>vzdiffFD_S2_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==3)   & ((df_dvcs.vzdiff<vzdiffFD_S3_lb) | (df_dvcs.vzdiff>vzdiffFD_S3_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==4)   & ((df_dvcs.vzdiff<vzdiffFD_S4_lb) | (df_dvcs.vzdiff>vzdiffFD_S4_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==5)   & ((df_dvcs.vzdiff<vzdiffFD_S5_lb) | (df_dvcs.vzdiff>vzdiffFD_S5_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==6)   & ((df_dvcs.vzdiff<vzdiffFD_S6_lb) | (df_dvcs.vzdiff>vzdiffFD_S6_ub)), "PFid"] = 0

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
            vzdiffCD_lb,    vzdiffCD_ub    = -2.737, 2.096
            vzdiffFD_S1_lb, vzdiffFD_S1_ub = -4.435, 3.429
            vzdiffFD_S2_lb, vzdiffFD_S2_ub = -4.646, 2.978
            vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.922, 3.040
            vzdiffFD_S4_lb, vzdiffFD_S4_ub = -4.646, 3.493
            vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.901, 3.750
            vzdiffFD_S6_lb, vzdiffFD_S6_ub = -3.846, 3.623

            df_dvcs.loc[ (df_dvcs.Psector>4000) & ((df_dvcs.vzdiff<vzdiffCD_lb)   | (df_dvcs.vzdiff>vzdiffCD_ub)  ), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==1)   & ((df_dvcs.vzdiff<vzdiffFD_S1_lb) | (df_dvcs.vzdiff>vzdiffFD_S1_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==2)   & ((df_dvcs.vzdiff<vzdiffFD_S2_lb) | (df_dvcs.vzdiff>vzdiffFD_S2_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==3)   & ((df_dvcs.vzdiff<vzdiffFD_S3_lb) | (df_dvcs.vzdiff>vzdiffFD_S3_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==4)   & ((df_dvcs.vzdiff<vzdiffFD_S4_lb) | (df_dvcs.vzdiff>vzdiffFD_S4_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==5)   & ((df_dvcs.vzdiff<vzdiffFD_S5_lb) | (df_dvcs.vzdiff>vzdiffFD_S5_ub)), "PFid"] = 0
            df_dvcs.loc[ (df_dvcs.Psector==6)   & ((df_dvcs.vzdiff<vzdiffFD_S6_lb) | (df_dvcs.vzdiff>vzdiffFD_S6_ub)), "PFid"] = 0
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

        # dealing with duplicates
        df_dvcs = df_dvcs.sort_values(by=['reconGam', 'Ge', 'Pe'], ascending = [True, False, False])
        df_dvcs = df_dvcs.loc[~df_dvcs.event.duplicated(), :]
        df_dvcs = df_dvcs.sort_values(by='event')

        self.df_dvcs = df_dvcs                        

    def pi02gSubtraction(self):
        #exclude dvpi0 from dvcs. use only when both set up.
        df_epg = self.df_epg
        pi0to2gammas = df_epg["event"].isin(self.df_dvpi0p["event"])
        df_epg = df_epg[~pi0to2gammas]
        self.df_epg = df_epg

    def save(self):
        # df_x = self.df_dvcs
        # df_protonDet = self.df_protonDet
        # df = pd.merge(df_x, df_protonDet, how = 'inner', on ='event')
        self.df = self.df_dvcs

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
    parser.add_argument("-nf","--nofid", help="no additional fiducial cuts", action = "store_true")


    args = parser.parse_args()

    if args.entry_start:
        args.entry_start = int(args.entry_start)
    if args.entry_stop:
        args.entry_stop = int(args.entry_stop)

    converter = root2pickle(args.fname, entry_start = args.entry_start, entry_stop = args.entry_stop, pol = args.polarity, detRes = args.detRes, logistics = args.logistics, width = args.width, nofid = args.nofid)
    df = converter.df

    df.to_pickle(args.out)