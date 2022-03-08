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
from scipy.stats import skewnorm

class root2pickle():
    #class to read root to make epg pairs, inherited from epg
    def __init__(self, fname, entry_start = None, entry_stop = None, pol = "inbending", detRes = False, logistics = False):
        self.fname = fname
        self.readEPGG(entry_start = entry_start, entry_stop = entry_stop, pol = pol, detRes = detRes, logistics = logistics)
        self.saveDVpi0vars()
        self.makeDVpi0P(pol = pol)
        self.save()

    def readFile(self):
        #read root using uproot
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        #close file for saving memory
        self.file = None
        self.tree = None

    def readEPGG(self, entry_start = None, entry_stop = None, pol = "inbending", detRes = False, logistics = False):
        #save data into df_epg, df_epgg for parent class epg
        self.readFile()

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()
        eleKeysRec = ["Epx", "Epy", "Epz", "Eedep", "Esector"]
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pstat", "Psector"]
        proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz"])
        proKeysRec.extend(["Pchi2pid", "Pchi2track", "PNDFtrack"])
        gamKeysRec = ["Gpx", "Gpy", "Gpz", "Gedep", "GcX", "GcY", "Gsector"]

        if detRes:
            eleKeysRec.extend(["Evx", "Evy", "Evz"])
            eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
            eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
            gamKeysRec.extend(["Gedep1", "Gedep2", "Gedep3"])
            proKeysRec.extend(["Pvz"])
            proKeysRec.extend(["PCvt1Hitx", "PCvt1Hity", "PCvt1Hitz", "PCvt3Hitx", "PCvt3Hity", "PCvt3Hitz", "PCvt5Hitx", "PCvt5Hity", "PCvt5Hitz", "PCvt7Hitx", "PCvt7Hity", "PCvt7Hitz", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"])
            proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PDc3Hitx", "PDc3Hity", "PDc3Hitz"])
            eleKeysRec.extend(["startTime"])
            proKeysRec.extend(["PFtof1aTime", "PFtof1bTime", "PFtof2Time", "PCtofTime"])
            # proKeysRec.extend(["Pchi2pid", "Pchi2track", "PNDFtrack"])
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
        # df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float, "Evx": float, "Evy": float, "Evz": float})
        # df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float, "Pvz": float})
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        # df_protonRec = df_protonRec.astype({"PDc1Hitx": float, "PDc1Hity": float, "PDc1Hitz": float})
        # df_protonRec = df_protonRec.astype({"PDc3Hitx": float, "PDc3Hity": float, "PDc3Hitz": float})
        # df_protonRec = df_protonRec.astype({"PCvt1Hitx": float, "PCvt1Hity": float, "PCvt1Hitz": float})
        # df_protonRec = df_protonRec.astype({"PCvt3Hitx": float, "PCvt3Hity": float, "PCvt3Hitz": float})
        # df_protonRec = df_protonRec.astype({"PCvt5Hitx": float, "PCvt5Hity": float, "PCvt5Hitz": float})
        # df_protonRec = df_protonRec.astype({"PCvt7Hitx": float, "PCvt7Hity": float, "PCvt7Hitz": float})
        # df_protonRec = df_protonRec.astype({"PCvt12Hitx": float, "PCvt12Hity": float, "PCvt12Hitz": float})
        df_gammaRec = df_gammaRec.astype({"Gpx": float, "Gpy": float, "Gpz": float, "Gedep": float, "GcX": float, "GcY": float})

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

        df_protonRecFD = df_protonRec.loc[df_protonRec.Psector<7, :]
        df_protonRecCD = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), :]
        df_protonRecOthers = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta>=75), :]

        correction = False

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

        #inbending
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
            #proton correction
            print("correction applied for " + pol)

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

        ele = [df_electronRec['Epx'], df_electronRec['Epy'], df_electronRec['Epz']]
        df_electronRec.loc[:, 'Ep'] = mag(ele)
        df_electronRec.loc[:,'ESamplFrac'] = df_electronRec.Eedep/ df_electronRec.Ep
        gam = [df_gammaRec['Gpx'], df_gammaRec['Gpy'], df_gammaRec['Gpz']]
        df_gammaRec.loc[:, 'Gp'] = mag(gam)
        df_gammaRec.loc[:,'GSamplFrac'] = df_gammaRec.Gedep/ df_gammaRec.Gp

        df_gg = pd.merge(df_gammaRec, df_gammaRec,
                         how='outer', on='event', suffixes=("", "2"))
        df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
        df_gg = df_gg.drop(['GIndex', 'GIndex2'], axis = 1)

        if correction:
            #photon correction of df_gg
            gam = [df_gg['Gpx'], df_gg['Gpy'], df_gg['Gpz']]
            df_gg.loc[:, 'Gp'] = mag(gam)
            df_gg.loc[:, 'Gtheta'] = getTheta(gam)
            df_gg.loc[:, 'Gphi'] = getPhi(gam)
            #photon correction of df_gamma
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

        if detRes:
            df_protonRec.loc[:, "PAngleDiff"] = df_protonRec.loc[:, "PDc3theta"] - df_protonRec.loc[:, "PDc1theta"]
            df_gg = df_gg.loc[:, ~df_gg.columns.duplicated()]
            df_gg.loc[:, "Gedep2_tot"] = df_gg.Gedep12 + df_gg.Gedep22 + df_gg.Gedep32
            df_electronRec = df_electronRec.drop(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz", "EDc1theta", "EDc3theta"], axis = 1)
            df_protonRec = df_protonRec.drop(["PCvt1Hitx", "PCvt1Hity", "PCvt1Hitz", "PCvt3Hitx", "PCvt3Hity", "PCvt3Hitz", "PCvt5Hitx", "PCvt5Hity", "PCvt5Hitz", "PCvt7Hitx", "PCvt7Hity", "PCvt7Hitz", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"], axis = 1)
            df_protonRec = df_protonRec.drop(["PDc3Hitx", "PDc3Hity", "PDc3Hitz", "PDc3theta"], axis = 1)
        else:
            df_protonRec = df_protonRec.drop(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PDc1theta"], axis = 1)

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
        df_epgg.loc[:,'y'] = df_epgg['nu']/ebeam
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
        df_epgg.loc[:,'openingAngle'] = angle(gam, gam2)
        
        df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)

        # encode other binning
        for Q2bin in range(len(Q2bin_i)):
            #square Q2 binning
            df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]), "Q2bin"] = Q2bin
            #adaptive xB binning
            for xBbin in range(len(xBbin_i[Q2bin])):
                if Q2bin < len(Q2bin_i) -1:
                    if xBbin == 0:
                        df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin #0
                        df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin #0
                    elif xBbin < len(xBbin_i[Q2bin])-1:
                        df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin
                        df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin
                    else:
                        df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
                        df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin
                else:
                    df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
                    df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin #0

                Q2xBbin = Q2xBbin + 1
        for tbin in range(len(tbin_i)):
            #square t binning
            df_epgg.loc[(df_epgg.t1>=tbin_i[tbin]) & (df_epgg.t1<tbin_f[tbin]), "tbin"] = tbin
            # df_epgg.loc[(df_epgg.t2>=tbin_i[tbin]) & (df_epgg.t2<tbin_f[tbin]), "tbin2"] = tbin
        for phibin in range(len(phibin_i)):
            #square phi binning
            df_epgg.loc[(df_epgg.phi1>=phibin_i[phibin]) & (df_epgg.phi1<phibin_f[phibin]), "phibin"] = phibin
            # df_epgg.loc[(df_epgg.phi2>=phibin_i[phibin]) & (df_epgg.phi2<phibin_f[phibin]), "phibin2"] = phibin

        df_epgg.loc[(df_epgg.Q2xBbin>0)&(df_epgg.tbin>0), "Q2xBtbin"] = len(tbin_i) * df_epgg.loc[(df_epgg.Q2xBbin>0)&(df_epgg.tbin>0), "Q2xBbin"] + df_epgg.loc[(df_epgg.Q2xBbin>0)&(df_epgg.tbin>0), "tbin"]
        # df_epgg.loc[(df_epgg.Q2bin>0)&(df_epgg.xBbin>0)&(df_epgg.tbin2>0), "Q2xBtbin2"] = df_epgg.Q2bin.astype(str) + df_epgg.xBbin.astype(str) + df_epgg.tbin2.astype(str)
        df_epgg.loc[(df_epgg.Q2xBbin>0)&(df_epgg.tbin>0), "Q2xBtphibin"] = len(phibin_i) * df_epgg.loc[(df_epgg.Q2xBbin>0)&(df_epgg.tbin>0), "Q2xBtbin"] + df_epgg.loc[(df_epgg.Q2xBbin>0)&(df_epgg.tbin>0), "phibin"]

        self.df_epgg = df_epgg

    def makeDVpi0P(self, pol = "inbending"):
        #make dvpi0 pairs
        df_dvpi0p = self.df_epgg

        #common cuts
        cut_xBupper = df_dvpi0p.loc[:, "xB"] < 1  # xB
        cut_xBlower = df_dvpi0p.loc[:, "xB"] > 0  # xB
        cut_Q2 = df_dvpi0p.loc[:, "Q2"] > 1  # Q2
        cut_W = df_dvpi0p.loc[:, "W"] > 2  # W
        cut_Ee = df_dvpi0p["Ee"] > 2  # Ee
        cut_Ge2 = df_dvpi0p["Ge2"] > 0.6  # Ge cut. Ge>3 for DVCS module.
        cut_Esector = (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector"]) & (df_dvpi0p["Esector"]!=df_dvpi0p["Gsector2"])
        cut_Psector = ~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector"]) ) & ~( ((df_dvpi0p["Pstat"]//10)%10>0) & (df_dvpi0p["Psector"]==df_dvpi0p["Gsector2"]) )
        cut_Ppmax = df_dvpi0p.Pp < 1.6  # Pp
        cut_Pthetamin = df_dvpi0p.Ptheta > 0  # Ptheta
        # cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])
        cut_common = cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge2 & cut_Esector & cut_Psector & cut_Ppmax & cut_Pthetamin

        df_dvpi0p = df_dvpi0p[cut_common]

        # proton reconstruction quality
        # cut_FD_proton = (df_epgg.loc[:, "Psector"]<7) & (df_epgg.loc[:, "Ptheta"]<35)
        # cut_CD_proton = (df_epgg.loc[:, "Psector"]>7) & (df_epgg.loc[:, "Ptheta"]>45) & (df_epgg.loc[:, "Ptheta"]<65)
        # cut_proton = (cut_FD_proton)|(cut_CD_proton)
        cut_proton = 1

        df_dvpi0p.loc[:, "config"] = 0

        if pol == "inbending":
            #CDFT
            cut_Pp1_CDFT = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CDFT = df_dvpi0p.Psector>7
            cut_Ptheta1_CDFT = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CDFT = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CDFT = df_dvpi0p.Gsector>7
            cut_GFid_CDFT = df_dvpi0p.GFid==1
            cut_GFid2_CDFT = df_dvpi0p.GFid2==1
            cut_mpi01_CDFT = df_dvpi0p["Mpi0"] < 0.149#0.157  # mpi0
            cut_mpi02_CDFT = df_dvpi0p["Mpi0"] > 0.126#0.118  # mpi0
            cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] < 0.610#0.914  # mmep
            cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] > -0.384#-0.715  # mmep
            cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] < 1.641#2.155  # mmegg
            cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] > 0.0974#-0.417  # mmegg
            cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] < 0.481#0.799  # meepgg
            cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] > -0.474#-0.792  # meepgg
            cut_mpt_CDFT = df_dvpi0p["MPt"] < 0.1272#0.189  # mpt
            cut_recon_CDFT = df_dvpi0p["reconPi"] < 0.955  # recon gam angle
            cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] < 9.259#15.431  # coplanarity angle
            cut_mmepgg1_CDFT = df_dvpi0p["MM2_epgg"] < 0.02564#0.0440  # mmepgg
            cut_mmepgg2_CDFT = df_dvpi0p["MM2_epgg"] > -0.02944#-0.0478  # mmepgg

            cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & 
                        cut_GFid_CDFT & cut_GFid2_CDFT &
                        cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
                        cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
                        cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


            #CD
            cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CD = df_dvpi0p.Psector>7
            cut_Ptheta1_CD = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CD = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CD = df_dvpi0p.Gsector<7
            cut_Gsector2_CD = df_dvpi0p.Gsector2<7
            cut_GFid_CD = df_dvpi0p.GFid==1
            cut_GFid2_CD = df_dvpi0p.GFid2==1
            cut_mpi01_CD = df_dvpi0p["Mpi0"] < 0.162  # mpi0
            cut_mpi02_CD = df_dvpi0p["Mpi0"] > 0.107  # mpi0
            cut_mmep1_CD = df_dvpi0p["MM2_ep"] < 0.354  # mmep
            cut_mmep2_CD = df_dvpi0p["MM2_ep"] > -0.283  # mmep
            cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < 1.922  # mmegg
            cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > 0.007  # mmegg
            cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < 0.822  # meepgg
            cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > -0.677  # meepgg
            cut_mpt_CD = df_dvpi0p["MPt"] < 0.176  # mpt
            cut_recon_CD = df_dvpi0p["reconPi"] < 1.476  # recon gam angle
            cut_coplanarity_CD = df_dvpi0p["coplanarity"] < 10.203  # coplanarity angle
            cut_mmepgg1_CD = df_dvpi0p["MM2_epgg"] < 0.0208  # mmepgg
            cut_mmepgg2_CD = df_dvpi0p["MM2_epgg"] > -0.0250  # mmepgg

            cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_Gsector2_CD & 
                        cut_GFid_CD & cut_GFid2_CD &
                        cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
                        cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
                        cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

            #FD
            cut_Pp1_FD = df_dvpi0p.Pp > 0.42  # Pp
            cut_Psector_FD = df_dvpi0p.Psector<7
            cut_Ptheta1_FD = df_dvpi0p.Ptheta<FD_Ptheta_ub
            cut_Ptheta2_FD = df_dvpi0p.Ptheta>FD_Ptheta_lb
            cut_Gsector_FD = df_dvpi0p.Gsector<7
            cut_Gsector2_FD = df_dvpi0p.Gsector2<7
            cut_GFid_FD = df_dvpi0p.GFid==1
            cut_GFid2_FD = df_dvpi0p.GFid2==1
            cut_mpi01_FD = df_dvpi0p["Mpi0"] < 0.178  # mpi0
            cut_mpi02_FD = df_dvpi0p["Mpi0"] > 0.0910  # mpi0
            cut_mmep1_FD = df_dvpi0p["MM2_ep"] < 0.335  # mmep
            cut_mmep2_FD = df_dvpi0p["MM2_ep"] > -0.271  # mmep
            cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < 1.762  # mmegg
            cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > 0.117  # mmegg
            cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < 0.816 # meepgg
            cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > -0.685  # meepgg
            cut_mpt_FD = df_dvpi0p["MPt"] < 0.180  # mpt
            cut_recon_FD = df_dvpi0p["reconPi"] < 1.363  # recon gam angle
            cut_coplanarity_FD = df_dvpi0p["coplanarity"] < 9.190  # coplanarity angle
            cut_mmepgg1_FD = df_dvpi0p["MM2_epgg"] < 0.0189  # mmepgg
            cut_mmepgg2_FD = df_dvpi0p["MM2_epgg"] > -0.0224  # mmepgg

            cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_Gsector2_FD &
                        cut_GFid_FD & cut_GFid2_FD &
                        cut_mmep1_FD & cut_mmep2_FD & cut_mpi01_FD & cut_mpi02_FD & 
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
            cut_GFid2_CDFT = df_dvpi0p.GFid2==1
            cut_mpi01_CDFT = df_dvpi0p["Mpi0"] < 0.151#0.160  # mpi0
            cut_mpi02_CDFT = df_dvpi0p["Mpi0"] > 0.124#0.115  # mpi0
            cut_mmep1_CDFT = df_dvpi0p["MM2_ep"] < 0.575#0.892  # mmep
            cut_mmep2_CDFT = df_dvpi0p["MM2_ep"] > -0.378#-0.694  # mmep
            cut_mmegg1_CDFT = df_dvpi0p["MM2_egg"] < 1.665#2.184  # mmegg
            cut_mmegg2_CDFT = df_dvpi0p["MM2_egg"] > 0.107#-0.412  # mmegg
            cut_meepgg1_CDFT = df_dvpi0p["ME_epgg"] < 0.514#0.844  # meepgg
            cut_meepgg2_CDFT = df_dvpi0p["ME_epgg"] > -0.476#-0.806  # meepgg
            cut_mpt_CDFT = df_dvpi0p["MPt"] < 0.146#0.210  # mpt
            cut_recon_CDFT = df_dvpi0p["reconPi"] < 1.114#1.630  # recon gam angle
            cut_coplanarity_CDFT = df_dvpi0p["coplanarity"] < 10.69#17.817  # coplanarity angle
            cut_mmepgg1_CDFT = df_dvpi0p["MM2_epgg"] < 0.0324#0.0549  # mmepgg
            cut_mmepgg2_CDFT = df_dvpi0p["MM2_epgg"] > -0.035#-0.0575  # mmepgg

            cut_CDFT = (cut_Pp1_CDFT & cut_Psector_CDFT & cut_Ptheta1_CDFT & cut_Ptheta2_CDFT & cut_Gsector_CDFT & 
                        cut_GFid_CDFT & cut_GFid2_CDFT &
                        cut_mmep1_CDFT & cut_mmep2_CDFT & cut_mpi01_CDFT & cut_mpi02_CDFT & 
                        cut_mmegg1_CDFT & cut_mmegg2_CDFT & cut_meepgg1_CDFT & cut_meepgg2_CDFT &
                        cut_mpt_CDFT & cut_recon_CDFT & cut_coplanarity_CDFT & cut_mmepgg1_CDFT & cut_mmepgg2_CDFT)


            #CD
            cut_Pp1_CD = df_dvpi0p.Pp > 0.3  # Pp
            cut_Psector_CD = df_dvpi0p.Psector>7
            cut_Ptheta1_CD = df_dvpi0p.Ptheta<CD_Ptheta_ub
            cut_Ptheta2_CD = df_dvpi0p.Ptheta>CD_Ptheta_lb
            cut_Gsector_CD = df_dvpi0p.Gsector<7
            cut_Gsector2_CD = df_dvpi0p.Gsector2<7
            cut_GFid_CD = df_dvpi0p.GFid==1
            cut_GFid2_CD = df_dvpi0p.GFid2==1
            cut_mpi01_CD = df_dvpi0p["Mpi0"] < 0.163  # mpi0
            cut_mpi02_CD = df_dvpi0p["Mpi0"] > 0.106  # mpi0
            cut_mmep1_CD = df_dvpi0p["MM2_ep"] < 0.294  # mmep
            cut_mmep2_CD = df_dvpi0p["MM2_ep"] > -0.218  # mmep
            cut_mmegg1_CD = df_dvpi0p["MM2_egg"] < 1.876  # mmegg
            cut_mmegg2_CD = df_dvpi0p["MM2_egg"] > -0.0142  # mmegg
            cut_meepgg1_CD = df_dvpi0p["ME_epgg"] < 0.700  # meepgg
            cut_meepgg2_CD = df_dvpi0p["ME_epgg"] > -0.597  # meepgg
            cut_mpt_CD = df_dvpi0p["MPt"] < 0.194  # mpt
            cut_recon_CD = df_dvpi0p["reconPi"] < 1.761  # recon gam angle
            cut_coplanarity_CD = df_dvpi0p["coplanarity"] < 9.530  # coplanarity angle
            cut_mmepgg1_CD = df_dvpi0p["MM2_epgg"] < 0.0182  # mmepgg
            cut_mmepgg2_CD = df_dvpi0p["MM2_epgg"] > -0.0219  # mmepgg

            cut_CD = (cut_Pp1_CD & cut_Psector_CD & cut_Ptheta1_CD & cut_Ptheta2_CD & cut_Gsector_CD & cut_Gsector2_CD & 
                        cut_GFid_CD & cut_GFid2_CD &
                        cut_mmep1_CD & cut_mmep2_CD & cut_mpi01_CD & cut_mpi02_CD & 
                        cut_mmegg1_CD & cut_mmegg2_CD & cut_meepgg1_CD & cut_meepgg2_CD &
                        cut_mpt_CD & cut_recon_CD & cut_coplanarity_CD & cut_mmepgg1_CD & cut_mmepgg2_CD)

            #FD
            cut_Pp1_FD = df_dvpi0p.Pp > 0.5  # Pp
            cut_Psector_FD = df_dvpi0p.Psector<7
            cut_Ptheta1_FD = df_dvpi0p.Ptheta<FD_Ptheta_ub
            cut_Ptheta2_FD = df_dvpi0p.Ptheta>FD_Ptheta_lb
            cut_Gsector_FD = df_dvpi0p.Gsector<7
            cut_Gsector2_FD = df_dvpi0p.Gsector2<7
            cut_GFid_FD = df_dvpi0p.GFid==1
            cut_GFid2_FD = df_dvpi0p.GFid2==1
            cut_mpi01_FD = df_dvpi0p["Mpi0"] < 0.164  # mpi0
            cut_mpi02_FD = df_dvpi0p["Mpi0"] > 0.105  # mpi0
            cut_mmep1_FD = df_dvpi0p["MM2_ep"] < 0.323  # mmep
            cut_mmep2_FD = df_dvpi0p["MM2_ep"] > -0.256  # mmep
            cut_mmegg1_FD = df_dvpi0p["MM2_egg"] < 1.828  # mmegg
            cut_mmegg2_FD = df_dvpi0p["MM2_egg"] > 0.0491  # mmegg
            cut_meepgg1_FD = df_dvpi0p["ME_epgg"] < 0.754  # meepgg
            cut_meepgg2_FD = df_dvpi0p["ME_epgg"] > -0.583  # meepgg
            cut_mpt_FD = df_dvpi0p["MPt"] < 0.177  # mpt
            cut_recon_FD = df_dvpi0p["reconPi"] < 1.940  # recon gam angle
            cut_coplanarity_FD = df_dvpi0p["coplanarity"] < 7.498  # coplanarity angle
            cut_mmepgg1_FD = df_dvpi0p["MM2_epgg"] < 0.0195  # mmepgg
            cut_mmepgg2_FD = df_dvpi0p["MM2_epgg"] > -0.0240  # mmepgg

            cut_FD = (cut_Pp1_FD & cut_Psector_FD & cut_Ptheta1_FD & cut_Ptheta2_FD & cut_Gsector_FD & cut_Gsector2_FD &
                        cut_GFid_FD & cut_GFid2_FD &
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
        df_dvpi0p = df_dvpi0p.sort_values(by=['closeness', 'Psector', 'Gsector'], ascending = [True, True, True])
        df_dvpi0p = df_dvpi0p.loc[~df_dvpi0p.event.duplicated(), :]
        df_dvpi0p = df_dvpi0p.sort_values(by='event')        
        self.df_dvpi0p = df_dvpi0p #done with saving x

    def save(self):
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
    parser.add_argument("-l","--logistics", help="include logistics", action = "store_true")
    
    args = parser.parse_args()

    if args.entry_start:
        args.entry_start = int(args.entry_start)
    if args.entry_stop:
        args.entry_stop = int(args.entry_stop)

    converter = root2pickle(args.fname, entry_start = args.entry_start, entry_stop = args.entry_stop, pol = args.polarity, detRes = args.detRes, logistics = args.logistics)
    df = converter.df

    df.to_pickle(args.out)