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
        self.save(raw = raw, pol = pol)


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
        '''save data into df_epg, df_epgg for parent class epg'''
        self.readFile()

        # data frames and their keys to read Z part
        df_electronGen = pd.DataFrame()
        df_protonGen = pd.DataFrame()

        eleKeysGen = ["GenEpx", "GenEpy", "GenEpz"]
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
        eleKeysRec = ["Epx", "Epy", "Epz", "Eedep", "Evz", "Esector"]
        eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
        eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
        eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc2Hitx", "EDc2Hity", "EDc2Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
        eleKeysRec.extend(["Enphe"])
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "Psector", "Pchi2pid"]
        proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz", "PCvt12Hitx", "PCvt12Hity", "PCvt12Hitz"])
        proKeysRec.extend(["PDc2Hitx", "PDc2Hity", "PDc2Hitz", "PDc3Hitx", "PDc3Hity", "PDc3Hitz"])
        proKeysRec.extend(["PFtof1aTime", "PFtof1bTime", "PFtof2Time", "PCtofTime"])
        # proKeysRec.extend(["Pchi2pid", "Pchi2track", "PNDFtrack"])

        if detRes:
            eleKeysRec.extend(["Evx", "Evy"])
            # eleKeysRec.extend(["EDc1Hitx", "EDc1Hity", "EDc1Hitz", "EDc3Hitx", "EDc3Hity", "EDc3Hitz"])
            # eleKeysRec.extend(["Eedep1", "Eedep2", "Eedep3"])
            # eleKeysRec.extend(["EcalU1", "EcalV1", "EcalW1"])
            eleKeysRec.extend(["EcalU2", "EcalV2", "EcalW2"])
            eleKeysRec.extend(["EcalU3", "EcalV3", "EcalW3"])
            # eleKeysRec.extend(["Enphe"])
            eleKeysRec.extend(["EhtccX", "EhtccY", "EhtccZ"])
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
            df_electronRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        for key in proKeysRec:
            df_protonRec[key] = ak.to_dataframe(self.tree[key].array(library="ak", entry_start=entry_start, entry_stop=entry_stop))
        self.closeFile()

        #convert data type to standard double
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float})
        ele = [df_electronRec['Epx'], df_electronRec['Epy'], df_electronRec['Epz']]
        df_electronRec.loc[:, 'Ep'] = mag(ele)
        df_electronRec.loc[:,'ESamplFrac'] = df_electronRec.Eedep/ df_electronRec.Ep

        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')
        #apply fiducial cuts
        print(len(df_electronRec), len(df_protonRec))
        if nofid:
            df_electronRec.loc[:, "EFid"] = 1
            df_protonRec.loc[:, "PFid"] = 1
        else:
            df_electronRec = electronFiducial(df_electronRec, pol = pol, mc = False, fidlevel = fidlevel)
            df_protonRec = protonFiducial(df_protonRec, pol = pol)
            print(len(df_electronRec), len(df_protonRec))
            coincidence = reduce(np.intersect1d, (df_electronRec.event, df_protonRec.event))
            df_electronRec = df_electronRec.loc[df_electronRec.event.isin(coincidence), :]
            df_protonRec = df_protonRec.loc[df_protonRec.event.isin(coincidence), :]
            print(len(df_electronRec), len(df_protonRec))

            #electron pcal dead wires
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

        #prepare for proton energy loss corrections correction
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)
        df_protonRec.loc[:, "PpOrig"] = df_protonRec.Pp
        df_protonRec.loc[:, "PeOrig"] = df_protonRec.Pe
        df_protonRec.loc[:, "PthetaOrig"] = df_protonRec.Ptheta
        df_protonRec.loc[:, "PphiOrig"] = df_protonRec.Pphi

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

        if nocorr:
            print("no correction applied")
            df_protonRec = pd.concat([df_protonRecFD, df_protonRecCD, df_protonRecOthers])
        else:
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


            df_protonRecFD_1.loc[:, "PpEloss"] = CorrectedPp_FD_1
            df_protonRecFD_1.loc[:, "PthetaEloss"] = CorrectedPtheta_FD_1
            df_protonRecFD_1.loc[:, "PphiEloss"] = CorrectedPphi_FD_1
            df_protonRecFD_1.loc[:, "Pband"] = "lower"

            df_protonRecFD_2.loc[:, "PpEloss"] = CorrectedPp_FD_2
            df_protonRecFD_2.loc[:, "PthetaEloss"] = CorrectedPtheta_FD_2
            df_protonRecFD_2.loc[:, "Pband"] = "upper"

            df_protonRecCD.loc[:, "PpEloss"] = CorrectedPp_CD
            df_protonRecCD.loc[:, "PthetaEloss"] = CorrectedPtheta_CD
            df_protonRecCD.loc[:, "PphiEloss"] = CorrectedPphi_CD

            if noeloss or nopcorr:
                print("no energy loss correction applied.")
                pass
            else:
                print("energy loss correction applied for " + pol)
                df_protonRecCD.loc[:, "Pp"] = CorrectedPp_CD
                df_protonRecCD.loc[:, "Ptheta"] = CorrectedPtheta_CD
                df_protonRecCD.loc[:, "Pphi"] = CorrectedPphi_CD

                df_protonRecFD_1.loc[:, "Pp"] = CorrectedPp_FD_1
                df_protonRecFD_1.loc[:, "Ptheta"] = CorrectedPtheta_FD_1
                df_protonRecFD_1.loc[:, "Pphi"] = CorrectedPphi_FD_1

                df_protonRecFD_2.loc[:, "Pp"] = CorrectedPp_FD_2
                df_protonRecFD_2.loc[:, "Ptheta"] = CorrectedPtheta_FD_2
                df_protonRecFD_2.loc[:, "Pphi"] = CorrectedPphi_FD_2


            if nopcorr:
                print("no proton correction applied.")
                pass
            else:
                #correcting reconstruction bias after proton energy loss correction
                df_protonRecFD = pd.concat([df_protonRecFD_1, df_protonRecFD_2])
                print("applying the proton kinematic corrections for " + pol)

                df_protonRecCD.loc[:, "Pp"] = df_protonRecCD.Pp + 0.01
                df_protonRecCD.loc[:, "Ptheta"] = df_protonRecCD.Ptheta - 0.002129*df_protonRecCD.Ptheta**2 + 0.198*df_protonRecCD.Ptheta - 4.762 -0.2/(1+np.exp((df_protonRecCD.Pp-0.55)/(-0.05)))
                df_protonRecCD.loc[:, "Pphi"] = df_protonRecCD.Pphi
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

        # proton fiducial cuts
        if nofid:
            df_protonRec.loc[:, "PFid"] = 1
        else:
            df_protonRec.loc[:, "PFid"] = 0

            df_protonRec.loc[df_protonRec.Psector<7, "PFid"] = 1 #FD fid done by previous pipeline

            cut_CD = df_protonRec.Psector > 7
            if nocorr:
                cut_right = cut_CD
            elif fidlevel == 'mid':
                cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta)
            elif fidlevel == 'tight':
                cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta-5)
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
            
            # if pol == "inbending":
            #     pchi2CD_lb,   pchi2CD_ub   = -5.000, 6.345
            #     pchi2FD_S1_lb, pchi2FD_S1_ub = -3.296, 3.508
            #     pchi2FD_S2_lb, pchi2FD_S2_ub = -3.552, 4.000
            #     pchi2FD_S3_lb, pchi2FD_S3_ub = -3.446, 3.937
            #     pchi2FD_S4_lb, pchi2FD_S4_ub = -2.747, 3.190
            #     pchi2FD_S5_lb, pchi2FD_S5_ub = -2.851, 3.418
            #     pchi2FD_S6_lb, pchi2FD_S6_ub = -3.174, 3.514
            # elif pol == "outbending":
            #     pchi2CD_lb,   pchi2CD_ub   = -5.592,  6.785
            #     pchi2FD_S1_lb, pchi2FD_S1_ub = -3.905, 4.088
            #     pchi2FD_S2_lb, pchi2FD_S2_ub = -3.411, 3.939
            #     pchi2FD_S3_lb, pchi2FD_S3_ub = -4.042, 5.954
            #     pchi2FD_S4_lb, pchi2FD_S4_ub = -3.820, 5.065
            #     pchi2FD_S5_lb, pchi2FD_S5_ub = -3.384, 4.232
            #     pchi2FD_S6_lb, pchi2FD_S6_ub = -5.077, 5.100

            # df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.Pchi2pid<pchi2CD_lb)   | (df_protonRec.Pchi2pid>pchi2CD_ub)  ), "PFid"] = 0
            # df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.Pchi2pid<pchi2FD_S1_lb) | (df_protonRec.Pchi2pid>pchi2FD_S1_ub)), "PFid"] = 0
            # df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.Pchi2pid<pchi2FD_S2_lb) | (df_protonRec.Pchi2pid>pchi2FD_S2_ub)), "PFid"] = 0
            # df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.Pchi2pid<pchi2FD_S3_lb) | (df_protonRec.Pchi2pid>pchi2FD_S3_ub)), "PFid"] = 0
            # df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.Pchi2pid<pchi2FD_S4_lb) | (df_protonRec.Pchi2pid>pchi2FD_S4_ub)), "PFid"] = 0
            # df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.Pchi2pid<pchi2FD_S5_lb) | (df_protonRec.Pchi2pid>pchi2FD_S5_ub)), "PFid"] = 0
            # df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.Pchi2pid<pchi2FD_S6_lb) | (df_protonRec.Pchi2pid>pchi2FD_S6_ub)), "PFid"] = 0
        
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

        # binning kinematics
        df_ep.loc[:,'Q2'] = -((self.ebeam - df_ep['Ee'])**2 - mag2(VGS))
        df_ep.loc[:,'nu'] = (self.ebeam - df_ep['Ee'])
        df_ep.loc[:,'y'] = df_ep['nu']/self.ebeam
        df_ep.loc[:,'xB'] = df_ep['Q2'] / 2.0 / M / df_ep['nu']
        df_ep.loc[:,'t1'] = 2 * M * (df_ep['Pe'] - M)
        df_ep.loc[:,'t1Orig'] = 2 * M * (df_ep['PeOrig'] - M)
        df_ep.loc[:,'W'] = np.sqrt(np.maximum(0, (self.ebeam + M - df_ep['Ee'])**2 - mag2(VGS)))

        # trento angles
        df_ep.loc[:,'phi1'] = angle(v3l, v3h)
        df_ep.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_ep['phi1'], df_ep['phi1'])

        # exclusivity variables
        df_ep.loc[:,'MM2_ep'] = (-M - self.ebeam + df_ep["Ee"] + df_ep["Pe"])**2 - mag2(VmissG)

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