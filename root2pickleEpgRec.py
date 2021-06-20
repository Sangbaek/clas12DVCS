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
    def __init__(self, fname, entry_stop = None, gen = "dvcs", pol = "inbending"):
        self.fname = fname

        self.readEPGG(entry_stop, gen = gen, pol = pol)
        self.saveDVCSvars()
        self.makeDVCS()
        self.saveDVpi0vars()
        self.makeDVpi0()
        self.pi02gSubtraction()
        self.saveRaw()

    def readFile(self):
        #read root using uproot
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        #close file for saving memory
        self.file = None
        self.tree = None

    def readEPGG(self, entry_stop = None, gen = "dvcs", pol = "inbending"):
        #save data into df_epg, df_epgg for parent class epg
        self.readFile()

        # data frames and their keys to read Z part
        df_electronGen = pd.DataFrame()
        df_protonGen = pd.DataFrame()
        df_gammaGen = pd.DataFrame()
        # eleKeysGen = ["GenEpx", "GenEpy", "GenEpz", "GenEvx", "GenEvy", "GenEvz"]
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
        # df_electronGen = df_electronGen.astype({"GenEpx": float, "GenEpy": float, "GenEpz": float, "GenEvx": float, "GenEvy": float, "GenEvz": float})
        df_electronGen = df_electronGen.astype({"GenEpx": float, "GenEpy": float, "GenEpz": float})
        df_protonGen = df_protonGen.astype({"GenPpx": float, "GenPpy": float, "GenPpz": float})
        df_gammaGen = df_gammaGen.astype({"GenGpx": float, "GenGpy": float, "GenGpz": float})

        #set up a dummy index for merging
        df_electronGen.loc[:,'event'] = df_electronGen.index
        df_protonGen.loc[:,'event'] = df_protonGen.index
        df_gammaGen.loc[:,'event'] = df_gammaGen.index.get_level_values('entry')

        #sort columns for readability
        # df_electronGen = df_electronGen.loc[:, ["event", "GenEpx", "GenEpy", "GenEpz", "GenEvz"]]
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

        if gen == "dvcs":
            df_gammaGen = df_gammaGen[df_gammaGen.index.get_level_values('subentry')==0]
            gamGen = [df_gammaGen["GenGpx"], df_gammaGen["GenGpy"], df_gammaGen["GenGpz"]]
            df_gammaGen.loc[:, 'GenGp'] = mag(gamGen)
            df_gammaGen.loc[:, 'GenGtheta'] = getTheta(gamGen)
            df_gammaGen.loc[:, 'GenGphi'] = getPhi(gamGen)

            df_z = pd.merge(df_z, df_gammaGen, how='inner', on='event')
            self.df_z = df_z    #done with saving z

        elif gen == "pi0":
            #two g's to one gg.
            gamGen = [df_gammaGen["GenGpx"], df_gammaGen["GenGpy"], df_gammaGen["GenGpz"]]
            df_gammaGen.loc[:, 'GenGp'] = mag(gamGen)

            gam1 = df_gammaGen[df_gammaGen.index.get_level_values('subentry')==0]
            gam1 = gam1.reset_index(drop=True)
            gam2 = df_gammaGen[df_gammaGen.index.get_level_values('subentry')==1]
            gam2 = gam2.reset_index(drop=True)

            gam1.loc[:,"GenGp2"] = gam2.loc[:,"GenGp"]
            gam1.loc[:,"GenGpx2"] = gam2.loc[:,"GenGpx"]
            gam1.loc[:,"GenGpy2"] = gam2.loc[:,"GenGpy"]
            gam1.loc[:,"GenGpz2"] = gam2.loc[:,"GenGpz"]
            df_gammaGen = gam1

            #sort GenG indices so that GenGp > GenGp2. This is because Gp > Gp2 at reconstruction level.
            df_gammaGencopy = copy(df_gammaGen)
            df_gammaGencopy.loc[:, "GenGp"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGp"], df_gammaGen.loc[:, "GenGp2"])
            df_gammaGencopy.loc[:, "GenGpx"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGpx"], df_gammaGen.loc[:, "GenGpx2"])
            df_gammaGencopy.loc[:, "GenGpy"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGpy"], df_gammaGen.loc[:, "GenGpy2"])
            df_gammaGencopy.loc[:, "GenGpz"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGpz"], df_gammaGen.loc[:, "GenGpz2"])
            df_gammaGencopy.loc[:, "GenGp2"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGp2"], df_gammaGen.loc[:, "GenGp"])
            df_gammaGencopy.loc[:, "GenGpx2"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGpx2"], df_gammaGen.loc[:, "GenGpx"])
            df_gammaGencopy.loc[:, "GenGpy2"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGpy2"], df_gammaGen.loc[:, "GenGpy"])
            df_gammaGencopy.loc[:, "GenGpz2"] = np.where(df_gammaGen["GenGp"]>df_gammaGen["GenGp2"], df_gammaGen.loc[:, "GenGpz2"], df_gammaGen.loc[:, "GenGpz"])
            df_gammaGen = df_gammaGencopy

            gamGen = [df_gammaGen["GenGpx"], df_gammaGen["GenGpy"], df_gammaGen["GenGpz"]]
            # df_gammaGen.loc[:, 'GenGp'] = mag(gamGen)
            df_gammaGen.loc[:, 'GenGtheta'] = getTheta(gamGen)
            df_gammaGen.loc[:, 'GenGphi'] = getPhi(gamGen)

            gamGen2 = [df_gammaGen["GenGpx2"], df_gammaGen["GenGpy2"], df_gammaGen["GenGpz2"]]
            debug = df_gammaGen.loc[:, 'GenGp2'] == mag(gamGen2)
            df_gammaGen.loc[:, 'GenGtheta2'] = getTheta(gamGen2)
            df_gammaGen.loc[:, 'GenGphi2'] = getPhi(gamGen2)

            df_z = pd.merge(df_z, df_gammaGen, how='inner', on='event')
            self.df_z = df_z    #done with saving z

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        df_gammaRec = pd.DataFrame()
        # eleKeysRec = ["Epx", "Epy", "Epz", "Evx", "Evy", "Evz", "Esector"]
        # proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Psector"]
        eleKeysRec = ["Epx", "Epy", "Epz", "Evz", "Esector"]
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
        # df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float, "Evx": float, "Evy": float, "Evz": float})
        # df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float, "Pvz": float})
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float, "Evz": float})
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

        df_protonRecFD = df_protonRec.loc[df_protonRec.Psector<7, :]
        df_protonRecCD = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), :]

        #inbending
        # if pol == "inbending":
        #     #FD part
        #     const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [-0.0123049 + 0.00028887*df_protonRecFD.Ptheta, -0.138227479 + 8.07557430*0.001*df_protonRecFD.Ptheta -1.34807927*0.0001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -0.0275235])
        #     coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [0.01528006 - 0.00024079*df_protonRecFD.Ptheta, 5.65817597*0.01 -2.36903348*0.001*df_protonRecFD.Ptheta + 4.93780046*0.00001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 0.03998975])    

        #     CorrectedPp_FD = const_FD + coeff_FD/df_protonRecFD.loc[:, "Pp"] + df_protonRecFD.loc[:, "Pp"]

        #     const_FD = np.select([df_protonRecFD.Ptheta<19.5, (df_protonRecFD.Ptheta>=19.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<39), (df_protonRecFD.Ptheta>=39) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [2.63643690*0.01, 0.50047232 -0.03834672 *df_protonRecFD.Ptheta + 0.00071967*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 6.91308654 - 0.439839300*df_protonRecFD.Ptheta +6.83075548*0.001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 1.59424606, 1.47198581*10])
        #     coeff_FD = np.select([df_protonRecFD.Ptheta<19.5, (df_protonRecFD.Ptheta>=19.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<39), (df_protonRecFD.Ptheta>=39) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [-1.46440415, 74.99891704  -6.1576777*df_protonRecFD.Ptheta + 0.11469137*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 682.909471 - 43.9551177 * df_protonRecFD.Ptheta + 0.682383790 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, -8.19627119, -23.55701865])    
        #     coeff2_FD = np.select([df_protonRecFD.Ptheta<19.5, (df_protonRecFD.Ptheta>=19.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<39), (df_protonRecFD.Ptheta>=39) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [-3.47690993, 47.71351973 -4.34918241*df_protonRecFD.Ptheta + 0.08841191*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 100.33995753 - 6.96600416*df_protonRecFD.Ptheta + 0.11223046*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -1.25261927, -0.40113733])    

        #     CorrectedPtheta_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Ptheta"]

        #     const_FD = np.select([df_protonRecFD.Ptheta<16.5, (df_protonRecFD.Ptheta>=16.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [-0.190662844, -0.20725736 -0.00675627 *df_protonRecFD.Ptheta + 0.0007863*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 12.1881698 - 0.78906294*df_protonRecFD.Ptheta +0.01297898*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -4.59743066*10])
        #     coeff_FD = np.select([df_protonRecFD.Ptheta<16.5, (df_protonRecFD.Ptheta>=16.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [6.48745941, 142.96379788  -16.66339055*df_protonRecFD.Ptheta + 0.51311212*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 2.1853046 + 5.78521226 * df_protonRecFD.Ptheta - 0.09727796 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, 7.46969457*10])    
        #     coeff2_FD = np.select([df_protonRecFD.Ptheta<16.5, (df_protonRecFD.Ptheta>=16.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
        #                       [-3.14646608, 17.39529095 -1.78403359*df_protonRecFD.Ptheta + 0.0335692*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -1.03655317*10 + 0.161333213*df_protonRecFD.Ptheta -1.29625675*0.001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -4.41246899*0.1])    

        #     CorrectedPphi_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Pphi"]

        #     #CD part
        #     const_CD = 1.93686914 - 0.116288824*df_protonRecCD.Ptheta + 0.00223685833*df_protonRecCD.Ptheta**2 - 1.40771969 * 10**(-5)*df_protonRecCD.Ptheta**3
        #     coeff_CD = -0.738047800 + 0.0443343685*df_protonRecCD.Ptheta - 8.50985972*10**(-4)*df_protonRecCD.Ptheta*df_protonRecCD.Ptheta + 5.36810280 * 10**(-6) * df_protonRecCD.Ptheta**3

        #     CorrectedPp_CD = const_CD + coeff_CD/df_protonRecCD.loc[:, "Pp"] + df_protonRecCD.loc[:, "Pp"]

        #     const_CD = -1.09849291*100 + 8.86664014 * df_protonRecCD.Ptheta - 0.26643881 * df_protonRecCD.Ptheta**2 + 3.53814210 * 10**(-3) * df_protonRecCD.Ptheta**3 - 1.75297107 * 10**(-5) * df_protonRecCD.Ptheta**4
        #     coeff_CD = 9.52034523*100 -5.74808292 * 10 * df_protonRecCD.Ptheta + 1.15386949 * df_protonRecCD.Ptheta**2 - 7.57970373 * 0.001 * df_protonRecCD.Ptheta**3
        #     coeff2_CD = -2.00387313*100 + 1.18979079 * 10 * df_protonRecCD.Ptheta - 2.37730217*0.1 * df_protonRecCD.Ptheta**2 + 1.55153003*0.001*df_protonRecCD.Ptheta**3

        #     CorrectedPtheta_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Ptheta"]

        #     const_CD = 4.94546178 -3.26662886*0.1 * df_protonRecCD.Ptheta +  7.39069603 * 0.001 * df_protonRecCD.Ptheta**2 -6.83599356*10**(-5) * df_protonRecCD.Ptheta**3 + 2.12303103*10**(-7) * df_protonRecCD.Ptheta**4
        #     coeff_CD = 1.72181613*10**(5) -1.36827111*10**(4) * df_protonRecCD.Ptheta + 4.00923146*10**(2) * df_protonRecCD.Ptheta**2 - 5.12792347 * df_protonRecCD.Ptheta**3 + 2.41793167*10**(-2) * df_protonRecCD.Ptheta**4
        #     coeff2_CD =  1.20477219*10**(2) -5.86630228 * df_protonRecCD.Ptheta + 7.44007875*10**(-2) * df_protonRecCD.Ptheta**2 -2.42652473*10**(-4) * df_protonRecCD.Ptheta**3
        #     CorrectedPphi_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Pphi"]

        # elif pol == "outbending":
        #     #FD part
        #     const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27)],
        #                       [0.02067157-0.0009827*df_protonRecFD.Ptheta, -0.11216694 + 0.0069912*df_protonRecFD.Ptheta - 0.00011733 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta])
        #     coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27)],
        #                       [-0.03334437+0.00177781*df_protonRecFD.Ptheta, 0.0402797945 - 0.00197220505*df_protonRecFD.Ptheta + 4.50918200*10**(-5) * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta])

        #     CorrectedPp_FD = const_FD + coeff_FD/df_protonRecFD.loc[:, "Pp"] + df_protonRecFD.loc[:, "Pp"]

        #     const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
        #                       [0, -1.79343987 +0.105559096 *df_protonRecFD.Ptheta + -0.00157174358*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -0.123044632])
        #     coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
        #                       [0, -27.4344526 + 1.61037587* df_protonRecFD.Ptheta - 0.0242300381* df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, -7.52117236])    
        #     coeff2_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
        #                       [0, -45.2983842 +2.51745350*df_protonRecFD.Ptheta - 0.0365942178*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -3.52825441])    

        #     CorrectedPtheta_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Ptheta"]

        #     const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
        #                       [0, 5.37967179 -0.324630795 *df_protonRecFD.Ptheta + 0.00476947696*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -0.0224918574])
        #     coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
        #                       [0, 7.25038499*1000 + -413.586911* df_protonRecFD.Ptheta + 5.91815405 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, 55.6319490])    
        #     coeff2_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
        #                       [0, -124.626261 + 6.77668728*df_protonRecFD.Ptheta - 0.0960045129*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -5.12646023])    

        #     CorrectedPphi_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Pphi"]

        #     #CD part
        #     const_CD = 1.92657376 - 0.113836734*df_protonRecCD.Ptheta + 0.00215038526*df_protonRecCD.Ptheta**2 - 1.32525053 * 10**(-5)*df_protonRecCD.Ptheta**3
        #     coeff_CD = -0.755650043 + 0.0445538936*df_protonRecCD.Ptheta - 8.38241864*10**(-4)*df_protonRecCD.Ptheta*df_protonRecCD.Ptheta + 5.16887255 * 10**(-6) * df_protonRecCD.Ptheta**3

        #     CorrectedPp_CD = const_CD + coeff_CD/df_protonRecCD.loc[:, "Pp"] + df_protonRecCD.loc[:, "Pp"]

        #     const_CD = -5.79024055*10 + 4.67197531 * df_protonRecCD.Ptheta - 0.140156897 * df_protonRecCD.Ptheta**2 + 1.85853057 * 10**(-3) * df_protonRecCD.Ptheta**3 - 9.19989908 * 10**(-6) * df_protonRecCD.Ptheta**4
        #     coeff_CD = 2.99700765*1000 - 2.18027982 * 10**2 * df_protonRecCD.Ptheta + 5.84757503 * df_protonRecCD.Ptheta**2 - 6.80409195 * 0.01 * df_protonRecCD.Ptheta**3 + 2.89244618 * 0.0001 * df_protonRecCD.Ptheta**4
        #     coeff2_CD = -1.82237904*100 + 1.10153549 * 10 * df_protonRecCD.Ptheta - 2.24699931*0.1 * df_protonRecCD.Ptheta**2 + 1.49390960*0.001*df_protonRecCD.Ptheta**3

        #     CorrectedPtheta_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Ptheta"]

        #     const_CD = 7.58761670 - 5.28224578*0.1 * df_protonRecCD.Ptheta +  1.31580117 * 0.01 * df_protonRecCD.Ptheta**2 -1.41738951*10**(-4) * df_protonRecCD.Ptheta**3 + 5.62884363*10**(-7) * df_protonRecCD.Ptheta**4
        #     coeff_CD = 1.07644097*10**(5) - 8.67994639*10**(3) * df_protonRecCD.Ptheta + 2.57187193*10**(2) * df_protonRecCD.Ptheta**2 - 3.31379317 * df_protonRecCD.Ptheta**3 + 1.56896621*10**(-2) * df_protonRecCD.Ptheta**4
        #     coeff2_CD =  1.92263184*10**(2) -1.00870704 * 10 * df_protonRecCD.Ptheta + 1.56575252*10**(-1) * df_protonRecCD.Ptheta**2 -7.71489734*10**(-4) * df_protonRecCD.Ptheta**3
        #     CorrectedPphi_CD = const_CD + coeff_CD*np.exp(coeff2_CD*df_protonRecCD.loc[:, "Pp"]) + df_protonRecCD.loc[:, "Pphi"]

        df_protonRec.loc[df_protonRec["Psector"]<7, "Pp"] = CorrectedPp_FD
        df_protonRec.loc[df_protonRec["Psector"]<7, "Ptheta"] = CorrectedPtheta_FD
        df_protonRec.loc[df_protonRec["Psector"]<7, "Pphi"] = CorrectedPphi_FD

        df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), "Pp"] = CorrectedPp_CD
        df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), "Ptheta"] = CorrectedPtheta_CD
        df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), "Pphi"] = CorrectedPphi_CD

        df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
        df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
        df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

        # df_gammaRec = df_gammaRec[df_gammaRec["Gsector"]<7]

        df_gg = pd.merge(df_gammaRec, df_gammaRec,
                         how='outer', on='event', suffixes=("", "2"))
        df_gg = df_gg[df_gg["GIndex"] < df_gg["GIndex2"]]
        df_ep = pd.merge(df_electronRec, df_protonRec, how='outer', on='event')

        df_epgg = pd.merge(df_ep, df_gg, how='outer', on='event')
        df_epgg = df_epgg[~np.isnan(df_epgg["Ppx"])]
        df_epgg = df_epgg[~np.isnan(df_epgg["Gpx"])]
        df_epgg = df_epgg[~np.isnan(df_epgg["Gpx2"])]

        self.df_epgg = df_epgg #temporarily save df_epgg

        df_epg = pd.merge(df_ep, df_gammaRec, how='outer', on='event')
        df_epg = df_epg[~np.isnan(df_epg["Ppx"])]
        df_epg = df_epg[~np.isnan(df_epg["Gpx"])]

        self.df_epg = df_epg #temporarily save df_epgg

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

        # proton reconstruction quality
        cut_FD_proton = (df_epgg.loc[:, "Psector"]<7) & (df_epgg.loc[:, "Ptheta"]<37)
        cut_CD_proton = (df_epgg.loc[:, "Psector"]>7) & (df_epgg.loc[:, "Ptheta"]>40) & (df_epgg.loc[:, "Ptheta"]<66)
        cut_proton = (cut_FD_proton)|(cut_CD_proton)

        # Exclusivity cuts
        cut_mmep = df_epgg.loc[:, "MM2_ep"] < 0.7  # mmep
        cut_meepgg = df_epgg.loc[:, "ME_epgg"] < 0.7  # meepgg
        cut_mpt = df_epgg.loc[:, "MPt"] < 0.2  # mpt
        cut_recon = df_epgg.loc[:, "reconPi"] < 2  # recon gam angle
        cut_pi0upper = df_epgg.loc[:, "Mpi0"] < 0.2
        cut_pi0lower = df_epgg.loc[:, "Mpi0"] > 0.07
        cut_sector = (df_epgg.loc[:, "Esector"]!=df_epgg.loc[:, "Gsector"]) & (df_epgg.loc[:, "Esector"]!=df_epgg.loc[:, "Gsector2"])
        cut_Vz = np.abs(df_epgg["Evz"] - df_epgg["Pvz"]) < 2.5 + 2.5 / mag([df_epgg["Ppx"], df_epgg["Ppy"], df_epgg["Ppz"]])

        df_dvpi0 = df_epgg.loc[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_proton & cut_mmep & cut_meepgg & cut_Vz &
                           cut_mpt & cut_recon & cut_pi0upper & cut_pi0lower & cut_sector, :]

        #For an event, there can be two gg's passed conditions above.
        #Take only one gg's that makes pi0 invariant mass
        #This case is very rare.
        #For now, duplicated proton is not considered.
        df_dvpi0 = df_dvpi0.sort_values(by=['closeness', 'Psector', 'Gsector'], ascending = [True, True, True])
        df_dvpi0 = df_dvpi0.loc[~df_dvpi0.event.duplicated(), :]
        df_dvpi0 = df_dvpi0.sort_values(by='event')        
        self.df_dvpi0 = df_dvpi0 #done with saving x

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
        self.df_epg = df_epg

    def makeDVCS(self):
        #make dvcs pairs
        df_dvcs = self.df_epg
        df_dvcs = df_dvcs[df_dvcs["MM2_eg"] > 0]  # mmeg

        cut_xBupper = df_dvcs["xB"] < 1  # xB
        cut_xBlower = df_dvcs["xB"] > 0  # xB
        cut_Q2 = df_dvcs["Q2"] > 1  # Q2
        cut_W = df_dvcs["W"] > 2  # W
        cut_Ee = df_dvcs["Ee"] > 2  # Ee
        cut_Ge = df_dvcs["Ge"] > 3  # Ge
        cut_Pp = mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]]) > 0.12  # Pp
        cut_Vz = np.abs(df_dvcs["Evz"] - df_dvcs["Pvz"]) < 2.5 + 2.5 / mag([df_dvcs["Ppx"], df_dvcs["Ppy"], df_dvcs["Ppz"]])

        # proton reconstruction quality
        # cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) & (df_dvcs.loc[:, "Ptheta"]<37)
        # cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) & (df_dvcs.loc[:, "Ptheta"]>40) & (df_dvcs.loc[:, "Ptheta"]<66)
        cut_FD_proton = (df_dvcs.loc[:, "Psector"]<7) #& (df_dvcs.loc[:, "Ptheta"]<37)
        cut_CD_proton = (df_dvcs.loc[:, "Psector"]>7) & (df_dvcs.loc[:, "Ptheta"]<66) #& (df_dvcs.loc[:, "Ptheta"]>40) 
        cut_proton = (cut_FD_proton)|(cut_CD_proton)

        #   Exclusivity cuts
        cut_mmepg = np.abs(df_dvcs["MM2_epg"]) < 0.1  # mmepg
        cut_mmep = np.abs(df_dvcs["MM2_ep"]) < 0.6  # mmep
        cut_mmegupper = df_dvcs["MM2_eg"] < 3  # mmeg
        cut_mmeglower = df_dvcs["MM2_eg"] > 0  # mmeg
        cut_meepgupper = df_dvcs["ME_epg"] < 1.5  # meepg
        cut_meepglower = df_dvcs["ME_epg"] > -0.5  # meepg
        cut_mpt = df_dvcs["MPt"] < 0.25  # mpt
        cut_cone = df_dvcs["coneAngle"] > 5  # coneangle
        cut_recon = df_dvcs["reconGam"] < 2.5  # recon gam angle
        cut_coplanarity = df_dvcs["coplanarity"] < 25  # coplanarity angle
        if "Esector" in df_dvcs:
            cut_sector = df_dvcs["Esector"]!=df_dvcs["Gsector"]
        else:
            cut_sector = 1

        df_dvcs = df_dvcs[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_proton & cut_Ee & cut_Ge & cut_Pp & cut_Vz & cut_mmepg & cut_mmep &
                         cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon & cut_sector]

        #dealing with duplicates
        # df_dvcs = df_dvcs.sort_values(by=['Ge', 'Psector', 'Gsector'], ascending = [False, True, True])
        # df_dvcs = df_dvcs.loc[~df_dvcs.event.duplicated(), :]
        # df_dvcs = df_dvcs.sort_values(by='event')
        self.df_dvcs = df_dvcs               
         

    def pi02gSubtraction(self):
        #exclude dvpi0 from dvcs. use only when both set up.
        df_dvcs = self.df_dvcs
        pi0to2gammas = df_dvcs["event"].isin(self.df_dvpi0["event"])
        df_dvcs = df_dvcs[~pi0to2gammas]
        self.df_dvcs = df_dvcs

    def saveRaw(self):
        df_x = self.df_dvcs
        df_z = self.df_z
        df = pd.merge(df_x, df_z, how = 'inner', on='event')
        self.df = df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-g","--generator", help="choose dvcs or pi0", default = "dvcs")
    parser.add_argument("-p","--polarity", help="polarity", default = "inbending")
    
    args = parser.parse_args()

    converter = root2pickle(args.fname, entry_stop = args.entry_stop, gen = args.generator, pol = args.polarity)
    df = converter.df

    df.to_pickle(args.out)