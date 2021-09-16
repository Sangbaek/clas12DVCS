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

        self.readEP(entry_stop, pol = pol)
        self.saveRaw()

    def readFile(self):
        #read root using uproot
        self.file = uproot.open(self.fname)
        self.tree = self.file["T"]

    def closeFile(self):
        #close file for saving memory
        self.file = None
        self.tree = None

    def readEP(self, entry_stop = None, pol = "inbending"):
        #save data into df_epg, df_ep for parent class epg
        self.readFile()

        # data frames and their keys to read Z part
        df_electronGen = pd.DataFrame()
        df_protonGen = pd.DataFrame()
        # eleKeysGen = ["GenEpx", "GenEpy", "GenEpz", "GenEvx", "GenEvy", "GenEvz"]
        eleKeysGen = ["GenEpx", "GenEpy", "GenEpz"]
        proKeysGen = ["GenPpx", "GenPpy", "GenPpz"]
        # read keys
        for key in eleKeysGen:
            df_electronGen[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        for key in proKeysGen:
            df_protonGen[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)

        #convert data type to standard double
        # df_electronGen = df_electronGen.astype({"GenEpx": float, "GenEpy": float, "GenEpz": float, "GenEvx": float, "GenEvy": float, "GenEvz": float})
        df_electronGen = df_electronGen.astype({"GenEpx": float, "GenEpy": float, "GenEpz": float})
        df_protonGen = df_protonGen.astype({"GenPpx": float, "GenPpy": float, "GenPpz": float})

        #set up a dummy index for merging
        df_electronGen.loc[:,'event'] = df_electronGen.index
        df_protonGen.loc[:,'event'] = df_protonGen.index

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

        self.df_z = pd.merge(df_electronGen, df_protonGen, how='inner', on='event')

        # data frames and their keys to read X part
        df_electronRec = pd.DataFrame()
        df_protonRec = pd.DataFrame()
        # eleKeysRec = ["Epx", "Epy", "Epz", "Evx", "Evy", "Evz", "Esector"]
        # proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Psector"]
        eleKeysRec = ["Epx", "Epy", "Epz", "Evz", "Esector"]
        proKeysRec = ["Ppx", "Ppy", "Ppz", "Pvz", "Psector"]
        proKeysRec.extend(["PDc1Hitx", "PDc1Hity", "PDc1Hitz"])
        # read them
        for key in eleKeysRec:
            df_electronRec[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        for key in proKeysRec:
            df_protonRec[key] = self.tree[key].array(library="pd", entry_stop=entry_stop)
        self.closeFile()

        #convert data type to standard double
        # df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float, "Evx": float, "Evy": float, "Evz": float})
        # df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float, "Pvz": float})
        df_electronRec = df_electronRec.astype({"Epx": float, "Epy": float, "Epz": float, "Evz": float})
        df_protonRec = df_protonRec.astype({"Ppx": float, "Ppy": float, "Ppz": float, "Pvz": float})

        #set up a dummy index for merging
        df_electronRec.loc[:,'event'] = df_electronRec.index
        df_protonRec.loc[:,'event'] = df_protonRec.index.get_level_values('entry')

        #save only FD protons and photons
        # df_protonRec = df_protonRec[df_protonRec["Psector"]<7]
        #proton momentum correction
        pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]
        df_protonRec.loc[:, 'Pp'] = mag(pro)
        df_protonRec.loc[:, 'Ptheta'] = getTheta(pro)
        df_protonRec.loc[:, 'Pphi'] = getPhi(pro)

        df_protonRecFD = df_protonRec.loc[df_protonRec.Psector<7, :]
        df_protonRecCD = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), :]
        df_protonRecOthers = df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta>=75), :]

        correction = False

        #inbending
        if pol == "inbending":
            def corr(x, t):
                x0, x1, x2, x3 = x
                return x0 + x1*np.power(t-np.ones(len(t))*0.3, x3)

            df_protonRecFD = df_protonRecFD.loc[df_protonRec.Pp > 0.3, :]
            df_protonRecFD.loc[:, "DC1theta"] = getTheta([df_protonRecFD.PDc1Hitx, df_protonRecFD.PDc1Hity, df_protonRecFD.PDc1Hitz])
            best_params = [-53.14680163254601, 79.61307254040804, 0.3, 0.05739232362022314]
            df_protonRecFD_check1 = df_protonRecFD.loc[df_protonRecFD.DC1theta < corr(best_params, df_protonRecFD.Pp), :]
            df_protonRecFD_check2 = df_protonRecFD.loc[df_protonRecFD.DC1theta >= corr(best_params, df_protonRecFD.Pp), :]
            # FD part
            # const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [-0.0123049 + 0.00028887*df_protonRecFD.Ptheta, -0.138227479 + 8.07557430*0.001*df_protonRecFD.Ptheta -1.34807927*0.0001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -0.0275235])
            # coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [0.01528006 - 0.00024079*df_protonRecFD.Ptheta, 5.65817597*0.01 -2.36903348*0.001*df_protonRecFD.Ptheta + 4.93780046*0.00001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 0.03998975])    

            # CorrectedPp_FD = const_FD + coeff_FD/df_protonRecFD.loc[:, "Pp"] + df_protonRecFD.loc[:, "Pp"]

            # const_FD = np.select([df_protonRecFD.Ptheta<19.5, (df_protonRecFD.Ptheta>=19.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<39), (df_protonRecFD.Ptheta>=39) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [2.63643690*0.01, 0.50047232 -0.03834672 *df_protonRecFD.Ptheta + 0.00071967*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 6.91308654 - 0.439839300*df_protonRecFD.Ptheta +6.83075548*0.001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 1.59424606, 1.47198581*10])
            # coeff_FD = np.select([df_protonRecFD.Ptheta<19.5, (df_protonRecFD.Ptheta>=19.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<39), (df_protonRecFD.Ptheta>=39) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [-1.46440415, 74.99891704  -6.1576777*df_protonRecFD.Ptheta + 0.11469137*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 682.909471 - 43.9551177 * df_protonRecFD.Ptheta + 0.682383790 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, -8.19627119, -23.55701865])    
            # coeff2_FD = np.select([df_protonRecFD.Ptheta<19.5, (df_protonRecFD.Ptheta>=19.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<39), (df_protonRecFD.Ptheta>=39) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [-3.47690993, 47.71351973 -4.34918241*df_protonRecFD.Ptheta + 0.08841191*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 100.33995753 - 6.96600416*df_protonRecFD.Ptheta + 0.11223046*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -1.25261927, -0.40113733])    

            # CorrectedPtheta_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Ptheta"]

            # const_FD = np.select([df_protonRecFD.Ptheta<16.5, (df_protonRecFD.Ptheta>=16.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [-0.190662844, -0.20725736 -0.00675627 *df_protonRecFD.Ptheta + 0.0007863*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 12.1881698 - 0.78906294*df_protonRecFD.Ptheta +0.01297898*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -4.59743066*10])
            # coeff_FD = np.select([df_protonRecFD.Ptheta<16.5, (df_protonRecFD.Ptheta>=16.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [6.48745941, 142.96379788  -16.66339055*df_protonRecFD.Ptheta + 0.51311212*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, 2.1853046 + 5.78521226 * df_protonRecFD.Ptheta - 0.09727796 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, 7.46969457*10])    
            # coeff2_FD = np.select([df_protonRecFD.Ptheta<16.5, (df_protonRecFD.Ptheta>=16.5) & (df_protonRecFD.Ptheta<27), (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<42), df_protonRecFD.Ptheta>=42],
            #                   [-3.14646608, 17.39529095 -1.78403359*df_protonRecFD.Ptheta + 0.0335692*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -1.03655317*10 + 0.161333213*df_protonRecFD.Ptheta -1.29625675*0.001*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -4.41246899*0.1])    

            # CorrectedPphi_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Pphi"]

            const_FD = -0.00051894 - 0.00018104 * df_protonRecFD_check1.Ptheta
            coeff_FD = 3.29466917*10**(-3) +  5.73663160*10**(-4) * df_protonRecFD_check1.Ptheta - 1.40807209 * 10**(-5) * df_protonRecFD_check1.Ptheta * df_protonRecFD_check1.Ptheta
            CorrectedPp_FD_1 = const_FD + coeff_FD/df_protonRecFD_check1.loc[:, "Pp"] + df_protonRecFD_check1.loc[:, "Pp"]

            const_FD = np.select([df_protonRecFD_check1.Ptheta<8, (df_protonRecFD_check1.Ptheta>=8) & (df_protonRecFD_check1.Ptheta<20), (df_protonRecFD_check1.Ptheta>=20)],
                              [-0.3828741 + 0.05081733 * df_protonRecFD_check1.Ptheta, 0, -0.2425042 + 0.01725155 * df_protonRecFD_check1.Ptheta])
            coeff_FD = np.select([df_protonRecFD_check1.Ptheta<8, (df_protonRecFD_check1.Ptheta>=8) & (df_protonRecFD_check1.Ptheta<20), (df_protonRecFD_check1.Ptheta>=20)],
                              [-29.52213218 + 6.87345942 * df_protonRecFD_check1.Ptheta, 0,  2.08762216 - 0.11464056 * df_protonRecFD_check1.Ptheta])    
            coeff2_FD = np.select([df_protonRecFD_check1.Ptheta<8, (df_protonRecFD_check1.Ptheta>=8) & (df_protonRecFD_check1.Ptheta<20), (df_protonRecFD_check1.Ptheta>=20)],
                              [6.63357631 - 1.66750966 * df_protonRecFD_check1.Ptheta, 0, 1.37736738 - 0.09848285 * df_protonRecFD_check1.Ptheta])    

            CorrectedPtheta_FD_1 = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD_check1.loc[:, "Pp"]) + df_protonRecFD_check1.loc[:, "Ptheta"]

            const_FD = np.select([df_protonRecFD_check1.Ptheta<15, df_protonRecFD_check1.Ptheta>=15],
                              [0, 4.53922278 - 0.72676857*df_protonRecFD_check1.Ptheta + 0.02197507*df_protonRecFD_check1.Ptheta*df_protonRecFD_check1.Ptheta])
            coeff_FD = np.select([df_protonRecFD_check1.Ptheta<15, df_protonRecFD_check1.Ptheta>=15],
                              [0, -13.23900015  + 1.70084877*df_protonRecFD_check1.Ptheta  -0.04636718*df_protonRecFD_check1.Ptheta*df_protonRecFD_check1.Ptheta])
            coeff2_FD = np.select([df_protonRecFD_check1.Ptheta<15, df_protonRecFD_check1.Ptheta>=15],
                              [0, -3.25471872 + 0.38338761*df_protonRecFD_check1.Ptheta  -0.0122272*df_protonRecFD_check1.Ptheta*df_protonRecFD_check1.Ptheta])
            CorrectedPphi_FD_1 = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD_check1.loc[:, "Pp"]) + df_protonRecFD_check1.loc[:, "Pphi"]

            const_FD = -3.03346359*10**(-1) + 1.83368163*10**(-2)*df_protonRecFD_check2.Ptheta - 2.86486404*10**(-4)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta
            coeff_FD =  2.01023276*10**(-1) - 1.13312215*10**(-2)*df_protonRecFD_check2.Ptheta + 1.82487916*10**(-4)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta
            CorrectedPp_FD_2 = const_FD + coeff_FD/df_protonRecFD_check2.loc[:, "Pp"] + df_protonRecFD_check2.loc[:, "Pp"]

            const_FD = np.select([df_protonRecFD_check2.Ptheta<30, (df_protonRecFD_check2.Ptheta>=30)],
                              [2.98140652*10**(2) - 2.75743096*10**(1)*df_protonRecFD_check2.Ptheta + 8.35818159*10**(-1)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta - 8.26147941*10**(-3)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta, 1.34032831*10**(2) - 1.12452272*10**(1)*df_protonRecFD_check2.Ptheta + 3.12125195*10**(-1)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta -2.86545130*10**(-3)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta])
            coeff_FD = np.select([df_protonRecFD_check2.Ptheta<30, (df_protonRecFD_check2.Ptheta>=30)],
                              [1.82790220*10**(6) - 1.98950907*10**(5)*df_protonRecFD_check2.Ptheta + 7.20647164*10**(3)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta - 8.68597440*10**(1)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta, 9.37101732*10**(2) - 7.92142664*10**(1)*df_protonRecFD_check2.Ptheta + 2.22452861*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta -2.08191298*10**(-2)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta])    
            coeff2_FD = np.select([df_protonRecFD_check2.Ptheta<30, (df_protonRecFD_check2.Ptheta>=30)],
                              [1.09369187*10**(4) - 1.29345623*10**(3)*df_protonRecFD_check2.Ptheta + 5.07524944*10**(1)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta - 6.61520585*10**(-1)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta, 1.02129004*10**(3) - 8.43398593*10**(1)*df_protonRecFD_check2.Ptheta + 2.29760161*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta  -2.07085860*10**(-2)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta])
            CorrectedPtheta_FD_2 = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD_check2.loc[:, "Pp"]) + df_protonRecFD_check2.loc[:, "Ptheta"]

            const_FD = 0.54697831 -0.04896981*df_protonRecFD_check2.Ptheta +  0.00111376*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta
            coeff_FD = -4.06733541*10**2 + 2.43696202*10*df_protonRecFD_check2.Ptheta -3.36144736*10**(-1)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta
            coeff2_FD = 2.06378660*10 - 1.42866062*df_protonRecFD_check2.Ptheta + 2.01085440*10**(-2)*df_protonRecFD_check2.Ptheta*df_protonRecFD_check2.Ptheta
            CorrectedPphi_FD_2 = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD_check2.loc[:, "Pp"]) + df_protonRecFD_check2.loc[:, "Pphi"]



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
            const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27)],
                              [0.02067157-0.0009827*df_protonRecFD.Ptheta, -0.11216694 + 0.0069912*df_protonRecFD.Ptheta - 0.00011733 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta])
            coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27)],
                              [-0.03334437+0.00177781*df_protonRecFD.Ptheta, 0.0402797945 - 0.00197220505*df_protonRecFD.Ptheta + 4.50918200*10**(-5) * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta])

            CorrectedPp_FD = const_FD + coeff_FD/df_protonRecFD.loc[:, "Pp"] + df_protonRecFD.loc[:, "Pp"]

            const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
                              [0, -1.79343987 +0.105559096 *df_protonRecFD.Ptheta + -0.00157174358*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -0.123044632])
            coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
                              [0, -27.4344526 + 1.61037587* df_protonRecFD.Ptheta - 0.0242300381* df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, -7.52117236])    
            coeff2_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
                              [0, -45.2983842 +2.51745350*df_protonRecFD.Ptheta - 0.0365942178*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -3.52825441])    

            CorrectedPtheta_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Ptheta"]

            const_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
                              [0, 5.37967179 -0.324630795 *df_protonRecFD.Ptheta + 0.00476947696*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -0.0224918574])
            coeff_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
                              [0, 7.25038499*1000 + -413.586911* df_protonRecFD.Ptheta + 5.91815405 * df_protonRecFD.Ptheta * df_protonRecFD.Ptheta, 55.6319490])    
            coeff2_FD = np.select([df_protonRecFD.Ptheta<27, (df_protonRecFD.Ptheta>=27) & (df_protonRecFD.Ptheta<38), df_protonRecFD.Ptheta>=38],
                              [0, -124.626261 + 6.77668728*df_protonRecFD.Ptheta - 0.0960045129*df_protonRecFD.Ptheta*df_protonRecFD.Ptheta, -5.12646023])    

            CorrectedPphi_FD = const_FD + coeff_FD*np.exp(coeff2_FD*df_protonRecFD.loc[:, "Pp"]) + df_protonRecFD.loc[:, "Pphi"]

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
            print("correction applied for " + pol)

            df_protonRecCD.loc[:, "Pp"] = CorrectedPp_CD
            df_protonRecCD.loc[:, "Ptheta"] = CorrectedPtheta_CD
            df_protonRecCD.loc[:, "Pphi"] = CorrectedPphi_CD

            if pol == "inbending":
                df_protonRecFD_check1.loc[:, "Pp"] = CorrectedPp_FD_1
                df_protonRecFD_check1.loc[:, "Ptheta"] = CorrectedPtheta_FD_1
                df_protonRecFD_check1.loc[:, "Pphi"] = CorrectedPphi_FD_1

                df_protonRecFD_check2.loc[:, "Pp"] = CorrectedPp_FD_2
                df_protonRecFD_check2.loc[:, "Ptheta"] = CorrectedPtheta_FD_2
                df_protonRecFD_check2.loc[:, "Pphi"] = CorrectedPphi_FD_2

                df_protonRecFD = pd.concat([df_protonRecFD_check1, df_protonRecFD_check2])
                df_protonRecFD = df_protonRecFD.drop("DC1theta", axis = 1)

            if pol == "outbending":
                df_protonRecFD.loc[:, "Pp"] = CorrectedPp_FD
                df_protonRecFD.loc[:, "Ptheta"] = CorrectedPtheta_FD
                df_protonRecFD.loc[:, "Pphi"] = CorrectedPphi_FD

            df_protonRec = pd.concat([df_protonRecFD, df_protonRecCD, df_protonRecOthers])

            # df_protonRec.loc[df_protonRec.Psector<7, :] = df_protonRecFD
            # df_protonRec.loc[(df_protonRec.Psector>7) & (df_protonRec.Ptheta<75), :] = df_protonRecCD

            df_protonRec.loc[:, "Ppx"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.cos(np.radians(df_protonRec.loc[:, "Pphi"]))
            df_protonRec.loc[:, "Ppy"] = df_protonRec.loc[:, "Pp"]*np.sin(np.radians(df_protonRec.loc[:, "Ptheta"]))*np.sin(np.radians(df_protonRec.loc[:, "Pphi"]))
            df_protonRec.loc[:, "Ppz"] = df_protonRec.loc[:, "Pp"]*np.cos(np.radians(df_protonRec.loc[:, "Ptheta"]))
            pro = [df_protonRec['Ppx'], df_protonRec['Ppy'], df_protonRec['Ppz']]

        df_protonRec.loc[:, 'Pe'] = getEnergy(pro, M)

        # df_gammaRec = df_gammaRec[df_gammaRec["Gsector"]<7]

        df_ep = pd.merge(df_electronRec, df_protonRec, how='outer', on='event')

        df_ep = df_ep[~np.isnan(df_ep["Ppx"])]

        self.df_ep = df_ep #temporarily save df_ep

    def saveRaw(self):
        df_x = self.df_ep
        df_z = self.df_z
        df = pd.merge(df_x, df_z, how = 'inner', on='event')
        self.df = df


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