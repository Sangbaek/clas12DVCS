#!/usr/bin/env python3
"""
A simple script to save data in pickle.
"""

import pandas as pd
import numpy as np
import argparse
from copy import copy
from utils.const import *
from utils.physics import *

class lund2pickle():
    #class to read lund format to make epg pairs, inherited from epg
    def __init__(self, fname, entry_stop = None):
        self.fname = fname
        self.readEPG(entry_stop)

    def readFile(self):
        #read tsv file using python built-in functions
        with open(self.fname,"r") as file:
            self.data = file.read()

    def closeFile(self):
        #close file for saving memory
        self.data = None

    def readEPG(self, entry_stop = None):
        #save data into df_epg, df_epgg for parent class epg
        self.readFile()
        partArray = []

        txtlst = self.data.split("\n")
        txts = iter(txtlst)
        skip = 0
        for ind, line in enumerate(txtlst[:-1]):
            num_particles = line.split()[0]
            weight        = line.split()[-1]
            if skip > 0:
                skip = skip - 1
                continue
            if num_particles == "3":
                eleLine = txtlst[ind+1]
                proLine = txtlst[ind+2]
                gamLine = txtlst[ind+3]
                skip = 3 
            elif num_particles == "4":
                eleLine = txtlst[ind+1]
                proLine = txtlst[ind+2]
                gamLine = txtlst[ind+3]
                radLine = txtlst[ind+4]
                skip = 4
            elif num_particles == "5":
                eleLine = txtlst[ind+1]
                proLine = txtlst[ind+2]
                gamLine = txtlst[ind+3]
                radLine = txtlst[ind+4]
                radLine2 = txtlst[ind+5]
                skip = 5
            eleQuantities = eleLine.split()
            Epx = eleQuantities[6]
            Epy = eleQuantities[7]
            Epz = eleQuantities[8]
            proLine = txtlst[ind+2]
            proQuantities = proLine.split()
            Ppx = proQuantities[6]
            Ppy = proQuantities[7]
            Ppz = proQuantities[8]
            gamLine = txtlst[ind+3]
            gamQuantities = gamLine.split()
            Gpx = gamQuantities[6]
            Gpy = gamQuantities[7]
            Gpz = gamQuantities[8]
            if num_particles == "3":
                Gpx2 = 0
                Gpy2 = 0
                Gpz2 = 0
            elif num_particles == "4":
                radQuantities = radLine.split()
                Gpx2 = radQuantities[6]
                Gpy2 = radQuantities[7]
                Gpz2 = radQuantities[8]
            elif num_particles == "5":
                radQuantities = radLine.split()
                Gpx2 = radQuantities[6]
                Gpy2 = radQuantities[7]
                Gpz2 = radQuantities[8]
                radQuantities2 = radLine2.split()
                Gpx3 = radQuantities2[6]
                Gpy3 = radQuantities2[7]
                Gpz3 = radQuantities2[8]

            partArray.append([float(Epx), float(Epy), float(Epz), float(Ppx), float(Ppy), float(Ppz), float(Gpx), float(Gpy), float(Gpz), float(Gpx2), float(Gpy2), float(Gpz2), float(Gpx3), float(Gpy3), float(Gpz3), float(weight)])

        df_epgg = pd.DataFrame(partArray, columns = ["Epx", "Epy", "Epz", "Ppx", "Ppy", "Ppz", "Gpx", "Gpy", "Gpz", "Gpx2", "Gpy2", "Gpz2", "Gpx3", "Gpy3", "Gpz3", "weight"])
        self.ebeam = 10.604
        self.pbeam = np.sqrt(self.ebeam * self.ebeam - me * me)

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

        VGS = [-df_epgg['Epx'], -df_epgg['Epy'], 10.604 - df_epgg['Epz']]
        self.beam = [0, 0, self.pbeam]
        v3l = cross(self.beam, ele)
        v3h = cross(pro, VGS)
        v3g = cross(VGS, gam)

        df_epgg.loc[:,'Q2'] = -((10.604 - df_epgg['Ee'])**2 - mag2(VGS))
        df_epgg.loc[:,'nu'] = (10.604 - df_epgg['Ee'])
        df_epgg.loc[:,'xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
        df_epgg.loc[:,'t1'] = 2 * M * (df_epgg['Pe'] - M)

        VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
                    df_epgg["Ppy"], self.pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]
        VmissP = [-df_epgg["Epx"] - df_epgg["Gpx"] - df_epgg["Gpx2"], -df_epgg["Epy"] -
                    df_epgg["Gpy"] - df_epgg["Gpy2"], self.pbeam - df_epgg["Epz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        Vmiss = [-df_epgg["Epx"] - df_epgg["Ppx"] - df_epgg["Gpx"] - df_epgg["Gpx2"],
                    -df_epgg["Epy"] - df_epgg["Ppy"] - df_epgg["Gpy"] - df_epgg["Gpy2"],
                    self.pbeam - df_epgg["Epz"] - df_epgg["Ppz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
        costheta = cosTheta(VGS, gam)

        df_epgg.loc[:,'phi1'] = angle(v3l, v3h)
        df_epgg.loc[:,'phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                  df_epgg['phi1'], df_epgg['phi1'])

        self.df = df_epgg

        self.closeFile()

if __name__ == "__main__":
    print("run lund2pickle for radiative DVCS generator")

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)

    args = parser.parse_args()

    converter = lund2pickle(args.fname, entry_stop = args.entry_stop)
    df = converter.df

    df.to_pickle(args.out)