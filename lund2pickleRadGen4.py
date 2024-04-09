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
        kinArray = []

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
                logLine = txtlst[ind]
                eleLine = txtlst[ind+1]
                proLine = txtlst[ind+2]
                gamLine = txtlst[ind+3]
                skip = 3 
            elif num_particles == "4":
                logLine = txtlst[ind]
                eleLine = txtlst[ind+1]
                proLine = txtlst[ind+2]
                gamLine = txtlst[ind+3]
                radLine = txtlst[ind+4]
                skip = 4
            elif num_particles == "5":
                logLine = txtlst[ind]
                eleLine = txtlst[ind+1]
                proLine = txtlst[ind+2]
                gamLine = txtlst[ind+3]
                radLine = txtlst[ind+4]
                radLine2 = txtlst[ind+5]
                skip = 5
            logQuantities = logLine.split()
            eleQuantities = eleLine.split()
            proQuantities = proLine.split()
            gamQuantities = gamLine.split()
            dsigma = logQuantities[-1]
            helicity = logQuantities[4]
            xB = eleQuantities[1]
            radMode = eleQuantities[5]
            Epx = float(eleQuantities[6])
            Epy = float(eleQuantities[7])
            Epz = float(eleQuantities[8])
            Q2 = eleQuantities[9]
            t1 = eleQuantities[10]
            phi1 = proQuantities[1]
            Ppx = float(proQuantities[6])
            Ppy = float(proQuantities[7])
            Ppz = float(proQuantities[8])
            beamE = gamQuantities[1]
            dsigma_born = gamQuantities[10]
            Gpx = float(gamQuantities[6])
            Gpy = float(gamQuantities[7])
            Gpz = float(gamQuantities[8])
            if num_particles == "3":
                Gpx2 = 0
                Gpy2 = 0
                Gpz2 = 0
                Gpx3 = 0
                Gpy3 = 0
                Gpz3 = 0
            elif num_particles == "4":
                radQuantities = radLine.split()
                Gpx2 = radQuantities[6]
                Gpy2 = radQuantities[7]
                Gpz2 = radQuantities[8]
                Gpx3 = 0
                Gpy3 = 0
                Gpz3 = 0
            elif num_particles == "5":
                radQuantities = radLine.split()
                Gpx2 = radQuantities[6]
                Gpy2 = radQuantities[7]
                Gpz2 = radQuantities[8]
                radQuantities2 = radLine2.split()
                Gpx3 = radQuantities2[6]
                Gpy3 = radQuantities2[7]
                Gpz3 = radQuantities2[8]
            Etheta = (180.0/np.pi)*np.arctan2(np.sqrt(Epx*Epx+Epy*Epy),Epz)
            Ptheta = (180.0/np.pi)*np.arctan2(np.sqrt(Ppx*Ppx+Ppy*Ppy),Ppz)
            Gtheta = (180.0/np.pi)*np.arctan2(np.sqrt(Gpx*Gpx+Gpy*Gpy),Gpz)

            if (Ptheta<40) and (Gtheta<5):
                config = 0; # FDFT
            elif (Ptheta<40) and (Gtheta>=5):
                config = 1; #FDFD
            elif (Ptheta>=40) and (Gtheta>=5):
                config = 2; #CDFD
            elif (Ptheta>=40) and (Gtheta<5):
                config = 3; #CDFT
            else:
                config = -1;#errorneous bit

            kinArray.append([float(xB), float(Q2), float(t1), np.degrees(float(phi1)), float(dsigma), float(dsigma_born), int(radMode), int(float(helicity)), config, float(beamE)])

        df_epgg = pd.DataFrame(kinArray, columns = ["xB", "Q2", "t1", "phi1", "GenWeight", "BornWeight", "radMode", "helicity", "config", "beamE"])
        self.df = df_epgg

        self.closeFile()

if __name__ == "__main__":
    print("run lund2pickle for radiative DVCS generator")

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-model","--model", help="cross section model", default=None)
    parser.add_argument("-bin","--bin", help="bin number", default=None)
    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)

    args = parser.parse_args()
    if args.model:
        dfs = []
        for filenum in range(1, 68+1):
            file = "/volatile/clas12/sangbaek/dvcs_related/sim_rad_gen/{0}/fall2018_inb3/{0}_{1}_{2}.dat".format(args.model, args.bin, filenum)
            converter = lund2pickle(file, entry_stop = args.entry_stop)
            df = converter.df
            dfs.append(df)
        dfs = pd.concat(dfs)
        dfs = dfs.reset_index()
        dfs = dfs.loc[:, dfs.columns[1:]]
        outfile = "/volatile/clas12/sangbaek/dvcs_related/sim_rad_gen/{0}/root/fall2018_inb3/{1}.pkl".format(args.model, args.bin)
        dfs.to_pickle(outfile)

        dfs = []
        for filenum in range(1+68, 68+1+68):
            file = "/volatile/clas12/sangbaek/dvcs_related/sim_rad_gen/{0}/fall2018_outb3/{0}_{1}_{2}.dat".format(args.model, args.bin, filenum)
            converter = lund2pickle(file, entry_stop = args.entry_stop)
            df = converter.df
            dfs.append(df)
        dfs = pd.concat(dfs)
        dfs = dfs.reset_index()
        dfs = dfs.loc[:, dfs.columns[1:]]
        outfile = "/volatile/clas12/sangbaek/dvcs_related/sim_rad_gen/{0}/root/fall2018_outb3/{1}.pkl".format(args.model, args.bin)
        dfs.to_pickle(outfile)

    else:
        converter = lund2pickle(args.fname, entry_stop = args.entry_stop)
        df = converter.df

        df.to_pickle(args.out)