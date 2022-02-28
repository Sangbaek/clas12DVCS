#!/usr/bin/env python3
"""
A script for skimming
"""

import uproot
import pandas as pd
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i","--infile", help="input", default="/volatile/clas12/sangbaek/nov2021/convPkl_full/inb/dvcs/3987.pkl")
    parser.add_argument("-o","--output", help="output directory", default="/SimtoDat/v2")
    
    args = parser.parse_args()

    infile = pd.read_pickle(args.infile)

    ParticleVars1 = ["Epx", "Epy", "Epz", "Ep", "Etheta", "Ephi", "Esector"]
    ParticleVars2 = ["Ppx", "Ppy", "Ppz", "Pp", "Ptheta", "Pphi", "Psector", "Pstat"]
    ParticleVars3 = ["Gpx", "Gpy", "Gpz", "Gp", "Gtheta", "Gphi", "Gsector", "GFid"]
    ParticleVars4 = ["Gpx2", "Gpy2", "Gpz2", "Gp2", "Gtheta2", "Gphi2", "Gsector2", "GFid2"]
    Exclvars1 = ['event', 'Mpx', 'Mpy', 'Mpz','Q2','nu','y','xB','t1','t2','W','phi1','phi2','MM2_epg','ME_epg','MM2_ep','MM2_eg','MPt','coneAngle','reconGam','coplanarity', 'config']
    Exclvars2 = ['event', 'Mpx', 'Mpy', 'Mpz','Q2','nu','xB','t1','t2','W','phi1','phi2','MPt', 'MM2_ep', 'MM2_egg', 'MM2_epgg', 'ME_epgg', 'Mpi0', 'reconPi', "Pie", 'coplanarity', 'coneAngle1', 'coneAngle2', 'closeness', 'config']
    DVCSvars = ParticleVars1 + ParticleVars2+ ParticleVars3 + Exclvars1
    DVpi0Pvars = ParticleVars1 + ParticleVars2 + ParticleVars3 + ParticleVars4 + Exclvars2

    if "MM2_epgg" in infile.columns():
    	skims = DVpi0Pvars
    else:
    	skims = DVCSvars

    outfile = infile.loc[:, skims]

    outfile.to_pickle(args.outfile)
