import pandas as pd
import numpy as np
import argparse
from utils.const import *
from utils.physics import *
import subprocess, os
from gepard.fits import th_KM15
import gepard as g

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-file", "--filename", type = str)
    parser.add_argument("-Nind", "--Nind", type = int)
    parser.add_argument("-Nmax", "--Nmax", type = int)

    args = parser.parse_args()
    filename = args.filename
    Nind = args.Nind
    Nmax = args.Nmax

    df = pd.read_pickle(args.filename)

    if not df.index.is_unique:
      df.index = np.linspace(0, len(df)-1, len(df), dtype = int)

    indices = np.array_split(df.index.to_numpy(), Nmax)[Nind]

    df.loc[:, "BornWeight_KM"] = 0
    df.loc[:, "RadWeight_KM"] = 0

    for index in indices:
      this_row = df.loc[df.index == index, :]

      ebeam = 10.604
      pbeam = np.sqrt(ebeam * ebeam - me * me)
      GenxB = this_row.GenxB[index]
      GenQ2 = this_row.GenQ2[index]
      Gent = this_row.Gent[index]
      Genphi = this_row.Genphi[index]

      ds_born = printKM(GenxB, GenQ2, Gent, np.radians(Genphi))

      if this_row.radMode[index] == 1:
        ds_rad = ds_born
      elif this_row.radMode[index] == 2:
        ebeam = 10.604 - this_row.GenGp2[index]
        pbeam = np.sqrt(ebeam * ebeam - me * me)
        GenEe = np.sqrt(this_row.GenEp[index]**2 + me**2)
        VGS = [-this_row.GenEpx[index], -this_row.GenEpy[index], pbeam - this_row.GenEpz[index]]
        GenQ2 = -((ebeam - GenEe)**2 - mag2(VGS))
        Gennu = (ebeam - GenEe)
        GenxB = GenQ2/2.0/M/Gennu
        ds_rad = printKM(GenxB, GenQ2, Gent, np.radians(Genphi))
      elif this_row.radMode[index] == 3:
        GenEe = np.sqrt((this_row.GenEp[index] + this_row.GenGp2[index])**2 + me**2)
        VGS = [-(this_row.GenEpx + this_row.GenGpx2)[index], -(this_row.GenEpy + this_row.GenGpy2)[index], pbeam - (this_row.GenEpz + this_row.GenGpz2)[index]]
        GenQ2 = -((ebeam - GenEe)**2 - mag2(VGS))
        Gennu = (ebeam - GenEe)
        GenxB = GenQ2/2.0/M/Gennu
        ds_rad = printKM(GenxB, GenQ2, Gent, np.radians(Genphi))
      df.loc[df.index == index, "BornWeight_KM"] = ds_born/(2*np.pi/1000)
      df.loc[df.index == index, "RadWeight_KM"] = ds_rad/(2*np.pi/1000)

      new_filename = filename.replace('.pkl', '_{}_{}.pkl'.format(Nind, Nmax))
      df.loc[df.index.isin(indices), ["BornWeight_KM", "RadWeight_KM"]].to_pickle(new_filename)