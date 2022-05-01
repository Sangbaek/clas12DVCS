import pandas as pd
import numpy as np
import argparse
from utils.const import *
from utils.physics import *
import subprocess, os
from gepard.fits import th_KM15
import gepard as g

def printDVCS(xB, Q2, t, phi, heli):
    my_env = os.environ.copy()
    my_env["CLASDVCS_PDF"] = "/volatile/clas12/sangbaek/rad2/"
    dstot = subprocess.check_output(['/volatile/clas12/sangbaek/rad2/dvcsgen', '--printdstot', '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '1', '--gpd', '101', '--phi', str(phi), '--heli', str(heli), '--vv2cut', '0.6', '--delta', '0.1', '--w', '3.61'], env = my_env)
    lines = dstot.decode("utf-8").split("\n")
    dstotline = lines[-2]
    if "xsec_obs" in dstotline:
        return float(dstotline.split()[1]), float(dstotline.split()[2])
    else:
        return 0, 0

def printDVCSarray(xBarray, Q2array, tarray, phiarray, heliarray):
    XsecObs = []
    XsecBorn = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        heliarray = heliarray.to_numpy()
        
    for xB, Q2, t, phi, heli in zip(xBarray, Q2array, tarray, phiarray, heliarray):
        XsecObs.append(printDVCS(xB, Q2, t, phi, heli)[0])
        XsecBorn.append(printDVCS(xB, Q2, t, phi, heli)[1])
    return XsecObs, XsecBorn

def printDVCSKMarray(xBarray, Q2array, tarray, phiarray, heliarray):
    KMarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi, heli in zip(xBarray, Q2array, tarray, phiarray, heliarray):
        KMarray.append(printDVCSKM(xB, Q2, t, phi, heli))
    return KMarray

def printDVCSKM(xB, Q2, t, phi, heli):
    phi = np.pi - phi
    pt1 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phi,
                  process='ep2epgamma', exptype='fixed target', frame ='trento',
                  in1energy=10.604, in1charge=-1, in1polarization=heli)
#     pt2 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phi,
#                    process='ep2epgamma', exptype='fixed target', frame = 'trento',
#                    in1energy=10.604, in1charge=-1, in1polarization=+1)
    return th_KM15.XS(pt1)/(2*np.pi/1000)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="file to be reweighted", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-S","--entry_start", help="entry_start", default="0")
    parser.add_argument("-s","--entry_stop", help="entry_stop", default="100")
    parser.add_argument("-o","--ofname", help="file to be saved", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    
    args = parser.parse_args()

    df = pd.read_pickle(args.fname)
    entry_start = int(args.entry_start)
    if args.entry_stop == "inf":
        df = df.loc[(df.index>=entry_start), :]
    else:
        entry_stop = int(args.entry_stop)
        df = df.loc[(df.index>=entry_start) & (df.index<entry_stop), :]
    XsecObs, XsecBorn = printDVCSarray(df.xB.to_numpy(), df.Q2.to_numpy(), df.t1.to_numpy(), np.radians(df.phi1.to_numpy()), df.helicity)
    XsecKMBorn = printDVCSKMarray(df.xB.to_numpy(), df.Q2.to_numpy(), df.t1.to_numpy(), np.radians(df.phi1.to_numpy()), df.helicity)

    df.loc[:, "BHrad"] = XsecObs
    df.loc[:, "BHborn"] = XsecBorn
    df.loc[:, "KMborn"] = XsecKMBorn
    df.to_pickle(args.ofname)