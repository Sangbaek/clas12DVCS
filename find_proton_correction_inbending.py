import uproot
import pandas as pd
import numpy as np
import argparse
from copy import copy
from utils.const import *
from utils.physics import *
from utils.fiducial import *
from utils.kinCorrection import *
import itertools
import os
import awkward as ak
from scipy.stats import entropy

ebeam = 10.604
pbeam = np.sqrt(ebeam * ebeam - me * me)
beam  = [0, 0, pbeam]

def saveDVCSvars(df_epg):
    df_epg = copy(df_epg)
    #set up dvcs variables

    ele = [df_epg['Epx'], df_epg['Epy'], df_epg['Epz']]
    df_epg.loc[:, 'Ep'] = mag(ele)
    df_epg.loc[:, 'Ee'] = getEnergy(ele, me)
    df_epg.loc[:, 'Etheta'] = getTheta(ele)
    df_epg.loc[:, 'Ephi'] = getPhi(ele)

    pro = [df_epg['Ppx'], df_epg['Ppy'], df_epg['Ppz']]
    df_epg.loc[:, 'Pp'] = mag(pro)
    df_epg.loc[:, 'Pe'] = getEnergy(pro, M)
    df_epg.loc[:, 'Ptheta'] = getTheta(pro)
    df_epg.loc[:, 'Pphi'] = getPhi(pro)

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
    # df_epg.loc[:,'t1Orig'] = 2 * M * (df_epg['PeOrig'] - M)
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

    df_epg.loc[:,'closeness2'] = np.abs(df_epg.MM2_ep)

    eps = 2*M*df_epg.xB / np.sqrt(df_epg.Q2)
    df_epg.loc[:,'ycol1'] = (df_epg.Q2-df_epg.t2)/(df_epg.Q2-df_epg.xB*df_epg.t2)
    df_epg.loc[:,'ycol2'] = 1 - (1-df_epg.xB)*df_epg.t2/df_epg.Q2
    df_epg.loc[:,'ymax1'] = 2*(np.sqrt(1+eps**2)-1)/(eps**2)
    df_epg.loc[:,'ymax2'] = 1 - (M**2)*(df_epg.xB**2)/df_epg.Q2
    df_epg.loc[:,'tmin1'] = df_epg.Q2*(2*(1-df_epg.xB)*(1-np.sqrt(1+eps**2))+eps**2)/(4*df_epg.xB*(1-df_epg.xB) + eps**2)
    df_epg.loc[:,'tmin2'] = M*M*(df_epg.xB**2)/(1-df_epg.xB+df_epg.xB*M*M/df_epg.Q2)
    df_epg.loc[:,'tcol'] = df_epg.Q2*(df_epg.Q2-2*df_epg.xB*M*ebeam)/df_epg.xB/(df_epg.Q2-2*M*ebeam)

    df_epg.loc[:, 'vzdiff'] = df_epg.Evz - df_epg.Pvz

    df_epg.loc[:, 'xBbin'] = np.zeros(len(df_epg.xB), dtype = 'int') - 1
    df_epg.loc[:, 'Q2bin'] = np.zeros(len(df_epg.Q2), dtype = 'int') - 1
    df_epg.loc[:, 'tbin'] = np.zeros(len(df_epg.t1), dtype = 'int') - 1
    df_epg.loc[:, 'phibin'] = np.zeros(len(df_epg.phi1), dtype = 'int') - 1
    for xB in newxBbins2:
        df_epg.xBbin = df_epg.xBbin + (df_epg.xB>xB).astype("int").to_numpy(dtype = 'int')

    for Q2 in newQ2bins2:
        df_epg.Q2bin = df_epg.Q2bin + (df_epg.Q2>Q2).astype("int").to_numpy(dtype = 'int')

    for t1 in newtbins:
        df_epg.tbin = df_epg.tbin + (df_epg.t1>t1).astype("int").to_numpy(dtype = 'int')

    for phi1 in phibins:
        df_epg.phibin = df_epg.phibin + (df_epg.phi1>phi1).astype("int").to_numpy(dtype = 'int')
    return df_epg

def saveDVpi0vars(df_epgg):
    df_epgg = copy(df_epgg)
    #set up pi0 variables
    # useful objects
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
    # df_epgg.loc[:,'t1Orig'] = 2 * M * (df_epgg['PeOrig'] - M)
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
    df_epgg.loc[:, "closeness2"] = np.abs(df_epgg.loc[:, "MM2_ep"] - .1349766**2)

    df_epgg.loc[:, 'vzdiff'] = df_epgg.Evz - df_epgg.Pvz

    df_epgg.loc[:, 'xBbin'] = np.zeros(len(df_epgg.xB), dtype = 'int') - 1
    df_epgg.loc[:, 'Q2bin'] = np.zeros(len(df_epgg.Q2), dtype = 'int') - 1
    df_epgg.loc[:, 'tbin'] = np.zeros(len(df_epgg.t1), dtype = 'int') - 1
    df_epgg.loc[:, 'phibin'] = np.zeros(len(df_epgg.phi1), dtype = 'int') - 1
    for xB in newxBbins2:
        df_epgg.xBbin = df_epgg.xBbin + (df_epgg.xB>xB).astype("int").to_numpy(dtype = 'int')

    for Q2 in newQ2bins2:
        df_epgg.Q2bin = df_epgg.Q2bin + (df_epgg.Q2>Q2).astype("int").to_numpy(dtype = 'int')

    for t1 in newtbins:
        df_epgg.tbin = df_epgg.tbin + (df_epgg.t1>t1).astype("int").to_numpy(dtype = 'int')

    for phi1 in phibins:
        df_epgg.phibin = df_epgg.phibin + (df_epgg.phi1>phi1).astype("int").to_numpy(dtype = 'int')

    # # encode unassigned bin as -1
    # df_epgg.loc[:, "Q2bin"] = -1
    # df_epgg.loc[:, "xBbin"] = -1
    # df_epgg.loc[:, "tbin"] = -1
    # # df_epgg.loc[:, "tbin2"] = -1
    # df_epgg.loc[:, "phibin"] = -1
    # # df_epgg.loc[:, "phibin2"] = -1
    # df_epgg.loc[:, "Q2xBbin"] = -1
    # df_epgg.loc[:, "Q2xBtbin"] = -1
    # # df_epgg.loc[:, "Q2xBtbin2"] = -1
    # df_epgg.loc[:, "Q2xBtphibin"] = -1
    # Q2xBbin = 0

    # # encode all binning
    # for Q2bin in range(len(Q2bin_i)):
    #     #square Q2 binning
    #     df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]), "Q2bin"] = Q2bin
    #     #adaptive xB binning
    #     for xBbin in range(len(xBbin_i[Q2bin])):
    #         if Q2bin < len(Q2bin_i) -1:
    #             if xBbin == 0:
    #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin #0
    #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin #0
    #             elif xBbin < len(xBbin_i[Q2bin])-1:
    #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "xBbin"] = xBbin
    #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.xB<xBbin_f[Q2bin][xBbin]), "Q2xBbin"] = Q2xBbin
    #             else:
    #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
    #                 df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.xB>=xBbin_i[Q2bin][xBbin]) & (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin
    #         else:
    #             df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "xBbin"] = xBbin
    #             df_epgg.loc[(df_epgg.Q2>=Q2bin_i[Q2bin]) & (df_epgg.Q2<Q2bin_f[Q2bin]) & (df_epgg.Q2<=2*M*(10.604-2)*df_epgg.xB)& (df_epgg.Q2>=(4-M*M)*df_epgg.xB/(1-df_epgg.xB)), "Q2xBbin"] = Q2xBbin #0

    #         Q2xBbin = Q2xBbin + 1
    # for tbin in range(len(tbin_i)):
    #     #square t binning
    #     df_epgg.loc[(df_epgg.t1>=tbin_i[tbin]) & (df_epgg.t1<tbin_f[tbin]), "tbin"] = tbin
    #     # df_epgg.loc[(df_epgg.t2>=tbin_i[tbin]) & (df_epgg.t2<tbin_f[tbin]), "tbin2"] = tbin
    # for phibin in range(len(phibin_i)):
    #     #square phi binning
    #     df_epgg.loc[(df_epgg.phi1>=phibin_i[phibin]) & (df_epgg.phi1<phibin_f[phibin]), "phibin"] = phibin
    #     # df_epgg.loc[(df_epgg.phi2>=phibin_i[phibin]) & (df_epgg.phi2<phibin_f[phibin]), "phibin2"] = phibin

    # df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtbin"] = len(tbin_i) * df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBbin"] + df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "tbin"]
    # # df_epgg.loc[(df_epgg.Q2bin>0)&(df_epgg.xBbin>0)&(df_epgg.tbin2>0), "Q2xBtbin2"] = df_epgg.Q2bin.astype(str) + df_epgg.xBbin.astype(str) + df_epgg.tbin2.astype(str)
    # df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtphibin"] = len(phibin_i) * df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "Q2xBtbin"] + df_epgg.loc[(df_epgg.Q2xBbin>=0)&(df_epgg.tbin>=0), "phibin"]

    # df_epgg = df_epgg.astype({"Q2bin": int, "xBbin": int, "tbin": int, "phibin": int, "Q2xBbin": int, "Q2xBtbin": int, "Q2xBtphibin": int})

    return df_epgg

charge_inb, charge_outb = (30046676.501082145, 32024144.64472983)
target_thickness = 6.022e23 * 0.07151 * 5 / 1.00794#1.079e23 # from inclusive note
electric_charge = 1.602176634e-19
luminosity_inb =  charge_inb * 10**(-9) * target_thickness / electric_charge * 10**(-24) * 10**(-9)
luminosity_outb =  charge_outb * 10**(-9) * target_thickness / electric_charge * 10**(-24) * 10**(-9)

bin_scheme = np.loadtxt('/work/clas12/sangbaek/km15gen/bin_scheme.csv', delimiter = ',')
fringe_bin_scheme = np.loadtxt('/work/clas12/sangbaek/km15gen/fringe_bin_scheme.csv', delimiter = ',')

N = int(1e6)
print("Calculate volume")
bin_volume_bulk = []
bin_volume_fringe = []
for bin in range(1, 147+1):
    this_xBmin, this_xBmax, this_Q2min, this_Q2max, this_t1min, this_t1max = bin_scheme[bin-1]
    xB = np.random.uniform(this_xBmin, this_xBmax, N)
    Q2 = np.random.uniform(this_Q2min, this_Q2max, N)
    t = np.random.uniform(this_t1min, this_t1max, N)
    tmaxs = tmax(xB, Q2, t, 0)
    tmins = tmin(xB, Q2, t, 0)
    tcond = (np.abs(t)<np.abs(tmaxs)) & (np.abs(t) > np.abs(tmins))
    yd = y(xB, Q2, t, 0)
    ycond = (yd>0.19) & (yd<0.85)
    cond = tcond & ycond
    bin_volume_bulk.append( np.sum(cond)/N *(this_xBmax-this_xBmin)*(this_Q2max-this_Q2min)*(this_t1max-this_t1min)*2*np.pi )

for bin in range(1, 159+1):
    this_xBmin, this_xBmax, this_Q2min, this_Q2max, this_t1min, this_t1max = fringe_bin_scheme[bin-1]
    xB = np.random.uniform(this_xBmin, this_xBmax, N)
    Q2 = np.random.uniform(this_Q2min, this_Q2max, N)
    t = np.random.uniform(this_t1min, this_t1max, N)
    tmaxs = tmax(xB, Q2, t, 0)
    tmins = tmin(xB, Q2, t, 0)
    tcond = (np.abs(t)<np.abs(tmaxs)) & (np.abs(t) > np.abs(tmins))
    yd = y(xB, Q2, t, 0)
    ycond = (yd>0.19) & (yd<0.85)
    cond = tcond & ycond
    bin_volume_fringe.append( np.sum(cond)/N *(this_xBmax-this_xBmin)*(this_Q2max-this_Q2min)*(this_t1max-this_t1min)*2*np.pi )

inb_list = np.loadtxt('/work/clas12/sangbaek/inb_run_list_pass1', dtype = int)
df_exp_dvcs_inb = []
df_exp_pi0_inb = []
for run in inb_list:
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_1/pkl/{:d}.pkl".format(run))
    df = df.loc[:, ["MM2_ep", "reconGam", "Pp", "Ptheta", "Psector"]]
    df_exp_dvcs_inb.append(df)
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/pi0/excl_level_1/pkl/{:d}.pkl".format(run))
    df = df.loc[:, ["MM2_ep", "reconPi", "Pp", "Ptheta", "Psector"]]
    df_exp_pi0_inb.append(df)
df_exp_dvcs_inb = pd.concat(df_exp_dvcs_inb).reset_index()
df_exp_dvcs_inb = df_exp_dvcs_inb.loc[:, df_exp_dvcs_inb.columns[1:]]
df_exp_pi0_inb = pd.concat(df_exp_pi0_inb).reset_index()
df_exp_pi0_inb = df_exp_pi0_inb.loc[:, df_exp_pi0_inb.columns[1:]]

df_sim_dvcs_inb = []
for bin in range(147):
    try:
        df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl/{}.pkl".format(bin+1))
        df = df.loc[:, ["MM2_ep", "reconGam", "Pp", "Ptheta", "Psector", "weights"]]
        df.loc[:, "weights"] = df.GenWeight * bin_volume_bulk[bin] * luminosity_inb / (680000*9999/10000)
        df_sim_dvcs_inb.append(df)
    except:
        continue
df_sim_dvcs_inb = pd.concat(df_sim_dvcs_inb).reset_index()
df_sim_dvcs_inb = df_sim_dvcs_inb.loc[:, df_sim_dvcs_inb.columns[1:]]

df_sim_pi0_1gamma_inb = []
for filenum in range(20):
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/pkl/1/{}.pkl".format(filenum+1))
    df = df.loc[:, ["MM2_ep", "reconGam", "Pp", "Ptheta", "Psector"]]
    df_sim_pi0_1gamma_inb.append(df)
df_sim_pi0_1gamma_inb = pd.concat(df_sim_pi0_1gamma_inb).reset_index()
df_sim_pi0_1gamma_inb = df_sim_pi0_1gamma_inb.loc[:, df_sim_pi0_1gamma_inb.columns[1:]]


df_sim_pi0_2gamma_inb = []
for filenum in range(20):
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/pkl/1/{}.pkl".format(filenum+1))
    df = df.loc[:, ["MM2_ep", "reconPi", "Pp", "Ptheta", "Psector"]]
    df_sim_pi0_2gamma_inb.append(df)
df_sim_pi0_2gamma_inb = pd.concat(df_sim_pi0_2gamma_inb).reset_index()
df_sim_pi0_2gamma_inb = df_sim_pi0_2gamma_inb.loc[:, df_sim_pi0_2gamma_inb.columns[1:]]

df_correction = pd.DataFrame()
Pbins   = [0.4, 0.6, 0.8, 1, 1.2, 1.4]

for Psector in [1, 2, 3, 4, 5, 6, "CD"]:
    for i in range(len(Pbins)-1):
        Pmin = Pbins[i]
        Pmax = Pbins[i+1]
        print ("{} < p < {}, Sector {}, {} Polarity".format(Pmin, Pmax, Psector, "Inbending"))
        Pcenter = (Pmin + Pmax)/2.
        
        scores  = []
        dps     = []
        dthetas = []
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_exp_pi0 = copy(df_exp_pi0_inb                 .loc[ (df_exp_pi0_inb.Pp             >Pmin) & (df_exp_pi0_inb.Pp         < Pmax) &  (df_exp_pi0_inb.Psector         == Psector), :])
        elif Psector == "CD":
            df_exp_pi0 = copy(df_exp_pi0_inb                 .loc[ (df_exp_pi0_inb.Pp             >Pmin) & (df_exp_pi0_inb.Pp         < Pmax) &  (df_exp_pi0_inb.Psector         >4000), :])
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_sim_pi0 = copy(df_sim_pi0_2gamma_inb          .loc[ (df_sim_pi0_2gamma_inb.Pp      >Pmin) & (df_sim_pi0_2gamma_inb.Pp  < Pmax) &  (df_sim_pi0_2gamma_inb.Psector  == Psector), :])
        elif Psector == "CD":
            df_sim_pi0 = copy(df_sim_pi0_2gamma_inb          .loc[ (df_sim_pi0_2gamma_inb.Pp      >Pmin) & (df_sim_pi0_2gamma_inb.Pp  < Pmax) &  (df_sim_pi0_2gamma_inb.Psector  >4000), :])
        sim_dist_pi0_mm2ep  , bins_pi0_mm2ep      = np.histogram(df_sim_pi0.loc[:, "MM2_ep"] , bins = 100, density = True)#, weights = df_sim_pi0.weights)
        sim_dist_pi0_reconPi, bins_pi0_reconPi    = np.histogram(df_sim_pi0.loc[:, "reconPi"], bins = 100, density = True)#, weights = df_sim_pi0.weights)
    
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_exp_dvcs = copy(df_exp_dvcs_inb               .loc[ (df_exp_dvcs_inb.Pp       >Pmin)       & (df_exp_dvcs_inb.Pp        < Pmax)  &  (df_exp_dvcs_inb.Psector       == Psector), :])
        elif Psector == "CD":
            df_exp_dvcs = copy(df_exp_dvcs_inb               .loc[ (df_exp_dvcs_inb.Pp       >Pmin)       & (df_exp_dvcs_inb.Pp        < Pmax)  &  (df_exp_dvcs_inb.Psector       >4000), :])
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_sim_dvcs = copy(df_sim_dvcs_inb               .loc[ (df_sim_dvcs_inb.Pp       >Pmin)       & (df_sim_dvcs_inb.Pp        < Pmax)  &  (df_sim_dvcs_inb.Psector       == Psector), :])
        elif Psector == "CD":
            df_sim_dvcs = copy(df_sim_dvcs_inb               .loc[ (df_sim_dvcs_inb.Pp       >Pmin)       & (df_sim_dvcs_inb.Pp        < Pmax)  &  (df_sim_dvcs_inb.Psector       >4000), :])
        sim_dist_dvcs_mm2ep   , bins_dvcs_mm2ep       = np.histogram(df_sim_dvcs.loc[:, "MM2_ep"]  , bins = 100, density = True, weights = df_sim_dvcs.weights)
        sim_dist_dvcs_reconGam, bins_dvcs_reconGam    = np.histogram(df_sim_dvcs.loc[:, "reconGam"], bins = 100, density = True, weights = df_sim_dvcs.weights)
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_bkg_dvcs = copy(df_sim_pi0_1gamma_inb         .loc[ (df_sim_pi0_1gamma_inb.Pp >Pmin)       & (df_sim_pi0_1gamma_inb.Pp < Pmax)   &  (df_sim_pi0_1gamma_inb.Psector == Psector), :])
        elif Psector == "CD":
            df_bkg_dvcs = copy(df_sim_pi0_1gamma_inb         .loc[ (df_sim_pi0_1gamma_inb.Pp >Pmin)       & (df_sim_pi0_1gamma_inb.Pp < Pmax)   &  (df_sim_pi0_1gamma_inb.Psector >4000), :])
        bkg_dist_dvcs_mm2ep  , _     = np.histogram(df_bkg_dvcs.loc[:, "MM2_ep"]  , bins = bins_dvcs_mm2ep, density = True)#, weights = df_bkg_dvcs.weights)
        bkg_dist_dvcs_reconGam, _    = np.histogram(df_bkg_dvcs.loc[:, "reconGam"], bins = bins_dvcs_reconGam, density = True)#, weights = df_bkg_dvcs.weights)
    
        try:
            cont  =  len(df_bkg_dvcs)/ len(df_sim_pi0) * len(df_exp_pi0) / len(df_exp_dvcs)
        except:
            cont  = 0
    
        sim_dist_epg_mm2ep    = (1 - cont)* sim_dist_dvcs_mm2ep    + cont* bkg_dist_dvcs_mm2ep
        sim_dist_epg_reconGam = (1 - cont)* sim_dist_dvcs_reconGam + cont* bkg_dist_dvcs_reconGam
    
        for trial in range(10000):
    
            dp     = np.random.uniform(-0.05, 0.05)
            dtheta = np.random.uniform(-2, 2)
            dphi   = 0
        
            df_exp_pi0 = copy(df_exp_pi0_inb)
            df_exp_pi0.loc[:, "Pp"]      = df_exp_pi0.Pp     + dp
            df_exp_pi0.loc[:, "Ptheta"]  = df_exp_pi0.Ptheta + dtheta
            df_exp_pi0.loc[:, "Pphi"]    = df_exp_pi0.Pphi   + dphi
            df_exp_pi0.loc[:, "Ppx"]     = df_exp_pi0.loc[:, "Pp"]*np.sin(np.radians(df_exp_pi0.loc[:, "Ptheta"]))*np.cos(np.radians(df_exp_pi0.loc[:, "Pphi"]))
            df_exp_pi0.loc[:, "Ppy"]     = df_exp_pi0.loc[:, "Pp"]*np.sin(np.radians(df_exp_pi0.loc[:, "Ptheta"]))*np.sin(np.radians(df_exp_pi0.loc[:, "Pphi"]))
            df_exp_pi0.loc[:, "Ppz"]     = df_exp_pi0.loc[:, "Pp"]*np.cos(np.radians(df_exp_pi0.loc[:, "Ptheta"]))
            
            if Psector in [1, 2, 3, 4, 5, 6]:
                df_exp_pi0 = copy(df_exp_pi0                 .loc[ (df_exp_pi0.Pp             >Pmin) & (df_exp_pi0.Pp         < Pmax) &  (df_exp_pi0.Psector         == Psector), :])
            elif Psector == "CD":
                df_exp_pi0 = copy(df_exp_pi0                 .loc[ (df_exp_pi0.Pp             >Pmin) & (df_exp_pi0.Pp         < Pmax) &  (df_exp_pi0.Psector         >4000), :])
            df_exp_pi0                   = saveDVpi0vars(df_exp_pi0)
            exp_dist_pi0_mm2ep  , _      = np.histogram(df_exp_pi0.loc[:, "MM2_ep" ]  , bins = bins_pi0_mm2ep, density = True  )#, weights = df_exp_pi0.weights)
            exp_dist_pi0_reconPi, _      = np.histogram(df_exp_pi0.loc[:, "reconPi"]  , bins = bins_pi0_reconPi, density = True)#, weights = df_exp_pi0.weights)
    
            df_exp_dvcs = copy(df_exp_dvcs_inb)
            df_exp_dvcs.loc[:, "Pp"]      = df_exp_dvcs.Pp     + dp
            df_exp_dvcs.loc[:, "Ptheta"]  = df_exp_dvcs.Ptheta + dtheta
            df_exp_dvcs.loc[:, "Pphi"]    = df_exp_dvcs.Pphi   + dphi
            df_exp_dvcs.loc[:, "Ppx"]     = df_exp_dvcs.loc[:, "Pp"]*np.sin(np.radians(df_exp_dvcs.loc[:, "Ptheta"]))*np.cos(np.radians(df_exp_dvcs.loc[:, "Pphi"]))
            df_exp_dvcs.loc[:, "Ppy"]     = df_exp_dvcs.loc[:, "Pp"]*np.sin(np.radians(df_exp_dvcs.loc[:, "Ptheta"]))*np.sin(np.radians(df_exp_dvcs.loc[:, "Pphi"]))
            df_exp_dvcs.loc[:, "Ppz"]     = df_exp_dvcs.loc[:, "Pp"]*np.cos(np.radians(df_exp_dvcs.loc[:, "Ptheta"]))
            if Psector in [1, 2, 3, 4, 5, 6]:
                df_exp_dvcs = copy(df_exp_dvcs                .loc[ (df_exp_dvcs.Pp       >Pmin) & (df_exp_dvcs.Pp       < Pmax) &  (df_exp_dvcs.Psector       == Psector), :])
            elif Psector == "CD":
                df_exp_dvcs = copy(df_exp_dvcs                .loc[ (df_exp_dvcs.Pp       >Pmin) & (df_exp_dvcs.Pp       < Pmax) &  (df_exp_dvcs.Psector       >4000), :])
    
            df_exp_dvcs                    = saveDVCSvars(df_exp_dvcs)
            exp_dist_epg_mm2ep   , _       = np.histogram(df_exp_dvcs.loc[:, "MM2_ep" ] , bins = bins_dvcs_mm2ep, density = True)#, weights = df_exp_dvcs.weights)
            exp_dist_epg_reconGam, _       = np.histogram(df_exp_dvcs.loc[:, "reconGam"], bins = bins_dvcs_reconGam, density = True)#, weights = df_exp_dvcs.weights)
    
            score = entropy(exp_dist_pi0_mm2ep, sim_dist_pi0_mm2ep) * entropy(exp_dist_pi0_reconPi, sim_dist_pi0_reconPi) * entropy(exp_dist_epg_mm2ep, sim_dist_epg_mm2ep) * entropy(exp_dist_epg_reconGam, sim_dist_epg_reconGam)
            scores  .append(score)
            dps     .append(dp)
            dthetas .append(dtheta)
            # print(score)
            # plt.show()
        dps             = np.array(dps)
        dthetas         = np.array(dthetas)
        optimal_dps     = dps    [np.argsort(scores)[:20]]
        optimal_dthetas = dthetas[np.argsort(scores)[:20]]
        optimal_dp      = np.mean(optimal_dps)
        optmial_dtheta  = np.mean(optimal_dthetas)
        optimal_dp_std     = np.std(optimal_dps)
        optimal_dtheta_std = np.std(optimal_dthetas)
    
        this_row = pd.DataFrame.from_dict({"Psector": [Psector], "p": [Pcenter], "dp": [optimal_dp], "dtheta": [optimal_dtheta], "dp_std": [optimal_dp_std], "dtheta_std": [optimal_dtheta_std], "polarity": ["inbending"]})
    
        df_correction = pd.concat([df_correction, this_row])

for Psector in [1, 2, 3, 4, 5, 6, "CD"]:
    for i in range(len(Pbins)-1):
        Pmin = Pbins[i]
        Pmax = Pbins[i+1]
        print ("{} < p < {}, Sector {}, {} Polarity".format(Pmin, Pmax, Psector, "Outbending"))
        Pcenter = (Pmin + Pmax)/2.
        
        scores  = []
        dps     = []
        dthetas = []
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_exp_pi0 = copy(df_exp_pi0_outb                 .loc[ (df_exp_pi0_outb.Pp             >Pmin) & (df_exp_pi0_outb.Pp         < Pmax) &  (df_exp_pi0_outb.Psector         == Psector), :])
        elif Psector == "CD":
            df_exp_pi0 = copy(df_exp_pi0_outb                 .loc[ (df_exp_pi0_outb.Pp             >Pmin) & (df_exp_pi0_outb.Pp         < Pmax) &  (df_exp_pi0_outb.Psector         >4000), :])
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_sim_pi0 = copy(df_sim_pi0_2gamma_outb          .loc[ (df_sim_pi0_2gamma_outb.Pp      >Pmin) & (df_sim_pi0_2gamma_outb.Pp  < Pmax) &  (df_sim_pi0_2gamma_outb.Psector  == Psector), :])
        elif Psector == "CD":
            df_sim_pi0 = copy(df_sim_pi0_2gamma_outb          .loc[ (df_sim_pi0_2gamma_outb.Pp      >Pmin) & (df_sim_pi0_2gamma_outb.Pp  < Pmax) &  (df_sim_pi0_2gamma_outb.Psector  >4000), :])
        sim_dist_pi0_mm2ep  , bins_pi0_mm2ep      = np.histogram(df_sim_pi0.loc[:, "MM2_ep"] , bins = 100, density = True)#, weights = df_sim_pi0.weights)
        sim_dist_pi0_reconPi, bins_pi0_reconPi    = np.histogram(df_sim_pi0.loc[:, "reconPi"], bins = 100, density = True)#, weights = df_sim_pi0.weights)
    
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_exp_dvcs = copy(df_exp_dvcs_outb               .loc[ (df_exp_dvcs_outb.Pp       >Pmin)       & (df_exp_dvcs_outb.Pp        < Pmax)  &  (df_exp_dvcs_outb.Psector       == Psector), :])
        elif Psector == "CD":
            df_exp_dvcs = copy(df_exp_dvcs_outb               .loc[ (df_exp_dvcs_outb.Pp       >Pmin)       & (df_exp_dvcs_outb.Pp        < Pmax)  &  (df_exp_dvcs_outb.Psector       >4000), :])
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_sim_dvcs = copy(df_sim_dvcs_outb               .loc[ (df_sim_dvcs_outb.Pp       >Pmin)       & (df_sim_dvcs_outb.Pp        < Pmax)  &  (df_sim_dvcs_outb.Psector       == Psector), :])
        elif Psector == "CD":
            df_sim_dvcs = copy(df_sim_dvcs_outb               .loc[ (df_sim_dvcs_outb.Pp       >Pmin)       & (df_sim_dvcs_outb.Pp        < Pmax)  &  (df_sim_dvcs_outb.Psector       >4000), :])
        sim_dist_dvcs_mm2ep   , bins_dvcs_mm2ep       = np.histogram(df_sim_dvcs.loc[:, "MM2_ep"]  , bins = 100, density = True, weights = df_sim_dvcs.weights)
        sim_dist_dvcs_reconGam, bins_dvcs_reconGam    = np.histogram(df_sim_dvcs.loc[:, "reconGam"], bins = 100, density = True, weights = df_sim_dvcs.weights)
        if Psector in [1, 2, 3, 4, 5, 6]:
            df_bkg_dvcs = copy(df_sim_pi0_1gamma_outb         .loc[ (df_sim_pi0_1gamma_outb.Pp >Pmin)       & (df_sim_pi0_1gamma_outb.Pp < Pmax)   &  (df_sim_pi0_1gamma_outb.Psector == Psector), :])
        elif Psector == "CD":
            df_bkg_dvcs = copy(df_sim_pi0_1gamma_outb         .loc[ (df_sim_pi0_1gamma_outb.Pp >Pmin)       & (df_sim_pi0_1gamma_outb.Pp < Pmax)   &  (df_sim_pi0_1gamma_outb.Psector >4000), :])
        bkg_dist_dvcs_mm2ep  , _     = np.histogram(df_bkg_dvcs.loc[:, "MM2_ep"]  , bins = bins_dvcs_mm2ep, density = True)#, weights = df_bkg_dvcs.weights)
        bkg_dist_dvcs_reconGam, _    = np.histogram(df_bkg_dvcs.loc[:, "reconGam"], bins = bins_dvcs_reconGam, density = True)#, weights = df_bkg_dvcs.weights)
    
        try:
            cont  =  len(df_bkg_dvcs)/ len(df_sim_pi0) * len(df_exp_pi0) / len(df_exp_dvcs)
        except:
            cont  = 0
    
        sim_dist_epg_mm2ep    = (1 - cont)* sim_dist_dvcs_mm2ep    + cont* bkg_dist_dvcs_mm2ep
        sim_dist_epg_reconGam = (1 - cont)* sim_dist_dvcs_reconGam + cont* bkg_dist_dvcs_reconGam
    
        for trial in range(10000):
    
            dp     = np.random.uniform(-0.05, 0.05)
            dtheta = np.random.uniform(-2, 2)
            dphi   = 0
        
            df_exp_pi0 = copy(df_exp_pi0_outb)
            df_exp_pi0.loc[:, "Pp"]      = df_exp_pi0.Pp     + dp
            df_exp_pi0.loc[:, "Ptheta"]  = df_exp_pi0.Ptheta + dtheta
            df_exp_pi0.loc[:, "Pphi"]    = df_exp_pi0.Pphi   + dphi
            df_exp_pi0.loc[:, "Ppx"]     = df_exp_pi0.loc[:, "Pp"]*np.sin(np.radians(df_exp_pi0.loc[:, "Ptheta"]))*np.cos(np.radians(df_exp_pi0.loc[:, "Pphi"]))
            df_exp_pi0.loc[:, "Ppy"]     = df_exp_pi0.loc[:, "Pp"]*np.sin(np.radians(df_exp_pi0.loc[:, "Ptheta"]))*np.sin(np.radians(df_exp_pi0.loc[:, "Pphi"]))
            df_exp_pi0.loc[:, "Ppz"]     = df_exp_pi0.loc[:, "Pp"]*np.cos(np.radians(df_exp_pi0.loc[:, "Ptheta"]))
            
            if Psector in [1, 2, 3, 4, 5, 6]:
                df_exp_pi0 = copy(df_exp_pi0                 .loc[ (df_exp_pi0.Pp             >Pmin) & (df_exp_pi0.Pp         < Pmax) &  (df_exp_pi0.Psector         == Psector), :])
            elif Psector == "CD":
                df_exp_pi0 = copy(df_exp_pi0                 .loc[ (df_exp_pi0.Pp             >Pmin) & (df_exp_pi0.Pp         < Pmax) &  (df_exp_pi0.Psector         >4000), :])
            df_exp_pi0                   = saveDVpi0vars(df_exp_pi0)
            exp_dist_pi0_mm2ep  , _      = np.histogram(df_exp_pi0.loc[:, "MM2_ep" ]  , bins = bins_pi0_mm2ep, density = True  )#, weights = df_exp_pi0.weights)
            exp_dist_pi0_reconPi, _      = np.histogram(df_exp_pi0.loc[:, "reconPi"]  , bins = bins_pi0_reconPi, density = True)#, weights = df_exp_pi0.weights)
    
            df_exp_dvcs = copy(df_exp_dvcs_outb)
            df_exp_dvcs.loc[:, "Pp"]      = df_exp_dvcs.Pp     + dp
            df_exp_dvcs.loc[:, "Ptheta"]  = df_exp_dvcs.Ptheta + dtheta
            df_exp_dvcs.loc[:, "Pphi"]    = df_exp_dvcs.Pphi   + dphi
            df_exp_dvcs.loc[:, "Ppx"]     = df_exp_dvcs.loc[:, "Pp"]*np.sin(np.radians(df_exp_dvcs.loc[:, "Ptheta"]))*np.cos(np.radians(df_exp_dvcs.loc[:, "Pphi"]))
            df_exp_dvcs.loc[:, "Ppy"]     = df_exp_dvcs.loc[:, "Pp"]*np.sin(np.radians(df_exp_dvcs.loc[:, "Ptheta"]))*np.sin(np.radians(df_exp_dvcs.loc[:, "Pphi"]))
            df_exp_dvcs.loc[:, "Ppz"]     = df_exp_dvcs.loc[:, "Pp"]*np.cos(np.radians(df_exp_dvcs.loc[:, "Ptheta"]))
            if Psector in [1, 2, 3, 4, 5, 6]:
                df_exp_dvcs = copy(df_exp_dvcs                .loc[ (df_exp_dvcs.Pp       >Pmin) & (df_exp_dvcs.Pp       < Pmax) &  (df_exp_dvcs.Psector       == Psector), :])
            elif Psector == "CD":
                df_exp_dvcs = copy(df_exp_dvcs                .loc[ (df_exp_dvcs.Pp       >Pmin) & (df_exp_dvcs.Pp       < Pmax) &  (df_exp_dvcs.Psector       >4000), :])
    
            df_exp_dvcs                    = saveDVCSvars(df_exp_dvcs)
            exp_dist_epg_mm2ep   , _       = np.histogram(df_exp_dvcs.loc[:, "MM2_ep" ] , bins = bins_dvcs_mm2ep, density = True)#, weights = df_exp_dvcs.weights)
            exp_dist_epg_reconGam, _       = np.histogram(df_exp_dvcs.loc[:, "reconGam"], bins = bins_dvcs_reconGam, density = True)#, weights = df_exp_dvcs.weights)
    
            score = entropy(exp_dist_pi0_mm2ep, sim_dist_pi0_mm2ep) * entropy(exp_dist_pi0_reconPi, sim_dist_pi0_reconPi) * entropy(exp_dist_epg_mm2ep, sim_dist_epg_mm2ep) * entropy(exp_dist_epg_reconGam, sim_dist_epg_reconGam)
            scores  .append(score)
            dps     .append(dp)
            dthetas .append(dtheta)
            # print(score)
            # plt.show()
        dps             = np.array(dps)
        dthetas         = np.array(dthetas)
        optimal_dps     = dps    [np.argsort(scores)[:20]]
        optimal_dthetas = dthetas[np.argsort(scores)[:20]]
        optimal_dp      = np.mean(optimal_dps)
        optmial_dtheta  = np.mean(optimal_dthetas)
        optimal_dp_std     = np.std(optimal_dps)
        optimal_dtheta_std = np.std(optimal_dthetas)
    
        this_row = pd.DataFrame.from_dict({"Psector": [Psector], "p": [Pcenter], "dp": [optimal_dp], "dtheta": [optimal_dtheta], "dp_std": [optimal_dp_std], "dtheta_std": [optimal_dtheta_std], "polarity": ["outbending"]})
    
        df_correction = pd.concat([df_correction, this_row])

df_correction.to_pickle("df_correction_info_inb.pkl")