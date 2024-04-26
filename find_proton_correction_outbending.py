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

outb_list = np.loadtxt('/work/clas12/sangbaek/outb_run_list_pass1', dtype = int)
df_exp_dvcs_outb = []
df_exp_pi0_outb = []
for run in outb_list:
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_1/pkl/{:d}.pkl".format(run))
    df_exp_dvcs_outb.append(df)
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/pi0/excl_level_1/pkl/{:d}.pkl".format(run))
    df_exp_pi0_outb.append(df)
df_exp_dvcs_outb = pd.concat(df_exp_dvcs_outb).reset_index()
df_exp_dvcs_outb = df_exp_dvcs_outb.loc[:, df_exp_dvcs_outb.columns[1:]]
df_exp_pi0_outb = pd.concat(df_exp_pi0_outb).reset_index()
df_exp_pi0_outb = df_exp_pi0_outb.loc[:, df_exp_pi0_outb.columns[1:]]

df_sim_dvcs_outb = []
for bin in range(147):
    try:
        df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl/{}.pkl".format(bin+1))
        df.loc[:, "weights"] = df.GenWeight * bin_volume_bulk[bin] * luminosity_outb / (680000*9999/10000)
        df_sim_dvcs_outb.append(df)
    except:
        continue
df_sim_dvcs_outb = pd.concat(df_sim_dvcs_outb).reset_index()
df_sim_dvcs_outb = df_sim_dvcs_outb.loc[:, df_sim_dvcs_outb.columns[1:]]

df_sim_pi0_1gamma_outb = []
for filenum in range(20):
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/pkl/1/{}.pkl".format(filenum+1))
    df_sim_pi0_1gamma_outb.append(df)
df_sim_pi0_1gamma_outb = pd.concat(df_sim_pi0_1gamma_outb).reset_index()
df_sim_pi0_1gamma_outb = df_sim_pi0_1gamma_outb.loc[:, df_sim_pi0_1gamma_outb.columns[1:]]


df_sim_pi0_2gamma_outb = []
for filenum in range(20):
    df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/pkl/1/{}.pkl".format(filenum+1))
    df_sim_pi0_2gamma_outb.append(df)
df_sim_pi0_2gamma_outb = pd.concat(df_sim_pi0_2gamma_outb).reset_index()
df_sim_pi0_2gamma_outb = df_sim_pi0_2gamma_outb.loc[:, df_sim_pi0_2gamma_outb.columns[1:]]

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

df_correction.to_pickle("df_correction_info_outb.pkl")