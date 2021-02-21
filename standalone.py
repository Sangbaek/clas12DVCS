import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from lmfit import Model, Parameter, report_fit
from scipy.optimize import curve_fit
from copy import copy
from utils.physics import *
import icecream as ic


M = 0.938272081
me = 0.5109989461 * 0.001
ebeam = 10.604
pbeam = np.sqrt(ebeam * ebeam - me * me)
beam = [0, 0, pbeam]
target = [0, 0, 0]

file_in = uproot.open(
    "~/Dropbox (MIT)/data/dvcs_inb.root")
tree_in = file_in["T"]
file_out = uproot.open(
    "~/Dropbox (MIT)/data/dvcs_outb.root")
tree_out = file_out["T"]

# initialize df
df_electron_in = pd.DataFrame()
df_proton_in = pd.DataFrame()
df_gamma_in = pd.DataFrame()

df_electron_out = pd.DataFrame()
df_proton_out = pd.DataFrame()
df_gamma_out = pd.DataFrame()

for key in ["Epx", "Epy", "Epz", "Evz", "Esector"]:
    df_electron_in[key] = tree_in[key].array(library="pd")#, entry_stop=10000)

for key in ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "PorigIndx", "Psector"]:
    df_proton_in[key] = tree_in[key].array(library="pd")#, entry_stop=10000)

for key in ["Gpx", "Gpy", "Gpz", "Gstat", "GorigIndx", "Gsector"]:
    df_gamma_in[key] = tree_in[key].array(library="pd")#, entry_stop=10000)

# for key in ["Epx", "Epy", "Epz", "Evz"]:
#     df_electron_out[key] = tree_out[key].array(library="pd")#, entry_stop=100)

# for key in ["Ppx", "Ppy", "Ppz", "Pvz", "Pstat", "PorigIndx"]:
#     df_proton_out[key] = tree_out[key].array(library="pd")#, entry_stop=100)

# for key in ["Gpx", "Gpy", "Gpz", "Gstat", "GorigIndx"]:
#     df_gamma_out[key] = tree_out[key].array(library="pd")#, entry_stop=100)

# df_electron = pd.concat([df_electron_in, df_electron_out])
# df_proton = pd.concat([df_proton_in, df_proton_out])
# df_gamma = pd.concat([df_gamma_in, df_gamma_out])

df_electron = df_electron_in
df_proton = df_proton_in
df_gamma = df_gamma_in

# change the data type to double
df_electron = df_electron.astype({"Epx": float, "Epy": float, "Epz": float})
df_proton = df_proton.astype({"Ppx": float, "Ppy": float, "Ppz": float})
df_gamma = df_gamma.astype({"Gpx": float, "Gpy": float, "Gpz": float})

df_electron['event'] = df_electron.index.get_level_values('entry')
df_proton['event'] = df_proton.index.get_level_values('entry')
df_gamma['event'] = df_gamma.index.get_level_values('entry')

df_gg = pd.merge(df_gamma, df_gamma,
                 how='outer', on='event', suffixes=("", "2"))
df_gg = df_gg[df_gg["GorigIndx"] < df_gg["GorigIndx2"]]

df_ep = pd.merge(df_electron, df_proton, how='outer', on='event')
df_epg = pd.merge(df_ep, df_gamma, how='outer', on='event')
df_epgg = pd.merge(df_ep, df_gg, how='outer', on='event')
df_epgg = df_epgg[~np.isnan(df_epgg["Gpx"])]

df_epg = df_epg[np.abs(df_epg["Evz"] - df_epg["Pvz"]) < 2.5 +
                2.5 / mag([df_epg["Epx"], df_epg["Epy"], df_epg["Epz"]])]
df_epgg = df_epgg[np.abs(df_epgg["Evz"] - df_epgg["Pvz"]) < 2.5 +
                  2.5 / mag([df_epgg["Epx"], df_epgg["Epy"], df_epgg["Epz"]])]

# dvpip
# basic particle energies
df_epgg['Ee'] = np.sqrt(me**2 + df_epgg["Epx"]**2 + df_epgg["Epy"]**2 + df_epgg["Epz"]**2)
df_epgg['Pe'] = np.sqrt(M**2 + df_epgg["Ppx"]**2 + df_epgg["Ppy"]**2 + df_epgg["Ppz"]**2)
df_epgg['Ge'] = np.sqrt(df_epgg["Gpx"]**2 + df_epgg["Gpy"]**2 + df_epgg["Gpz"]**2)
df_epgg['Ge2'] = np.sqrt(df_epgg["Gpx2"]**2 + df_epgg["Gpy2"]**2 + df_epgg["Gpz2"]**2)

# useful objects
ele = [df_epgg['Epx'], df_epgg['Epy'], df_epgg['Epz']]
pro = [df_epgg['Ppx'], df_epgg['Ppy'], df_epgg['Ppz']]
gam = [df_epgg['Gpx'], df_epgg['Gpy'], df_epgg['Gpz']]
gam2 = [df_epgg['Gpx2'], df_epgg['Gpy2'], df_epgg['Gpz2']]
pi0 = vecAdd(gam, gam2)
VGS = [-df_epgg['Epx'], -df_epgg['Epy'], pbeam - df_epgg['Epz']]
v3l = cross(beam, ele)
v3h = cross(pro, VGS)
v3g = cross(VGS, gam)
VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
            df_epgg["Ppy"], pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]

# binning kinematics
df_epgg['Q2'] = -((ebeam - df_epgg['Ee'])**2 - mag2(VGS))
df_epgg['nu'] = (ebeam - df_epgg['Ee'])
df_epgg['xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
df_epgg['W'] = np.sqrt((ebeam + M - df_epgg['Ee'])**2 - mag2(VGS))
df_epgg['MPt'] = np.sqrt((df_epgg["Epx"] + df_epgg["Ppx"] + df_epgg["Gpx"] + df_epgg["Gpx2"])**2 +
                         (df_epgg["Epy"] + df_epgg["Ppy"] + df_epgg["Gpy"] + df_epgg["Gpy2"])**2)

# exclusivity variables
df_epgg['MM2_ep'] = (-M - ebeam + df_epgg["Ee"] +
                     df_epgg["Pe"])**2 - mag2(VmissPi0)
df_epgg['ME_epgg'] = (M + ebeam - df_epgg["Ee"] - df_epgg["Pe"] - df_epgg["Ge"] - df_epgg["Ge2"])
df_epgg['Mpi0'] = pi0InvMass(gam, gam2)
df_epgg['reconPi'] = angle(VmissPi0, pi0)
df_epgg["Pie"] = df_epgg['Ge'] + df_epgg['Ge2']

# Kinematic cuts
#   xB>0 && xB<1 && W>2 && Q2>1 && VE.e()>2 && VG1.e() > 2 && VP.p() > 0.12
cut_pFD = df_epgg["Pstat"] < 4000  # FD
cut_pCD = df_epgg["Pstat"] > 4000  # CD
cut_gFD = df_epgg["Gstat"] > 2000  # FD
cut_gFT = df_epgg["Gstat"] < 2000  # FT
cut_g2FD = df_epgg["Gstat2"] > 2000  # FD
cut_g2FT = df_epgg["Gstat2"] < 2000  # FT
cut_FTFT = cut_gFT & cut_g2FT
cut_FDFT = (cut_gFD & cut_g2FD) | (cut_gFT & cut_g2FT)
cut_FDFD = cut_gFD & cut_g2FD
cut_xBupper = df_epgg["xB"] < 1  # xB
cut_xBlower = df_epgg["xB"] > 0  # xB
cut_Q2 = df_epgg["Q2"] > 1  # Q2
cut_W = df_epgg["W"] > 2  # W

# Exclusivity cuts
cut_mmep = df_epgg["MM2_ep"] < 0.7  # mmep
cut_meepgg = df_epgg["ME_epgg"] < 0.7  # meepgg
cut_mpt = df_epgg["MPt"] < 0.2  # mpt
cut_recon = df_epgg["reconPi"] < 2  # recon gam angle
cut_pi0upper = df_epgg["Mpi0"] < 0.2
cut_pi0lower = df_epgg["Mpi0"] > 0.07
# cut_pi0energy = df_epgg["Pie"] > 3
cut_maxE = 1 #df_epgg[["Ge", "Ge2"]].max(1) > 2
cut_minE = 1 #df_epgg[["Ge", "Ge2"]].min(1) > 0.8

mpi0_edges = np.linspace(0.5, 1.5, 101)
df_dvpi0 = df_epgg[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_mmep & cut_meepgg &
                   cut_mpt & cut_recon & cut_pi0upper & cut_pi0lower]

# dvcs
# basic particle energies
df_epg['Ee'] = np.sqrt(me**2 + df_epg["Epx"]**2 +
                       df_epg["Epy"]**2 + df_epg["Epz"]**2)
df_epg['Pe'] = np.sqrt(M**2 + df_epg["Ppx"]**2 +
                       df_epg["Ppy"]**2 + df_epg["Ppz"]**2)
df_epg['Ge'] = np.sqrt(df_epg["Gpx"]**2 + df_epg["Gpy"]**2 + df_epg["Gpz"]**2)

# useful objects
ele = [df_epg['Epx'], df_epg['Epy'], df_epg['Epz']]
pro = [df_epg['Ppx'], df_epg['Ppy'], df_epg['Ppz']]
gam = [df_epg['Gpx'], df_epg['Gpy'], df_epg['Gpz']]
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

# binning kinematics
df_epg['Q2'] = -((ebeam - df_epg['Ee'])**2 - mag2(VGS))
df_epg['nu'] = (ebeam - df_epg['Ee'])
df_epg['xB'] = df_epg['Q2'] / 2.0 / M / df_epg['nu']
df_epg['t1'] = 2 * M * (df_epg['Pe'] - M)
df_epg['t2'] = (M * df_epg['Q2'] + 2 * M * df_epg['nu'] * (df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta))\
/ (M + df_epg['nu'] - np.sqrt(df_epg['nu'] * df_epg['nu'] + df_epg['Q2']) * costheta)
df_epg['W'] = np.sqrt((ebeam + M - df_epg['Ee'])**2 - mag2(VGS))

# trento angles
df_epg['phi1'] = angle(v3l, v3h)
df_epg['phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                          df_epg['phi1'], df_epg['phi1'])
df_epg['phi2'] = angle(v3l, v3g)
df_epg['phi2'] = np.where(dot(VGS, cross(v3l, v3g)) <
                          0, 360.0 - df_epg['phi2'], df_epg['phi2'])

# exclusivity variables
df_epg['MM2_epg'] = (-M - ebeam + df_epg["Ee"] +
                     df_epg["Pe"] + df_epg["Ge"])**2 - mag2(Vmiss)
df_epg['ME_epg'] = (M + ebeam - df_epg["Ee"] - df_epg["Pe"] - df_epg["Ge"])
df_epg['MM2_ep'] = (-M - ebeam + df_epg["Ee"] + df_epg["Pe"])**2 - mag2(VmissG)
df_epg['MM2_eg'] = (-M - ebeam + df_epg["Ee"] + df_epg["Ge"])**2 - mag2(VmissP)
df_epg = df_epg[df_epg["MM2_eg"] > 0]  # mmeg
df_epg['MPt'] = np.sqrt((df_epg["Epx"] + df_epg["Ppx"] + df_epg["Gpx"])**2 +
                        (df_epg["Epy"] + df_epg["Ppy"] + df_epg["Gpy"])**2)
# 180.*/np.pi*np.arccos((df_epg["Epx"]*df_epg["Gpx"]+df_epg["Epy"]*df_epg["Gpy"]+df_epg["Epz"]*df_epg["Gpz"])/np.sqrt(df_epg["Epx"]*df_epg["Epx"]+df_epg["Epy"]*df_epg["Epy"]+df_epg["Epz"]*df_epg["Epz"])/np.sqrt(df_epg["Gpx"]*df_epg["Gpx"]+df_epg["Gpy"]*df_epg["Gpy"]+df_epg["Gpz"]*df_epg["Gpz"]))
df_epg['coneAngle'] = angle(ele, gam)
df_epg['reconGam'] = angle(gam, VmissG)
df_epg['coplanarity'] = angle(v3h, v3g)

#	Kinematic cuts
#   xB>0 && xB<1 && W>2 && Q2>1 && VE.e()>2 && VG1.e() > 2 && VP.p() > 0.12
cut_forward = df_epg["Pstat"] < 4000  # FD
cut_central = df_epg["Pstat"] > 4000  # CD
cut_xBupper = df_epg["xB"] < 1  # xB
cut_xBlower = df_epg["xB"] > 0  # xB
cut_Q2 = df_epg["Q2"] > 1  # Q2
cut_W = df_epg["W"] > 2  # W
cut_Ee = df_epg["Ee"] > 2  # Ee
cut_Ge = df_epg["Ge"] > 2  # Ge
cut_Pp = mag([df_epg["Ppx"], df_epg["Ppy"], df_epg["Ppz"]]) > 0.12  # Pp

#	Exclusivity cuts
cut_mmepg = np.abs(df_epg["MM2_epg"]) < 0.04  # mmepg
cut_mmegupper = np.sqrt(df_epg["MM2_eg"]) < 1.7  # mmeg
cut_mmeglower = np.sqrt(df_epg["MM2_eg"]) > 0.1  # mmeg
cut_meepgupper = df_epg["ME_epg"] < 1.2  # meepg
cut_meepglower = df_epg["ME_epg"] > -0.5  # meepg
cut_mpt = df_epg["MPt"] < 0.12  # mpt
cut_cone = df_epg["coneAngle"] > 10  # coneangle
cut_recon = df_epg["reconGam"] < 1.1  # recon gam angle
# cut on coplanarity
# cut on mmep

df_dvcs = df_epg[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Pp & cut_mmepg &
                 cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon]

print(df_dvcs)
exit()

df_dvcsForward = df_epg[cut_forward & cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Pp &
                        cut_mmepg & cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon]
df_dvcsCentral = df_epg[cut_central & cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_Ee & cut_Ge & cut_Pp &
                        cut_mmepg & cut_mmegupper & cut_mmeglower & cut_meepgupper & cut_meepglower & cut_mpt & cut_cone & cut_recon]
# dvcs end

pi0to2gammas = df_dvcs["event"].isin(df_dvpi0["event"])
df_dvcs = df_dvcs[~pi0to2gammas]

xB_edges = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,\
			0.45, 0.5, 0.55, 0.6, 0.7, 0.85, 1]
Q2_edges = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,\
			4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 12.]
t_edges = [0.09, 0.15, 0.20, 0.3, 0.4, 0.60, 1.00, 1.5, 2.0]
phi_edges = np.linspace(0, 360, 31)

fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
for xBind in range(0,len(xB_edges)-1):
	for Q2ind in range(0, len(Q2_edges)-1):
		binQ2lower = df_dvcsCentral["Q2"]>Q2_edges[Q2ind]
		binQ2upper = df_dvcsCentral["Q2"]<Q2_edges[Q2ind+1]
		binxBlower = df_dvcsCentral["xB"]>xB_edges[xBind]
		binxBupper = df_dvcsCentral["xB"]<xB_edges[xBind+1]
		plotter = df_dvcsCentral[binQ2lower&binQ2upper&binxBlower&binxBupper]
		axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, histtype='stepfilled', facecolor='none', edgecolor='k')
		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
		axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,10000])
		axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
plt.subplots_adjust(wspace=0, hspace=0)
plt.tick_params(axis='both', labelsize=0, length = 0)
plt.savefig("CD_phi.pdf")
plt.clf()

fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
for xBind in range(0,len(xB_edges)-1):
	for Q2ind in range(0, len(Q2_edges)-1):
		binQ2lower = df_dvcsForward["Q2"]>Q2_edges[Q2ind]
		binQ2upper = df_dvcsForward["Q2"]<Q2_edges[Q2ind+1]
		binxBlower = df_dvcsForward["xB"]>xB_edges[xBind]
		binxBupper = df_dvcsForward["xB"]<xB_edges[xBind+1]
		plotter = df_dvcsForward[binQ2lower&binQ2upper&binxBlower&binxBupper]
		axs[len(Q2_edges)-2-Q2ind,xBind].hist(plotter["phi2"],bins=phi_edges, histtype='stepfilled', facecolor='none', edgecolor='k')
		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
		axs[len(Q2_edges)-2-Q2ind,xBind].set_ylim([1,10000])
		axs[len(Q2_edges)-2-Q2ind,xBind].set_yscale('log')
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("FD_phi.pdf")
plt.clf()

cmap = copy(plt.cm.get_cmap("jet"))

fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
for xBind in range(0,len(xB_edges)-1):
	for Q2ind in range(0, len(Q2_edges)-1):
		binQ2lower = df_dvcsCentral["Q2"]>Q2_edges[Q2ind]
		binQ2upper = df_dvcsCentral["Q2"]<Q2_edges[Q2ind+1]
		binxBlower = df_dvcsCentral["xB"]>xB_edges[xBind]
		binxBupper = df_dvcsCentral["xB"]<xB_edges[xBind+1]
		plotter = df_dvcsCentral[binQ2lower&binQ2upper&binxBlower&binxBupper]
		axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=cmap, cmin=1)
		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("CD_tphi.pdf")
plt.clf()

fig, axs = plt.subplots(len(Q2_edges)-1, len(xB_edges)-1)
for xBind in range(0,len(xB_edges)-1):
	for Q2ind in range(0, len(Q2_edges)-1):
		binQ2lower = df_dvcsForward["Q2"]>Q2_edges[Q2ind]
		binQ2upper = df_dvcsForward["Q2"]<Q2_edges[Q2ind+1]
		binxBlower = df_dvcsForward["xB"]>xB_edges[xBind]
		binxBupper = df_dvcsForward["xB"]<xB_edges[xBind+1]
		plotter = df_dvcsForward[binQ2lower&binQ2upper&binxBlower&binxBupper]
		axs[len(Q2_edges)-2-Q2ind,xBind].hist2d(plotter["phi2"],plotter["t2"],bins=[phi_edges, t_edges], cmap=cmap, cmin=1)
		axs[len(Q2_edges)-2-Q2ind,xBind].set_title("")
		axs[len(Q2_edges)-2-Q2ind,xBind].xaxis.set_visible(False)
		axs[len(Q2_edges)-2-Q2ind,xBind].yaxis.set_visible(False)
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_xticklabels(), visible=False)
		plt.setp( axs[len(Q2_edges)-2-Q2ind,xBind].get_yticklabels(), visible=False)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("FD_tphi.pdf")
plt.clf()


def dsigma(phi, A, B, C, D, E):
    """
    rudimentary fitting function, should be in 7 parameters.
    """
    return A + B * np.cos(phi * np.pi / 180.) + C * np.sin(phi * np.pi / 180.) + \
        D * np.cos(2 * phi * np.pi / 180.) + E * np.cos(2 * phi * np.pi / 180.)#/ (1 + F * np.cos(phi * np.pi / 180.)) / (1 + G * np.cos(phi * np.pi / 180.))


for tind in range(0, len(t_edges) - 1):
    fig, axs = plt.subplots(len(Q2_edges) - 1, len(xB_edges) - 1)
    for xBind in range(0, len(xB_edges) - 1):
        for Q2ind in range(0, len(Q2_edges) - 1):
            binQ2lower = df_dvcsCentral["Q2"] > Q2_edges[Q2ind]
            binQ2upper = df_dvcsCentral["Q2"] < Q2_edges[Q2ind + 1]
            binxBlower = df_dvcsCentral["xB"] > xB_edges[xBind]
            binxBupper = df_dvcsCentral["xB"] < xB_edges[xBind + 1]
            bintlower = df_dvcsCentral["t2"] > t_edges[tind]
            bintupper = df_dvcsCentral["t2"] < t_edges[tind + 1]
            plotter = df_dvcsCentral[binQ2lower & binQ2upper &
                                     binxBlower & binxBupper & bintlower & bintupper]
            axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
            plt.setp(axs[len(Q2_edges) - 2 - Q2ind,
                         xBind].get_xticklabels(), visible=False)
            plt.setp(axs[len(Q2_edges) - 2 - Q2ind,
                         xBind].get_yticklabels(), visible=False)
            if (len(plotter)==0):
              continue
            axs[len(Q2_edges) - 2 - Q2ind,
                xBind].hist(plotter["phi2"],
                            bins=phi_edges,
                            linewidth=0.5,
                            histtype='stepfilled',
                            facecolor='none',
                            edgecolor='k')
            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_title("")
            axs[len(Q2_edges) - 2 - Q2ind, xBind].xaxis.set_visible(False)
            axs[len(Q2_edges) - 2 - Q2ind, xBind].yaxis.set_visible(False)

            data_entries, bins = np.histogram(plotter["phi2"], bins=phi_edges)
            binscenters = np.array([0.5 * (bins[i] + bins[i + 1])
                                    for i in range(len(bins) - 1)])
            popt, pcov = curve_fit(
                dsigma, xdata=binscenters, ydata=data_entries, p0=[
                    100.0, 1.0, 1.0, 1.0, 1.0])
            print(popt)
            phispace = np.arange(0, 360.1, 1)
            axs[len(Q2_edges) - 2 - Q2ind, xBind].plot(phispace, dsigma(phispace,
                                                                        *popt), color='red', linewidth=1.0, label=r'Fitted function')

            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_ylim([1, 5000])
            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_yscale('log')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tick_params(axis='both', labelsize=0, length=0)
    plt.savefig("central_{}.pdf".format(tind))
    plt.clf()

for tind in range(0, len(t_edges) - 1):
    fig, axs = plt.subplots(len(Q2_edges) - 1, len(xB_edges) - 1)
    for xBind in range(0, len(xB_edges) - 1):
        for Q2ind in range(0, len(Q2_edges) - 1):
            binQ2lower = df_dvcsForward["Q2"] > Q2_edges[Q2ind]
            binQ2upper = df_dvcsForward["Q2"] < Q2_edges[Q2ind + 1]
            binxBlower = df_dvcsForward["xB"] > xB_edges[xBind]
            binxBupper = df_dvcsForward["xB"] < xB_edges[xBind + 1]
            bintlower = df_dvcsForward["t2"] > t_edges[tind]
            bintupper = df_dvcsForward["t2"] < t_edges[tind + 1]
            plotter = df_dvcsForward[binQ2lower & binQ2upper &
                                     binxBlower & binxBupper & bintlower & bintupper]
            axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
            plt.setp(axs[len(Q2_edges) - 2 - Q2ind,
                         xBind].get_xticklabels(), visible=False)
            plt.setp(axs[len(Q2_edges) - 2 - Q2ind,
                         xBind].get_yticklabels(), visible=False)
            if (len(plotter)==0):
              continue
            axs[len(Q2_edges) - 2 - Q2ind,
                xBind].hist(plotter["phi2"],
                            bins=phi_edges,
                            linewidth=0.5,
                            histtype='stepfilled',
                            facecolor='none',
                            edgecolor='k')
            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_title("")
            axs[len(Q2_edges) - 2 - Q2ind, xBind].xaxis.set_visible(False)
            axs[len(Q2_edges) - 2 - Q2ind, xBind].yaxis.set_visible(False)

            data_entries, bins = np.histogram(plotter["phi2"], bins=phi_edges)
            binscenters = np.array([0.5 * (bins[i] + bins[i + 1])
                                    for i in range(len(bins) - 1)])
            popt, pcov = curve_fit(
                dsigma, xdata=binscenters, ydata=data_entries, p0=[
                    100.0, 1.0, 1.0, 1.0, 1.0])
            print(popt)
            phispace = np.arange(0, 360.1, 1)
            axs[len(Q2_edges) - 2 - Q2ind, xBind].plot(phispace, dsigma(phispace,
                                                                        *popt), color='red', linewidth=1.0, label=r'Fitted function')

            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_ylim([1, 5000])
            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_yscale('log')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tick_params(axis='both', labelsize=0, length=0)
    plt.savefig("forward_{}.pdf".format(tind))
    plt.clf()

for tind in range(0, len(t_edges) - 1):
    fig, axs = plt.subplots(len(Q2_edges) - 1, len(xB_edges) - 1)
    for xBind in range(0, len(xB_edges) - 1):
        for Q2ind in range(0, len(Q2_edges) - 1):
            binQ2lower = df_dvcs["Q2"] > Q2_edges[Q2ind]
            binQ2upper = df_dvcs["Q2"] < Q2_edges[Q2ind + 1]
            binxBlower = df_dvcs["xB"] > xB_edges[xBind]
            binxBupper = df_dvcs["xB"] < xB_edges[xBind + 1]
            bintlower = df_dvcs["t2"] > t_edges[tind]
            bintupper = df_dvcs["t2"] < t_edges[tind + 1]
            plotter = df_dvcs[binQ2lower & binQ2upper &
                              binxBlower & binxBupper & bintlower & bintupper]
            axs[len(Q2_edges)-2-Q2ind,xBind].set_xticks([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_yticks([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_xticklabels([])
            axs[len(Q2_edges)-2-Q2ind,xBind].set_yticklabels([])
            plt.setp(axs[len(Q2_edges) - 2 - Q2ind,
                         xBind].get_xticklabels(), visible=False)
            plt.setp(axs[len(Q2_edges) - 2 - Q2ind,
                         xBind].get_yticklabels(), visible=False)
            if (len(plotter)==0):
              continue
            axs[len(Q2_edges) - 2 - Q2ind,
                xBind].hist(plotter["phi2"],
                            bins=phi_edges,
                            linewidth=0.5,
                            histtype='stepfilled',
                            facecolor='none',
                            edgecolor='k')
            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_title("")
            axs[len(Q2_edges) - 2 - Q2ind, xBind].xaxis.set_visible(False)
            axs[len(Q2_edges) - 2 - Q2ind, xBind].yaxis.set_visible(False)

            data_entries, bins = np.histogram(plotter["phi2"], bins=phi_edges)
            binscenters = np.array([0.5 * (bins[i] + bins[i + 1])
                                    for i in range(len(bins) - 1)])
            popt, pcov = curve_fit(
                dsigma, xdata=binscenters, ydata=data_entries, p0=[
                    100.0, 1.0, 1.0, 1.0, 1.0])
            print(popt)
            phispace = np.arange(0, 360.1, 1)
            axs[len(Q2_edges) - 2 - Q2ind, xBind].plot(phispace, dsigma(phispace,
                                                                        *popt), color='red', linewidth=1.0, label=r'Fitted function')

            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_ylim([1, 5000])
            axs[len(Q2_edges) - 2 - Q2ind, xBind].set_yscale('log')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tick_params(axis='both', labelsize=0, length=0)
    plt.savefig("dvcs_{}.pdf".format(tind))
    plt.clf()