#!/usr/bin/env python3
"""
Script to reproduce some of plots used in Sangbaek's thesis.
The original plots were produced in the jupyter notebook.
Some of them had problems in that
(1) they are not clearly visible,
(2) they do not have color bar scale for the 2d histogram.
"""
import gc
import matplotlib.pyplot as plt
from copy import copy
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from utils.const import *
from utils.physics import *
from utils.fiducial import *
from matplotlib.colors import LogNorm
import argparse
from glob import glob
import itertools

degree = r"${}^{\circ}$"
GeV = "GeV"
GeV2 = "GeV"+r"${}^{2}$"
GeVc = "GeV/c"
GeVc2 = "(GeV/c)"+r"${}^{2}$"

import matplotlib
# initial settings
# cmap = copy(matplotlib.colormaps["jet"])
# cmap.set_under('w',0)
# cmap.set_bad('w',0)
pgf_with_latex = {
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,			# use LaTeX to write all text
    "font.family": "sans-serif",		
    "font.sans-serif": "Helvetica",
    "font.size": 25,				# default font size
    "axes.titlepad": 20,			# x and y label size
    "axes.labelsize": 24,			# x and y label size
    "axes.titlesize": 24,		  # subfigure title size, i.e. title size when one figure
    "legend.fontsize": 22,			# legend size
    "xtick.labelsize": 23,			# x axis tick label size
    "ytick.labelsize": 23,			# y axis tick label 
    "figure.titlesize": 25,         # Figure title size, useful when you have multiple plots in one canvas.
    "pgf.preamble": r"\usepackage{xcolor}",     # xcolor for colours
    "figure.autolayout": False
}
matplotlib.rcParams.update(pgf_with_latex)
from matplotlib.colors import LinearSegmentedColormap
cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
 [0.0589714286, 0.6837571429, 0.7253857143], 
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
 [0.7184095238, 0.7411333333, 0.3904761905], 
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
 [0.9763, 0.9831, 0.0538]]

parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
cmap = parula_map

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

'''
ch2
1. exp sf precut (inb)
2. mc  sf precut (inb)
3. exp e' dc pre + post (inb)
4. exp e' dc pre + post (outb)
5. exp e' pcal pre + post (inb)
6. exp e' pcal pre + post (outb)
7. exp SF pre (inb) V, W
8. exp edep vs edep (pre)
9. exp E/p vs E/p (inb) (pre)
10. exp chi vs p (inb) (pre)
11. exp p' dc pre + post (inb)
13. exp p' dc pre + post (outb)
14. exp gamma calorimeter pre + post (inb + outb)
ch5
15. exp proton CVT vs theta, phi(pre)
16. exp + mc theta distribution (pre)
17. exp + mc CVT theta distribution (pre)
18. exp + mc phi cvt distribution (pre)
19. exp photon FT (pre + post)
'''
'''
dvcsSample = []
for bin in range(1, 147+1):
  print(bin)
  df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_1_nocorr_nofid/{}.pkl".format(bin))
  df = df.sample(len(df)//10)
  dvcsSample.append(df)
dvcsSample = pd.concat(dvcsSample)
dvcsSample = dvcsSample.reset_index()
dvcsSample = dvcsSample.loc[:, dvcsSample.columns[1:]]

expSample = []
for bin in range(1, 147+1):
  print(bin)
  df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_inb/dvcs/excl_level_2/restructured_1_nocorr_nofid/{}.pkl".format(bin))
  # df = df.sample(len(df)//10)
  expSample.append(df)
expSample = pd.concat(expSample)
expSample = expSample.reset_index()
expSample = expSample.loc[:, expSample.columns[1:]]

expSampleOutb = []
for bin in range(1, 147+1):
  print(bin)
  df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_outb/dvcs/excl_level_2/restructured_1_nocorr_nofid/{}.pkl".format(bin))
  # df = df.sample(len(df)//10)
  expSampleOutb.append(df)
expSampleOutb = pd.concat(expSampleOutb)
expSampleOutb = expSampleOutb.reset_index()
expSampleOutb = expSampleOutb.loc[:, expSampleOutb.columns[1:]]

#1. exp sf precut (inb)
A = [0.286, 0.280, 0.275, 0.273, 0.271, 0.276]
B = [-0.040, -0.038, -0.034, -0.033, -0.032, -0.034]
C = [-0.0030, -0.0012, -0.0014, -0.0007, 0.0005, -0.0014]
D = [0.017, 0.019, 0.017, 0.0157, 0.016, 0.017]
E = [-0.0012, -0.003, -0.002, 0.0003, -0.00135, -0.002]
F = [-0.0012, -0.00135, -0.00129, -0.0013, -0.001, -0.001]

pcal_sf_mu    = [A, B, C]
pcal_sf_sigma = [D, E, F]

fig, axs = plt.subplots(2, 3, figsize = (15,10))
dummy = np.linspace(0.2, 2.5, 101)
for xind in range(0, 2):
  for yind in range(0, 3):
    Esector = 3*(xind) + yind + 1
    h = axs[xind, yind].hist2d(expSample.loc[expSample.Esector==Esector, "Eedep"],  expSample.loc[expSample.Esector==Esector, "ESamplFrac"], cmin =1, bins = [np.linspace(0, 2.5, 100), np.linspace(0, 0.35, 100)], cmap = parula_map, rasterized = True, norm = LogNorm())
    ticklabels = [one, ten, hundred, thousand]
    ticks = [1, 10, 100, 1000]
    cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
    # cbar.ax.set_yticklabels(ticklabels)

    ecal_e_sampl_mu_0 = pcal_sf_mu[0][Esector - 1]
    ecal_e_sampl_mu_1 = pcal_sf_mu[1][Esector - 1]
    ecal_e_sampl_mu_2 = pcal_sf_mu[2][Esector - 1]
    ecal_e_sampl_sigm_0 = pcal_sf_sigma[0][Esector - 1]
    ecal_e_sampl_sigm_1 = pcal_sf_sigma[1][Esector - 1]
    ecal_e_sampl_sigm_2 = pcal_sf_sigma[2][Esector - 1]
    mean =  ecal_e_sampl_mu_0   + ecal_e_sampl_mu_1  /dummy + ecal_e_sampl_mu_2  /dummy/dummy
    sigma = ecal_e_sampl_sigm_0 + ecal_e_sampl_sigm_1/dummy + ecal_e_sampl_sigm_2/dummy/dummy
    # axs[xind, yind].plot(dummy, mean+3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
    axs[xind, yind].plot(dummy, mean-3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
    axs[xind, yind].set_title(r"$e'$" + " Samp. Frac. Sector {}".format(Esector))
    axs[xind, yind].set_xlim([0, 2.5])
    axs[xind, yind].set_ylim([0, 0.37])
    #         axs[xind, yind].set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
    #         axs[xind, yind].set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
    axs[xind, yind].set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
    axs[xind, yind].set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
    axs[xind, yind].set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
    axs[xind, yind].set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
    axs[xind, yind].set_xlabel(r"$E_{dep.}$" + " ["+GeV+"]")
    axs[xind, yind].set_ylabel(r"$E_{dep.}/p_{e'}$")
plt.tight_layout()
plt.savefig("plots/ch2/precutSampling.pdf")
plt.clf()


#2. mc  sf precut (inb)
A_sim = [0.29]*6
B_sim = [-0.040]*6
C_sim = [-0.0029]*6
D_sim = [0.015]*6
E_sim = [-0.00053]*6
F_sim = [-0.0014]*6


pcal_sf_mu    = [A_sim, B_sim, C_sim]
pcal_sf_sigma = [D_sim, E_sim, F_sim]

fig, axs = plt.subplots(2, 3, figsize = (15,10))
dummy = np.linspace(0.2, 2.5, 101)
for xind in range(0, 2):
  for yind in range(0, 3):
    Esector = 3*xind + yind  +1
    h = axs[xind, yind].hist2d(dvcsSample.loc[dvcsSample.Esector==Esector, "Eedep"],  dvcsSample.loc[dvcsSample.Esector==Esector, "ESamplFrac"], cmin =1, bins = [np.linspace(0, 2.5, 100), np.linspace(0, 0.35, 100)], cmap = parula_map, rasterized = True, norm = LogNorm(), weights = dvcsSample.loc[dvcsSample.Esector==Esector, "weights"])
    ticklabels = [one, ten, hundred, thousand, tenthousands]
    ticks = [1, 10, 100, 1000, 10000]
    cbar = plt.colorbar(h[3], ax = axs[xind, yind], ticks = ticks)
    # cbar.ax.set_yticklabels(ticklabels)
    ecal_e_sampl_mu_0 = pcal_sf_mu[0][Esector - 1]
    ecal_e_sampl_mu_1 = pcal_sf_mu[1][Esector - 1]
    ecal_e_sampl_mu_2 = pcal_sf_mu[2][Esector - 1]
    ecal_e_sampl_sigm_0 = pcal_sf_sigma[0][Esector - 1]
    ecal_e_sampl_sigm_1 = pcal_sf_sigma[1][Esector - 1]
    ecal_e_sampl_sigm_2 = pcal_sf_sigma[2][Esector - 1]
    mean =  ecal_e_sampl_mu_0   + ecal_e_sampl_mu_1  /dummy + ecal_e_sampl_mu_2  /dummy/dummy
    sigma = ecal_e_sampl_sigm_0 + ecal_e_sampl_sigm_1/dummy + ecal_e_sampl_sigm_2/dummy/dummy
    # axs[xind, yind].plot(dummy, mean+3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
    axs[xind, yind].plot(dummy, mean-3.5*sigma, color = 'k', linestyle = '--', linewidth = 5)
    axs[xind, yind].set_title(r"$e'$" + " Samp. Frac. Sector {}".format(Esector))
    axs[xind, yind].set_xlim([0, 2.5])
    axs[xind, yind].set_ylim([0, 0.37])
    #         axs[xind, yind].set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
    #         axs[xind, yind].set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
    axs[xind, yind].set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
    axs[xind, yind].set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
    axs[xind, yind].set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
    axs[xind, yind].set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
    axs[xind, yind].set_xlabel(r"$E_{dep.}$" + " ["+GeV+"]")
    axs[xind, yind].set_ylabel(r"$E_{dep.}/p_{e'}$")
plt.tight_layout()
plt.savefig("plots/ch2/precutSamplingMC.pdf")
plt.clf()

# 3. exp e' dc pre + post (inb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample.EDc3Hitx, expSample.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, thousand, tenthousands]
ticks = [1, 10, 100, 1000, 10000]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-300, 300])
axs.set_ylim([-300, 300])
axs.set_xticks([-300, -150, 0, 150, 300])
axs.set_yticks([-300, -150, 0, 150, 300])
axs.set_title("(a) "+r"$e'$"+" DC Outmost Layer Hits, Pre-fiducial (Inb.)")
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidDC.pdf")
plt.clf()

expSample_post = electronFiducial(expSample)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample_post.EDc3Hitx, expSample_post.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, thousand, tenthousands]
ticks = [1, 10, 100, 1000, 10000]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-300, 300])
axs.set_ylim([-300, 300])
axs.set_xticks([-300, -150, 0, 150, 300])
axs.set_yticks([-300, -150, 0, 150, 300])
axs.set_title("(b) " + r"$e'$"+" DC Outmost Layer Hits, Post-fiducial (Inb.)")
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/postcut_efidDC.pdf")
plt.clf()
# 4. exp e' dc pre + post (outb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSampleOutb.EDc3Hitx, expSampleOutb.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-300, 300])
axs.set_ylim([-300, 300])
axs.set_xticks([-300, -150, 0, 150, 300])
axs.set_yticks([-300, -150, 0, 150, 300])
axs.set_title("(a) " + r"$e'$"+" DC Outmost Layer Hits, Pre-fiducial (Outb.)")
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidDCOutb.pdf")
plt.clf()

expSampleOutb_post = electronFiducial(expSampleOutb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSampleOutb_post.EDc3Hitx, expSampleOutb_post.EDc3Hity, bins = np.linspace(-300, 300, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-300, 300])
axs.set_ylim([-300, 300])
axs.set_xticks([-300, -150, 0, 150, 300])
axs.set_yticks([-300, -150, 0, 150, 300])
axs.set_title("(b) " + r"$e'$"+" DC Outmost Layer Hits, Post-fiducial (Outb.)")
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/postcut_efidDCOutb.pdf")
plt.clf()


# 5. exp e' pcal pre + post (inb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample.EcalHx1, expSample.EcalHy1, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, thousand, tenthousands]
ticks = [1, 10, 100, 1000, 10000]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-400, 400])
axs.set_ylim([-400, 400])
axs.set_xticks([-400, -200, 0, 200, 400])
axs.set_yticks([-400, -200, 0, 200, 400])
axs.set_title("(a) "+r"$e'$"+" PCAL Hits, Pre-fiducial (Inb.)")
axs.set_xlabel(r"$Hx_{\mathrm{PCAL}}$"+ " [cm]")
axs.set_ylabel(r"$Hy_{\mathrm{PCAL}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidPCAL.pdf")
plt.clf()

expSample_post = electronFiducial(expSample)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample_post.EcalHx1, expSample_post.EcalHy1, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, thousand, tenthousands]
ticks = [1, 10, 100, 1000, 10000]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-400, 400])
axs.set_ylim([-400, 400])
axs.set_xticks([-400, -200, 0, 200, 400])
axs.set_yticks([-400, -200, 0, 200, 400])
axs.set_title("(b) " + r"$e'$"+" PCAL Hits, Post-fiducial (Inb.)")
axs.set_xlabel(r"$Hx_{\mathrm{PCAL}}$"+ " [cm]")
axs.set_ylabel(r"$Hy_{\mathrm{PCAL}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/postcut_efidPCAL.pdf")
plt.clf()
# 6. exp e' pcal pre + post (outb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSampleOutb.EcalHx1, expSampleOutb.EcalHy1, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-400, 400])
axs.set_ylim([-400, 400])
axs.set_xticks([-400, -200, 0, 200, 400])
axs.set_yticks([-400, -200, 0, 200, 400])
axs.set_title("(a) " + r"$e'$"+" PCAL Hits, Pre-fiducial (Outb.)")
axs.set_xlabel(r"$Hx_{\mathrm{PCAL}}$"+ " [cm]")
axs.set_ylabel(r"$Hy_{\mathrm{PCAL}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidPCALOutb.pdf")
plt.clf()

expSampleOutb_post = electronFiducial(expSampleOutb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSampleOutb_post.EcalHx1, expSampleOutb_post.EcalHy1, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlim([-400, 400])
axs.set_ylim([-400, 400])
axs.set_xticks([-400, -200, 0, 200, 400])
axs.set_yticks([-400, -200, 0, 200, 400])
axs.set_title("(b) " + r"$e'$"+" PCAL Hits, Post-fiducial (Outb.)")
axs.set_xlabel(r"$Hx_{\mathrm{PCAL}}$"+ " [cm]")
axs.set_ylabel(r"$Hy_{\mathrm{PCAL}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/postcut_efidPCALOutb.pdf")
plt.clf()



#7. exp SF pre (inb) V, W
fig, axs = plt.subplots(1, 1, figsize = (6, 5))
h = axs.hist2d(expSample.EcalV1, expSample.ESamplFrac, bins = [np.linspace(0, 30, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, five+times+hundred]
ticks = [1, 10, 100, 500]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.axvline(19, color = 'k', linestyle = '--', linewidth = 5)
axs.set_xlabel(r"$l_{V}$"+ " [cm]")
axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
axs.set_xticks([0, 5, 10, 15, 19, 20, 25, 30])
axs.set_xticklabels([0, 5, 10, 15, 19, "", 25, 30])
axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
axs.set_title("(a) "+r"$e'$"+" PCAL "+r"$l_{V}$" +" Pre-fiducial (Inb.)")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidV.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (6, 5))
h = axs.hist2d(expSample.EcalW1, expSample.ESamplFrac, bins = [np.linspace(0, 30, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, five+times+hundred]
ticks = [1, 10, 100, 500]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.axvline(19, color = 'k', linestyle = '--', linewidth = 5)
axs.set_xlabel(r"$l_{W}$"+ " [cm]")
axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
axs.set_xticks([0, 5, 10, 15, 19, 20, 25, 30])
axs.set_xticklabels([0, 5, 10, 15, 19, "", 25, 30])
axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
axs.set_title("(b) "+r"$e'$"+" PCAL "+r"$l_{W}$" +" Pre-fiducial (Inb.)")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidW.pdf")
plt.clf()


fig, axs = plt.subplots(1, 1, figsize = (6, 5))
h = axs.hist2d(expSample.EcalU1, expSample.ESamplFrac, bins = [np.linspace(0, 400, 101), np.linspace(0.16, 0.34, 101)], cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, five+times+hundred]
ticks = [1, 10, 100, 500]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.axvline(395, color = 'k', linestyle = '--', linewidth = 5)
axs.set_xlabel(r"$l_{W}$"+ " [cm]")
axs.set_ylabel(r"$E_{dep.}/p_{e'}$")
axs.set_xticks([0, 100, 200, 300, 395, 400])
axs.set_xticklabels([0, 100, 200, 300, 395, ''])
axs.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
axs.set_yticklabels([0.15, 0.2, 0.25, 0.3, 0.35])
axs.set_title("(c) "+r"$e'$"+" PCAL "+r"$l_{U}$" +" Pre-fiducial (Inb.)")
plt.tight_layout()
plt.savefig("plots/ch2/precut_efidU.pdf")
plt.clf()

#8. exp edep vs edep (pre)
#pass (not used)
#9. exp + mc vz (inb) (pre)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))

axs.hist(expSample.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'k', density = True, label = 'Experimental Data')
axs.hist(dvcsSample.Evz, bins = np.linspace(-20, 15, 35*4+1), histtype = 'step', color = 'r', density = True, label = 'Simulation', weights = dvcsSample.weights)
axs.axvline(-8, color = 'k', linestyle = '--', linewidth = 5)
axs.axvline(2, color = 'k', linestyle = '--', linewidth = 5)
axs.set_xlim([-20, 20])
axs.set_xticks([-20, -10, -8, 0, 2, 10, 20 ])
axs.set_xticklabels([-20, "", -8, "", 2, 10, 20])
axs.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
axs.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2])
# plt.hist(dvcsSample.Enphe, bins = np.linspace(0, 50, 51), density = True, histtype = 'step')
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 1)
axs.set_xlabel("$vz_{e'}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precutVz_inb.pdf")
plt.clf()

# 10. exp E/p vs E/p (inb) (pre)
#pass (use the same plots w/v2)
# 11. exp chi vs p (inb) (pre)
#pass (use the same plots w/v2)
# 12. exp p' dc pre + post (inb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample.PDc3Hitx, expSample.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, three+times+hundred]
ticks = [1, 10, 100, 300]
# cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(a) " + r"$p'$"+" DC Outmost Layer Hits, Pre-fiducial (Inb.)")
# axs.set_xticks([-400, -200, 0, 200, 400])
# axs.set_yticks([-400, -200, 0, 200, 400])
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfidDC.pdf")
plt.clf()

expSample_post = protonFiducial(expSample)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample_post.PDc3Hitx, expSample_post.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticklabels = [one, ten, hundred, three+times+hundred]
ticks = [1, 10, 100, 300]
# cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(b) " + r"$p'$"+" DC Outmost Layer Hits, Post-fiducial (Inb.)")
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/postcut_pfidDC.pdf")
plt.clf()

# 13. exp p' dc pre + post (outb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSampleOutb.PDc3Hitx, expSampleOutb.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 300]
ticklabels = [one, ten, hundred, three + times + hundred]
# cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(a) " + r"$p'$"+" DC Outmost Layer Hits, Pre-fiducial (Outb.)")
# axs.set_xticks([-400, -200, 0, 200, 400])
# axs.set_yticks([-400, -200, 0, 200, 400])
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfidDCOutb.pdf")
plt.clf()

expSampleOutb_post = protonFiducial(expSampleOutb)
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSampleOutb_post.PDc3Hitx, expSampleOutb_post.PDc3Hity, bins = np.linspace(-400, 400, 100), cmin =1 , cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 300]
ticklabels = [one, ten, hundred, three + times + hundred]
# cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(b) " + r"$p'$"+" DC Outmost Layer Hits, Post-fiducial (Outb.)")
axs.set_xlabel(r"$x_{\mathrm{DC}}$"+ " [cm]")
axs.set_ylabel(r"$y_{\mathrm{DC}}$"+ " [cm]")
plt.tight_layout()
plt.savefig("plots/ch2/postcut_pfidDCOutb.pdf")
plt.clf()




# 14. exp gamma calorimeter pre + post (inb + outb)
expSample_post     = gammaFiducial(expSample)
expSampleOutb_post = gammaFiducial(expSampleOutb)

expSample_all = pd.concat([expSample.loc[:, ["GcX", "GcY", "Gsector", "config"]], expSampleOutb.loc[:, ["GcX", "GcY", "Gsector", "config"]]])
expSample_all_post = pd.concat([expSample_post.loc[:, ["GcX", "GcY", "Gsector", "config"]], expSampleOutb_post.loc[:, ["GcX", "GcY", "Gsector", "config"]]])


fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample_all.loc[expSample_all.config<3, "GcX"], expSample_all.loc[expSample_all.config<3, "GcY"], bins = np.linspace(-450, 450, 100), cmin = 1, cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_xlabel("$x_{\mathrm{PCAL}}$"+" (cm)")
axs.set_ylabel("$y_{\mathrm{PCAL}}$"+" (cm)")
axs.set_xlim([-450, 450])
axs.set_ylim([-450, 450])
axs.set_title("(a) " + r"$\gamma$"+" PCAL Hits, Pre-fiducial")
axs.set_xticks([-450, -300, -150, 0, 150, 300, 450])
axs.set_yticks([-450, -300, -150, 0, 150, 300, 450])
plt.tight_layout()
plt.savefig("plots/ch2/precut_gfidPCAL.pdf")
plt.clf()


fig, axs = plt.subplots(1, 1, figsize = (8, 5))
ang = -np.radians((expSample_all_post.loc[expSample_all_post.Gsector<7, "Gsector"]-1) * 60)
GcX_rot = expSample_all_post.loc[expSample_all_post.Gsector<7, "GcY"] * np.sin(ang) + expSample_all_post.loc[expSample_all_post.Gsector<7, "GcX"] * np.cos(ang)
GcY_rot = expSample_all_post.loc[expSample_all_post.Gsector<7, "GcY"] * np.cos(ang) - expSample_all_post.loc[expSample_all_post.Gsector<7, "GcX"] * np.sin(ang)
h = axs.hist2d(GcX_rot, GcY_rot, bins = np.linspace(-450, 450, 100), cmin = 1, cmap = cmap, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(b) " + r"$\gamma$"+" PCAL Hits, Post-fiducial")
axs.set_xlabel("$x_{\mathrm{PCAL}}$"+" (cm)")
axs.set_ylabel("$y_{\mathrm{PCAL}}$"+" (cm)")
axs.set_xlim([-450, 450])
axs.set_ylim([-450, 450])
plt.tight_layout()
plt.savefig("plots/ch2/postcut_gfidPCAL.pdf")
plt.clf()

# 15. exp proton CVT vs theta, phi(pre)
# 16. exp + mc theta distribution (pre)
# 17. exp + mc CVT theta distribution (pre)
# 18. exp + mc phi cvt distribution (pre)



expSample_post     = protonFiducial(expSample)
expSampleOutb_post = protonFiducial(expSampleOutb)

expSample_all = pd.concat([expSample.loc[:, ["Ptheta", "PCvt12theta", "config"]], expSampleOutb.loc[:, ["Ptheta", "PCvt12theta", "config"]]])
expSample_all_post = pd.concat([expSample_post.loc[:, ["Ptheta", "PCvt12theta", "config"]], expSampleOutb_post.loc[:, ["Ptheta", "PCvt12theta", "config"]]])


def linearfit(args, x):
  x = np.array(x)
  a, b = args
  return a + b * x
fig, axs = plt.subplots(1, 1, figsize = (8, 5))
expSample.loc[(expSample.config>1), "Ptheta"].hist(bins = np.linspace(30, 80, 101), histtype = 'step', density = True, ax = axs, color = 'k', label = 'Experimental Data')
dvcsSample.loc[(dvcsSample.config>1), "Ptheta"].hist(bins = np.linspace(30, 80, 101), histtype = 'step', density = True, ax = axs, color = 'r', label = 'Simulation', weights= dvcsSample.loc[dvcsSample.config>1, "weights"])
plt.axvline(64.23, color = 'k', linestyle = '--', linewidth = 2)
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 0.5)
plt.xlabel(r"$\theta_{p'}$" + " ["+degree+"]")
axs.set_title("(a)", loc = 'left')
axs.set_xlim([30, 85])
axs.set_xticks([30, 40, 50, 60, 64.23, 70, 80])
axs.set_xticklabels([30, 40, 50, 60, "", 70, 80])
axs.set_yticks([0, 0.05, 0.1])
axs.set_ylim([0, 0.1])
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd.pdf")
plt.clf()


fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample.loc[expSample.config>1, "Ptheta"], expSample.loc[expSample.config>1, "PCvt12theta"], bins = [np.linspace(30,70, 101), np.linspace(40, 80, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 3000]
ticklabels = [one, ten, hundred, thousand, three + times + thousand]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
plt.ylabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
plt.xlabel(r"$\theta_{~p'}$" + " ["+degree+"]")
axs.set_xticks([30, 35,40, 45,50,55, 60, 65, 70])
axs.set_xticklabels([30, 35,40, 45,50,55, 60, 65, 70])
axs.set_yticks([40, 50, 60, 70, 80])
axs.set_yticklabels([40, 50, 60, 70, 80])
x1 = np.linspace((2.924+46.5)/1.274, 64.23, 101)
axs.plot(x1, linearfit([-2.924, 1.274], x1), color = 'k', linewidth = 4, linestyle = '--')
x1 = np.linspace((3.523+46.5)/(1.046), 64.23, 101)
axs.plot(x1, linearfit([-3.523, 1.046], x1), color = 'k', linewidth = 4, linestyle = '--')
x1 = np.linspace((2.924+46.5)/1.274, (3.523+46.5)/(1.046), 101)
axs.plot(x1, 46.5 + 0*x1, color = 'k', linewidth = 4, linestyle = '--')
y1 = np.linspace(-2.924+1.274*64.23, -3.523+1.046*64.23, 101)
axs.plot(y1*0 + 64.23, y1, color = 'k', linewidth = 4, linestyle = '--')
axs.set_title("(a)", loc = 'left')
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd2.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample.loc[expSample.config>1, "PCvt12phi"], expSample.loc[expSample.config>1, "PCvt12theta"], bins = [np.linspace(-180,180, 361), np.linspace(40, 80, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100]
ticklabels = [one, ten, hundred]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
# plt.axvline(62)
plt.ylabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
plt.xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
axs.set_xticks([-180, -90, 0, 90, 180])
axs.set_xticklabels([-180, -90, 0, 90, 180])
axs.set_yticks([40, 50, 60, 70, 80])
axs.set_yticklabels([40, 50, 60, 70, 80])
plt.axvline(-95,  color = 'k', linestyle = '--', linewidth = 4)
plt.axvline(-80,  color = 'k', linestyle = '--', linewidth = 4)
plt.axvline(25,  color = 'k', linestyle = '--', linewidth = 4)
plt.axvline(40,  color = 'k', linestyle = '--', linewidth = 4)
plt.axvline(143,  color = 'k', linestyle = '--', linewidth = 4)
plt.axvline(158,  color = 'k', linestyle = '--', linewidth = 4)
# plt.axhline(46.5,  color = 'k', linestyle = '--', linewidth = 4)
axs.set_title("(b)", loc = 'left')
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd3.pdf")
plt.clf()


fig, axs = plt.subplots(1, 1, figsize = (8, 5))
expSample.loc[(expSample.config>1), "PCvt12theta"].hist(bins = np.linspace(40, 90, 101), histtype = 'step', density = True, ax = axs, color = 'k', label = 'Experimental Data')
dvcsSample.loc[(dvcsSample.config>1), "PCvt12theta"].hist(bins = np.linspace(40, 90, 101), histtype = 'step', density = True, ax = axs, color = 'r', label = 'Simulation', weights = dvcsSample.loc[(dvcsSample.config>1), "weights"])
plt.axvline(46.5, color = 'k', linestyle = '--', linewidth = 2)
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 0.5)
plt.xlabel(r"$\theta_{\mathrm{CVT}}$" + " ["+degree+"]")
axs.set_title("(a)", loc = 'left')
axs.set_xlim([40, 90])
axs.set_xticks([40, 46.5, 50, 60, 70, 80, 90])
axs.set_xticklabels([40, "", 50, 60, 70, 80, 90])
axs.set_yticks([0, 0.03, 0.06])
axs.set_ylim([0, 0.06])
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd_cvt.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (8, 5))
expSample.loc[(expSample.config>1), "PCvt12phi"].hist(bins = np.linspace(-180, 180, 181), histtype = 'step', density = True, ax = axs, color = 'k', label = 'Experimental Data')
dvcsSample.loc[(dvcsSample.config>1), "PCvt12phi"].hist(bins = np.linspace(-180, 180, 181), histtype = 'step', density = True, ax = axs, color = 'r', label = 'Simulation', weights = dvcsSample.loc[(dvcsSample.config>1), "weights"])
plt.xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
axs.set_xticks([-180, -90, 0, 90, 180])
axs.set_xticklabels([-180, -90, 0, 90, 180])
plt.legend(loc='lower left', bbox_to_anchor = (0.6, 0.6), framealpha = 1)
axs.set_xlim([-180, 180])
axs.set_ylim([0, 0.006])
axs.set_yticks([0, 0.002, 0.004, 0.006])
axs.set_title("(a)", loc = 'left')
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd6.pdf")
plt.clf()


fig, axs = plt.subplots(1, 1, figsize = (8, 5))
dvcsPthetahist, bins = np.histogram(dvcsSample.loc[(dvcsSample.config>1), "Ptheta"], bins = np.linspace(60, 70, 101), weights = dvcsSample.loc[(dvcsSample.config>1), "weights"])
expPthetahist, bins = np.histogram(expSample.loc[(expSample.config>1), "Ptheta"], bins = np.linspace(60, 70, 101))
bincenters = (bins[:-1] + bins[1:])/2
axs.errorbar(bincenters, expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist)), xerr = 0.03, yerr =expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))*np.sqrt(1/dvcsPthetahist + 1/expPthetahist), linestyle = '', color = 'k')
axs.axvline(64.23, color = 'k', linestyle = '--')
fom = expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))
for j in range(60, 100):
  moving_average = np.sum(fom[0:j])/j
  current = fom[j]
  if current > moving_average*1.05:
    print(j, bincenters[j], moving_average, fom[j])
    break
axs.set_xlabel(r"$\theta_{p'}$"+ " ["+degree+"]")
axs.set_ylabel(r"$n_{exp.}/n_{sim.}$")
axs.set_xlim([60, 70])
axs.set_xticks([60, 62, 64.23, 66, 68, 70])
axs.set_xticklabels([60, 62, 64.23, 66, 68, 70])
axs.set_title("(b)", loc = 'left')
axs.set_yticks([0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5])
axs.set_ylim([0.75, 2.5])
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd4.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (8, 5))
bins = np.linspace(40, 80, 81)
dvcsPthetahist, bins = np.histogram(dvcsSample.loc[(dvcsSample.config>1), "PCvt12theta"], bins = bins, weights = dvcsSample.loc[(dvcsSample.config>1), "weights"])
expPthetahist, bins = np.histogram(expSample.loc[(expSample.config>1), "PCvt12theta"], bins = bins)
bincenters = (bins[:-1] + bins[1:])/2
axs.errorbar(bincenters, expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist)), xerr = 0.25, yerr = expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))*np.sqrt(1/dvcsPthetahist + 1/expPthetahist), linestyle = '', color = 'k')
for j in range(79, 0, -1):
  moving_average = np.sum((expPthetahist/ dvcsPthetahist)[j:])/(len(bins)-1-j)
  current = (expPthetahist/ dvcsPthetahist)[j]
  if current > moving_average*1.05:
    print(j, bincenters[j], moving_average,(expPthetahist/ dvcsPthetahist)[j])
    break
axs.axvline(46.5, color = 'k', linestyle = '--')
axs.set_xlim([40, 80])
axs.set_xticks([40, 46.5, 50, 60, 70, 80])
axs.set_xticklabels([40, 46.5, 50, 60, 70, 80])
axs.set_xlabel(r"$\theta_{\mathrm{CVT}}$"+ " ["+degree+"]")
axs.set_ylabel(r"$n_{exp.}/n_{sim.}$")
axs.set_yticks([0.5, 1, 1.5, 2, 2.5])
axs.set_ylim([0.5, 2.7])
axs.set_title("(b)", loc = 'left')
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd5.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (8, 5))
dvcsPthetahist, bins = np.histogram(dvcsSample.loc[(dvcsSample.config>1), "PCvt12phi"], bins = np.linspace(-180, 180, 181), weights = dvcsSample.loc[(dvcsSample.config>1), "weights"])
expPthetahist, bins = np.histogram(expSample.loc[(expSample.config>1), "PCvt12phi"], bins = np.linspace(-180, 180, 181))
bincenters = (bins[:-1] + bins[1:])/2
axs.errorbar(bincenters, expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist)), xerr = 0.25, yerr = expPthetahist/np.sum(expPthetahist)/ (dvcsPthetahist/np.sum(dvcsPthetahist))*np.sqrt(1/dvcsPthetahist + 1/expPthetahist), linestyle = '', color = 'k')
axs.set_yticks([0, 0.5, 1, 1.5])
axs.set_ylim([0, 1.5])
axs.set_xlim([-180, 180])
axs.set_xticks([-180, -90, 0, 90, 180])
axs.set_xticklabels([-180, -90, 0, 90, 180])
axs.axvline(-95,  color ='k', linestyle = '--')
axs.axvline(-80,  color ='k', linestyle = '--')
axs.axvline(25,  color ='k', linestyle = '--')
axs.axvline(40,  color ='k', linestyle = '--')
axs.axvline(143,  color ='k', linestyle = '--')
axs.axvline(158,  color ='k', linestyle = '--')
axs.set_xlabel(r"$\phi_{\mathrm{CVT}}$" + " ["+degree+"]")
axs.set_ylabel(r"$n_{exp.}/n_{sim.}$")
axs.set_title("(b)", loc = 'left')
plt.tight_layout()
plt.savefig("plots/ch2/precut_pfid_cd7.pdf")
plt.clf()

# 19. exp photon FT (inb+outb, pre + post)
expSample_post     = gammaFiducial(expSample)
expSampleOutb_post = gammaFiducial(expSampleOutb)

expSample_all = pd.concat([expSample.loc[:, ["GcX", "GcY", "Gsector", "config"]], expSampleOutb.loc[:, ["GcX", "GcY", "Gsector", "config"]]])
expSample_all_post = pd.concat([expSample_post.loc[:, ["GcX", "GcY", "Gsector", "config"]], expSampleOutb_post.loc[:, ["GcX", "GcY", "Gsector", "config"]]])

fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample_all.loc[expSample_all.config==3, "GcX"], expSample_all.loc[expSample_all.config==3, "GcY"], bins = [np.linspace(-20, 20, 101), np.linspace(-20, 20, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(a) " + r"$\gamma$"+" FT-Cal Hits, Pre-fiducial")
axs.set_ylabel(r"$y_{\mathrm{FT}}$" + " ["+degree+"]")
axs.set_xlabel(r"$x_{\mathrm{FT}}$" + " ["+degree+"]")
theta = np.linspace(0, 2*np.pi, 101)
circleCenterX1 = -8.419
circleCenterY1 = 9.889
circleRadius1 = 1.6

circleCenterX2 = -9.89
circleCenterY2 = -5.327
circleRadius2 = 1.6

circleCenterX3 = -6.15
circleCenterY3 = -13
circleRadius3 = 2.3
    
circleCenterX4 = 3.7
circleCenterY4 = -6.5
circleRadius4 = 2

plt.plot(circleRadius1*np.cos(theta) + circleCenterX1, circleRadius1*np.sin(theta) + circleCenterY1, color = 'k', linewidth = 1, linestyle = '--')
plt.plot(circleRadius2*np.cos(theta) + circleCenterX2, circleRadius2*np.sin(theta) + circleCenterY2, color = 'k', linewidth = 1, linestyle = '--')
plt.plot(circleRadius3*np.cos(theta) + circleCenterX3, circleRadius3*np.sin(theta) + circleCenterY3, color = 'k', linewidth = 1, linestyle = '--')
plt.plot(circleRadius4*np.cos(theta) + circleCenterX4, circleRadius4*np.sin(theta) + circleCenterY4, color = 'k', linewidth = 1, linestyle = '--')
plt.tight_layout()
plt.savefig("plots/ch2/precut_gfidFT.pdf")
plt.clf()

fig, axs = plt.subplots(1, 1, figsize = (8, 5))
h = axs.hist2d(expSample_all_post.loc[expSample_all_post.config==3, "GcX"], expSample_all_post.loc[expSample_all_post.config==3, "GcY"], bins = [np.linspace(-20, 20, 101), np.linspace(-20, 20, 101)], cmap = cmap, cmin = 1, rasterized = True, norm = LogNorm())
ticks = [1, 10, 100, 1000, 10000]
ticklabels = [one, ten, hundred, thousand, tenthousands]
cbar = plt.colorbar(h[3], ax = axs, ticks = ticks)
# cbar.ax.set_yticklabels(ticklabels)
axs.set_title("(b) " + r"$\gamma$"+" FT-Cal Hits, Post-fiducial")
axs.set_ylabel(r"$y_{\mathrm{FT}}$" + " ["+degree+"]")
axs.set_xlabel(r"$x_{\mathrm{FT}}$" + " ["+degree+"]")
theta = np.linspace(0, 2*np.pi, 101)
circleCenterX1 = -8.419
circleCenterY1 = 9.889
circleRadius1 = 1.6

circleCenterX2 = -9.89
circleCenterY2 = -5.327
circleRadius2 = 1.6

circleCenterX3 = -6.15
circleCenterY3 = -13
circleRadius3 = 2.3
    
circleCenterX4 = 3.7
circleCenterY4 = -6.5
circleRadius4 = 2

plt.plot(circleRadius1*np.cos(theta) + circleCenterX1, circleRadius1*np.sin(theta) + circleCenterY1, color = 'k', linewidth = 1, linestyle = '--')
plt.plot(circleRadius2*np.cos(theta) + circleCenterX2, circleRadius2*np.sin(theta) + circleCenterY2, color = 'k', linewidth = 1, linestyle = '--')
plt.plot(circleRadius3*np.cos(theta) + circleCenterX3, circleRadius3*np.sin(theta) + circleCenterY3, color = 'k', linewidth = 1, linestyle = '--')
plt.plot(circleRadius4*np.cos(theta) + circleCenterX4, circleRadius4*np.sin(theta) + circleCenterY4, color = 'k', linewidth = 1, linestyle = '--')
plt.tight_layout()
plt.savefig("plots/ch2/postcut_gfidFT.pdf")
plt.clf()
'''
# '''
#draw nominal
topo   = {1: "FD", 2:"CD", 3: "CDFT"}

#inbending
yticklabel_1 = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
yticklabel_2 = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
yticklabel_3 = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
yticklabel_1[-1]  = "{}\n{}".format(yticklabel_1[-1], yticklabel_1[0])
yticklabel_1[0]   = ''
yticklabel_2[0]   = ''
yticklabel_3[0]  = "0\n"
yticklabel_3[-1]  = "{}\n{}".format(yticklabel_3[-1], yticklabel_3[0])

for polarity in ["inb", "outb"]:
  for config in [1, 2, 3]:
    if polarity == "inb" and config == 1:
      cuts_dvcs_3sigma = cuts_dvcs_FD_Inb_3sigma
      cuts_dvpi0p_3sigma = cuts_dvpi0p_FD_Inb_3sigma
    if polarity == "inb" and config == 2:
      cuts_dvcs_3sigma = cuts_dvcs_CD_Inb_3sigma
      cuts_dvpi0p_3sigma = cuts_dvpi0p_CD_Inb_3sigma
    if polarity == "inb" and config == 3:
      cuts_dvcs_3sigma = cuts_dvcs_CDFT_Inb_3sigma
      cuts_dvpi0p_3sigma = cuts_dvpi0p_CDFT_Inb_3sigma
    if polarity == "outb" and config == 1:
      cuts_dvcs_3sigma = cuts_dvcs_FD_Outb_3sigma
      cuts_dvpi0p_3sigma = cuts_dvpi0p_FD_Outb_3sigma
    if polarity == "outb" and config == 2:
      cuts_dvcs_3sigma = cuts_dvcs_CD_Outb_3sigma
      cuts_dvpi0p_3sigma = cuts_dvpi0p_CD_Outb_3sigma
    if polarity == "outb" and config == 3:
      cuts_dvcs_3sigma = cuts_dvcs_CDFT_Outb_3sigma
      cuts_dvpi0p_3sigma = cuts_dvpi0p_CDFT_Outb_3sigma
    for varind in range(8):

      # DVCS plots
      var    = dvcsvars[varind]
      label  = dvcstitles[varind]
      unit   = dvcsunits[varind]
      exp_epg   = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{0}/dvcs/excl_level_2/pkl_7_nominal/fall2018_{0}.pkl".format(polarity))
      if var == "coneAngle":
          xub = np.poly1d(cuts_dvcs_3sigma["coneAngle_ub"])(exp_epg.loc[exp_epg.config==config].Etheta).max()
          xlb = np.poly1d(cuts_dvcs_3sigma["coneAngle_lb"])(exp_epg.loc[exp_epg.config==config].Etheta).min()
      else:
          xub = cuts_dvcs_3sigma["{}_ub".format(var)]
          if "{}_lb".format(var) in cuts_dvcs_3sigma.keys():
              xlb = cuts_dvcs_3sigma["{}_lb".format(var)]
          else:
              xlb = 0


      xtick        = [xlb, xub]
      xticklabel_1 = ["{:.3f}".format(i) for i in xtick]
      xticklabel_2 = ["{:.3f}".format(i) for i in xtick]
      xticklabel_3 = ["{:.3f}".format(i) for i in xtick]
      xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
      xticklabel_1[-1] = ''
      xticklabel_2[-1] = ''
      xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])

      ind = 0 
      filenum = 0
      
      fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))

      yub    = np.zeros((5, 5))
      filled = np.zeros((5, 5))
      for integrated_binnum in range(1, 147+1):
        exp_epg   = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/dvcs/excl_level_2/restructured_7_nominal/{}.pkl".format(polarity, integrated_binnum))
        exp_epgg  = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(polarity, integrated_binnum))
        sim_epg    = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/dvcs_km15/excl_level_1/restructured_7_nominal/{}.pkl".format(polarity, integrated_binnum))
        sim_bkg_1g = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/pi0_1gamma/excl_level_1/restructured_7_nominal/2/{}.pkl".format(polarity, integrated_binnum))
        sim_bkg_2g = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/pi0_2gamma/excl_level_1/restructured_7_nominal/2/{}.pkl".format(polarity, integrated_binnum))
        df3 = exp_epg   .loc[(exp_epg   .config == config )  , :]
        df1 = sim_epg   .loc[(sim_epg   .config == config )  , :]
        df2 = sim_bkg_1g.loc[(sim_bkg_1g.config == config )  , :]
        df4 = exp_epgg  .loc[(exp_epgg  .config == config )  , :]
        df5 = sim_bkg_2g.loc[(sim_bkg_2g.config == config )  , :]
      
        if len(df3)>50:
          pass
        else:
          continue
        if var == "coneAngle":
          bins = np.linspace(xlb, xub, 100+1)
        else:
          bins = np.linspace(xlb, xub, 20+1)

        expDist     , _  = np.histogram(df3.loc[:, var], bins, density = True)
        simDist_dvcs, _     = np.histogram(df1.loc[:, var], bins, density = True, weights = df1.weights)
        try:
          cont_thisbin = len(df4)*len(df2)/len(df5)/len(df3)
          simDist_dvpi0, _ = np.histogram(df2.loc[:, var], bins, density = True)
          simDist = (1-cont)*simDist_dvcs + cont*simDist_dvpi0
        except:
          cont_thisbin = 0
          simDist      = simDist_dvcs
        
        xind = ind//5
        yind = ind%5
        ind  = ind + 1
        yub[xind, yind] = np.max([expDist, simDist])
        filled[xind, yind] = 1
                                    
        axs[xind, yind].hist(bins[:-1], bins, weights = expDist, histtype = 'step', color='b', linewidth=3)
        axs[xind, yind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=3)

        annotation = "{}".format(integrated_binnum)
        axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)

        axs[xind, yind].set_ylim([0, 1.2*yub[xind, yind]])
        axs[xind, yind].set_xlim([xlb, xub])

        axs[xind, yind].get_xaxis().set_visible(False)
        axs[xind, yind].get_yaxis().set_visible(False)

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_2)

        elif (xind == 5 -1) and (yind < 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_1)

        if (xind == 0) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_2)

        elif (xind < 5 -1) and (yind == 0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_1)

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_3)

        if (xind == 5 -1) and (yind == 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_3)

          plt.subplots_adjust(wspace=0, hspace=0)
          plt.savefig("plots/q16/dset_e/{}/dvcs_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
          filenum = filenum + 1
          plt.clf()
          ind = 0
          fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
          yub    = np.zeros((5, 5))
          filled = np.zeros((5, 5))

      for xind, yind in itertools.product(range(5), range(5)):      

        if filled[xind, yind]:
          continue

        axs[xind, yind].get_xaxis().set_visible(False)
        axs[xind, yind].get_yaxis().set_visible(False)
        axs[xind, yind].set_xlim([xlb, xub])

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_2)

        elif (xind == 5 -1) and (yind < 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_1)

        if (xind == 0) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_2)

        elif (xind < 5 -1) and (yind == 0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_1)

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_3)

        if (xind == 5 -1) and (yind == 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_3)

          plt.subplots_adjust(wspace=0, hspace=0)
          plt.savefig("plots/q16/dset_e/{}/dvcs_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
          plt.clf()

      # DVπ0P plots
      var   = pi0vars[varind]
      label = pi0titles[varind]
      unit  = pi0units[varind]


      xub = cuts_dvpi0p_3sigma["{}_ub".format(var)]
      if "{}_lb".format(var) in cuts_dvpi0p_3sigma.keys():
          xlb = cuts_dvpi0p_3sigma["{}_lb".format(var)]
      else:
          xlb = 0

      xtick        = [xlb, xub]
      xticklabel_1 = ["{:.3f}".format(i) for i in xtick]
      xticklabel_2 = ["{:.3f}".format(i) for i in xtick]
      xticklabel_3 = ["{:.3f}".format(i) for i in xtick]
      xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
      xticklabel_1[-1] = ''
      xticklabel_2[-1] = ''
      xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])

      ind = 0 
      filenum = 0
      
      fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))

      yub    = np.zeros((5, 5))
      filled = np.zeros((5, 5))
      for integrated_binnum in range(1, 147+1):
        exp_epgg  = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/exp_fall2018_{}/pi0/excl_level_2/restructured_7_nominal/{}.pkl".format(polarity, integrated_binnum))
        sim_bkg_2g = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_{}/pi0_2gamma/excl_level_1/restructured_7_nominal/2/{}.pkl".format(polarity, integrated_binnum))
        df4 = exp_epgg  .loc[(exp_epgg  .config == config )  , :]
        df5 = sim_bkg_2g.loc[(sim_bkg_2g.config == config )  , :]
      
        if len(df4)>20:
          pass
        else:
          continue

        bins = np.linspace(xlb, xub, 20+1)
        expDist     , bins  = np.histogram(df4.loc[:, var],   20, density = True)
        simDist     , _     = np.histogram(df5.loc[:, var], bins, density = True)
        
        xind = ind//5
        yind = ind%5
        ind  = ind + 1
        yub[xind, yind] = np.max([expDist, simDist])
        filled[xind, yind] = 1
                                    
        axs[xind, yind].hist(bins[:-1], bins, weights = expDist, histtype = 'step', color='b', linewidth=3)
        axs[xind, yind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=3)

        annotation = "{}".format(integrated_binnum)
        axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)

        axs[xind, yind].set_ylim([0, 1.2*yub[xind, yind]])
        axs[xind, yind].set_xlim([xlb, xub])

        axs[xind, yind].get_xaxis().set_visible(False)
        axs[xind, yind].get_yaxis().set_visible(False)

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_2)

        elif (xind == 5 -1) and (yind < 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_1)

        if (xind == 0) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_2)

        elif (xind < 5 -1) and (yind == 0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_1)

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_3)

        if (xind == 5 -1) and (yind == 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_3)
            
          plt.subplots_adjust(wspace=0, hspace=0)
          plt.savefig("plots/q16/dset_e/{}/pi0_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
          plt.clf()
          filenum = filenum + 1
          ind = 0
          fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
          yub    = np.zeros((5, 5))
          filled = np.zeros((5, 5))
      for xind, yind in itertools.product(range(5), range(5)):      

        if filled[xind, yind]:
          continue

        axs[xind, yind].get_xaxis().set_visible(False)
        axs[xind, yind].get_yaxis().set_visible(False)

        axs[xind, yind].set_xlim([xlb, xub])

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_2)

        elif (xind == 5 -1) and (yind < 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_1)

        if (xind == 0) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_2)

        elif (xind < 5 -1) and (yind == 0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_1)

        if (xind == 5 -1) and (yind==0):
          axs[xind, yind].get_yaxis().set_visible(True)
          ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
          axs[xind, yind].set_yticks(ytick)
          axs[xind, yind].set_yticklabels(yticklabel_3)

        if (xind == 5 -1) and (yind == 5 -1):
          axs[xind, yind].get_xaxis().set_visible(True)
          axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
          axs[xind, yind].set_xticks(xtick)
          axs[xind, yind].set_xticklabels(xticklabel_3)

          plt.subplots_adjust(wspace=0, hspace=0)
          plt.savefig("plots/q16/dset_e/{}/pi0_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
          plt.clf()
# '''