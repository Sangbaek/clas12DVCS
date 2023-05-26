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
cmap = copy(plt.cm.get_cmap("jet"))
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from utils.const import *
from utils.physics import *
from utils.fiducial import *
from matplotlib.colors import LogNorm
import argparse
from glob import glob

cmap.set_under('w',0)
cmap.set_bad('w',0)

degree = r"${}^{\circ}$"
GeV = "GeV"
GeV2 = "GeV"+r"${}^{2}$"
GeVc = "GeV/c"
GeVc2 = "(GeV/c)"+r"${}^{2}$"

cmap.set_under('w',0)
cmap.set_bad('w',0)

import matplotlib
# initial settings
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

parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-q","--question", help="question no.", type = int)
args = parser.parse_args()


point = r'$.$'
zero = r'$0$'
one = r"$1$"
two = r"$2$"
three = r"$3$"
four = r"$4$"
five = r"$5$"
six = r"$6$"
seven = r"$7$"
eight = r"$8$"
nine = r"$9$"
times = r"$\times$"
ten = r"$10$"
hundred = r"$10^2$"
thousand = r"$10^3$"
tenthousands = r"$10^4$"
hundredthousands = r'$10^5$'
tenth = r"$10^{-1}$"
hundredth = r"$10^{-2}$"
thousandth = r"$10^{-3}$"


if args.question == 16:
  dfs = []
  for file in glob("../convPkl/vgg/inb/*.pkl"):
    dfs.append(pd.read_pickle(file))
  vgg_inb    = pd.concat(dfs)
  
  dfs = []
  for file in glob("../convPkl/bh/inb/*.pkl"):
    dfs.append(pd.read_pickle(file))
  bh_inb    = pd.concat(dfs)
  
  dfs = []
  for file in glob("../convPkl/bkg_1g/inb/*.pkl"):
    dfs.append(pd.read_pickle(file))
  bkg_1g_inb    = pd.concat(dfs)
  
  dfs = []
  for file in glob("../convPkl/bkg_2g/inb/*.pkl"):
    dfs.append(pd.read_pickle(file))
  bkg_2g_inb    = pd.concat(dfs)
  
  dfs = []
  for file in glob("../convPkl/exp/epg/inb/*.pkl"):
    dfs.append(pd.read_pickle(file))
  epg_inb    = pd.concat(dfs)
  
  dfs = []
  for file in glob("../convPkl/exp/epgg/inb/*.pkl"):
    dfs.append(pd.read_pickle(file))
  epgg_inb    = pd.concat(dfs)

  dfs = []
  topo   = {1: "FD", 2:"CD", 3: "CDFT"}

  for config in [1, 2, 3]:
    for varind in range(8):
      var    = dvcsvars[varind]
      label  = dvcstitles[varind]
      unit   = dvcsunits[varind]
      if var in ["MPt", "reconGam", "coplanarity"]:
        lb        = 0
      if config == 1:
        ub        = cuts_dvcs_FD_Inb_3sigma["{}_ub".format(var)]
      if config == 2:
        ub        = cuts_dvcs_CD_Inb_3sigma["{}_ub".format(var)]
      if config == 3:
        ub        = cuts_dvcs_CDFT_Inb_3sigma["{}_ub".format(var)]
      elif var == "coneAngle":
        if config == 1:
          lb        = 25
          ub        = 45
        if config == 2:
          lb        = 10
          ub        = 30
        if config == 3:
          lb        = 10
          ub        = 30
      else:
        if config == 1:
          lb        = cuts_dvcs_FD_Inb_3sigma["{}_lb".format(var)]
          ub        = cuts_dvcs_FD_Inb_3sigma["{}_ub".format(var)]
        if config == 2:
          lb        = cuts_dvcs_CD_Inb_3sigma["{}_lb".format(var)]
          ub        = cuts_dvcs_CD_Inb_3sigma["{}_ub".format(var)]
        if config == 3:
          lb        = cuts_dvcs_CDFT_Inb_3sigma["{}_lb".format(var)]
          ub        = cuts_dvcs_CDFT_Inb_3sigma["{}_ub".format(var)]
      bins          = np.linspace(lb, ub, 21)
      ind = 0 
      filenum = 0
      fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
      for tbin in range(len(tbins)-1):
        for xBbin in range(len(newxBbins2)-1):
          for Q2bin in range(len(newxBbins2)-1):
            df_bh = bh_inb.loc[(bh_inb.t1>newtbins[tbin]) & (bh_inb.t1<newtbins[tbin+1]), :]
            df_vgg = vgg_inb.loc[(vgg_inb.t1>newtbins[tbin]) & (vgg_inb.t1<newtbins[tbin+1]), :]
            df_bkg_1g = bkg_1g_inb.loc[(bkg_1g_inb.t1>newtbins[tbin]) & (bkg_1g_inb.t1<newtbins[tbin+1]), :]
            df_bkg_2g = bkg_2g_inb.loc[(bkg_2g_inb.t1>newtbins[tbin]) & (bkg_2g_inb.t1<newtbins[tbin+1]), :]
            df_epg = epg_inb.loc[(epg_inb.t1>newtbins[tbin]) & (epg_inb.t1<newtbins[tbin+1]), :]
            df_epgg = epgg_inb.loc[(epgg_inb.t1>newtbins[tbin]) & (epgg_inb.t1<newtbins[tbin+1]), :]
            df_bh = df_bh.loc[(df_bh.xB>newxBbins2[xBbin]) & (df_bh.xB<newxBbins2[xBbin+1]), :]
            df_bh = df_bh.loc[(df_bh.Q2>newQ2bins2[Q2bin]) & (df_bh.Q2<newQ2bins2[Q2bin+1]), :]
            df_vgg = df_vgg.loc[(df_vgg.xB>newxBbins2[xBbin]) & (df_vgg.xB<newxBbins2[xBbin+1]), :]
            df_vgg = df_vgg.loc[(df_vgg.Q2>newQ2bins2[Q2bin]) & (df_vgg.Q2<newQ2bins2[Q2bin+1]), :]
            df_bkg_1g = df_bkg_1g.loc[(df_bkg_1g.xB>newxBbins2[xBbin]) & (df_bkg_1g.xB<newxBbins2[xBbin+1]), :]
            df_bkg_1g = df_bkg_1g.loc[(df_bkg_1g.Q2>newQ2bins2[Q2bin]) & (df_bkg_1g.Q2<newQ2bins2[Q2bin+1]), :]
            df_bkg_2g = df_bkg_2g.loc[(df_bkg_2g.xB>newxBbins2[xBbin]) & (df_bkg_2g.xB<newxBbins2[xBbin+1]), :]
            df_bkg_2g = df_bkg_2g.loc[(df_bkg_2g.Q2>newQ2bins2[Q2bin]) & (df_bkg_2g.Q2<newQ2bins2[Q2bin+1]), :]
            df_epg = df_epg.loc[(df_epg.xB>newxBbins2[xBbin]) & (df_epg.xB<newxBbins2[xBbin+1]), :]
            df_epg = df_epg.loc[(df_epg.Q2>newQ2bins2[Q2bin]) & (df_epg.Q2<newQ2bins2[Q2bin+1]), :]
            df_epgg = df_epgg.loc[(df_epgg.xB>newxBbins2[xBbin]) & (df_epgg.xB<newxBbins2[xBbin+1]), :]
            df_epgg = df_epgg.loc[(df_epgg.Q2>newQ2bins2[Q2bin]) & (df_epgg.Q2<newQ2bins2[Q2bin+1]), :]
            bh_var     = df_bh.loc[df_bh.config == config, var]
            bkg_1g_var = df_bkg_1g.loc[df_bkg_1g.config == config, var]
            vgg_var    = df_vgg.loc[df_vgg.config == config, var]
            epg_var    = df_epg.loc[df_epg.config == config, var]
            
            if len(epg_var)>20:
              pass
            else:
              continue
            
            try:
              cont_thisbin = len(bkg_1g_var)/len(df_bkg_2g.loc[df_bkg_2g.config==config, :])*len(df_epgg.loc[df_epgg.config==config, :])/len(epg_var)
            except:
              cont_thisbin = 0
            
            xind = ind//5
            yind = ind%5
            ind = ind + 1
            
            bh_hist, _    = np.histogram(bh_var, bins = bins, density = True)
            vgg_hist, _   = np.histogram(vgg_var, bins = bins, density = True)
            epg_hist, _   = np.histogram(epg_var, bins = bins, density = True)
            if cont_thisbin > 0:
              bkg_1g_hist, _= np.histogram(bkg_1g_var, bins = bins, density = True)
            if np.isnan(bkg_1g_hist).any():
              bkg_1g_hist = np.zeros(len(bh_hist))
            else:
              bkg_1g_hist = np.zeros(len(bh_hist))
              bincenters = (bins[1:]+bins[:-1])/2.
              sim_hist   = cont_thisbin*bkg_1g_hist + (1-cont_thisbin)*bh_hist
              sim_hist2   = cont_thisbin*bkg_1g_hist + (1-cont_thisbin)*vgg_hist
                       
            
            axs[xind, yind].step(bincenters, epg_hist, where = 'mid', color = 'tab:blue')
            axs[xind, yind].step(bincenters, sim_hist, where = 'mid', color = 'tab:red')
            axs[xind, yind].step(bincenters, sim_hist2, where = 'mid', color = 'tab:green')
            annotation = "({}, {}, {})".format(xBbin, Q2bin, tbin)
            axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)
            axs[xind, yind].set_ylim(bottom = 0)
            
            
            if ind >24:
              for ax in axs.flatten():
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
              for ax in axs[-1, :]:
                ax.get_xaxis().set_visible(True)
                xtick        = np.linspace(lb, ub, 5)
                xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
                xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
                xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
                xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
                xticklabel_1[-1] = ''
                xticklabel_2[-1] = ''
                xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
                ax.set_xlabel( label + " [" + unit + "]")
                ax.set_xticks(xtick)
                ax.set_xticklabels(xticklabel_1)
                axs[-1, 0].set_xticklabels(xticklabel_2)
                axs[-1, -1].set_xticklabels(xticklabel_3)
                
                plt.subplots_adjust(wspace=0, hspace=0)
                plt.savefig("plots/q16/inb_{}_{}_{}.pdf".format(topo[config], var, filenum))
                filenum = filenum + 1
                plt.clf()
                ind = 0
              fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
            
      for ax in axs.flatten():
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
      for ax in axs[-1, :]:
        ax.get_xaxis().set_visible(True)
        xtick        = np.linspace(lb, ub, 5)
        xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
        xticklabel_1[-1]  = ''
        xticklabel_2[-1] = ''
        xticklabel_3[0] = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
        ax.set_xlabel( label + " [" + unit + "]")
        ax.set_xticks(xtick)
        ax.set_xticklabels(xticklabel_1)
        axs[-1, 0].set_xticklabels(xticklabel_2)
        axs[-1, -1].set_xticklabels(xticklabel_3)
      
      plt.subplots_adjust(wspace=0, hspace=0)
     ()
      plt.savefig("plots/q16/inb_{}_{}_{}.pdf".format(topo[config], var, filenum))
      plt.clf()
