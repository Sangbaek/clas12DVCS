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
cmap = copy(matplotlib.colormaps["jet"])
cmap.set_under('w',0)
cmap.set_bad('w',0)
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
parser.add_argument("-pl","--plot", help="plot no.", type = str)
parser.add_argument("-p","--polarity", help="polarity: inbending or outbending", type = str, default = 'inbending')
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

  if args.polarity == 'inbending':
    polarity = 'inb'
    pol_cap  = 'Inb'
  elif args.polarity == 'outbending':
    polarity = 'outb'
    pol_cap  = 'Outb'
  else:
    print("Polarity undefined")
    exit()

  # key = ('inb', 2, 'a')
  # for key in [, ('inb', 2, 'a'), ('inb', 3, 'a')]:
  #   xlbs_epg[key] = [5,  -10, 0, 0, -1, -0.1, -1, 0]
  #   xubs_epg[key] = [40,  10, 5, 10, 1, 0.1,  1, 1]
  #   yubs_epg[key] = [0.15, 2, 3, 0.4, 3, 250, 4, 25]
    # xticks_epg[key, 0] = [10, 20, 30]
    # xticks_epg[key, 1] = [0.5, 1, 1.5]
    # xticks_epg[key, 2] = [0, 0.2, 0.4, 0.6]
    # xticks_epg[key, 3] = [0, 2, 4, 6]
    # xticks_epg[key, 4] = [-0.3, 0, 0.3]
    # xticks_epg[key, 5] = [-0.01, 0, 0.01]
    # xticks_epg[key, 6] = [-0.3, 0, 0.3]
    # xticks_epg[key, 7] = [0, 0.04, 0.08]
    # yticks_epg[key, 0] = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12]
    # yticks_epg[key, 1] = [0.5, 1, 1.5, 2]
    # yticks_epg[key, 2] = [1, 2, 3]
    # yticks_epg[key, 3] = [0.2, 0.4]
    # yticks_epg[key, 4] = [1, 2, 3]
    # yticks_epg[key, 5] = [50, 100, 150, 200, 250]
    # yticks_epg[key, 6] = [1,2,3]
    # yticks_epg[key, 7] = [5, 10, 15, 20, 25]


  # key = ('inb', 1, 'a')
  # for key in [('inb', 1, 'a'), ('inb', 2, 'a'), ('inb', 3, 'a')]:
  #   xlbs_epgg[key] = [0.1,    0, 0, 0, -.4, -0.02, -0.5, 0]
  #   xubs_epgg[key] = [0.16,  1.6, 1, 8,  .4,  0.02, 0.5, 0.15]
  #   yubs_epgg[key] = [50, 2, 2.5, 0.5, 3, 250, 4, 25]

    # xticks_epgg[key, 1] = [0, 0.4, 0.8, 1.2, 1.6]
    # xticks_epgg[key, 2] = [0, 0.25, 0.5 ,0.75, 1]
    # xticks_epgg[key, 3] = [0, 2, 4, 6, 8]
    # xticks_epgg[key, 4] = [-0.4, -0.2, 0, 0.2, 0.4]
    # xticks_epgg[key, 5] = [-0.02, 0, 0.02]
    # xticks_epgg[key, 6] = [-0.5, 0, 0.5]
    # xticks_epgg[key, 7] = [0, 0.05, 0.1, 0.15]
    # yticks_epgg[key, 0] = [50, 100]
    # yticks_epgg[key, 1] = [0.5, 1, 1.5, 2]
    # yticks_epgg[key, 2] = [.5, 1, 1.5, 2, 2.5]
    # yticks_epgg[key, 3] = [0.25, 0.5]
    # yticks_epgg[key, 4] = [1, 2, 3]
    # yticks_epgg[key, 5] = [50, 100, 150, 200]
    # yticks_epgg[key, 6] = [2, 4]
    # yticks_epgg[key, 7] = [5, 10, 15, 20, 25]

  # path_exp_epg  = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/exp/epg/{}/'.format(args.plot, polarity)
  # path_exp_epgg = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/exp/epgg/{}/'.format(args.plot, polarity)
  # path_sim_epg  = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/sim/epg/{}/'.format(args.plot, polarity)
  # path_sim_epgg = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/sim/epgg/{}/'.format(args.plot, polarity)
  topo   = {1: "FD", 2:"CD", 3: "CDFT"}
  
  path_exp_epg  = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/exp/epg/{}/'.format(args.plot, polarity)
  path_exp_epgg = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/exp/epgg/{}/'.format(args.plot, polarity)
  path_sim_epg  = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/sim/epg/{}/'.format(args.plot, polarity)
  path_sim_bkg_1g = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/sim/bkg_1g/{}/'.format(args.plot, polarity)
  path_sim_bkg_2g = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/sim/bkg_2g/{}/'.format(args.plot, polarity)

  dsets = []
  for path in [path_exp_epg, path_exp_epgg, path_sim_epg, path_sim_bkg_1g, path_sim_bkg_2g]:
    dfs = []
    for file in glob('{}/*.pkl'.format(path)):
      dfs.append(pd.read_pickle(file))
    dfs = pd.concat(dfs)
    dsets.append(dfs)

  exp_epg, exp_epgg, sim_epg, sim_bkg_1g, sim_bkg_2g = dsets

  exp_label = 'Experimental Data'
  sim_label = 'Simulation'
  
  if args.plot == 'a':
    sim_label = 'Simulation\n (dvcsgen only)'

  for config in [1, 2, 3]:
    df3 = exp_epg.loc[exp_epg.config == config, :]
    df1 = sim_epg.loc[sim_epg.config == config, :]
    df2 = sim_bkg_1g.loc[sim_bkg_1g.config == config, :]
    df4 = exp_epgg.loc[exp_epgg.config == config, :]
    df5 = sim_bkg_2g.loc[sim_bkg_2g.config == config, :]

    # DVCS plots
    varstoplot = dvcsvars
    title = dvcstitles
    unit = dvcsunits
    key  = (polarity, config, args.plot)

    fig, axs = plt.subplots(2, 4, figsize = (16, 8))
    for yind, xind in itertools.product(range(2), range(4)):
      ind = 4*yind + xind

      xlb = xlbs_epg[key][ind]
      xub = xubs_epg[key][ind]
      yub = yubs_epg[key][ind]

      bins = np.linspace(xlb, xub, 101)

      simDist_dvcs, _  = np.histogram(df1.loc[:, varstoplot[ind]], bins, density = True)
      simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)

      axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = exp_label)
      if args.plot == 'a':
        simDist = simDist_dvcs
      else:
        cont             = len(df4)*len(df2)/len(df5)/len(df3)
        simDist = (1-cont)*simDist_dvcs + cont*simDist_dvpi0
      axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = sim_label)
      axs[yind, xind].set_title(title[ind])
      if (unit[ind]):
        axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
      else:
        axs[yind, xind].set_xlabel(title[ind])
      axs[yind, xind].set_xlim([xlb, xub])
      axs[yind, xind].set_ylim([0, yub])
      if (key, ind) in xticks_epg.keys():
        xtick = xticks_epg[key, ind]
        axs[yind, xind].set_xticks(xtick)
      if (key, ind) in yticks_epg.keys():
        ytick = yticks_epg[key, ind]
        axs[yind, xind].set_yticks(ytick)

    handles, labels = axs[0, 0].get_legend_handles_labels()
    lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
    plt.tight_layout()
    plt.savefig("plots/q16/dset_{}/{}/dvcs{}{}excl.pdf".format(args.plot, polarity, pol_cap, topo[config]), bbox_extra_artists=[lgd], bbox_inches = 'tight')
    plt.clf()

    # DVπ0P plots
    varstoplot = pi0vars
    title = pi0titles
    unit = pi0units
    key  = (polarity, config, args.plot)

    fig, axs = plt.subplots(2, 4, figsize = (16, 8))
    for yind, xind in itertools.product(range(2), range(4)):
      ind = 4*yind + xind

      xlb = xlbs_epgg[key][ind]
      xub = xubs_epgg[key][ind]
      yub = yubs_epgg[key][ind]
      if args.plot == 'a':
        bins = np.linspace(xlb, xub, 41)
      else:
        bins = np.linspace(xlb, xub, 101)

      simDist, _ = np.histogram(df5[varstoplot[ind]], bins, density = True)

      axs[yind, xind].hist(df4[varstoplot[ind]], bins = bins, histtype = 'step', edgecolor='b', density=True, linewidth=1, label = exp_label)
      axs[yind, xind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=1, label = sim_label)
      axs[yind, xind].set_title(title[ind])
      if (unit[ind]):
        axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
      else:
        axs[yind, xind].set_xlabel(title[ind])
      axs[yind, xind].set_xlim([xlb, xub])
      axs[yind, xind].set_ylim([0, yub])

      if (key, ind) in xticks_epgg.keys():
        xtick = xticks_epgg[key, ind]
        axs[yind, xind].set_xticks(xtick)
      if (key, ind) in yticks_epgg.keys():
        ytick = yticks_epgg[key, ind]
        axs[yind, xind].set_yticks(ytick)

    handles, labels = axs[0, 0].get_legend_handles_labels()
    lgd = plt.figlegend(handles,labels, loc='center left', fontsize= 30, title_fontsize = 30, bbox_to_anchor = (1.0, 0.5))
    plt.tight_layout()
    plt.savefig("plots/q16/dset_{}/{}/pi0{}{}excl.pdf".format(args.plot, polarity, pol_cap, topo[config]), bbox_extra_artists=[lgd], bbox_inches = 'tight')
    plt.clf()

  # dfs = []
  # topo   = {1: "FD", 2:"CD", 3: "CDFT"}

  # for config in [1, 2, 3]:
  #   for varind in range(8):
  #     var    = dvcsvars[varind]
  #     label  = dvcstitles[varind]
  #     unit   = dvcsunits[varind]
  #     if var in ["MPt", "reconGam", "coplanarity"]:
  #       lb        = 0
  #       if config == 1:
  #         ub        = cuts_dvcs_FD_Inb_3sigma["{}_ub".format(var)]
  #       if config == 2:
  #         ub        = cuts_dvcs_CD_Inb_3sigma["{}_ub".format(var)]
  #       if config == 3:
  #         ub        = cuts_dvcs_CDFT_Inb_3sigma["{}_ub".format(var)]
  #     elif var == "coneAngle":
  #       if config == 1:
  #         lb        = 25
  #         ub        = 45
  #       if config == 2:
  #         lb        = 10
  #         ub        = 30
  #       if config == 3:
  #         lb        = 10
  #         ub        = 30
  #     else:
  #       if config == 1:
  #         lb        = cuts_dvcs_FD_Inb_3sigma["{}_lb".format(var)]
  #         ub        = cuts_dvcs_FD_Inb_3sigma["{}_ub".format(var)]
  #       if config == 2:
  #         lb        = cuts_dvcs_CD_Inb_3sigma["{}_lb".format(var)]
  #         ub        = cuts_dvcs_CD_Inb_3sigma["{}_ub".format(var)]
  #       if config == 3:
  #         lb        = cuts_dvcs_CDFT_Inb_3sigma["{}_lb".format(var)]
  #         ub        = cuts_dvcs_CDFT_Inb_3sigma["{}_ub".format(var)]
  #     bins          = np.linspace(lb, ub, 21)
  #     bincenters    = (bins[1:]+bins[:-1])/2.
  #     ind = 0 
  #     filenum = 0
  #     fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
  #     for tbin in range(len(tbins)-1):
  #       for xBbin in range(len(newxBbins2)-1):
  #         for Q2bin in range(len(newxBbins2)-1):
  #           df_bh = bh_inb.loc[(bh_inb.t1>newtbins[tbin]) & (bh_inb.t1<newtbins[tbin+1]), :]
  #           df_vgg = vgg_inb.loc[(vgg_inb.t1>newtbins[tbin]) & (vgg_inb.t1<newtbins[tbin+1]), :]
  #           df_bkg_1g = bkg_1g_inb.loc[(bkg_1g_inb.t1>newtbins[tbin]) & (bkg_1g_inb.t1<newtbins[tbin+1]), :]
  #           df_bkg_2g = bkg_2g_inb.loc[(bkg_2g_inb.t1>newtbins[tbin]) & (bkg_2g_inb.t1<newtbins[tbin+1]), :]
  #           df_epg = epg_inb.loc[(epg_inb.t1>newtbins[tbin]) & (epg_inb.t1<newtbins[tbin+1]), :]
  #           df_epgg = epgg_inb.loc[(epgg_inb.t1>newtbins[tbin]) & (epgg_inb.t1<newtbins[tbin+1]), :]
  #           df_bh = df_bh.loc[(df_bh.xB>newxBbins2[xBbin]) & (df_bh.xB<newxBbins2[xBbin+1]), :]
  #           df_bh = df_bh.loc[(df_bh.Q2>newQ2bins2[Q2bin]) & (df_bh.Q2<newQ2bins2[Q2bin+1]), :]
  #           df_vgg = df_vgg.loc[(df_vgg.xB>newxBbins2[xBbin]) & (df_vgg.xB<newxBbins2[xBbin+1]), :]
  #           df_vgg = df_vgg.loc[(df_vgg.Q2>newQ2bins2[Q2bin]) & (df_vgg.Q2<newQ2bins2[Q2bin+1]), :]
  #           df_bkg_1g = df_bkg_1g.loc[(df_bkg_1g.xB>newxBbins2[xBbin]) & (df_bkg_1g.xB<newxBbins2[xBbin+1]), :]
  #           df_bkg_1g = df_bkg_1g.loc[(df_bkg_1g.Q2>newQ2bins2[Q2bin]) & (df_bkg_1g.Q2<newQ2bins2[Q2bin+1]), :]
  #           df_bkg_2g = df_bkg_2g.loc[(df_bkg_2g.xB>newxBbins2[xBbin]) & (df_bkg_2g.xB<newxBbins2[xBbin+1]), :]
  #           df_bkg_2g = df_bkg_2g.loc[(df_bkg_2g.Q2>newQ2bins2[Q2bin]) & (df_bkg_2g.Q2<newQ2bins2[Q2bin+1]), :]
  #           df_epg = df_epg.loc[(df_epg.xB>newxBbins2[xBbin]) & (df_epg.xB<newxBbins2[xBbin+1]), :]
  #           df_epg = df_epg.loc[(df_epg.Q2>newQ2bins2[Q2bin]) & (df_epg.Q2<newQ2bins2[Q2bin+1]), :]
  #           df_epgg = df_epgg.loc[(df_epgg.xB>newxBbins2[xBbin]) & (df_epgg.xB<newxBbins2[xBbin+1]), :]
  #           df_epgg = df_epgg.loc[(df_epgg.Q2>newQ2bins2[Q2bin]) & (df_epgg.Q2<newQ2bins2[Q2bin+1]), :]
  #           bh_var     = df_bh.loc[df_bh.config == config, var]
  #           bkg_1g_var = df_bkg_1g.loc[df_bkg_1g.config == config, var]
  #           vgg_var    = df_vgg.loc[df_vgg.config == config, var]
  #           epg_var    = df_epg.loc[df_epg.config == config, var]
            
  #           if len(epg_var)>20:
  #             pass
  #           else:
  #             continue
            
  #           try:
  #             cont_thisbin = len(bkg_1g_var)/len(df_bkg_2g.loc[df_bkg_2g.config==config, :])*len(df_epgg.loc[df_epgg.config==config, :])/len(epg_var)
  #           except:
  #             cont_thisbin = 0
            
  #           xind = ind//5
  #           yind = ind%5
  #           ind = ind + 1
            
  #           bh_hist, _    = np.histogram(bh_var, bins = bins, density = True)
  #           vgg_hist, _   = np.histogram(vgg_var, bins = bins, density = True)
  #           epg_hist, _   = np.histogram(epg_var, bins = bins, density = True)
  #           if cont_thisbin > 0:
  #             bkg_1g_hist, _= np.histogram(bkg_1g_var, bins = bins, density = True)
  #             if np.isnan(bkg_1g_hist).any():
  #               bkg_1g_hist = np.zeros(len(bh_hist))
  #           else:
  #             bkg_1g_hist = np.zeros(len(bh_hist))
  #           sim_hist   = cont_thisbin*bkg_1g_hist + (1-cont_thisbin)*bh_hist
  #           sim_hist2   = cont_thisbin*bkg_1g_hist + (1-cont_thisbin)*vgg_hist
                       
            
  #           axs[xind, yind].step(bincenters, epg_hist, where = 'mid', color = 'tab:blue')
  #           axs[xind, yind].step(bincenters, sim_hist, where = 'mid', color = 'tab:red')
  #           axs[xind, yind].step(bincenters, sim_hist2, where = 'mid', color = 'tab:green')
  #           annotation = "({}, {}, {})".format(xBbin, Q2bin, tbin)
  #           axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)
  #           axs[xind, yind].set_ylim(bottom = 0)
            
            
  #           if ind >24:
  #             for ax in axs.flatten():
  #               ax.get_xaxis().set_visible(False)
  #               ax.get_yaxis().set_visible(False)
  #             for ax in axs[-1, :]:
  #               ax.get_xaxis().set_visible(True)
  #               xtick        = np.linspace(lb, ub, 5)
  #               xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
  #               xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
  #               xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
  #               xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
  #               xticklabel_1[-1] = ''
  #               xticklabel_2[-1] = ''
  #               xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
  #               ax.set_xlabel( label + " [" + unit + "]")
  #               ax.set_xticks(xtick)
  #               ax.set_xticklabels(xticklabel_1)
  #             axs[-1, 0] .set_xticklabels(xticklabel_2)
  #             axs[-1, -1].set_xticklabels(xticklabel_3)
                
  #             plt.subplots_adjust(wspace=0, hspace=0)
  #             plt.savefig("plots/q16/inb_{}_{}_{}.pdf".format(topo[config], var, filenum))
  #             filenum = filenum + 1
  #             plt.clf()
  #             ind = 0
  #             fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
            
  #     for ax in axs.flatten():
  #       ax.get_xaxis().set_visible(False)
  #       ax.get_yaxis().set_visible(False)
  #     for ax in axs[-1, :]:
  #       ax.get_xaxis().set_visible(True)
  #       xtick        = np.linspace(lb, ub, 5)
  #       xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
  #       xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
  #       xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
  #       xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
  #       xticklabel_1[-1]  = ''
  #       xticklabel_2[-1] = ''
  #       xticklabel_3[0] = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
  #       ax.set_xlabel( label + " [" + unit + "]")
  #       ax.set_xticks(xtick)
  #       ax.set_xticklabels(xticklabel_1)
  #       ax.set_xlim([lb, ub])
  #     axs[-1, 0].set_xticklabels(xticklabel_2)
  #     axs[-1, -1].set_xticklabels(xticklabel_3)
      
  #     plt.subplots_adjust(wspace=0, hspace=0)
  #     plt.savefig("plots/q16/inb_{}_{}_{}.pdf".format(topo[config], var, filenum))
  #     plt.clf()
