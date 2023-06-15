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

  if args.plot == 'a':
    print("Plot at ifarm.")
    exit()

  if args.polarity == 'inbending':
    polarity = 'inb'
    pol_cap  = 'Inb'
  elif args.polarity == 'outbending':
    polarity = 'outb'
    pol_cap  = 'Outb'
  else:
    print("Polarity undefined")
    exit()

  # path_exp_epg  = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/exp/epg/{}/'.format(args.plot, polarity)
  # path_exp_epgg = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/exp/epgg/{}/'.format(args.plot, polarity)
  # path_sim_epg  = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/sim/epg/{}/'.format(args.plot, polarity)
  # path_sim_epgg = '/volatile/clas12/sangbaek/jan2023/q16_18/dset_{}/sim/epgg/{}/'.format(args.plot, polarity)
  topo   = {1: "FD", 2:"CD", 3: "CDFT"}
  
  plot_ver = args.plot
  if plot_ver == 'e':
    plot_ver = 'd'
  path_exp_epg  = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/exp/epg/{}/'.format(plot_ver, polarity)
  path_exp_epgg = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/exp/epgg/{}/'.format(plot_ver, polarity)
  path_sim_epg  = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/sim/epg/{}/'.format(plot_ver, polarity)
  path_sim_bkg_1g = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/sim/bkg_1g/{}/'.format(plot_ver, polarity)
  path_sim_bkg_2g = '/Users/sangbaek.lee/CLAS12/q16_18/dset_{}/sim/bkg_2g/{}/'.format(plot_ver, polarity)

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

  # if args.plot in ['a', 'b', 'c', 'd']:
  if args.plot in ['b', 'c', 'd']:
    for config in [1, 2, 3]:
      key  = (polarity, config, args.plot)

      df3 = exp_epg   .loc[(exp_epg.config == config    ) & (exp_epg   .Ptheta < 64.23 )   , :]
      df1 = sim_epg   .loc[(sim_epg.config == config    ) & (sim_epg   .Ptheta < 64.23 )   , :]
      df2 = sim_bkg_1g.loc[(sim_bkg_1g.config == config ) & (sim_bkg_1g.Ptheta < 64.23 )   , :]
      df4 = exp_epgg  .loc[(exp_epgg.config == config   ) & (exp_epgg  .Ptheta < 64.23 )   , :]
      df5 = sim_bkg_2g.loc[(sim_bkg_2g.config == config ) & (sim_bkg_2g.Ptheta < 64.23 )   , :]

      # DVCS plots
      varstoplot = dvcsvars
      title = dvcstitles
      unit = dvcsunits

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

  elif args.plot == 'e':
    for config in [1, 2, 3]:
      key  = (polarity, config, args.plot)
      for varind in range(8):

        # DVCS plots
        var    = dvcsvars[varind]
        label  = dvcstitles[varind]
        unit   = dvcsunits[varind]

        xlb     = xlbs_epg[key][varind]
        xub     = xubs_epg[key][varind]

        ind = 0 
        filenum = 0
        
        fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
        for tbin in range(6):
          for xBbin in range(len(newxBbins2)-1):
            for Q2bin in range(len(newQ2bins2)-1):
              df3 = exp_epg   .loc[(exp_epg   .config == config ) & (exp_epg   .xBbin == xBbin ) & (exp_epg   .Q2bin == Q2bin ) & (exp_epg   .tbin == tbin ) & (exp_epg   .Ptheta < 64.23 )   , :]
              df1 = sim_epg   .loc[(sim_epg   .config == config ) & (sim_epg   .xBbin == xBbin ) & (sim_epg   .Q2bin == Q2bin ) & (sim_epg   .tbin == tbin ) & (sim_epg   .Ptheta < 64.23 )   , :]
              df2 = sim_bkg_1g.loc[(sim_bkg_1g.config == config ) & (sim_bkg_1g.xBbin == xBbin ) & (sim_bkg_1g.Q2bin == Q2bin ) & (sim_bkg_1g.tbin == tbin ) & (sim_bkg_1g.Ptheta < 64.23 )   , :]
              df4 = exp_epgg  .loc[(exp_epgg  .config == config ) & (exp_epgg  .xBbin == xBbin ) & (exp_epgg  .Q2bin == Q2bin ) & (exp_epgg  .tbin == tbin ) & (exp_epgg  .Ptheta < 64.23 )   , :]
              df5 = sim_bkg_2g.loc[(sim_bkg_2g.config == config ) & (sim_bkg_2g.xBbin == xBbin ) & (sim_bkg_2g.Q2bin == Q2bin ) & (sim_bkg_2g.tbin == tbin ) & (sim_bkg_2g.Ptheta < 64.23 )   , :]
            
              if len(df3)>50:
                pass
              else:
                continue

              expDist     , bins  = np.histogram(df3.loc[:, var],   20, density = True)
              simDist_dvcs, _     = np.histogram(df1.loc[:, var], bins, density = True)
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
              yub  = 1.2*np.max([expDist, simDist])
                                          
              axs[xind, yind].hist(bins[:-1], bins, weights = expDist, histtype = 'step', color='b', linewidth=3, label = exp_label)
              axs[xind, yind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=3, label = sim_label)

              annotation = "({}, {}, {})".format(xBbin, Q2bin, tbin)
              axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)
              axs[xind, yind].set_ylim(bottom = 0)
              axs[xind, yind].set_ylim(top    = yub)
              
              if ind >24:
                for ax in axs.flatten():
                  ax.get_xaxis().set_visible(False)
                  ax.get_yaxis().set_visible(False)
                for ax in axs[-1, :]:
                  ax.get_xaxis().set_visible(True)
                  xtick        = xticks_epg[key, varind]
                  xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
                  xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
                  xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
                  xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
                  xticklabel_1[-1] = ''
                  xticklabel_2[-1] = ''
                  xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
                  ax.set_xlabel(label + " [" + unit + "]", fontsize = 40)
                  ax.set_xticks(xtick)
                  ax.set_xticklabels(xticklabel_1)
                  ax.set_xlim([xlb, xub])
                axs[-1, 0] .set_xticklabels(xticklabel_2)
                axs[-1, 0] .set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
                axs[-1, -1].set_xticklabels(xticklabel_3)
                  
                plt.subplots_adjust(wspace=0, hspace=0)
                plt.savefig("plots/q16/dset_e/{}/dvcs_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
                filenum = filenum + 1
                plt.clf()
                ind = 0
                fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
              
        for ax in axs.flatten():
          ax.get_xaxis().set_visible(False)
          ax.get_yaxis().set_visible(False)
        for ax in axs[-1, :]:
          ax.get_xaxis().set_visible(True)
          xtick        = xticks_epg[key, varind]
          xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
          xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
          xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
          xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
          xticklabel_1[-1]  = ''
          xticklabel_2[-1] = ''
          xticklabel_3[0] = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
          ax.set_xlabel( label + " [" + unit + "]", fontsize = 40)
          ax.set_xticks(xtick)
          ax.set_xticklabels(xticklabel_1)
          ax.set_xlim([xlb, xub])
        axs[-1, 0] .set_xticklabels(xticklabel_2)
        axs[-1, 0] .set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
        axs[-1, -1].set_xticklabels(xticklabel_3)
        
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.savefig("plots/q16/dset_e/{}/dvcs_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
        plt.clf()


        # DVπ0P plots
        var   = pi0vars[varind]
        label = pi0titles[varind]
        unit  = pi0units[varind]

        xlb     = xlbs_epgg[key][varind]
        xub     = xubs_epgg[key][varind]

        ind = 0 
        filenum = 0
        
        fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
        for tbin in range(len(tbins)-1):
          for xBbin in range(len(newxBbins2)-1):
            for Q2bin in range(len(newQ2bins2)-1):
              df4 = exp_epgg  .loc[(exp_epgg  .config == config ) & (exp_epgg  .xBbin == xBbin ) & (exp_epgg  .Q2bin == Q2bin ) & (exp_epgg  .tbin == tbin ) & (exp_epgg  .Ptheta < 64.23 )   , :]
              df5 = sim_bkg_2g.loc[(sim_bkg_2g.config == config ) & (sim_bkg_2g.xBbin == xBbin ) & (sim_bkg_2g.Q2bin == Q2bin ) & (sim_bkg_2g.tbin == tbin ) & (sim_bkg_2g.Ptheta < 64.23 )   , :]
            
              if len(df4)>20:
                pass
              else:
                continue

              expDist     , bins  = np.histogram(df4.loc[:, var],   20, density = True)
              simDist     , _     = np.histogram(df5.loc[:, var], bins, density = True)
              
              xind = ind//5
              yind = ind%5
              ind  = ind + 1
              yub  = 1.2*np.max([expDist, simDist])
                                          
              axs[xind, yind].hist(bins[:-1], bins, weights = expDist, histtype = 'step', color='b', linewidth=3, label = exp_label)
              axs[xind, yind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=3, label = sim_label)

              annotation = "({}, {}, {})".format(xBbin, Q2bin, tbin)
              axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)
              axs[xind, yind].set_ylim(bottom = 0)
              axs[xind, yind].set_ylim(top    = yub)
              
              if ind >24:
                for ax in axs.flatten():
                  ax.get_xaxis().set_visible(False)
                  ax.get_yaxis().set_visible(False)
                for ax in axs[-1, :]:
                  ax.get_xaxis().set_visible(True)
                  xtick        = xticks_epgg[key, varind]
                  xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
                  xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
                  xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
                  xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
                  xticklabel_1[-1] = ''
                  xticklabel_2[-1] = ''
                  xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
                  ax.set_xlabel(label + " [" + unit + "]", fontsize = 40)
                  ax.set_xticks(xtick)
                  ax.set_xticklabels(xticklabel_1)
                  ax.set_xlim([xlb, xub])
                axs[-1, 0] .set_xticklabels(xticklabel_2)
                axs[-1, 0] .set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
                axs[-1, -1].set_xticklabels(xticklabel_3)
                  
                plt.subplots_adjust(wspace=0, hspace=0)
                plt.savefig("plots/q16/dset_e/{}/pi0_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
                filenum = filenum + 1
                plt.clf()
                ind = 0
                fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))
              
        for ax in axs.flatten():
          ax.get_xaxis().set_visible(False)
          ax.get_yaxis().set_visible(False)
        for ax in axs[-1, :]:
          ax.get_xaxis().set_visible(True)
          xtick        = xticks_epgg[key, varind]
          xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
          xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
          xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
          xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
          xticklabel_1[-1]  = ''
          xticklabel_2[-1] = ''
          xticklabel_3[0] = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])
          ax.set_xlabel( label + " [" + unit + "]", fontsize = 40)
          ax.set_xticks(xtick)
          ax.set_xticklabels(xticklabel_1)
          ax.set_xlim([xlb, xub])
        axs[-1, 0].set_xticklabels(xticklabel_2)
        axs[-1, 0] .set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
        axs[-1, -1].set_xticklabels(xticklabel_3)
        
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.savefig("plots/q16/dset_e/{}/pi0_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
        plt.clf()
