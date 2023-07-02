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

if args.question == 0:

  InbExp = []
  for file in glob("/volatile/clas12/sangbaek/jan2023/convPkl_nofid/exp/epg/inb/*.pkl"):
    df = pd.read_pickle(file)
    InbExp.append(df)
  InbExp = pd.concat(InbExp)
  InbSim = []
  for file in glob("/volatile/clas12/sangbaek/jan2023/convPkl_nofid/bh/inb/*.pkl"):
    df = pd.read_pickle(file)
    InbSim.append(df)
  InbSim = pd.concat(InbSim)

  OutbExp = []
  for file in glob("/volatile/clas12/sangbaek/jan2023/convPkl_nofid/exp/epg/outb/*.pkl"):
    df = pd.read_pickle(file)
    OutbExp.append(df)
  OutbExp = pd.concat(OutbExp)
  OutbSim = []
  for file in glob("/volatile/clas12/sangbaek/jan2023/convPkl_nofid/bh/outb/*.pkl"):
    df = pd.read_pickle(file)
    OutbSim.append(df)
  OutbSim = pd.concat(OutbSim)

  print("duplicated numbers Inb. Exp.: {} out of total {}".format(len(InbExp) - len(InbExp.event.unique()), len(InbExp)))
  print("duplicated numbers Inb. Sim.: {} out of total {}".format(len(InbSim) - len(InbSim.event.unique()), len(InbSim)))
  print("duplicated numbers Outb. Exp.: {} out of total {}".format(len(OutbExp) - len(OutbExp.event.unique()), len(OutbExp)))
  print("duplicated numbers Outb. Sim.: {} out of total {}".format(len(OutbSim) - len(OutbSim.event.unique()), len(OutbSim)))

  electronFiducialCounting(InbExp, pol = "inbending", mc = False, fidlevel = 'mid')
  protonFiducialCounting(InbExp, pol = "inbending")
  gammaFiducialCounting(InbExp)
  electronFiducialCounting(InbSim, pol = "inbending", mc = True, fidlevel = 'mid')
  protonFiducialCounting(InbSim, pol = "inbending")
  gammaFiducialCounting(InbSim)
  electronFiducialCounting(OutbExp, pol = "outbending", mc = False, fidlevel = 'mid')
  protonFiducialCounting(OutbExp, pol = "outbending")
  gammaFiducialCounting(OutbExp)
  electronFiducialCounting(OutbSim, pol = "outbending", mc = True, fidlevel = 'mid')
  protonFiducialCounting(OutbSim, pol = "outbending")
  gammaFiducialCounting(OutbSim)

  InbExp_effect_Epcal = len(InbExp.loc[InbExp.EFid_pcal>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Edw = len(InbExp.loc[InbExp.EFid_dw>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Edc = len(InbExp.loc[InbExp.EFid_dc>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Esf = len(InbExp.loc[InbExp.EFid_sf>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Evz = len(InbExp.loc[InbExp.EFid_vz>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Eedep = len(InbExp.loc[InbExp.EFid_edep>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Eap = len(InbExp.loc[InbExp.EFid_ap>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Pdc = len(InbExp.loc[(InbExp.Psector < 7) & (InbExp.PFid_dc>0), "event"].unique())/len(InbExp.loc[InbExp.Psector < 7, "event"].unique())
  InbExp_effect_Pcvt = len(InbExp.loc[(InbExp.Psector > 7) & (InbExp.PFid_cvt>0), "event"].unique())/len(InbExp.loc[InbExp.Psector > 7, "event"].unique())
  InbExp_effect_Pchi = len(InbExp.loc[InbExp.PFid_chi>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Pvz = len(InbExp.loc[InbExp.PFid_vz>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Gbeta = len(InbExp.loc[InbExp.GFid_beta>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Gpcal = len(InbExp.loc[(InbExp.Gsector < 7) & (InbExp.GFid_Pcal>0), "event"].unique())/len(InbExp.loc[InbExp.Gsector < 7, "event"].unique())
  InbExp_effect_Gpcal2 = len(InbExp.loc[(InbExp.Gsector < 7) & (InbExp.GFid_Pcal2>0), "event"].unique())/len(InbExp.loc[InbExp.Gsector < 7, "event"].unique())
  InbExp_effect_Gft = len(InbExp.loc[(InbExp.Gsector > 7) & (InbExp.GFid_FT>0), "event"].unique())/len(InbExp.loc[InbExp.Gsector > 7, "event"].unique())
  InbExp_effect_FDFD = len(InbExp.loc[(InbExp.EFid*InbExp.PFid*InbExp.GFid > 0) & (InbExp.config==1), "event"].unique())/sum(InbExp.config==1)
  InbExp_effect_CDFD = len(InbExp.loc[(InbExp.EFid*InbExp.PFid*InbExp.GFid > 0) & (InbExp.config==2), "event"].unique())/sum(InbExp.config==2)
  InbExp_effect_CDFT = len(InbExp.loc[(InbExp.EFid*InbExp.PFid*InbExp.GFid > 0) & (InbExp.config==3), "event"].unique())/sum(InbExp.config==3)

  InbSim_effect_Epcal = len(InbSim.loc[InbSim.EFid_pcal>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Edw = len(InbSim.loc[InbSim.EFid_dw>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Edc = len(InbSim.loc[InbSim.EFid_dc>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Esf = len(InbSim.loc[InbSim.EFid_sf>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Evz = len(InbSim.loc[InbSim.EFid_vz>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Eedep = len(InbSim.loc[InbSim.EFid_edep>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Eap = len(InbSim.loc[InbSim.EFid_ap>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Pdc = len(InbSim.loc[(InbSim.Psector < 7) & (InbSim.PFid_dc>0), "event"].unique())/len(InbSim.loc[InbSim.Psector < 7, "event"].unique())
  InbSim_effect_Pcvt = len(InbSim.loc[(InbSim.Psector > 7) & (InbSim.PFid_cvt>0), "event"].unique())/len(InbSim.loc[InbSim.Psector > 7, "event"].unique())
  InbSim_effect_Pchi = len(InbSim.loc[InbSim.PFid_chi>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Pvz = len(InbSim.loc[InbSim.PFid_vz>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Gbeta = len(InbSim.loc[InbSim.GFid_beta>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Gpcal = len(InbSim.loc[(InbSim.Gsector < 7) & (InbSim.GFid_Pcal>0), "event"].unique())/len(InbSim.loc[InbSim.Gsector < 7, "event"].unique())
  InbSim_effect_Gpcal2 = len(InbSim.loc[(InbSim.Gsector < 7) & (InbSim.GFid_Pcal2>0), "event"].unique())/len(InbSim.loc[InbSim.Gsector < 7, "event"].unique())
  InbSim_effect_Gft = len(InbSim.loc[(InbSim.Gsector > 7) & (InbSim.GFid_FT>0), "event"].unique())/len(InbSim.loc[InbSim.Gsector > 7, "event"].unique())
  InbSim_effect_FDFD = len(InbSim.loc[(InbSim.EFid*InbSim.PFid*InbSim.GFid > 0) & (InbSim.config==1), "event"].unique())/sum(InbSim.config==1)
  InbSim_effect_CDFD = len(InbSim.loc[(InbSim.EFid*InbSim.PFid*InbSim.GFid > 0) & (InbSim.config==2), "event"].unique())/sum(InbSim.config==2)
  InbSim_effect_CDFT = len(InbSim.loc[(InbSim.EFid*InbSim.PFid*InbSim.GFid > 0) & (InbSim.config==3), "event"].unique())/sum(InbSim.config==3)

  OutbExp_effect_Epcal = len(OutbExp.loc[OutbExp.EFid_pcal>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Edw = len(OutbExp.loc[OutbExp.EFid_dw>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Edc = len(OutbExp.loc[OutbExp.EFid_dc>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Esf = len(OutbExp.loc[OutbExp.EFid_sf>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Evz = len(OutbExp.loc[OutbExp.EFid_vz>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Eedep = len(OutbExp.loc[OutbExp.EFid_edep>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Eap = len(OutbExp.loc[OutbExp.EFid_ap>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Pdc = len(OutbExp.loc[(OutbExp.Psector < 7) & (OutbExp.PFid_dc>0), "event"].unique())/len(OutbExp.loc[OutbExp.Psector < 7, "event"].unique())
  OutbExp_effect_Pcvt = len(OutbExp.loc[(OutbExp.Psector > 7) & (OutbExp.PFid_cvt>0), "event"].unique())/len(OutbExp.loc[OutbExp.Psector > 7, "event"].unique())
  OutbExp_effect_Pchi = len(OutbExp.loc[OutbExp.PFid_chi>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Pvz = len(OutbExp.loc[OutbExp.PFid_vz>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Gbeta = len(OutbExp.loc[OutbExp.GFid_beta>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Gpcal = len(OutbExp.loc[(OutbExp.Gsector < 7) & (OutbExp.GFid_Pcal>0), "event"].unique())/len(OutbExp.loc[OutbExp.Gsector < 7, "event"].unique())
  OutbExp_effect_Gpcal2 = len(OutbExp.loc[(OutbExp.Gsector < 7) & (OutbExp.GFid_Pcal2>0), "event"].unique())/len(OutbExp.loc[OutbExp.Gsector < 7, "event"].unique())
  OutbExp_effect_Gft = len(OutbExp.loc[(OutbExp.Gsector > 7) & (OutbExp.GFid_FT>0), "event"].unique())/len(OutbExp.loc[OutbExp.Gsector > 7, "event"].unique())
  OutbExp_effect_FDFD = len(OutbExp.loc[(OutbExp.EFid*OutbExp.PFid*OutbExp.GFid > 0) & (OutbExp.config==1), "event"].unique())/sum(OutbExp.config==1)
  OutbExp_effect_CDFD = len(OutbExp.loc[(OutbExp.EFid*OutbExp.PFid*OutbExp.GFid > 0) & (OutbExp.config==2), "event"].unique())/sum(OutbExp.config==2)
  OutbExp_effect_CDFT = len(OutbExp.loc[(OutbExp.EFid*OutbExp.PFid*OutbExp.GFid > 0) & (OutbExp.config==3), "event"].unique())/sum(OutbExp.config==3)

  OutbSim_effect_Epcal = len(OutbSim.loc[OutbSim.EFid_pcal>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Edw = len(OutbSim.loc[OutbSim.EFid_dw>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Edc = len(OutbSim.loc[OutbSim.EFid_dc>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Esf = len(OutbSim.loc[OutbSim.EFid_sf>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Evz = len(OutbSim.loc[OutbSim.EFid_vz>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Eedep = len(OutbSim.loc[OutbSim.EFid_edep>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Eap = len(OutbSim.loc[OutbSim.EFid_ap>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Pdc = len(OutbSim.loc[(OutbSim.Psector < 7) & (OutbSim.PFid_dc>0), "event"].unique())/len(OutbSim.loc[OutbSim.Psector < 7, "event"].unique())
  OutbSim_effect_Pcvt = len(OutbSim.loc[(OutbSim.Psector > 7) & (OutbSim.PFid_cvt>0), "event"].unique())/len(OutbSim.loc[OutbSim.Psector > 7, "event"].unique())
  OutbSim_effect_Pchi = len(OutbSim.loc[OutbSim.PFid_chi>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Pvz = len(OutbSim.loc[OutbSim.PFid_vz>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Gbeta = len(OutbSim.loc[OutbSim.GFid_beta>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Gpcal = len(OutbSim.loc[(OutbSim.Gsector < 7) & (OutbSim.GFid_Pcal>0), "event"].unique())/len(OutbSim.loc[OutbSim.Gsector < 7, "event"].unique())
  OutbSim_effect_Gpcal2 = len(OutbSim.loc[(OutbSim.Gsector < 7) & (OutbSim.GFid_Pcal2>0), "event"].unique())/len(OutbSim.loc[OutbSim.Gsector < 7, "event"].unique())
  OutbSim_effect_Gft = len(OutbSim.loc[(OutbSim.Gsector > 7) & (OutbSim.GFid_FT>0), "event"].unique())/len(OutbSim.loc[OutbSim.Gsector > 7, "event"].unique())
  OutbSim_effect_FDFD = len(OutbSim.loc[(OutbSim.EFid*OutbSim.PFid*OutbSim.GFid > 0) & (OutbSim.config==1), "event"].unique())/sum(OutbSim.config==1)
  OutbSim_effect_CDFD = len(OutbSim.loc[(OutbSim.EFid*OutbSim.PFid*OutbSim.GFid > 0) & (OutbSim.config==2), "event"].unique())/sum(OutbSim.config==2)
  OutbSim_effect_CDFT = len(OutbSim.loc[(OutbSim.EFid*OutbSim.PFid*OutbSim.GFid > 0) & (OutbSim.config==3), "event"].unique())/sum(OutbSim.config==3)

  print("&Exp.&Sim.&Exp.:Sim.&Exp.&Sim.& Exp.:Sim.\\\\")
  print("&Inb.&Inb. &Inb. &Outb.&Outb. &Outb. \\\\\hline")
  print("$e'$ PCAL &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Epcal, 100*InbSim_effect_Epcal, 100*InbExp_effect_Epcal/InbSim_effect_Epcal, 100*OutbExp_effect_Epcal, 100*OutbSim_effect_Epcal, 100*OutbExp_effect_Epcal/OutbSim_effect_Epcal))
  print("$e'$ DC &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Edc, 100*InbSim_effect_Edc, 100*InbExp_effect_Edc/InbSim_effect_Edc, 100*OutbExp_effect_Edc, 100*OutbSim_effect_Edc, 100*OutbExp_effect_Edc/OutbSim_effect_Edc))
  print("$e'$ SF &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Esf, 100*InbSim_effect_Esf, 100*InbExp_effect_Esf/InbSim_effect_Esf, 100*OutbExp_effect_Esf, 100*OutbSim_effect_Esf, 100*OutbExp_effect_Esf/OutbSim_effect_Esf))
  print("$e'$ $vz$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Evz, 100*InbSim_effect_Evz, 100*InbExp_effect_Evz/InbSim_effect_Evz, 100*OutbExp_effect_Evz, 100*OutbSim_effect_Evz, 100*OutbExp_effect_Evz/OutbSim_effect_Evz))
  print("$e'$ $E_{dep.}$" + " &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Eedep, 100*InbSim_effect_Eedep, 100*InbExp_effect_Eedep/InbSim_effect_Eedep, 100*OutbExp_effect_Eedep, 100*OutbSim_effect_Eedep, 100*OutbExp_effect_Eedep/OutbSim_effect_Eedep))
  print("$e'$ anti-$\pi^-$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hline".format(100*InbExp_effect_Eap, 100*InbSim_effect_Eap, 100*InbExp_effect_Eap/InbSim_effect_Eap, 100*OutbExp_effect_Eap, 100*OutbSim_effect_Eap, 100*OutbExp_effect_Eap/OutbSim_effect_Eap))
  print("$p'$ DC &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Pdc, 100*InbSim_effect_Pdc, 100*InbExp_effect_Pdc/InbSim_effect_Pdc, 100*OutbExp_effect_Pdc, 100*OutbSim_effect_Pdc, 100*OutbExp_effect_Pdc/OutbSim_effect_Pdc))
  print("$p'$ CVT &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Pcvt, 100*InbSim_effect_Pcvt, 100*InbExp_effect_Pcvt/InbSim_effect_Pcvt, 100*OutbExp_effect_Pcvt, 100*OutbSim_effect_Pcvt, 100*OutbExp_effect_Pcvt/OutbSim_effect_Pcvt))
  print("$p'$ $\chi$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Pchi, 100*InbSim_effect_Pchi, 100*InbExp_effect_Pchi/InbSim_effect_Pchi, 100*OutbExp_effect_Pchi, 100*OutbSim_effect_Pchi, 100*OutbExp_effect_Pchi/OutbSim_effect_Pchi))
  print("$p'$ $vz_{e'}-vz_{p'}$" + " &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hline".format(100*InbExp_effect_Pvz, 100*InbSim_effect_Pvz, 100*InbExp_effect_Pvz/InbSim_effect_Pvz, 100*OutbExp_effect_Pvz, 100*OutbSim_effect_Pvz, 100*OutbExp_effect_Pvz/OutbSim_effect_Pvz))
  print("$\gamma$ $\\beta$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Gbeta, 100*InbSim_effect_Gbeta, 100*InbExp_effect_Gbeta/InbSim_effect_Gbeta, 100*OutbExp_effect_Gbeta, 100*OutbSim_effect_Gbeta, 100*OutbExp_effect_Gbeta/OutbSim_effect_Gbeta))
  print("$\gamma$ PCAL &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Gpcal, 100*InbSim_effect_Gpcal, 100*InbExp_effect_Gpcal/InbSim_effect_Gpcal, 100*OutbExp_effect_Gpcal, 100*OutbSim_effect_Gpcal, 100*OutbExp_effect_Gpcal/OutbSim_effect_Gpcal2))
  print("$\gamma$ FT &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hline".format(100*InbExp_effect_Gft, 100*InbSim_effect_Gft, 100*InbExp_effect_Gft/InbSim_effect_Gft, 100*OutbExp_effect_Gft, 100*OutbSim_effect_Gft, 100*OutbExp_effect_Gft/OutbSim_effect_Gft))
  print("(FD, FD) &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\".format(100*InbExp_effect_FDFD, 100*InbSim_effect_FDFD, 100*InbExp_effect_FDFD/InbSim_effect_FDFD, 100*OutbExp_effect_FDFD, 100*OutbSim_effect_FDFD, 100*OutbExp_effect_FDFD/OutbSim_effect_FDFD))
  print("(CD, FD) &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\".format(100*InbExp_effect_CDFD, 100*InbSim_effect_CDFD, 100*InbExp_effect_CDFD/InbSim_effect_CDFD, 100*OutbExp_effect_CDFD, 100*OutbSim_effect_CDFD, 100*OutbExp_effect_CDFD/OutbSim_effect_CDFD))
  print("(CD, FT) &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hhline".format(100*InbExp_effect_CDFT, 100*InbSim_effect_CDFT, 100*InbExp_effect_CDFT/InbSim_effect_CDFT, 100*OutbExp_effect_CDFT, 100*OutbSim_effect_CDFT, 100*OutbExp_effect_CDFT/OutbSim_effect_CDFT)+"{|=|=|=|=|=|=|=|}")
  print("\n")

  InbExp.loc[:, "EFid"] = 1
  InbExp.loc[:, "PFid"] = 1
  InbExp.loc[:, "GFid"] = 1
  InbSim.loc[:, "EFid"] = 1
  InbSim.loc[:, "PFid"] = 1
  InbSim.loc[:, "GFid"] = 1
  OutbExp.loc[:, "EFid"] = 1
  OutbExp.loc[:, "PFid"] = 1
  OutbExp.loc[:, "GFid"] = 1
  OutbSim.loc[:, "EFid"] = 1
  OutbSim.loc[:, "PFid"] = 1
  OutbSim.loc[:, "GFid"] = 1
  for df_gammaRec in [InbExp, InbSim, OutbExp, OutbSim]:
    ang = -np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
    GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
    GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)
    df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
    df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

  electronFiducialCounting(InbExp, pol = "inbending", mc = False, fidlevel = 'tight')
  protonFiducialCounting(InbExp, pol = "inbending", fidlevel = 'tight')
  gammaFiducialCounting(InbExp)
  electronFiducialCounting(InbSim, pol = "inbending", mc = True, fidlevel = 'tight')
  protonFiducialCounting(InbSim, pol = "inbending", fidlevel = 'tight')
  gammaFiducialCounting(InbSim)
  electronFiducialCounting(OutbExp, pol = "outbending", mc = False, fidlevel = 'tight')
  protonFiducialCounting(OutbExp, pol = "outbending", fidlevel = 'tight')
  gammaFiducialCounting(OutbExp)
  electronFiducialCounting(OutbSim, pol = "outbending", mc = True, fidlevel = 'tight')
  protonFiducialCounting(OutbSim, pol = "outbending", fidlevel = 'tight')
  gammaFiducialCounting(OutbSim)

  InbExp_effect_Epcal = len(InbExp.loc[InbExp.EFid_pcal>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Edc = len(InbExp.loc[InbExp.EFid_dc>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Esf = len(InbExp.loc[InbExp.EFid_sf>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Evz = len(InbExp.loc[InbExp.EFid_vz>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Eedep = len(InbExp.loc[InbExp.EFid_edep>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Eap = len(InbExp.loc[InbExp.EFid_ap>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Pdc = len(InbExp.loc[(InbExp.Psector < 7) & (InbExp.PFid_dc>0), "event"].unique())/len(InbExp.loc[InbExp.Psector < 7, "event"].unique())
  InbExp_effect_Pcvt = len(InbExp.loc[(InbExp.Psector > 7) & (InbExp.PFid_cvt>0), "event"].unique())/len(InbExp.loc[InbExp.Psector > 7, "event"].unique())
  InbExp_effect_Pchi = len(InbExp.loc[InbExp.PFid_chi>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Pvz = len(InbExp.loc[InbExp.PFid_vz>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Gbeta = len(InbExp.loc[InbExp.GFid_beta>0, "event"].unique())/len(InbExp.loc[:, "event"].unique())
  InbExp_effect_Gpcal = len(InbExp.loc[(InbExp.Gsector < 7) & (InbExp.GFid_Pcal>0), "event"].unique())/len(InbExp.loc[InbExp.Gsector < 7, "event"].unique())
  InbExp_effect_Gft = len(InbExp.loc[(InbExp.Gsector > 7) & (InbExp.GFid_FT>0), "event"].unique())/len(InbExp.loc[InbExp.Gsector > 7, "event"].unique())
  InbExp_effect_FDFD = len(InbExp.loc[(InbExp.EFid*InbExp.PFid*InbExp.GFid > 0) & (InbExp.config==1), "event"].unique())/sum(InbExp.config==1)
  InbExp_effect_CDFD = len(InbExp.loc[(InbExp.EFid*InbExp.PFid*InbExp.GFid > 0) & (InbExp.config==2), "event"].unique())/sum(InbExp.config==2)
  InbExp_effect_CDFT = len(InbExp.loc[(InbExp.EFid*InbExp.PFid*InbExp.GFid > 0) & (InbExp.config==3), "event"].unique())/sum(InbExp.config==3)

  InbSim_effect_Epcal = len(InbSim.loc[InbSim.EFid_pcal>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Edc = len(InbSim.loc[InbSim.EFid_dc>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Esf = len(InbSim.loc[InbSim.EFid_sf>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Evz = len(InbSim.loc[InbSim.EFid_vz>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Eedep = len(InbSim.loc[InbSim.EFid_edep>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Eap = len(InbSim.loc[InbSim.EFid_ap>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Pdc = len(InbSim.loc[(InbSim.Psector < 7) & (InbSim.PFid_dc>0), "event"].unique())/len(InbSim.loc[InbSim.Psector < 7, "event"].unique())
  InbSim_effect_Pcvt = len(InbSim.loc[(InbSim.Psector > 7) & (InbSim.PFid_cvt>0), "event"].unique())/len(InbSim.loc[InbSim.Psector > 7, "event"].unique())
  InbSim_effect_Pchi = len(InbSim.loc[InbSim.PFid_chi>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Pvz = len(InbSim.loc[InbSim.PFid_vz>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Gbeta = len(InbSim.loc[InbSim.GFid_beta>0, "event"].unique())/len(InbSim.loc[:, "event"].unique())
  InbSim_effect_Gpcal = len(InbSim.loc[(InbSim.Gsector < 7) & (InbSim.GFid_Pcal>0), "event"].unique())/len(InbSim.loc[InbSim.Gsector < 7, "event"].unique())
  InbSim_effect_Gft = len(InbSim.loc[(InbSim.Gsector > 7) & (InbSim.GFid_FT>0), "event"].unique())/len(InbSim.loc[InbSim.Gsector > 7, "event"].unique())
  InbSim_effect_FDFD = len(InbSim.loc[(InbSim.EFid*InbSim.PFid*InbSim.GFid > 0) & (InbSim.config==1), "event"].unique())/sum(InbSim.config==1)
  InbSim_effect_CDFD = len(InbSim.loc[(InbSim.EFid*InbSim.PFid*InbSim.GFid > 0) & (InbSim.config==2), "event"].unique())/sum(InbSim.config==2)
  InbSim_effect_CDFT = len(InbSim.loc[(InbSim.EFid*InbSim.PFid*InbSim.GFid > 0) & (InbSim.config==3), "event"].unique())/sum(InbSim.config==3)

  OutbExp_effect_Epcal = len(OutbExp.loc[OutbExp.EFid_pcal>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Edc = len(OutbExp.loc[OutbExp.EFid_dc>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Esf = len(OutbExp.loc[OutbExp.EFid_sf>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Evz = len(OutbExp.loc[OutbExp.EFid_vz>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Eedep = len(OutbExp.loc[OutbExp.EFid_edep>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Eap = len(OutbExp.loc[OutbExp.EFid_ap>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Pdc = len(OutbExp.loc[(OutbExp.Psector < 7) & (OutbExp.PFid_dc>0), "event"].unique())/len(OutbExp.loc[OutbExp.Psector < 7, "event"].unique())
  OutbExp_effect_Pcvt = len(OutbExp.loc[(OutbExp.Psector > 7) & (OutbExp.PFid_cvt>0), "event"].unique())/len(OutbExp.loc[OutbExp.Psector > 7, "event"].unique())
  OutbExp_effect_Pchi = len(OutbExp.loc[OutbExp.PFid_chi>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Pvz = len(OutbExp.loc[OutbExp.PFid_vz>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Gbeta = len(OutbExp.loc[OutbExp.GFid_beta>0, "event"].unique())/len(OutbExp.loc[:, "event"].unique())
  OutbExp_effect_Gpcal = len(OutbExp.loc[(OutbExp.Gsector < 7) & (OutbExp.GFid_Pcal>0), "event"].unique())/len(OutbExp.loc[OutbExp.Gsector < 7, "event"].unique())
  OutbExp_effect_Gft = len(OutbExp.loc[(OutbExp.Gsector > 7) & (OutbExp.GFid_FT>0), "event"].unique())/len(OutbExp.loc[OutbExp.Gsector > 7, "event"].unique())
  OutbExp_effect_FDFD = len(OutbExp.loc[(OutbExp.EFid*OutbExp.PFid*OutbExp.GFid > 0) & (OutbExp.config==1), "event"].unique())/sum(OutbExp.config==1)
  OutbExp_effect_CDFD = len(OutbExp.loc[(OutbExp.EFid*OutbExp.PFid*OutbExp.GFid > 0) & (OutbExp.config==2), "event"].unique())/sum(OutbExp.config==2)
  OutbExp_effect_CDFT = len(OutbExp.loc[(OutbExp.EFid*OutbExp.PFid*OutbExp.GFid > 0) & (OutbExp.config==3), "event"].unique())/sum(OutbExp.config==3)

  OutbSim_effect_Epcal = len(OutbSim.loc[OutbSim.EFid_pcal>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Edc = len(OutbSim.loc[OutbSim.EFid_dc>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Esf = len(OutbSim.loc[OutbSim.EFid_sf>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Evz = len(OutbSim.loc[OutbSim.EFid_vz>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Eedep = len(OutbSim.loc[OutbSim.EFid_edep>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Eap = len(OutbSim.loc[OutbSim.EFid_ap>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Pdc = len(OutbSim.loc[(OutbSim.Psector < 7) & (OutbSim.PFid_dc>0), "event"].unique())/len(OutbSim.loc[OutbSim.Psector < 7, "event"].unique())
  OutbSim_effect_Pcvt = len(OutbSim.loc[(OutbSim.Psector > 7) & (OutbSim.PFid_cvt>0), "event"].unique())/len(OutbSim.loc[OutbSim.Psector > 7, "event"].unique())
  OutbSim_effect_Pchi = len(OutbSim.loc[OutbSim.PFid_chi>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Pvz = len(OutbSim.loc[OutbSim.PFid_vz>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Gbeta = len(OutbSim.loc[OutbSim.GFid_beta>0, "event"].unique())/len(OutbSim.loc[:, "event"].unique())
  OutbSim_effect_Gpcal = len(OutbSim.loc[(OutbSim.Gsector < 7) & (OutbSim.GFid_Pcal>0), "event"].unique())/len(OutbSim.loc[OutbSim.Gsector < 7, "event"].unique())
  OutbSim_effect_Gft = len(OutbSim.loc[(OutbSim.Gsector > 7) & (OutbSim.GFid_FT>0), "event"].unique())/len(OutbSim.loc[OutbSim.Gsector > 7, "event"].unique())
  OutbSim_effect_FDFD = len(OutbSim.loc[(OutbSim.EFid*OutbSim.PFid*OutbSim.GFid > 0) & (OutbSim.config==1), "event"].unique())/sum(OutbSim.config==1)
  OutbSim_effect_CDFD = len(OutbSim.loc[(OutbSim.EFid*OutbSim.PFid*OutbSim.GFid > 0) & (OutbSim.config==2), "event"].unique())/sum(OutbSim.config==2)
  OutbSim_effect_CDFT = len(OutbSim.loc[(OutbSim.EFid*OutbSim.PFid*OutbSim.GFid > 0) & (OutbSim.config==3), "event"].unique())/sum(OutbSim.config==3)

  print("&Exp.&Sim.&Exp.:Sim.&Exp.&Sim.& Exp.:Sim.\\\\")
  print("&Inb.&Inb. &Inb. &Outb.&Outb. &Outb. \\\\\hline")
  print("$e'$ PCAL &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Epcal, 100*InbSim_effect_Epcal, 100*InbExp_effect_Epcal/InbSim_effect_Epcal, 100*OutbExp_effect_Epcal, 100*OutbSim_effect_Epcal, 100*OutbExp_effect_Epcal/OutbSim_effect_Epcal))
  print("$e'$ DC &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Edc, 100*InbSim_effect_Edc, 100*InbExp_effect_Edc/InbSim_effect_Edc, 100*OutbExp_effect_Edc, 100*OutbSim_effect_Edc, 100*OutbExp_effect_Edc/OutbSim_effect_Edc))
  print("$e'$ SF &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Esf, 100*InbSim_effect_Esf, 100*InbExp_effect_Esf/InbSim_effect_Esf, 100*OutbExp_effect_Esf, 100*OutbSim_effect_Esf, 100*OutbExp_effect_Esf/OutbSim_effect_Esf))
  print("$e'$ $vz$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Evz, 100*InbSim_effect_Evz, 100*InbExp_effect_Evz/InbSim_effect_Evz, 100*OutbExp_effect_Evz, 100*OutbSim_effect_Evz, 100*OutbExp_effect_Evz/OutbSim_effect_Evz))
  print("$e'$ $E_{dep.}$" + " &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Eedep, 100*InbSim_effect_Eedep, 100*InbExp_effect_Eedep/InbSim_effect_Eedep, 100*OutbExp_effect_Eedep, 100*OutbSim_effect_Eedep, 100*OutbExp_effect_Eedep/OutbSim_effect_Eedep))
  print("$e'$ anti-$\pi^-$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hline".format(100*InbExp_effect_Eap, 100*InbSim_effect_Eap, 100*InbExp_effect_Eap/InbSim_effect_Eap, 100*OutbExp_effect_Eap, 100*OutbSim_effect_Eap, 100*OutbExp_effect_Eap/OutbSim_effect_Eap))
  print("$p'$ DC &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Pdc, 100*InbSim_effect_Pdc, 100*InbExp_effect_Pdc/InbSim_effect_Pdc, 100*OutbExp_effect_Pdc, 100*OutbSim_effect_Pdc, 100*OutbExp_effect_Pdc/OutbSim_effect_Pdc))
  print("$p'$ CVT &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Pcvt, 100*InbSim_effect_Pcvt, 100*InbExp_effect_Pcvt/InbSim_effect_Pcvt, 100*OutbExp_effect_Pcvt, 100*OutbSim_effect_Pcvt, 100*OutbExp_effect_Pcvt/OutbSim_effect_Pcvt))
  print("$p'$ $\chi$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Pchi, 100*InbSim_effect_Pchi, 100*InbExp_effect_Pchi/InbSim_effect_Pchi, 100*OutbExp_effect_Pchi, 100*OutbSim_effect_Pchi, 100*OutbExp_effect_Pchi/OutbSim_effect_Pchi))
  print("$p'$ $vz_{e'}-vz_{p'}$" + " &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hline".format(100*InbExp_effect_Pvz, 100*InbSim_effect_Pvz, 100*InbExp_effect_Pvz/InbSim_effect_Pvz, 100*OutbExp_effect_Pvz, 100*OutbSim_effect_Pvz, 100*OutbExp_effect_Pvz/OutbSim_effect_Pvz))
  print("$\gamma$ $\\beta$ &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Gbeta, 100*InbSim_effect_Gbeta, 100*InbExp_effect_Gbeta/InbSim_effect_Gbeta, 100*OutbExp_effect_Gbeta, 100*OutbSim_effect_Gbeta, 100*OutbExp_effect_Gbeta/OutbSim_effect_Gbeta))
  print("$\gamma$ PCAL &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%\\\\".format(100*InbExp_effect_Gpcal, 100*InbSim_effect_Gpcal, 100*InbExp_effect_Gpcal/InbSim_effect_Gpcal, 100*OutbExp_effect_Gpcal, 100*OutbSim_effect_Gpcal, 100*OutbExp_effect_Gpcal/OutbSim_effect_Gpcal2))
  print("$\gamma$ FT &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hline".format(100*InbExp_effect_Gft, 100*InbSim_effect_Gft, 100*InbExp_effect_Gft/InbSim_effect_Gft, 100*OutbExp_effect_Gft, 100*OutbSim_effect_Gft, 100*OutbExp_effect_Gft/OutbSim_effect_Gft))
  print("(FD, FD) &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\".format(100*InbExp_effect_FDFD, 100*InbSim_effect_FDFD, 100*InbExp_effect_FDFD/InbSim_effect_FDFD, 100*OutbExp_effect_FDFD, 100*OutbSim_effect_FDFD, 100*OutbExp_effect_FDFD/OutbSim_effect_FDFD))
  print("(CD, FD) &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\".format(100*InbExp_effect_CDFD, 100*InbSim_effect_CDFD, 100*InbExp_effect_CDFD/InbSim_effect_CDFD, 100*OutbExp_effect_CDFD, 100*OutbSim_effect_CDFD, 100*OutbExp_effect_CDFD/OutbSim_effect_CDFD))
  print("(CD, FT) &{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\%&{:.1f}\% \\\\\hhline".format(100*InbExp_effect_CDFT, 100*InbSim_effect_CDFT, 100*InbExp_effect_CDFT/InbSim_effect_CDFT, 100*OutbExp_effect_CDFT, 100*OutbSim_effect_CDFT, 100*OutbExp_effect_CDFT/OutbSim_effect_CDFT)+"{|=|=|=|=|=|=|=|}")
  print("\n")

  exit()

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

    yticklabel_1 = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    yticklabel_2 = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    yticklabel_3 = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    yticklabel_1[-1]  = "{}\n{}".format(yticklabel_1[-1], yticklabel_1[0])
    yticklabel_1[0]   = ''
    yticklabel_2[0]   = ''
    yticklabel_3[0]  = "0\n"
    yticklabel_3[-1]  = "{}\n{}".format(yticklabel_3[-1], yticklabel_3[0])

    for config in [1, 2, 3]:
      key  = (polarity, config, args.plot)
      for varind in range(8):

        # DVCS plots
        var    = dvcsvars[varind]
        label  = dvcstitles[varind]
        unit   = dvcsunits[varind]

        xlb     = xlbs_epg[key][varind]
        xub     = xubs_epg[key][varind]

        xtick        = xticks_epg[key, varind]
        xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
        xticklabel_1[-1] = ''
        xticklabel_2[-1] = ''
        xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])

        ind = 0 
        filenum = 0
        
        fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))

        yub    = np.zeros((5, 5))
        filled = np.zeros((5, 5))
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
              yub[xind, yind] = np.max([expDist, simDist])
              filled[xind, yind] = 1
                                          
              axs[xind, yind].hist(bins[:-1], bins, weights = expDist, histtype = 'step', color='b', linewidth=3, label = exp_label)
              axs[xind, yind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=3, label = sim_label)

              annotation = "({}, {}, {})".format(xBbin, Q2bin, tbin)
              axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)

              axs[xind, yind].set_ylim([0, 1.2*yub[xind, yind]])
              axs[xind, yind].set_xlim([xlb, xub])

              axs[xind, yind].get_xaxis().set_visible(False)
              axs[xind, yind].get_yaxis().set_visible(False)

              if (xind == 5 -1) and (yind==0):
                axs[xind, yind].get_xaxis().set_visible(True)
                axs[xind, yind].set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
                axs[xind, yind].set_xticks(xtick, xticklabel_2)

              elif (xind == 5 -1) and (yind < 5 -1):
                axs[xind, yind].get_xaxis().set_visible(True)
                axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
                axs[xind, yind].set_xticks(xtick, xticklabel_1)

              if (xind == 0) and (yind==0):
                axs[xind, yind].get_yaxis().set_visible(True)
                ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
                axs[xind, yind].set_yticks(ytick, yticklabel_2)

              elif (xind < 5 -1) and (yind == 0):
                axs[xind, yind].get_yaxis().set_visible(True)
                ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
                axs[xind, yind].set_yticks(ytick, yticklabel_1)

              if (xind == 5 -1) and (yind==0):
                axs[xind, yind].get_yaxis().set_visible(True)
                ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
                axs[xind, yind].set_yticks(ytick, yticklabel_3)

              if (xind == 5 -1) and (yind == 5 -1):
                axs[xind, yind].get_xaxis().set_visible(True)
                axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
                axs[xind, yind].set_xticks(xtick, xticklabel_3)

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
            axs[xind, yind].set_xticks(xtick, xticklabel_2)

          elif (xind == 5 -1) and (yind < 5 -1):
            axs[xind, yind].get_xaxis().set_visible(True)
            axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
            axs[xind, yind].set_xticks(xtick, xticklabel_1)

          if (xind == 0) and (yind==0):
            axs[xind, yind].get_yaxis().set_visible(True)
            ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            axs[xind, yind].set_yticks(ytick, yticklabel_2)

          elif (xind < 5 -1) and (yind == 0):
            axs[xind, yind].get_yaxis().set_visible(True)
            ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            axs[xind, yind].set_yticks(ytick, yticklabel_1)

          if (xind == 5 -1) and (yind==0):
            axs[xind, yind].get_yaxis().set_visible(True)
            ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            axs[xind, yind].set_yticks(ytick, yticklabel_3)

          if (xind == 5 -1) and (yind == 5 -1):
            axs[xind, yind].get_xaxis().set_visible(True)
            axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
            axs[xind, yind].set_xticks(xtick, xticklabel_3)

            plt.subplots_adjust(wspace=0, hspace=0)
            plt.savefig("plots/q16/dset_e/{}/dvcs_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
            plt.clf()

        # DVπ0P plots
        var   = pi0vars[varind]
        label = pi0titles[varind]
        unit  = pi0units[varind]

        xlb     = xlbs_epgg[key][varind]
        xub     = xubs_epgg[key][varind]

        xtick        = xticks_epgg[key, varind]
        xticklabel_1 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_2 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_3 = ["{:.2f}".format(i) for i in xtick]
        xticklabel_1[0]  = "{}\n{}".format(xticklabel_1[-1], xticklabel_1[0])
        xticklabel_1[-1] = ''
        xticklabel_2[-1] = ''
        xticklabel_3[0]  = "{}\n{}".format(xticklabel_3[-1], xticklabel_3[0])

        ind = 0 
        filenum = 0
        
        fig, axs = plt.subplots(5, 5, figsize = (32.5, 45))

        yub    = np.zeros((5, 5))
        filled = np.zeros((5, 5))
        for tbin in range(6):
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
              yub[xind, yind] = np.max([expDist, simDist])
              filled[xind, yind] = 1
                                          
              axs[xind, yind].hist(bins[:-1], bins, weights = expDist, histtype = 'step', color='b', linewidth=3, label = exp_label)
              axs[xind, yind].hist(bins[:-1], bins, weights = simDist, histtype = 'step', color='r', linewidth=3, label = sim_label)

              annotation = "({}, {}, {})".format(xBbin, Q2bin, tbin)
              axs[xind, yind].annotate(annotation, xy = (0.1, 0.9), xytext = (0.02, 0.9), xycoords = 'axes fraction', fontsize = 30)

              axs[xind, yind].set_ylim([0, 1.2*yub[xind, yind]])
              axs[xind, yind].set_xlim([xlb, xub])

              axs[xind, yind].get_xaxis().set_visible(False)
              axs[xind, yind].get_yaxis().set_visible(False)

              if (xind == 5 -1) and (yind==0):
                axs[xind, yind].get_xaxis().set_visible(True)
                axs[xind, yind].set_xlabel("\n" + label + " [" + unit + "]", fontsize = 40)
                axs[xind, yind].set_xticks(xtick, xticklabel_2)

              elif (xind == 5 -1) and (yind < 5 -1):
                axs[xind, yind].get_xaxis().set_visible(True)
                axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
                axs[xind, yind].set_xticks(xtick, xticklabel_1)

              if (xind == 0) and (yind==0):
                axs[xind, yind].get_yaxis().set_visible(True)
                ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
                axs[xind, yind].set_yticks(ytick, yticklabel_2)

              elif (xind < 5 -1) and (yind == 0):
                axs[xind, yind].get_yaxis().set_visible(True)
                ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
                axs[xind, yind].set_yticks(ytick, yticklabel_1)

              if (xind == 5 -1) and (yind==0):
                axs[xind, yind].get_yaxis().set_visible(True)
                ytick        = [0 * yub[xind, yind], 0.2 * yub[xind, yind], 0.4 * yub[xind, yind], 0.6 * yub[xind, yind], 0.8 * yub[xind, yind], 1.0 * yub[xind, yind], 1.2 * yub[xind, yind]]
                axs[xind, yind].set_yticks(ytick, yticklabel_3)

              if (xind == 5 -1) and (yind == 5 -1):
                axs[xind, yind].get_xaxis().set_visible(True)
                axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
                axs[xind, yind].set_xticks(xtick, xticklabel_3)
                  
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
            axs[xind, yind].set_xticks(xtick, xticklabel_2)

          elif (xind == 5 -1) and (yind < 5 -1):
            axs[xind, yind].get_xaxis().set_visible(True)
            axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
            axs[xind, yind].set_xticks(xtick, xticklabel_1)

          if (xind == 0) and (yind==0):
            axs[xind, yind].get_yaxis().set_visible(True)
            ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            axs[xind, yind].set_yticks(ytick, yticklabel_2)

          elif (xind < 5 -1) and (yind == 0):
            axs[xind, yind].get_yaxis().set_visible(True)
            ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            axs[xind, yind].set_yticks(ytick, yticklabel_1)

          if (xind == 5 -1) and (yind==0):
            axs[xind, yind].get_yaxis().set_visible(True)
            ytick        = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            axs[xind, yind].set_yticks(ytick, yticklabel_3)

          if (xind == 5 -1) and (yind == 5 -1):
            axs[xind, yind].get_xaxis().set_visible(True)
            axs[xind, yind].set_xlabel(label + " [" + unit + "]", fontsize = 40)
            axs[xind, yind].set_xticks(xtick, xticklabel_3)

            plt.subplots_adjust(wspace=0, hspace=0)
            plt.savefig("plots/q16/dset_e/{}/pi0_{}_{}_{}.pdf".format(polarity, var, topo[config], filenum), bbox_inches='tight')
            plt.clf()
