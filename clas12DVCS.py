from utils.epg import *
import icecream as ic
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from copy import copy

# initial settings
pgf_with_latex = {
		"pgf.texsystem": "pdflatex",
		"text.usetex": True,            # use LaTeX to write all text
		"font.family": "sans-serif",         
		"font.sans-serif": "Helvetica",
		"font.size": 10,				# default font size
		"axes.labelsize": 10,			# x and y label size
		"axes.titlesize": 10,           # subfigure title size, i.e. title size when one figure
		"legend.fontsize": 12,			# legend size
		"xtick.labelsize": 10,			# x axis tick label size
		"ytick.labelsize": 10,			# y axis tick label 
		"figure.titlesize": 10,         # Figure title size, useful when you have multiple plots in one canvas.
		"pgf.preamble": r"\usepackage{xcolor}"     # xcolor for colours
}
matplotlib.rcParams.update(pgf_with_latex)


fname = "~/Dropbox (MIT)/data/dvcs_inb.root"
epg = epgFromROOT(fname)
dvpi0 = epg.getDVpi0()
dvcs = epg.getDVCS(sub2g=True)

fname_mc = "~/Dropbox (MIT)/data/dvcs_mc_inb.root"
epg_mc = epgFromROOT(fname_mc)
dvpi0_mc = epg_mc.getDVpi0()
dvcs_mc = epg_mc.getDVCS(sub2g=True)

varstoplot = ["nu", "W", "Q2", "xB", "t2", "phi2", "MM2_ep", "MM2_epg", "ME_epg"]
title = [r"$\nu$", r"$W$", r"$Q^{2}$", r"$x_{B}$", r"$-t$", r"$\phi$", "MM"+r"${}^{2}_{ep}$", "MM"+r"${}^{2}_{epg}$", "ME"+r"${}_{epg}$"]
binstarts = [0, 2, 1, 0, 0, 0, -0.5,-0.04, -0.5]
binends = [11, 10, 14, 1, 2.5, 360, 1.5, 0.04, 1.2]
fig, axs = plt.subplots(3, 3)
for xind in range(0,3):
	for yind in range(0, 3):
		ind = xind+3*(2-yind)
		axs[2-yind, xind].hist(dvcs[varstoplot[ind]], bins = np.linspace(binstarts[ind], binends[ind], 101), density = True, histtype='stepfilled', facecolor='none', edgecolor='b')
		axs[2-yind, xind].hist(dvcs_mc[varstoplot[ind]], bins = np.linspace(binstarts[ind], binends[ind], 101), density = True, histtype='stepfilled', facecolor='none', edgecolor='r')
		axs[2-yind, xind].set_title(title[ind])
plt.tight_layout()
plt.savefig("simComparison.pdf")
plt.clf()