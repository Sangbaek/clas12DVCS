import matplotlib
from utils.epg import *
import matplotlib.pyplot as plt
from copy import copy
cmap = copy(plt.cm.get_cmap("jet"))
from scipy.optimize import least_squares
from matplotlib.colors import LogNorm
cmap.set_under('w',0)
cmap.set_bad('w',0)

degree = r"${}^{\circ}$"
GeV = "GeV"
GeV2 = "GeV"+r"${}^{2}$"
GeVc = "GeV/c"
GeVc2 = "(GeV/c)"+r"${}^{2}$"

# initial settings
pgf_with_latex = {
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,            # use LaTeX to write all text
        "font.family": "sans-serif",         
        "font.sans-serif": "Helvetica",
        "font.size": 25,                # default font size
        "axes.labelsize": 24,           # x and y label size
        "axes.titlesize": 24,           # subfigure title size, i.e. title size when one figure
        "legend.fontsize": 22,          # legend size
        "xtick.labelsize": 23,          # x axis tick label size
        "ytick.labelsize": 23,          # y axis tick label 
        "figure.titlesize": 25,         # Figure title size, useful when you have multiple plots in one canvas.
        "pgf.preamble": r"\usepackage{xcolor}",     # xcolor for colours
        "figure.autolayout": True
}
matplotlib.rcParams.update(pgf_with_latex)

epgExpInb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/inb/exp/dvcs.pkl")
pi0ExpInb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/inb/exp/pi0.pkl")
dvcsSimInb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/inb/dvcs/4893.pkl")
bkgSimInb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/inb/bkg_1g/4076.pkl")
pi0SimInb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/inb/bkg_2g/4076.pkl")

epgExpInbCDFT = epgExpInb.loc[epgExpInb.config == 3]
dvcsSimInbCDFT = dvcsSimInb.loc[dvcsSimInb.config == 3]
bkgSimInbCDFT = bkgSimInb.loc[bkgSimInb.config == 3]
pi0ExpInbCDFT = pi0ExpInb.loc[(pi0ExpInb.config == 3)]
pi0SimInbCDFT = pi0SimInb.loc[(pi0SimInb.config == 3)]

epgExpInbCD = epgExpInb.loc[epgExpInb.config == 2]
dvcsSimInbCD = dvcsSimInb.loc[dvcsSimInb.config == 2]
bkgSimInbCD = bkgSimInb.loc[bkgSimInb.config == 2]
pi0ExpInbCD = pi0ExpInb.loc[(pi0ExpInb.config == 2)]
pi0SimInbCD = pi0SimInb.loc[(pi0SimInb.config == 2)]

epgExpInbFD = epgExpInb.loc[epgExpInb.config == 1]
dvcsSimInbFD = dvcsSimInb.loc[dvcsSimInb.config == 1]
bkgSimInbFD = bkgSimInb.loc[bkgSimInb.config == 1]
pi0ExpInbFD = pi0ExpInb.loc[(pi0ExpInb.config == 1)]
pi0SimInbFD = pi0SimInb.loc[(pi0SimInb.config == 1)]

epgExpOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/outb/exp/dvcs.pkl")
pi0ExpOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/outb/exp/pi0.pkl")
dvcsSimOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/outb/dvcs/4907.pkl")
bkgSimOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/outb/bkg_1g/4243.pkl")
pi0SimOutb = pd.read_pickle("/volatile/clas12/sangbaek/nov2021/convPkl_full_eb/outb/bkg_2g/4243.pkl")

epgExpOutbCDFT = epgExpOutb.loc[epgExpOutb.config == 3]
dvcsSimOutbCDFT = dvcsSimOutb.loc[dvcsSimOutb.config == 3]
bkgSimOutbCDFT = bkgSimOutb.loc[bkgSimOutb.config == 3]
pi0ExpOutbCDFT = pi0ExpOutb.loc[(pi0ExpOutb.config == 3)]
pi0SimOutbCDFT = pi0SimOutb.loc[(pi0SimOutb.config == 3)]

epgExpOutbCD = epgExpOutb.loc[epgExpOutb.config == 2]
dvcsSimOutbCD = dvcsSimOutb.loc[dvcsSimOutb.config == 2]
bkgSimOutbCD = bkgSimOutb.loc[bkgSimOutb.config == 2]
pi0ExpOutbCD = pi0ExpOutb.loc[(pi0ExpOutb.config == 2)]
pi0SimOutbCD = pi0SimOutb.loc[(pi0SimOutb.config == 2)]

epgExpOutbFD = epgExpOutb.loc[epgExpOutb.config == 1]
dvcsSimOutbFD = dvcsSimOutb.loc[dvcsSimOutb.config == 1]
bkgSimOutbFD = bkgSimOutb.loc[bkgSimOutb.config == 1]
pi0ExpOutbFD = pi0ExpOutb.loc[(pi0ExpOutb.config == 1)]
pi0SimOutbFD = pi0SimOutb.loc[(pi0SimOutb.config == 1)]

# varstoplot = ["Pp", "Ptheta", "Pphi", "Gp", "Gtheta", "Gphi",  "coneAngle", "MM2_eg", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
# title = [r"$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", r"$p_{\gamma}$", r"$\theta_{\gamma}$", r"$\phi_{\gamma}$", r"$\theta_{e'\gamma}$", r"$MM^2_{e'\gamma}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi_{\vec{L}\vec{\Gamma}}$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
# unit = [GeV, degree, degree, GeV, degree, degree, degree, GeV2,  degree, degree, GeV, GeV2, GeV2, GeVc]

varstoplot = ["coneAngle", "MM2_eg", "reconGam", "coplanarity", "ME_epg", "MM2_epg", "MM2_ep", "MPt"]
title = [r"$\theta_{e'\gamma}$", r"$MM^2_{e'\gamma}$", r"$\theta_{\gamma_{det.}\gamma_{rec.}}$", r"$\Delta\phi_{\vec{L}\vec{\Gamma}}$" , "ME"+r"${}_{epg}$", "MM"+r"${}^{2}_{epg}$", "MM"+r"${}^{2}_{ep}$", "MPt"+r"${}_{epg}$"]
unit = [degree, GeV2,  degree, degree, GeV, GeV2, GeV2, GeVc]

df3 = epgExpInbCDFT
df1 = dvcsSimInbCDFT
df2 = bkgSimInbCDFT
df4 = pi0ExpInbCDFT
df5 = pi0SimInbCDFT

fig, axs = plt.subplots(2, 4, figsize = (15, 10))
for yind in range(0, 2):
    for xind in range(0, 4):
        ind = 4*yind + xind
        if varstoplot[ind]:
            pass
        else:
            continue
        simDist_dvcs, bins = np.histogram(df1.loc[:, varstoplot[ind]], 100, density = True)
        simDist_dvpi0, _ = np.histogram(df2.loc[:, varstoplot[ind]], bins, density = True)
        simDist = (1-contInbCDFT)*simDist_dvcs + contInbCDFT*simDist_dvpi0
        bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        axs[yind, xind].step(bincenters, simDist, where='mid',color='b', linewidth=1)
        axs[yind, xind].hist(df3.loc[:,varstoplot[ind]], bins = bins, histtype='stepfilled', facecolor='none', edgecolor='k', density=True, linewidth=1)
        axs[yind, xind].set_title(title[ind])
        # axs[yind, xind].set_xlim([start, end])
        if (unit[ind]):
            axs[yind, xind].set_xlabel(title[ind]+" [" + unit[ind] +"]")
        else:
            axs[yind, xind].set_xlabel(title[ind])
plt.tight_layout()
# plt.savefig(outDir+"InbCDFT{}_{:.3f}_{:.3f}_{:.3f}.pdf".format(i, sigma1, sigma2, sigma3))
# plt.savefig(outDir+"InbCDFT{}_{:.3f}_{:.3f}.pdf".format(i, sigma1, sigma2))
# plt.savefig(outDir+"InbCDFT{}_{:.3f}.pdf".format(i, sigma2))
# plt.show()
plt.savefig("plots/dvcsInbCDFTexcl.pdf")
plt.clf()