import pandas as pd
import matplotlib.pyplot as plt
from utils.const import *
from utils.physics import *
from utils.fiducial import *
from glob import glob
from copy import copy
from matplotlib.colors import LogNorm
import uproot
import awkward as ak
from utils.fiducial import *
from utils.const import *
from utils.physics import *
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import matplotlib
cmap = matplotlib.colormaps["jet"]
from scipy.optimize import curve_fit
# matplotlib.use('Agg')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=SyntaxWarning)
pd.options.mode.chained_assignment = None
import pickle
import itertools
import matplotlib.image as mpimg
from matplotlib.ticker import EngFormatter
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
import argparse
parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--integrated_binnum", type = int)
args = parser.parse_args()

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

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# initial settings
pgf_with_latex = {
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,            # use LaTeX to write all text
        "font.family": "sans-serif",        
        "font.sans-serif": "Helvetica",
        "font.size": 25,                # default font size
        "axes.titlepad": 20,            # x and y label size
        "axes.labelsize": 24,           # x and y label size
        "axes.titlesize": 24,         # subfigure title size, i.e. title size when one figure
        "legend.fontsize": 22,          # legend size
        "xtick.labelsize": 23,          # x axis tick label size
        "ytick.labelsize": 23,          # y axis tick label 
        "figure.titlesize": 25,         # Figure title size, useful when you have multiple plots in one canvas.
        "pgf.preamble": r"\usepackage{xcolor}",     # xcolor for colours
        "figure.autolayout": False
}
matplotlib.rcParams.update(pgf_with_latex)


def tmind(xB, Q2, t, phi):
    return -Q2*(2*(1-xB)*(1-sqeps2(xB, Q2, t, phi))+eps2(xB, Q2, t, phi))/(4*xB*(1-xB)+eps2(xB, Q2, t, phi))

def draw_box(xmin, xmax, ymin, ymax, ax):
    dummy = np.linspace(0, 1, 101, dtype = float)
    ax.plot(xmin + 0*dummy, ymin + (ymax-ymin)*dummy, color = 'k', ls = '--')
    ax.plot(xmax + 0*dummy, ymin + (ymax-ymin)*dummy, color = 'k', ls = '--')
    ax.plot(xmin + (xmax-xmin)*dummy, ymin + 0*dummy, color = 'k', ls = '--')
    ax.plot(xmin + (xmax-xmin)*dummy, ymax + 0*dummy, color = 'k', ls = '--')
    return 0

def convert_log_to_tick_number(num):
    return np.maximum(np.ones_like(num), 10*(num-np.floor(num))) * 10**(np.floor(num))

chunks_inb  = [1, 2]
chunks_outb = [2, 3]

sim_current_inb = 45
sim_current_outb = 50

effective_current_inb = (40*charge_inb_40nA + 45* charge_inb_45nA + 50*charge_inb_50nA + 55*charge_inb_55nA)/charge_inb
effective_current_outb = (5*charge_outb_5nA + 40* charge_outb_40nA + 50*charge_outb_50nA)/charge_outb

luminosity = luminosity_inb + luminosity_outb

def quartic_efficiency(x, *par):
    a, b, c, d, e = par
    return a*x**4 + b*x**3 + c*x**2 + d*x + e

quartic = quartic_efficiency

def cubic(x, *par):
    a, b, c, d = par
    return a*x**3 + b*x**2 + c*x + d

def linear(x, *p):
    a, b = p
    return a * x  + b

def cosine_fitting(x, *args):
    a, b = args
    return a + b*np.cos(x)

def cosine_fitting_1(x, *args):
    a, b = args
    return a + b*np.cos(x)

def cosine_fitting_2(x, *args):
    a, b, c = args
    return a + b*np.cos(x) + c*np.cos(2*x)

def cosine_fitting_3(x, *args):
    a, b, c, d = args
    return a + b*np.cos(x) + c*np.cos(2*x)  + d*np.cos(3*x)

def engineering_to_latex(string):
    string = string.replace('e+00', '')
    for i in range(1, 10):
        string = string.replace('e+0{}'.format(i), '$\times 10^{}$'.format(i))
        string = string.replace('e-0{}'.format(i), '$\times 10^{{-{}}}$'.format(i))
    for i in range(10, 20):
        string = string.replace('e+{}'.format(i), '$\times 10^{}$'.format(i))
        string = string.replace('e-{}'.format(i), '$\times 10^{{-{}}}$'.format(i))
    return string

def engineering_to_latex2(string):
    string = string.replace('e+00', '')
    for i in range(1, 10):
        string = string.replace('e+0{}'.format(i), '\\times 10^{}'.format(i))
        string = string.replace('e-0{}'.format(i), '\\times 10^{{-{}}}'.format(i))
    for i in range(10, 20):
        string = string.replace('e+{}'.format(i), '\\times 10^{}'.format(i))
        string = string.replace('e-{}'.format(i), '\\times 10^{{-{}}}'.format(i))
    return string

def log_ticker(ticks):
    ticklabels = []
    for tick in ticks:
        if tick == 1:
            ticklabel = "$1$"
        else:
            ticklabel = r"$10^{{{:.0f}}}$".format(np.log10(tick))
        ticklabels.append(ticklabel)
    return ticklabels

def cosine_1_fitting_one_bin_th(df_display_this_bin, p0 = (1, 0)):

  popt_1_th_km15, pcov = curve_fit(cosine_fitting_1, df_display_this_bin.phi_display, df_display_this_bin.xsec_KM15_display_w, p0 = p0)
  cosine_1_th_km15 =  -popt_1_th_km15[1]
  cosine_1_th_km15_err =  np.sqrt(np.diag(pcov))[1]

  popt_1_th_vgg, pcov = curve_fit(cosine_fitting_1, df_display_this_bin.phi_display, df_display_this_bin.xsec_VGG_display_w, p0 = p0)
  cosine_1_th_vgg =  -popt_1_th_vgg[1]
  cosine_1_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
  
  popt_1_th_bh, pcov = curve_fit(cosine_fitting_1, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_display_w, p0 = p0)
  cosine_1_th_bh =  -popt_1_th_bh[1]
  cosine_1_th_bh_err =  np.sqrt(np.diag(pcov))[1]

  popt_1_th_bh_km15, pcov = curve_fit(cosine_fitting_1, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_KM15_display_w, p0 = p0)
  cosine_1_th_bh_km15 =  -popt_1_th_bh_km15[1]
  cosine_1_th_bh_km15_err =  np.sqrt(np.diag(pcov))[1]

  return cosine_1_th_km15, cosine_1_th_km15_err, cosine_1_th_vgg, cosine_1_th_vgg_err, cosine_1_th_bh, cosine_1_th_bh_err, cosine_1_th_bh_km15, cosine_1_th_bh_km15_err, popt_1_th_km15, popt_1_th_vgg, popt_1_th_bh, popt_1_th_bh_km15

def cosine_2_fitting_one_bin_th(df_display_this_bin, p0 = (1, 0, 0)):

  popt_2_th_km15, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_KM15_display_w, p0 = p0)
  cosine_2_th_km15 =  -popt_2_th_km15[1]
  cosine_2_th_km15_err =  np.sqrt(np.diag(pcov))[1]

  popt_2_th_vgg, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_VGG_display_w, p0 = p0)
  cosine_2_th_vgg =  -popt_2_th_vgg[1]
  cosine_2_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
  
  popt_2_th_bh, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_display_w, p0 = p0)
  cosine_2_th_bh =  -popt_2_th_bh[1]
  cosine_2_th_bh_err =  np.sqrt(np.diag(pcov))[1]

  popt_2_th_bh_km15, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_KM15_display_w, p0 = p0)
  cosine_2_th_bh_km15 =  -popt_2_th_bh_km15[1]
  cosine_2_th_bh_km15_err =  np.sqrt(np.diag(pcov))[1]

  return cosine_2_th_km15, cosine_2_th_km15_err, cosine_2_th_vgg, cosine_2_th_vgg_err, cosine_2_th_bh, cosine_2_th_bh_err, cosine_2_th_bh_km15, cosine_2_th_bh_km15_err, popt_2_th_km15, popt_2_th_vgg, popt_2_th_bh, popt_2_th_bh_km15

def cosine_3_fitting_one_bin_th(df_display_this_bin, p0 = (1, 0, 0, 0)):
    
  popt_3_th_km15, pcov = curve_fit(cosine_fitting_3, df_display_this_bin.phi_display, df_display_this_bin.xsec_KM15_display_w, p0 = p0)
  cosine_3_th_km15 =  -popt_3_th_km15[1]
  cosine_3_th_km15_err =  np.sqrt(np.diag(pcov))[1]

  popt_3_th_vgg, pcov = curve_fit(cosine_fitting_3, df_display_this_bin.phi_display, df_display_this_bin.xsec_VGG_display_w, p0 = p0)
  cosine_3_th_vgg =  -popt_3_th_vgg[1]
  cosine_3_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
  
  popt_3_th_bh, pcov = curve_fit(cosine_fitting_3, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_display_w, p0 = p0)
  cosine_3_th_bh =  -popt_3_th_bh[1]
  cosine_3_th_bh_err =  np.sqrt(np.diag(pcov))[1]

  popt_3_th_bh_km15, pcov = curve_fit(cosine_fitting_3, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_KM15_display_w, p0 = p0)
  cosine_3_th_bh_km15 =  -popt_3_th_bh_km15[1]
  cosine_3_th_bh_km15_err =  np.sqrt(np.diag(pcov))[1]

  return cosine_3_th_km15, cosine_3_th_km15_err, cosine_3_th_vgg, cosine_3_th_vgg_err, cosine_3_th_bh, cosine_3_th_bh_err, cosine_3_th_bh_km15, cosine_3_th_bh_km15_err, popt_3_th_km15, popt_3_th_vgg, popt_3_th_bh, popt_3_th_bh_km15

def cosine_1_fitting_one_bin_exp(df_this_bin, p0 = (1, 0)):

  xB_avg = df_this_bin.xB_avg_this_point.unique()[0]
  Q2_avg = df_this_bin.Q2_avg_this_point.unique()[0]
  t_avg = df_this_bin.t_avg_this_point.unique()[0]

  popt_1_pi0, pcov = curve_fit(cosine_fitting_1, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  cosine_1_exp_pi0 =  -popt_1_pi0[1]
  cosine_1_exp_pi0_stat_err =  np.sqrt(np.diag(pcov))[1]

  popt_1_pi0_min, _ = curve_fit(cosine_fitting_1, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  popt_1_pi0_max, _ = curve_fit(cosine_fitting_1, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)

  popt_1_pi0s = []

  for trial in range(10**3):
    popt, pcov = curve_fit(cosine_fitting_1, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + np.random.normal(np.ones_like(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w), absolute_sigma = True,  p0 = p0)
    popt_1_pi0s.append(popt)
  # popt, pcov = curve_fit(cosine_fitting_1, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, absolute_sigma = True,  p0 = p0)
  # popt_1_pi0s.append(popt)
  popt_1_pi0s = np.array(popt_1_pi0s)
  cosine_1_exp_pi0s = -popt_1_pi0s[:, 1]

  # popt_1_pi0_min             = popt_1_pi0s[np.argmax(cosine_1_exp_pi0s)]
  # popt_1_pi0_max             = popt_1_pi0s[np.argmin(cosine_1_exp_pi0s)]

  cosine_1_exp_pi0_syst_err   = np.sqrt( (0.5* (np.abs(np.max(cosine_1_exp_pi0s) - cosine_1_exp_pi0) + np.abs(np.min(cosine_1_exp_pi0s) - cosine_1_exp_pi0))) **2 + 0.3**2 + 0.0476**2)

  popt_1_bkgmerging_only, pcov = curve_fit(cosine_fitting_1, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_bkg_merging_w, sigma = df_this_bin.xsec_exp_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  cosine_1_exp_bkgmerging_only =  -popt_1_bkgmerging_only[1]

  return xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_bkgmerging_only, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_bkgmerging_only

def cosine_2_fitting_one_bin_exp(df_this_bin, p0 = (1, 0, 0)):

  xB_avg = df_this_bin.xB_avg_this_point.unique()[0]
  Q2_avg = df_this_bin.Q2_avg_this_point.unique()[0]
  t_avg = df_this_bin.t_avg_this_point.unique()[0]

  popt_2_pi0, pcov = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  cosine_2_exp_pi0 =  -popt_2_pi0[1]
  cosine_2_exp_pi0_stat_err =  np.sqrt(np.diag(pcov))[1]

  popt_2_pi0_min, _ = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  popt_2_pi0_max, _ = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)

  popt_2_pi0s = []

  for trial in range(10**3):
    popt, pcov = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + np.random.normal(np.ones_like(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w), absolute_sigma = True,  p0 = p0)
    popt_2_pi0s.append(popt)
  # popt, pcov = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, absolute_sigma = True,  p0 = p0)
  # popt_2_pi0s.append(popt)
  popt_2_pi0s = np.array(popt_2_pi0s)
  cosine_2_exp_pi0s = -popt_2_pi0s[:, 1]

  # popt_2_pi0_min             = popt_2_pi0s[np.argmax(cosine_2_exp_pi0s)]
  # popt_2_pi0_max             = popt_2_pi0s[np.argmin(cosine_2_exp_pi0s)]

  cosine_2_exp_pi0_syst_err   = np.sqrt( (0.5* (np.abs(np.max(cosine_2_exp_pi0s) - cosine_2_exp_pi0) + np.abs(np.min(cosine_2_exp_pi0s) - cosine_2_exp_pi0))) **2 + 0.3**2 + 0.0476**2)

  popt_2_bkgmerging_only, pcov = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_bkg_merging_w, sigma = df_this_bin.xsec_exp_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  cosine_2_exp_bkgmerging_only =  -popt_1_bkgmerging_only[1]

  return xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_bkgmerging_only, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_bkgmerging_only

def cosine_3_fitting_one_bin_exp(df_this_bin, p0 = (1, 0, 0, 0)):

  xB_avg = df_this_bin.xB_avg_this_point.unique()[0]
  Q2_avg = df_this_bin.Q2_avg_this_point.unique()[0]
  t_avg = df_this_bin.t_avg_this_point.unique()[0]

  popt_3_pi0, pcov = curve_fit(cosine_fitting_3, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  cosine_3_exp_pi0 =  -popt_3_pi0[1]
  cosine_3_exp_pi0_stat_err =  np.sqrt(np.diag(pcov))[1]

  popt_3_pi0_min, _ = curve_fit(cosine_fitting_3, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  popt_3_pi0_max, _ = curve_fit(cosine_fitting_3, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)

  popt_3_pi0s = []

  for trial in range(10**3):
    popt, pcov = curve_fit(cosine_fitting_3, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + np.random.normal(np.ones_like(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w), absolute_sigma = True,  p0 = p0)
    popt_3_pi0s.append(popt)
  # popt, pcov = curve_fit(cosine_fitting_3, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, absolute_sigma = True,  p0 = p0)
  # popt_3_pi0s.append(popt)
  popt_3_pi0s = np.array(popt_3_pi0s)
  cosine_3_exp_pi0s = -popt_3_pi0s[:, 1]

  # popt_3_pi0_min             = popt_3_pi0s[np.argmax(cosine_3_exp_pi0s)]
  # popt_3_pi0_max             = popt_3_pi0s[np.argmin(cosine_3_exp_pi0s)]

  cosine_3_exp_pi0_syst_err   = np.sqrt( (0.5* (np.abs(np.max(cosine_3_exp_pi0s) - cosine_3_exp_pi0) + np.abs(np.min(cosine_3_exp_pi0s) - cosine_3_exp_pi0))) **2 + 0.3**2 + 0.0476**2)

  popt_3_bkgmerging_only, pcov = curve_fit(cosine_fitting_3, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_bkg_merging_w, sigma = df_this_bin.xsec_exp_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = p0)
  cosine_3_exp_bkgmerging_only =  -popt_1_bkgmerging_only[1]

  return xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_bkgmerging_only, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_bkgmerging_only

exp_to_BH_inb_mean_mean      = 0.765
exp_to_BH_inb_mean_stat_err  = 0.149
exp_to_BH_outb_mean_mean     = 0.679
exp_to_BH_outb_mean_stat_err = 0.105

exp_to_BH_inb_mean2_mean      = 0.775
exp_to_BH_inb_mean2_stat_err  = 0.100
exp_to_BH_outb_mean2_mean     = 0.775
exp_to_BH_outb_mean2_stat_err = 0.100

# df_summary_table_rebinned_sig = pd.read_pickle("impact_study_dec2024/summary_table.sig.rebinned.pkl")
# df_summary_table_rebinned_bkg = pd.read_pickle("impact_study_dec2024/summary_table.bkg.rebinned.pkl")
# df_summary_table_rebinned_exp = pd.read_pickle("impact_study_dec2024/summary_table.exp.rebinned.pkl")
# df_summary_table_rebinned_gen = pd.read_pickle("summary_table.gen.rebinned.pkl")
# df_summary_table_rebinned     = pd.read_pickle("df_summary_table_rebinned.backup.pkl")


# pi0_sigma_inb_in_nb  = 2.7706591819052937
# survival_rate_inb    = 10430/83265 #100000/804311
# pi0_sigma_outb_in_nb = 5.9541707415266248 
# survival_rate_outb   = 9505.0/94878.0 #100000/985416

# df_exp_epg_inbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_inb/dvcs/excl_level_2/pkl_7_nominal/fall2018_inb.pkl")
# df_exp_pi0_inbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_inb/pi0/excl_level_2/pkl_7_nominal/fall2018_inb.pkl")

# df_sim_bkg_inbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/pkl_7_nominal/{}/fall2018_inb.pkl".format(i)) for i in [1, 2]])
# df_sim_pi0_inbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/pkl_7_nominal/{}/fall2018_inb.pkl".format(i)) for i in [1, 2]])
# df_sim_dvcs_inbs  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl_7_nominal/fall2018_inb.pkl")
# df_sim_dvcs_inbs_bkgmerging  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl_14_bkgmerging/fall2018_inb.pkl")

# df_exp_epg_outbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_outb/dvcs/excl_level_2/pkl_7_nominal/fall2018_outb.pkl")
# df_exp_pi0_outbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_outb/pi0/excl_level_2/pkl_7_nominal/fall2018_outb.pkl")

# df_sim_bkg_outbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/pkl_7_nominal/{}/fall2018_outb.pkl".format(i)) for i in [2, 3]])
# df_sim_pi0_outbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/pkl_7_nominal/{}/fall2018_outb.pkl".format(i)) for i in [2, 3]])
# df_sim_dvcs_outbs  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl_7_nominal/fall2018_outb.pkl")
# df_sim_dvcs_outbs_bkgmerging  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl_14_bkgmerging/fall2018_outb.pkl")

# df_sim_bkg_inbs.loc[:, "weights"]  = pi0_sigma_inb_in_nb * luminosity_inb * survival_rate_inb/10**8/2
# df_sim_bkg_outbs.loc[:, "weights"] = pi0_sigma_outb_in_nb * luminosity_outb * survival_rate_outb/10**8/2

# df_sim_bkg_inbs.loc[:, "eff_bkg_merging_inb"] = 0
# for integrated_binnum in df_summary_table_rebinned.integrated_binnum.unique():
#     for phi_binnum in range(24):
#         if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy()):
#             df_sim_dvcs_inbs.loc[(df_sim_dvcs_inbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_inbs.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
#             df_sim_dvcs_outbs.loc[(df_sim_dvcs_outbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_outbs.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
#             df_sim_bkg_inbs.loc[(df_sim_bkg_inbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_inbs.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
#             df_sim_bkg_outbs.loc[(df_sim_bkg_outbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_outbs.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
#         else:
#             df_sim_dvcs_inbs.loc[(df_sim_dvcs_inbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_inbs.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
#             df_sim_dvcs_outbs.loc[(df_sim_dvcs_outbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_outbs.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
#             df_sim_bkg_inbs.loc[(df_sim_bkg_inbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_inbs.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
#             df_sim_bkg_outbs.loc[(df_sim_bkg_outbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_outbs.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()

# df_sim_dvcs_inbs = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.eff_bkg_merging_inb>0, :]
# df_sim_dvcs_outbs = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.eff_bkg_merging_outb>0, :]
# df_sim_bkg_inbs = df_sim_bkg_inbs.loc[df_sim_bkg_inbs.eff_bkg_merging_inb>0, :]
# df_sim_bkg_outbs = df_sim_bkg_outbs.loc[df_sim_bkg_outbs.eff_bkg_merging_outb>0, :]


# df_exp_epg_inbs = df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum>=1) & (df_exp_epg_inbs.integrated_binnum<=147)]
# df_exp_epg_outbs = df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum>=1) & (df_exp_epg_outbs.integrated_binnum<=147)]

# fig, axs = plt.subplots(3, 1, figsize = (8, 12))

# # MM2_epg_exp_inb, bins = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.config == 1].MM2_epg, weights = 1/df_exp_epg_inbs.loc[df_exp_epg_inbs.config == 1].efficiency, bins = np.linspace(-0.02, 0.02, 81))
# # MM2_epg_sim_inb, bins = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == 1].MM2_epg, weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == 1].eff_bkg_merging_inb * df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == 1].efficiency_pi0 * df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == 1].weights, bins = np.linspace(-0.02, 0.02, 81))
# # MM2_epg_bkg_inb, bins = np.histogram(df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == 1].MM2_epg, weights = df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == 1].eff_bkg_merging_inb * df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == 1].efficiency_pi0 * df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == 1].weights, bins = np.linspace(-0.02, 0.02, 81))

# for config in [1, 2, 3]:
#   MM2_epg_exp_inb, bins = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.config == config].MM2_epg, weights = 1/df_exp_epg_inbs.loc[df_exp_epg_inbs.config == config].efficiency, bins = np.linspace(-0.02, 0.02, 81))
#   MM2_epg_sim_inb, bins = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == config].MM2_epg, weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == config].eff_bkg_merging_inb * df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == config].efficiency_pi0 * df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.config == config].weights, bins = np.linspace(-0.02, 0.02, 81))
#   MM2_epg_bkg_inb, bins = np.histogram(df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == config].MM2_epg, weights = df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == config].eff_bkg_merging_inb * df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == config].efficiency_pi0 * df_sim_bkg_inbs.loc[df_sim_bkg_inbs.config == config].weights, bins = np.linspace(-0.02, 0.02, 81))

#   MM2_epg_exp_outb, bins = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.config == config].MM2_epg, weights = 1/df_exp_epg_outbs.loc[df_exp_epg_outbs.config == config].efficiency, bins = np.linspace(-0.02, 0.02, 81))
#   MM2_epg_sim_outb, bins = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.config == config].MM2_epg, weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.config == config].eff_bkg_merging_outb * df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.config == config].efficiency_pi0 * df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.config == config].weights, bins = np.linspace(-0.02, 0.02, 81))
#   MM2_epg_bkg_outb, bins = np.histogram(df_sim_bkg_outbs.loc[df_sim_bkg_outbs.config == config].MM2_epg, weights = df_sim_bkg_outbs.loc[df_sim_bkg_outbs.config == config].eff_bkg_merging_outb * df_sim_bkg_outbs.loc[df_sim_bkg_outbs.config == config].efficiency_pi0 * df_sim_bkg_outbs.loc[df_sim_bkg_outbs.config == config].weights, bins = np.linspace(-0.02, 0.02, 81))

#   axs[config - 1].hist(bins[:-1], bins, weights = MM2_epg_exp_inb + MM2_epg_exp_outb, histtype = 'step', label = r'$\mathrm{Data}$', color = 'k')
#   # axs[config - 1].hist(bins[:-1], bins, weights = MM2_epg_sim_inb + MM2_epg_sim_outb, histtype = 'step', label = 'dvcs')
#   axs[config - 1].hist(bins[:-1], bins, weights = MM2_epg_sim_inb + MM2_epg_bkg_inb + MM2_epg_sim_outb + MM2_epg_bkg_outb, histtype = 'step', label = r'$\mathrm{MC}$', color = 'tab:orange')
#   axs[config - 1].hist(bins[:-1], bins, weights = MM2_epg_bkg_inb + MM2_epg_bkg_outb, histtype = 'step', label = r'$\mathrm{Bkg~MC}$', color = 'tab:blue')

#   axs[config -1].set_xlim([-.02, 0.02])
#   axs[config -1].axvline(.1349766**2, color = 'k', ls = '--', label = r"$m_{\pi^0}^2$")
# # axs[0].set_ylim([0, 3500])
#   axs[config -1].set_xticks([-0.02, -0.01, 0.000, 0.01, 0.02])
#   if config !=3:
#     axs[config -1].set_xticklabels(['']*5)
#   axs[config -1].set_xticks(np.linspace(-0.02, 0.02, 40+1), minor = True)
#   axs[config -1].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'major', length = 10)
#   axs[config -1].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'minor', length = 5)


# axs[0].set_yticks([0, 500, 1000, 1500], ['', r'$5\times10^2$', r'$\times10^3$', r'$1.5\times10^3$'])
# axs[0].set_yticks(np.linspace(0, 1500, 16), minor = True)

# axs[1].set_yticks([0, 1000, 2000, 3000, 4000], ['', r'$10^3$', r'$2\times10^3$', r'$3\times10^3$', r'$4\times10^3$'])
# axs[1].set_yticks(np.linspace(0, 4700, 48), minor = True)

# axs[2].set_yticks([0, 10000, 20000], ['', '$10^4$', r'$2\times10^4$'])
# axs[2].set_yticks(np.linspace(0, 23000, 24), minor = True)


# axs[0].legend(title = '', loc = 'upper right', bbox_to_anchor = (0.95, 0.95),  prop={'size': 20})

# axs[0].annotate(xy = (0.05, 0.85), xytext = (0.05, 0.85), text = r"$\mathrm{\mathbf{a.}}~p'~\mathrm{in~FD},~\gamma~\mathrm{in~FD}$", xycoords = 'axes fraction', fontsize = 20)
# axs[1].annotate(xy = (0.05, 0.85), xytext = (0.05, 0.85), text = r"$\mathrm{\mathbf{b.}}~p'~\mathrm{in~CD},~\gamma~\mathrm{in~FD}$", xycoords = 'axes fraction', fontsize = 20)
# axs[2].annotate(xy = (0.05, 0.85), xytext = (0.05, 0.85), text = r"$\mathrm{\mathbf{c.}}~p'~\mathrm{in~CD},~\gamma~\mathrm{in~FT}$", xycoords = 'axes fraction', fontsize = 20)

# axs[2].set_xlabel(r"$\mathrm{M}^2_{X}~(\mathrm{GeV}^2)$")
# plt.subplots_adjust(hspace=0)
# fig.text(-0.1 , 0.5, r"$\mathrm{Events}/(5\times 10^{-4}~\mathrm{GeV}^2)$", va='center', rotation='vertical')

# plt.savefig("addendum_v3/M2_X_distribution.new.pdf", bbox_inches = 'tight')
# plt.close()

# exp_weights_inb   = []
# sig_weights_inb   = []
# sig_weights_bh_inb   = []
# sig_weights_vgg_inb   = []
# bkg_weights_inb   = []
# exp_weights_stat_err_inb   = []
# sig_weights_stat_err_inb   = []
# bkg_weights_stat_err_inb   = []

# # sig_weights_syst_err_inb   = []
# sig_weights_syst_err_inb_up     = []
# sig_weights_syst_err_inb_down   = []
# bkg_weights_syst_err_inb   = []

# MM2_epg_exp_inb    = []
# MM2_epg_sig_inb    = []
# MM2_epg_bkg_inb    = []

# # for integrated_binnum, phi_binnum, phi_width in zip(*df_summary_table_rebinned.loc[:, ["integrated_binnum", "phi_binnum", "phi_width"]].to_numpy().T):
# for i in range(len(df_summary_table_rebinned)):
#     this_bin = df_summary_table_rebinned.iloc[i]
#     integrated_binnum = this_bin.integrated_binnum
#     phi_binnum = this_bin.phi_binnum
#     phi_width  = this_bin.phi_width
#     contamination_inb = this_bin.contamination_inb
#     eff_bkg_merging       = this_bin.eff_bkg_merging_inb
#     active_bin_inb        = this_bin.active_bin_inb
#     if active_bin_inb == 0:
#         continue
#     c_stat_ratio = this_bin.pi0_inb_exp_to_sim_stat_err_ratio
#     c_syst_ratio = this_bin.pi0_inb_exp_to_sim_syst_err_ratio
#     sig_syst_ratio_up   = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_inb_sim_syst_err_up_ratio
#     sig_syst_ratio_down = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_inb_sim_syst_err_down_ratio

#     df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "contamination"]   = contamination_inb
#     df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "phi_rebinnum"]    = phi_binnum
#     df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "phi_rebinwidth"]  = phi_width
#     df_dvcs_exp = df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), :]
#     MM2_epg_exp_inb.extend(df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
#     weight_exp = inverseHist(df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width)].efficiency.to_numpy())
#     exp_weights_inb.extend(weight_exp)
#     exp_weights_stat_err_inb.extend(weight_exp)

#     MM2_epg_bkg_inb.extend(df_sim_bkg_inbs.loc[(df_sim_bkg_inbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_inbs.phi_binnum >= phi_binnum) & (df_sim_bkg_inbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
#     # weight_bkg = np.array([c*np.sum(weight_exp)/len(df_bkg_sim)]*len(df_bkg_sim))
#     df_bkg_inb_efficiency = df_sim_bkg_inbs.loc[(df_sim_bkg_inbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_inbs.phi_binnum >= phi_binnum) & (df_sim_bkg_inbs.phi_binnum < phi_binnum + phi_width), "efficiency_pi0"]
#     weight_bkg = np.array(contamination_inb*np.sum(weight_exp)*df_bkg_inb_efficiency/np.sum(df_bkg_inb_efficiency))
#     bkg_weights_inb.extend(weight_bkg)
#     bkg_weights_syst_err_inb.extend(weight_bkg*c_syst_ratio)
#     try:
#         pi0_1gamma_stat_err_squared = 1/len(df_bkg_inb_efficiency)
#     except:
#         pi0_1gamma_stat_err_squared = 0
#     bkg_weights_stat_err_inb.extend(weight_bkg*np.sqrt(c_stat_ratio**2 + 1/len(df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), :]) + pi0_1gamma_stat_err_squared ) )

#     df_dvcs_sim = df_sim_dvcs_inbs.loc[(df_sim_dvcs_inbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_inbs.phi_binnum >= phi_binnum) & (df_sim_dvcs_inbs.phi_binnum < phi_binnum + phi_width), :]
#     MM2_epg_sig_inb.extend(df_dvcs_sim.MM2_epg.to_numpy())
#     weight_sig = ((1-contamination_inb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_pi0/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency)*np.sum(weight_exp)).to_numpy()
#     if contamination_inb <1 :
#         sig_weights_inb.extend(weight_sig)
#         sig_weights_bh_inb.extend(weight_sig_bh)
#         sig_weights_vgg_inb.extend(weight_sig_vgg)
#         # sig_weights_syst_err_inb.extend(weight_sig*c_syst_ratio)
#         sig_weights_syst_err_inb_up.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_up**2))
#         sig_weights_syst_err_inb_down.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_down**2))
#         if len(df_bkg_inb_efficiency):   
#             sig_weights_stat_err_inb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_inb)**2 + 1/len(df_dvcs_exp) + 1/len(df_bkg_inb_efficiency) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
#         else:
#             sig_weights_stat_err_inb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_inb)**2 + 1/len(df_dvcs_exp) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
#     else:
#         sig_weights_inb.extend([0]*len(weight_sig))
#         sig_weights_bh_inb.extend([0]*len(weight_sig_bh))
#         sig_weights_vgg_inb.extend([0]*len(weight_sig_vgg))
#         sig_weights_syst_err_inb.extend([0]*len(weight_sig))
#         sig_weights_stat_err_inb.extend([0]*len(weight_sig))


# exp_weights_outb   = []
# sig_weights_outb   = []
# sig_weights_bh_outb   = []
# sig_weights_vgg_outb   = []
# bkg_weights_outb   = []
# exp_weights_stat_err_outb   = []
# sig_weights_stat_err_outb   = []
# bkg_weights_stat_err_outb   = []

# # sig_weights_syst_err_outb   = []
# sig_weights_syst_err_outb_up     = []
# sig_weights_syst_err_outb_down   = []
# bkg_weights_syst_err_outb   = []
# bkg_weights_syst_err_outb   = []

# MM2_epg_exp_outb    = []
# MM2_epg_sig_outb    = []
# MM2_epg_bkg_outb    = []

# # for integrated_binnum, phi_binnum, phi_width in zip(*df_summary_table_rebinned.loc[:, ["integrated_binnum", "phi_binnum", "phi_width"]].to_numpy().T):
# for i in range(len(df_summary_table_rebinned)):
#     this_bin = df_summary_table_rebinned.iloc[i]
#     integrated_binnum = this_bin.integrated_binnum
#     phi_binnum = this_bin.phi_binnum
#     phi_width  = this_bin.phi_width
#     contamination_outb = this_bin.contamination_outb
#     eff_bkg_merging       = this_bin.eff_bkg_merging_outb
#     active_bin_outb        = this_bin.active_bin_outb
#     if active_bin_outb == 0:
#         continue
#     c_stat_ratio = this_bin.pi0_outb_exp_to_sim_stat_err_ratio
#     c_syst_ratio = this_bin.pi0_outb_exp_to_sim_syst_err_ratio
#     sig_syst_ratio_up   = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_outb_sim_syst_err_up_ratio
#     sig_syst_ratio_down = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_outb_sim_syst_err_down_ratio

#     df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "contamination"]   = contamination_outb
#     df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "phi_rebinnum"]    = phi_binnum
#     df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "phi_rebinwidth"]  = phi_width
#     df_dvcs_exp = df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), :]
#     MM2_epg_exp_outb.extend(df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
#     weight_exp = inverseHist(df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width)].efficiency.to_numpy())
#     exp_weights_outb.extend(weight_exp)
#     exp_weights_stat_err_outb.extend(weight_exp)

#     MM2_epg_bkg_outb.extend(df_sim_bkg_outbs.loc[(df_sim_bkg_outbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_outbs.phi_binnum >= phi_binnum) & (df_sim_bkg_outbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
#     # weight_bkg = np.array([c*np.sum(weight_exp)/len(df_bkg_sim)]*len(df_bkg_sim))
#     df_bkg_outb_efficiency = df_sim_bkg_outbs.loc[(df_sim_bkg_outbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_outbs.phi_binnum >= phi_binnum) & (df_sim_bkg_outbs.phi_binnum < phi_binnum + phi_width), "efficiency_pi0"]
#     weight_bkg = np.array(contamination_outb*np.sum(weight_exp)*df_bkg_outb_efficiency/np.sum(df_bkg_outb_efficiency))
#     bkg_weights_outb.extend(weight_bkg)
#     bkg_weights_syst_err_outb.extend(weight_bkg*c_syst_ratio)
#     try:
#         pi0_1gamma_stat_err_squared = 1/len(df_bkg_outb_efficiency)
#     except:
#         pi0_1gamma_stat_err_squared = 0
#     bkg_weights_stat_err_outb.extend(weight_bkg*np.sqrt(c_stat_ratio**2 + 1/len(df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), :]) + pi0_1gamma_stat_err_squared ) )

#     df_dvcs_sim = df_sim_dvcs_outbs.loc[(df_sim_dvcs_outbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_outbs.phi_binnum >= phi_binnum) & (df_sim_dvcs_outbs.phi_binnum < phi_binnum + phi_width), :]
#     MM2_epg_sig_outb.extend(df_dvcs_sim.MM2_epg.to_numpy())
#     # weight_sig = ((1-contamination_outb)*df_dvcs_sim.weights/np.sum(df_dvcs_sim.weights)*np.sum(weight_exp)).to_numpy()
#     weight_sig = ((1-contamination_outb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_pi0/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency)*np.sum(weight_exp)).to_numpy()
#     weight_sig_bh = ((1-contamination_outb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_bh/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency_bh)*np.sum(weight_exp)).to_numpy()
#     weight_sig_vgg = ((1-contamination_outb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_vgg/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency_vgg)*np.sum(weight_exp)).to_numpy()
#     # weight_sig = eff_bkg_merging*df_dvcs_sim.weights*df_dvcs_sim.efficiency
#     if contamination_outb <1 :
#         sig_weights_outb.extend(weight_sig)
#         sig_weights_bh_outb.extend(weight_sig_bh)
#         sig_weights_vgg_outb.extend(weight_sig_vgg)
#         # sig_weights_syst_err_outb.extend(weight_sig*c_syst_ratio)
#         sig_weights_syst_err_outb_up.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_up**2))
#         sig_weights_syst_err_outb_down.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_down**2))
#         if len(df_bkg_outb_efficiency):
#             sig_weights_stat_err_outb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_outb)**2 + 1/len(df_dvcs_exp) + 1/len(df_bkg_outb_efficiency) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
#         else:
#             sig_weights_stat_err_outb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_outb)**2 + 1/len(df_dvcs_exp) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
#     else:
#         sig_weights_outb.extend([0]*len(weight_sig))
#         sig_weights_bh_outb.extend([0]*len(weight_sig_bh))
#         sig_weights_vgg_outb.extend([0]*len(weight_sig_vgg))
#         sig_weights_syst_err_outb.extend([0]*len(weight_sig))
#         sig_weights_stat_err_outb.extend([0]*len(weight_sig))

# fig, axs = plt.subplots(2, 1, figsize = (8, 8), height_ratios=[2, 1])
# MM2_epg_exp_inb_hist, bins = np.histogram(np.array(MM2_epg_exp_inb).flatten(), weights = exp_weights_inb, bins = np.linspace(-0.02, 0.02, 81))
# MM2_epg_bkg_inb_hist, _ = np.histogram(np.array(MM2_epg_bkg_inb).flatten(), weights = bkg_weights_inb, bins = bins)
# MM2_epg_sig_inb_hist, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_inb, bins = bins)
# MM2_epg_sim_inb_hist  = MM2_epg_bkg_inb_hist + MM2_epg_sig_inb_hist
# MM2_epg_exp_outb_hist, _ = np.histogram(np.array(MM2_epg_exp_outb).flatten(), weights = exp_weights_outb, bins = bins)
# MM2_epg_bkg_outb_hist, _ = np.histogram(np.array(MM2_epg_bkg_outb).flatten(), weights = bkg_weights_outb, bins = bins)
# MM2_epg_sig_outb_hist, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_outb, bins = bins)
# MM2_epg_sim_outb_hist  = MM2_epg_bkg_outb_hist + MM2_epg_sig_outb_hist

# MM2_epg_exp_inb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_exp_inb).flatten(), weights = exp_weights_stat_err_inb, bins = bins)
# MM2_epg_bkg_inb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_bkg_inb).flatten(), weights = bkg_weights_stat_err_inb, bins = bins)
# MM2_epg_sig_inb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_stat_err_inb, bins = bins)
# MM2_epg_exp_outb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_exp_outb).flatten(), weights = exp_weights_stat_err_outb, bins = bins)
# MM2_epg_bkg_outb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_bkg_outb).flatten(), weights = bkg_weights_stat_err_outb, bins = bins)
# MM2_epg_sig_outb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_stat_err_outb, bins = bins)

# MM2_epg_exp_inb_hist_stat_err  = np.sqrt(MM2_epg_exp_inb_hist_stat_err)
# MM2_epg_bkg_inb_hist_stat_err  = np.sqrt(MM2_epg_bkg_inb_hist_stat_err)
# MM2_epg_sig_inb_hist_stat_err  = np.sqrt(MM2_epg_sig_inb_hist_stat_err)
# MM2_epg_exp_outb_hist_stat_err = np.sqrt(MM2_epg_exp_outb_hist_stat_err)
# MM2_epg_bkg_outb_hist_stat_err = np.sqrt(MM2_epg_bkg_outb_hist_stat_err)
# MM2_epg_sig_outb_hist_stat_err = np.sqrt(MM2_epg_sig_outb_hist_stat_err)

# MM2_epg_bkg_inb_hist_syst_err, _ = np.histogram(np.array(MM2_epg_bkg_inb).flatten(), weights = bkg_weights_syst_err_inb, bins = bins)
# MM2_epg_sig_inb_hist_syst_err_up, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_syst_err_inb_up, bins = bins)
# MM2_epg_sig_inb_hist_syst_err_down, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_syst_err_inb_down, bins = bins)
# MM2_epg_bkg_outb_hist_syst_err, _ = np.histogram(np.array(MM2_epg_bkg_outb).flatten(), weights = bkg_weights_syst_err_outb, bins = bins)
# MM2_epg_sig_outb_hist_syst_err_up, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_syst_err_outb_up, bins = bins)
# MM2_epg_sig_outb_hist_syst_err_down, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_syst_err_outb_down, bins = bins)

# MM2_epg_exp_hist         = MM2_epg_exp_inb_hist + MM2_epg_exp_outb_hist
# MM2_epg_bkg_hist         = MM2_epg_bkg_inb_hist + MM2_epg_bkg_outb_hist
# MM2_epg_sig_hist         = MM2_epg_sig_inb_hist + MM2_epg_sig_outb_hist
# MM2_epg_sim_hist         = MM2_epg_bkg_hist     + MM2_epg_sig_hist

# MM2_epg_exp_hist_stat_err     = np.sqrt(MM2_epg_exp_inb_hist_stat_err**2 + MM2_epg_exp_outb_hist_stat_err**2)
# MM2_epg_bkg_hist_stat_err     = np.sqrt(MM2_epg_bkg_inb_hist_stat_err**2 + MM2_epg_bkg_outb_hist_stat_err**2)
# MM2_epg_sig_hist_stat_err     = np.sqrt(MM2_epg_sig_inb_hist_stat_err**2 + MM2_epg_sig_outb_hist_stat_err**2)

# MM2_epg_bkg_hist_syst_err        = np.sqrt(MM2_epg_bkg_inb_hist_syst_err**2 + MM2_epg_bkg_outb_hist_syst_err**2)
# MM2_epg_sig_hist_syst_err_up     = np.sqrt(MM2_epg_sig_inb_hist_syst_err_up**2 + MM2_epg_sig_outb_hist_syst_err_up**2)
# MM2_epg_sig_hist_syst_err_down   = np.sqrt(MM2_epg_sig_inb_hist_syst_err_down**2 + MM2_epg_sig_outb_hist_syst_err_down**2)

# MM2_epg_sim_hist_stat_err        = np.sqrt(MM2_epg_bkg_hist_stat_err**2 + MM2_epg_sig_hist_stat_err**2)
# MM2_epg_sim_hist_syst_err_up     = np.sqrt(MM2_epg_bkg_hist_syst_err**2 + MM2_epg_sig_hist_syst_err_up**2)
# MM2_epg_sim_hist_syst_err_down   = np.sqrt(MM2_epg_bkg_hist_syst_err**2 + MM2_epg_sig_hist_syst_err_down**2)
# MM2_epg_sim_hist_syst_err        = (MM2_epg_sim_hist_syst_err_up + MM2_epg_sim_hist_syst_err_down)/2.

# MM2_epg_exp_hist_stat_err_ratio      = divideHist(MM2_epg_exp_hist_stat_err, MM2_epg_exp_hist)
# MM2_epg_sim_hist_stat_err_ratio      = divideHist(MM2_epg_sim_hist_stat_err, MM2_epg_sim_hist)
# MM2_epg_sim_hist_syst_err_ratio      = divideHist(MM2_epg_sim_hist_syst_err, MM2_epg_sim_hist)
# MM2_epg_sim_hist_syst_err_ratio_up   = divideHist(MM2_epg_sim_hist_syst_err_up, MM2_epg_sim_hist)
# MM2_epg_sim_hist_syst_err_ratio_down = divideHist(MM2_epg_sim_hist_syst_err_down, MM2_epg_sim_hist)

# MM2_epg_sim_hist_down    = MM2_epg_sim_hist - MM2_epg_sim_hist_syst_err
# MM2_epg_sim_hist_up      = MM2_epg_sim_hist + MM2_epg_sim_hist_syst_err

# bincenters = (bins[1:] + bins[:-1])/2.

# axs[0].hist(bins[:-1], bins = bins, weights = MM2_epg_exp_hist, histtype = 'step', color = 'k', label = "$\mathrm{Data}$")
# axs[0].hist(bins[:-1], bins = bins, weights = MM2_epg_sim_hist, histtype = 'step', color = 'tab:blue', label = '$S+B~\mathrm{(Sim.)}$')

# bin_fill_between = []
# MM2_epg_sim_hist_down_fill_between      = []
# MM2_epg_sim_hist_up_fill_between      = []
# for i in range(len(bins)):
#     if (i > 0):
#         bin_fill_between.append(bins[i])
#         MM2_epg_sim_hist_down_fill_between.append(MM2_epg_sim_hist_down[i-1])
#         MM2_epg_sim_hist_up_fill_between.append(MM2_epg_sim_hist_up[i-1])
#     if (i<len(bins)-1):
#         bin_fill_between.append(bins[i])
#         MM2_epg_sim_hist_down_fill_between.append(MM2_epg_sim_hist_down[i])
#         MM2_epg_sim_hist_up_fill_between.append(MM2_epg_sim_hist_up[i])

# # axs[0].fill_between(bin_fill_between, MM2_epg_sim_hist_down_fill_between, MM2_epg_sim_hist_up_fill_between, color = 'tab:blue', alpha = 0.5)

# axs[0].hist(bins[:-1], bins = bins, weights = 5*(MM2_epg_bkg_inb_hist + MM2_epg_bkg_outb_hist), histtype = 'step', color = 'tab:orange', label = "$5 \\times B~\mathrm{(Sim.)}$", zorder = -1)

# axs[1].errorbar(bincenters[8:-8], divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist)[8:-8], yerr = (divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist) * np.sqrt(MM2_epg_exp_hist_stat_err_ratio**2 + MM2_epg_sim_hist_stat_err_ratio**2))[8:-8]
#                , color = 'k', marker = 'o', ls = '')
# # axs[1].fill_between(bincenters, divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist) * (1 -  MM2_epg_sim_hist_syst_err_ratio), divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist)  * (1 +  MM2_epg_sim_hist_syst_err_ratio), color = 'r', alpha = 0.3)
# axs[1].fill_between(bincenters[8:-8], (1 -  MM2_epg_sim_hist_syst_err_ratio)[8:-8], (1 + MM2_epg_sim_hist_syst_err_ratio)[8:-8], color = 'k', alpha = 0.3)

# # axs[1].errorbar(bincenters, divideHist(MM2_epg_exp_inb_hist, MM2_epg_sim_inb_hist), color = 'k', marker = 'o', ls = '')
# # axs[1].errorbar(bincenters, divideHist(MM2_epg_exp_outb_hist, MM2_epg_sim_outb_hist), color = 'r', marker = 'o', ls = '')

# # axs[1].errorbar(bincenters, MM2_epg_exp_inb_hist/MM2_epg_sim_inb_hist#, yerr = (MM2_epg_exp_inb_hist/MM2_epg_sim_inb_hist) * np.sqrt(MM2_epg_exp_inb_hist_stat_err_ratio**2 + MM2_epg_sim_inb_hist_stat_err_ratio**2)
# #                , color = 'k', marker = 'o', ls = '')
# # axs[1].errorbar(bincenters, MM2_epg_exp_outb_hist/MM2_epg_sim_outb_hist#, yerr = (MM2_epg_exp_inb_hist/MM2_epg_sim_inb_hist) * np.sqrt(MM2_epg_exp_inb_hist_stat_err_ratio**2 + MM2_epg_sim_inb_hist_stat_err_ratio**2)
# #                , color = 'r', marker = 'o', ls = '')


# hist_pi0_inb, _ = np.histogram(df_exp_pi0_inbs.Mpi0**2, bins = bins)
# hist_pi0_outb, _ = np.histogram(df_exp_pi0_outbs.Mpi0**2, bins = bins)

# axs[0].axvline(.1349766**2, color = 'k', ls = '--', label = r"$m_{\pi^0}^2$")

# axs[0].set_xlim([-.02, 0.02])
# # axs[0].set_ylim([0, 3500])
# axs[0].set_xticks([-0.02, -0.01, 0.000, 0.01, 0.02])
# axs[0].set_xticklabels(['']*5)
# axs[0].set_xticks(np.linspace(-0.02, 0.02, 40+1), minor = True)

# axs[0].set_yticks([0, 10000, 20000, 30000], ['', '$1$', '$2$', '$3$'])
# # axs[0].set_yticklabels(['', r'$5\times 10^3$', r'$10^4$', r'$3000$'])
# axs[0].set_yticks(np.linspace(0, 30000, 31), minor = True)
# axs[0].set_ylabel(r"$\mathrm{Events}/(5\times 10^{-4}~\mathrm{GeV}^2)$")

# axs[0].annotate(xy = (-0.02, 1.03), xytext = (-0.02, 1.03), text = r'$\times 10^4$', xycoords = 'axes fraction', fontsize = 20)

# axs[1].set_xlim([-.02, 0.02])
# axs[1].set_ylim([0.5, 1.5])
# axs[1].set_xticks([-0.02, -0.01, 0, 0.01, 0.02])
# axs[1].set_xticklabels([r'${:.3f}$'.format(i) for i in [-0.02, -0.01, 0, 0.01, 0.02]])
# axs[1].set_xticks(np.linspace(-0.02, 0.02, 40+1), minor = True)
# axs[1].set_yticks(np.linspace(0.5, 1.5, 11), minor = True)

# axs[0].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'major', length = 10)
# axs[0].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'minor', length = 5)
# axs[1].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'major', length = 10)
# axs[1].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'minor', length = 5)
# axs[1].set_ylabel("$\mathrm{Data/Sim.}$", labelpad = 17)

# axs[0].legend(title = '', loc = 'upper right', bbox_to_anchor = (0.95, 0.95),  prop={'size': 15})

# axs[1].axhline(1, ls = '--', color = 'k')
# axs[1].axvline(.1349766**2, color = 'k', ls = '--', label = r"$m_{\pi^0}^2$")

# axs[1].set_xlabel(r"$\mathrm{M}^2_{X}~(\mathrm{GeV}^2)$")

# # axs[1].set_ylim([0.9, 1.1])

# # hist_exp_inb, bins = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector>9].MM2_epg, weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector>9].signal, bins = np.linspace(-0.005, 0.005, 101))
# # hist_sim_inb, _    = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector>9].MM2_epg, weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector>9].signal, bins = bins)

# # hist_exp_outb, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector>9].MM2_epg,  weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector>9].signal, bins = bins)
# # hist_sim_outb, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector>9].MM2_epg, weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector>9].signal, bins = bins)

# # hist_exp = hist_exp_inb + hist_exp_outb
# # hist_sim = hist_sim_inb + hist_sim_outb

# # plt.hist(bins[:-1], bins, weights =  hist_exp_inb, histtype = 'step')
# # plt.hist(bins[:-1], bins, weights =  hist_sim_inb, histtype = 'step')
# # plt.hist(bins[:-1], bins, weights =  hist_exp_outb, histtype = 'step')
# # plt.hist(bins[:-1], bins, weights =  hist_sim_outb, histtype = 'step')
# # plt.hist(bins[:-1], bins, weights =  (hist_exp_inb ) / (hist_sim_inb), histtype = 'step')
# # plt.hist(bins[:-1], bins, weights =  (hist_exp_outb ) / (hist_sim_outb), histtype = 'step')
# # axs[1].hist(bins[:-1], bins, weights =  (hist_exp_inb + hist_exp_outb) / (hist_sim_inb + hist_sim_outb), histtype = 'stepfilled')

# plt.subplots_adjust(hspace=0.05)
# # plt.savefig("MM2_epg_distribution.pdf", bbox_inches = 'tight')
# plt.savefig("addendum_v3/M2_X_distribution.new.3.pdf", bbox_inches = 'tight')
# plt.close()


'''
old figures for draft
df_exp_epg_inbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_inb/dvcs/excl_level_2/pkl_7_nominal/fall2018_inb.pkl")
df_exp_pi0_inbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_inb/pi0/excl_level_2/pkl_7_nominal/fall2018_inb.pkl")

df_sim_bkg_inbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/pi0_1gamma/excl_level_1/pkl_7_nominal/{}/fall2018_inb.pkl".format(i)) for i in chunks_inb])
df_sim_pi0_inbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/pi0_2gamma/excl_level_1/pkl_7_nominal/{}/fall2018_inb.pkl".format(i)) for i in chunks_inb])
df_sim_dvcs_inbs  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl_7_nominal/fall2018_inb.pkl")
df_sim_dvcs_inbs_bkgmerging  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl_14_bkgmerging/fall2018_inb.pkl")

df_exp_epg_outbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_outb/dvcs/excl_level_2/pkl_7_nominal/fall2018_outb.pkl")
df_exp_pi0_outbs   = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/exp_fall2018_outb/pi0/excl_level_2/pkl_7_nominal/fall2018_outb.pkl")

df_sim_bkg_outbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/pi0_1gamma/excl_level_1/pkl_7_nominal/{}/fall2018_outb.pkl".format(i)) for i in chunks_outb])
df_sim_pi0_outbs   = pd.concat([pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/pi0_2gamma/excl_level_1/pkl_7_nominal/{}/fall2018_outb.pkl".format(i)) for i in chunks_outb])
df_sim_dvcs_outbs  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl_7_nominal/fall2018_outb.pkl")
df_sim_dvcs_outbs_bkgmerging  = pd.read_pickle("/Users/sangbaek.lee/CLAS12/clas12DVCS/review_meeting/data/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl_14_bkgmerging/fall2018_outb.pkl")

exp_weights_inb   = []
sig_weights_inb   = []
sig_weights_bh_inb   = []
sig_weights_vgg_inb   = []
bkg_weights_inb   = []
exp_weights_stat_err_inb   = []
sig_weights_stat_err_inb   = []
bkg_weights_stat_err_inb   = []

# sig_weights_syst_err_inb   = []
sig_weights_syst_err_inb_up     = []
sig_weights_syst_err_inb_down   = []
bkg_weights_syst_err_inb   = []

MM2_epg_exp_inb    = []
MM2_epg_sig_inb    = []
MM2_epg_bkg_inb    = []

# for integrated_binnum, phi_binnum, phi_width in zip(*df_summary_table_rebinned.loc[:, ["integrated_binnum", "phi_binnum", "phi_width"]].to_numpy().T):
for i in range(len(df_summary_table_rebinned)):
    this_bin = df_summary_table_rebinned.iloc[i]
    integrated_binnum = this_bin.integrated_binnum
    phi_binnum = this_bin.phi_binnum
    phi_width  = this_bin.phi_width
    contamination_inb = this_bin.contamination_inb
    eff_bkg_merging       = this_bin.eff_bkg_merging_inb
    active_bin_inb        = this_bin.active_bin_inb
    if active_bin_inb == 0:
        continue
    c_stat_ratio = this_bin.pi0_inb_exp_to_sim_stat_err_ratio
    c_syst_ratio = this_bin.pi0_inb_exp_to_sim_syst_err_ratio
    sig_syst_ratio_up   = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_inb_sim_syst_err_up_ratio
    sig_syst_ratio_down = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_inb_sim_syst_err_down_ratio

    df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "contamination"]   = contamination_inb
    df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "phi_rebinnum"]    = phi_binnum
    df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "phi_rebinwidth"]  = phi_width
    df_dvcs_exp = df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), :]
    MM2_epg_exp_inb.extend(df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
    weight_exp = inverseHist(df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width)].efficiency.to_numpy())
    exp_weights_inb.extend(weight_exp)
    exp_weights_stat_err_inb.extend(weight_exp)

    MM2_epg_bkg_inb.extend(df_sim_bkg_inbs.loc[(df_sim_bkg_inbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_inbs.phi_binnum >= phi_binnum) & (df_sim_bkg_inbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
    # weight_bkg = np.array([c*np.sum(weight_exp)/len(df_bkg_sim)]*len(df_bkg_sim))
    df_bkg_inb_efficiency = df_sim_bkg_inbs.loc[(df_sim_bkg_inbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_inbs.phi_binnum >= phi_binnum) & (df_sim_bkg_inbs.phi_binnum < phi_binnum + phi_width), "efficiency"]
    weight_bkg = np.array(contamination_inb*np.sum(weight_exp)*df_bkg_inb_efficiency/np.sum(df_bkg_inb_efficiency))
    bkg_weights_inb.extend(weight_bkg)
    bkg_weights_syst_err_inb.extend(weight_bkg*c_syst_ratio)
    try:
        pi0_1gamma_stat_err_squared = 1/len(df_bkg_inb_efficiency)
    except:
        pi0_1gamma_stat_err_squared = 0
    bkg_weights_stat_err_inb.extend(weight_bkg*np.sqrt(c_stat_ratio**2 + 1/len(df_exp_epg_inbs.loc[(df_exp_epg_inbs.integrated_binnum == integrated_binnum) & (df_exp_epg_inbs.phi_binnum >= phi_binnum) & (df_exp_epg_inbs.phi_binnum < phi_binnum + phi_width), :]) + pi0_1gamma_stat_err_squared ) )

    df_dvcs_sim = df_sim_dvcs_inbs.loc[(df_sim_dvcs_inbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_inbs.phi_binnum >= phi_binnum) & (df_sim_dvcs_inbs.phi_binnum < phi_binnum + phi_width), :]
    MM2_epg_sig_inb.extend(df_dvcs_sim.MM2_epg.to_numpy())
    # weight_sig = ((1-contamination_inb)*df_dvcs_sim.weights/np.sum(df_dvcs_sim.weights)*np.sum(weight_exp)).to_numpy()
    weight_sig = ((1-contamination_inb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency)*np.sum(weight_exp)).to_numpy()
    weight_sig_bh = ((1-contamination_inb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_bh/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency_bh)*np.sum(weight_exp)).to_numpy()
    weight_sig_vgg = ((1-contamination_inb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_vgg/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency_vgg)*np.sum(weight_exp)).to_numpy()
    # weight_sig = eff_bkg_merging*df_dvcs_sim.weights*df_dvcs_sim.efficiency
    if contamination_inb <1 :
        sig_weights_inb.extend(weight_sig)
        sig_weights_bh_inb.extend(weight_sig_bh)
        sig_weights_vgg_inb.extend(weight_sig_vgg)
        # sig_weights_syst_err_inb.extend(weight_sig*c_syst_ratio)
        sig_weights_syst_err_inb_up.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_up**2))
        sig_weights_syst_err_inb_down.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_down**2))
        if len(df_bkg_inb_efficiency):   
            sig_weights_stat_err_inb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_inb)**2 + 1/len(df_dvcs_exp) + 1/len(df_bkg_inb_efficiency) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
        else:
            sig_weights_stat_err_inb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_inb)**2 + 1/len(df_dvcs_exp) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
    else:
        sig_weights_inb.extend([0]*len(weight_sig))
        sig_weights_bh_inb.extend([0]*len(weight_sig_bh))
        sig_weights_vgg_inb.extend([0]*len(weight_sig_vgg))
        sig_weights_syst_err_inb.extend([0]*len(weight_sig))
        sig_weights_stat_err_inb.extend([0]*len(weight_sig))


exp_weights_outb   = []
sig_weights_outb   = []
sig_weights_bh_outb   = []
sig_weights_vgg_outb   = []
bkg_weights_outb   = []
exp_weights_stat_err_outb   = []
sig_weights_stat_err_outb   = []
bkg_weights_stat_err_outb   = []

# sig_weights_syst_err_outb   = []
sig_weights_syst_err_outb_up     = []
sig_weights_syst_err_outb_down   = []
bkg_weights_syst_err_outb   = []
bkg_weights_syst_err_outb   = []

MM2_epg_exp_outb    = []
MM2_epg_sig_outb    = []
MM2_epg_bkg_outb    = []

# for integrated_binnum, phi_binnum, phi_width in zip(*df_summary_table_rebinned.loc[:, ["integrated_binnum", "phi_binnum", "phi_width"]].to_numpy().T):
for i in range(len(df_summary_table_rebinned)):
    this_bin = df_summary_table_rebinned.iloc[i]
    integrated_binnum = this_bin.integrated_binnum
    phi_binnum = this_bin.phi_binnum
    phi_width  = this_bin.phi_width
    contamination_outb = this_bin.contamination_outb
    eff_bkg_merging       = this_bin.eff_bkg_merging_outb
    active_bin_outb        = this_bin.active_bin_outb
    if active_bin_outb == 0:
        continue
    c_stat_ratio = this_bin.pi0_outb_exp_to_sim_stat_err_ratio
    c_syst_ratio = this_bin.pi0_outb_exp_to_sim_syst_err_ratio
    sig_syst_ratio_up   = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_outb_sim_syst_err_up_ratio
    sig_syst_ratio_down = np.sqrt(0.3**2 + 0.0476**2)#this_bin.dvcs_outb_sim_syst_err_down_ratio

    df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "contamination"]   = contamination_outb
    df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "phi_rebinnum"]    = phi_binnum
    df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "phi_rebinwidth"]  = phi_width
    df_dvcs_exp = df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), :]
    MM2_epg_exp_outb.extend(df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
    weight_exp = inverseHist(df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width)].efficiency.to_numpy())
    exp_weights_outb.extend(weight_exp)
    exp_weights_stat_err_outb.extend(weight_exp)

    MM2_epg_bkg_outb.extend(df_sim_bkg_outbs.loc[(df_sim_bkg_outbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_outbs.phi_binnum >= phi_binnum) & (df_sim_bkg_outbs.phi_binnum < phi_binnum + phi_width), "MM2_epg"].to_numpy())
    # weight_bkg = np.array([c*np.sum(weight_exp)/len(df_bkg_sim)]*len(df_bkg_sim))
    df_bkg_outb_efficiency = df_sim_bkg_outbs.loc[(df_sim_bkg_outbs.integrated_binnum == integrated_binnum) & (df_sim_bkg_outbs.phi_binnum >= phi_binnum) & (df_sim_bkg_outbs.phi_binnum < phi_binnum + phi_width), "efficiency"]
    weight_bkg = np.array(contamination_outb*np.sum(weight_exp)*df_bkg_outb_efficiency/np.sum(df_bkg_outb_efficiency))
    bkg_weights_outb.extend(weight_bkg)
    bkg_weights_syst_err_outb.extend(weight_bkg*c_syst_ratio)
    try:
        pi0_1gamma_stat_err_squared = 1/len(df_bkg_outb_efficiency)
    except:
        pi0_1gamma_stat_err_squared = 0
    bkg_weights_stat_err_outb.extend(weight_bkg*np.sqrt(c_stat_ratio**2 + 1/len(df_exp_epg_outbs.loc[(df_exp_epg_outbs.integrated_binnum == integrated_binnum) & (df_exp_epg_outbs.phi_binnum >= phi_binnum) & (df_exp_epg_outbs.phi_binnum < phi_binnum + phi_width), :]) + pi0_1gamma_stat_err_squared ) )

    df_dvcs_sim = df_sim_dvcs_outbs.loc[(df_sim_dvcs_outbs.integrated_binnum == integrated_binnum) & (df_sim_dvcs_outbs.phi_binnum >= phi_binnum) & (df_sim_dvcs_outbs.phi_binnum < phi_binnum + phi_width), :]
    MM2_epg_sig_outb.extend(df_dvcs_sim.MM2_epg.to_numpy())
    # weight_sig = ((1-contamination_outb)*df_dvcs_sim.weights/np.sum(df_dvcs_sim.weights)*np.sum(weight_exp)).to_numpy()
    weight_sig = ((1-contamination_outb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency)*np.sum(weight_exp)).to_numpy()
    weight_sig_bh = ((1-contamination_outb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_bh/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency_bh)*np.sum(weight_exp)).to_numpy()
    weight_sig_vgg = ((1-contamination_outb)*df_dvcs_sim.weights*df_dvcs_sim.efficiency_vgg/np.sum(df_dvcs_sim.weights*df_dvcs_sim.efficiency_vgg)*np.sum(weight_exp)).to_numpy()
    # weight_sig = eff_bkg_merging*df_dvcs_sim.weights*df_dvcs_sim.efficiency
    if contamination_outb <1 :
        sig_weights_outb.extend(weight_sig)
        sig_weights_bh_outb.extend(weight_sig_bh)
        sig_weights_vgg_outb.extend(weight_sig_vgg)
        # sig_weights_syst_err_outb.extend(weight_sig*c_syst_ratio)
        sig_weights_syst_err_outb_up.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_up**2))
        sig_weights_syst_err_outb_down.extend(weight_sig*np.sqrt(c_syst_ratio**2 + sig_syst_ratio_down**2))
        if len(df_bkg_outb_efficiency):
            sig_weights_stat_err_outb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_outb)**2 + 1/len(df_dvcs_exp) + 1/len(df_bkg_outb_efficiency) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
        else:
            sig_weights_stat_err_outb.extend(weight_sig**2 * (c_stat_ratio**2/(1-contamination_outb)**2 + 1/len(df_dvcs_exp) + np.sum(df_dvcs_sim.weights**2)/np.sum(df_dvcs_sim.weights)**2) )
    else:
        sig_weights_outb.extend([0]*len(weight_sig))
        sig_weights_bh_outb.extend([0]*len(weight_sig_bh))
        sig_weights_vgg_outb.extend([0]*len(weight_sig_vgg))
        sig_weights_syst_err_outb.extend([0]*len(weight_sig))
        sig_weights_stat_err_outb.extend([0]*len(weight_sig))

fig, axs = plt.subplots(2, 1, figsize = (8, 8), height_ratios=[2, 1])
MM2_epg_exp_inb_hist, bins = np.histogram(np.array(MM2_epg_exp_inb).flatten(), weights = exp_weights_inb, bins = np.linspace(-0.02, 0.02, 81))
MM2_epg_bkg_inb_hist, _ = np.histogram(np.array(MM2_epg_bkg_inb).flatten(), weights = bkg_weights_inb, bins = bins)
MM2_epg_sig_inb_hist, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_inb, bins = bins)
MM2_epg_sim_inb_hist  = MM2_epg_bkg_inb_hist + MM2_epg_sig_inb_hist
MM2_epg_exp_outb_hist, _ = np.histogram(np.array(MM2_epg_exp_outb).flatten(), weights = exp_weights_outb, bins = bins)
MM2_epg_bkg_outb_hist, _ = np.histogram(np.array(MM2_epg_bkg_outb).flatten(), weights = bkg_weights_outb, bins = bins)
MM2_epg_sig_outb_hist, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_outb, bins = bins)
MM2_epg_sim_outb_hist  = MM2_epg_bkg_outb_hist + MM2_epg_sig_outb_hist

MM2_epg_exp_inb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_exp_inb).flatten(), weights = exp_weights_stat_err_inb, bins = bins)
MM2_epg_bkg_inb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_bkg_inb).flatten(), weights = bkg_weights_stat_err_inb, bins = bins)
MM2_epg_sig_inb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_stat_err_inb, bins = bins)
MM2_epg_exp_outb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_exp_outb).flatten(), weights = exp_weights_stat_err_outb, bins = bins)
MM2_epg_bkg_outb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_bkg_outb).flatten(), weights = bkg_weights_stat_err_outb, bins = bins)
MM2_epg_sig_outb_hist_stat_err, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_stat_err_outb, bins = bins)

MM2_epg_exp_inb_hist_stat_err  = np.sqrt(MM2_epg_exp_inb_hist_stat_err)
MM2_epg_bkg_inb_hist_stat_err  = np.sqrt(MM2_epg_bkg_inb_hist_stat_err)
MM2_epg_sig_inb_hist_stat_err  = np.sqrt(MM2_epg_sig_inb_hist_stat_err)
MM2_epg_exp_outb_hist_stat_err = np.sqrt(MM2_epg_exp_outb_hist_stat_err)
MM2_epg_bkg_outb_hist_stat_err = np.sqrt(MM2_epg_bkg_outb_hist_stat_err)
MM2_epg_sig_outb_hist_stat_err = np.sqrt(MM2_epg_sig_outb_hist_stat_err)

MM2_epg_bkg_inb_hist_syst_err, _ = np.histogram(np.array(MM2_epg_bkg_inb).flatten(), weights = bkg_weights_syst_err_inb, bins = bins)
MM2_epg_sig_inb_hist_syst_err_up, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_syst_err_inb_up, bins = bins)
MM2_epg_sig_inb_hist_syst_err_down, _ = np.histogram(np.array(MM2_epg_sig_inb).flatten(), weights = sig_weights_syst_err_inb_down, bins = bins)
MM2_epg_bkg_outb_hist_syst_err, _ = np.histogram(np.array(MM2_epg_bkg_outb).flatten(), weights = bkg_weights_syst_err_outb, bins = bins)
MM2_epg_sig_outb_hist_syst_err_up, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_syst_err_outb_up, bins = bins)
MM2_epg_sig_outb_hist_syst_err_down, _ = np.histogram(np.array(MM2_epg_sig_outb).flatten(), weights = sig_weights_syst_err_outb_down, bins = bins)

MM2_epg_exp_hist         = MM2_epg_exp_inb_hist + MM2_epg_exp_outb_hist
MM2_epg_bkg_hist         = MM2_epg_bkg_inb_hist + MM2_epg_bkg_outb_hist
MM2_epg_sig_hist         = MM2_epg_sig_inb_hist + MM2_epg_sig_outb_hist
MM2_epg_sim_hist         = MM2_epg_bkg_hist     + MM2_epg_sig_hist

MM2_epg_exp_hist_stat_err     = np.sqrt(MM2_epg_exp_inb_hist_stat_err**2 + MM2_epg_exp_outb_hist_stat_err**2)
MM2_epg_bkg_hist_stat_err     = np.sqrt(MM2_epg_bkg_inb_hist_stat_err**2 + MM2_epg_bkg_outb_hist_stat_err**2)
MM2_epg_sig_hist_stat_err     = np.sqrt(MM2_epg_sig_inb_hist_stat_err**2 + MM2_epg_sig_outb_hist_stat_err**2)

MM2_epg_bkg_hist_syst_err        = np.sqrt(MM2_epg_bkg_inb_hist_syst_err**2 + MM2_epg_bkg_outb_hist_syst_err**2)
MM2_epg_sig_hist_syst_err_up     = np.sqrt(MM2_epg_sig_inb_hist_syst_err_up**2 + MM2_epg_sig_outb_hist_syst_err_up**2)
MM2_epg_sig_hist_syst_err_down   = np.sqrt(MM2_epg_sig_inb_hist_syst_err_down**2 + MM2_epg_sig_outb_hist_syst_err_down**2)

MM2_epg_sim_hist_stat_err        = np.sqrt(MM2_epg_bkg_hist_stat_err**2 + MM2_epg_sig_hist_stat_err**2)
MM2_epg_sim_hist_syst_err_up     = np.sqrt(MM2_epg_bkg_hist_syst_err**2 + MM2_epg_sig_hist_syst_err_up**2)
MM2_epg_sim_hist_syst_err_down   = np.sqrt(MM2_epg_bkg_hist_syst_err**2 + MM2_epg_sig_hist_syst_err_down**2)
MM2_epg_sim_hist_syst_err        = (MM2_epg_sim_hist_syst_err_up + MM2_epg_sim_hist_syst_err_down)/2.

MM2_epg_exp_hist_stat_err_ratio      = divideHist(MM2_epg_exp_hist_stat_err, MM2_epg_exp_hist)
MM2_epg_sim_hist_stat_err_ratio      = divideHist(MM2_epg_sim_hist_stat_err, MM2_epg_sim_hist)
MM2_epg_sim_hist_syst_err_ratio      = divideHist(MM2_epg_sim_hist_syst_err, MM2_epg_sim_hist)
MM2_epg_sim_hist_syst_err_ratio_up   = divideHist(MM2_epg_sim_hist_syst_err_up, MM2_epg_sim_hist)
MM2_epg_sim_hist_syst_err_ratio_down = divideHist(MM2_epg_sim_hist_syst_err_down, MM2_epg_sim_hist)

MM2_epg_sim_hist_down    = MM2_epg_sim_hist - MM2_epg_sim_hist_syst_err
MM2_epg_sim_hist_up      = MM2_epg_sim_hist + MM2_epg_sim_hist_syst_err

bincenters = (bins[1:] + bins[:-1])/2.

axs[0].hist(bins[:-1], bins = bins, weights = MM2_epg_exp_hist, histtype = 'step', color = 'k', label = r"$\mathrm{Data}$")
axs[0].hist(bins[:-1], bins = bins, weights = MM2_epg_sim_hist, histtype = 'step', color = 'tab:blue', label = r'$S+B~\mathrm{(Sim.)}$')

bin_fill_between = []
MM2_epg_sim_hist_down_fill_between      = []
MM2_epg_sim_hist_up_fill_between      = []
for i in range(len(bins)):
    if (i > 0):
        bin_fill_between.append(bins[i])
        MM2_epg_sim_hist_down_fill_between.append(MM2_epg_sim_hist_down[i-1])
        MM2_epg_sim_hist_up_fill_between.append(MM2_epg_sim_hist_up[i-1])
    if (i<len(bins)-1):
        bin_fill_between.append(bins[i])
        MM2_epg_sim_hist_down_fill_between.append(MM2_epg_sim_hist_down[i])
        MM2_epg_sim_hist_up_fill_between.append(MM2_epg_sim_hist_up[i])

# axs[0].fill_between(bin_fill_between, MM2_epg_sim_hist_down_fill_between, MM2_epg_sim_hist_up_fill_between, color = 'tab:blue', alpha = 0.5)

axs[0].hist(bins[:-1], bins = bins, weights = 5*(MM2_epg_bkg_inb_hist + MM2_epg_bkg_outb_hist), histtype = 'step', color = 'tab:orange', label = "$5 \\times B~\mathrm{(Sim.)}$", zorder = -1)

axs[1].errorbar(bincenters[8:-8], divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist)[8:-8], yerr = (divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist) * np.sqrt(MM2_epg_exp_hist_stat_err_ratio**2 + MM2_epg_sim_hist_stat_err_ratio**2))[8:-8]
               , color = 'k', marker = 'o', ls = '')
# axs[1].fill_between(bincenters, divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist) * (1 -  MM2_epg_sim_hist_syst_err_ratio), divideHist(MM2_epg_exp_hist, MM2_epg_sim_hist)  * (1 +  MM2_epg_sim_hist_syst_err_ratio), color = 'r', alpha = 0.3)
axs[1].fill_between(bincenters[8:-8], (1 -  MM2_epg_sim_hist_syst_err_ratio)[8:-8], (1 + MM2_epg_sim_hist_syst_err_ratio)[8:-8], color = 'k', alpha = 0.3)

# axs[1].errorbar(bincenters, divideHist(MM2_epg_exp_inb_hist, MM2_epg_sim_inb_hist), color = 'k', marker = 'o', ls = '')
# axs[1].errorbar(bincenters, divideHist(MM2_epg_exp_outb_hist, MM2_epg_sim_outb_hist), color = 'r', marker = 'o', ls = '')

# axs[1].errorbar(bincenters, MM2_epg_exp_inb_hist/MM2_epg_sim_inb_hist#, yerr = (MM2_epg_exp_inb_hist/MM2_epg_sim_inb_hist) * np.sqrt(MM2_epg_exp_inb_hist_stat_err_ratio**2 + MM2_epg_sim_inb_hist_stat_err_ratio**2)
#                , color = 'k', marker = 'o', ls = '')
# axs[1].errorbar(bincenters, MM2_epg_exp_outb_hist/MM2_epg_sim_outb_hist#, yerr = (MM2_epg_exp_inb_hist/MM2_epg_sim_inb_hist) * np.sqrt(MM2_epg_exp_inb_hist_stat_err_ratio**2 + MM2_epg_sim_inb_hist_stat_err_ratio**2)
#                , color = 'r', marker = 'o', ls = '')


hist_pi0_inb, _ = np.histogram(df_exp_pi0_inbs.Mpi0**2, bins = bins)
hist_pi0_outb, _ = np.histogram(df_exp_pi0_outbs.Mpi0**2, bins = bins)

axs[0].axvline(.1349766**2, color = 'k', ls = '--', label = r"$m_{\pi^0}^2$")

axs[0].set_xlim([-.02, 0.02])
# axs[0].set_ylim([0, 3500])
axs[0].set_xticks([-0.02, -0.01, 0.000, 0.01, 0.02])
axs[0].set_xticklabels(['']*5)
axs[0].set_xticks(np.linspace(-0.02, 0.02, 40+1), minor = True)

axs[0].set_yticks([0, 10000, 20000, 30000], ['', '$1$', '$2$', '$3$'])
# axs[0].set_yticklabels(['', r'$5\times 10^3$', r'$10^4$', r'$3000$'])
axs[0].set_yticks(np.linspace(0, 30000, 31), minor = True)
axs[0].set_ylabel(r"$\mathrm{Events}/(5\times 10^{-4}~\mathrm{GeV}^2)$")

axs[0].annotate(xy = (-0.02, 1.03), xytext = (-0.02, 1.03), text = r'$\times 10^4$', xycoords = 'axes fraction', fontsize = 20)

axs[1].set_xlim([-.02, 0.02])
axs[1].set_ylim([0.5, 1.5])
axs[1].set_xticks([-0.02, -0.01, 0, 0.01, 0.02])
axs[1].set_xticklabels([r'${:.3f}$'.format(i) for i in [-0.02, -0.01, 0, 0.01, 0.02]])
axs[1].set_xticks(np.linspace(-0.02, 0.02, 40+1), minor = True)
axs[1].set_yticks(np.linspace(0.5, 1.5, 11), minor = True)

axs[0].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'major', length = 10)
axs[0].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'minor', length = 5)
axs[1].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'major', length = 10)
axs[1].tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'minor', length = 5)
axs[1].set_ylabel("$\mathrm{Data/Sim.}$", labelpad = 17)

axs[0].legend(title = '', loc = 'upper right', bbox_to_anchor = (0.95, 0.95),  prop={'size': 15})

axs[1].axhline(1, ls = '--', color = 'k')
axs[1].axvline(.1349766**2, color = 'k', ls = '--', label = r"$m_{\pi^0}^2$")

axs[1].set_xlabel(r"$\mathrm{M}^2_{X}~(\mathrm{GeV}^2)$")

# axs[1].set_ylim([0.9, 1.1])

# hist_exp_inb, bins = np.histogram(df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector>9].MM2_epg, weights = df_exp_epg_inbs.loc[df_exp_epg_inbs.Psector>9].signal, bins = np.linspace(-0.005, 0.005, 101))
# hist_sim_inb, _    = np.histogram(df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector>9].MM2_epg, weights = df_sim_dvcs_inbs.loc[df_sim_dvcs_inbs.Psector>9].signal, bins = bins)

# hist_exp_outb, _ = np.histogram(df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector>9].MM2_epg,  weights = df_exp_epg_outbs.loc[df_exp_epg_outbs.Psector>9].signal, bins = bins)
# hist_sim_outb, _ = np.histogram(df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector>9].MM2_epg, weights = df_sim_dvcs_outbs.loc[df_sim_dvcs_outbs.Psector>9].signal, bins = bins)

# hist_exp = hist_exp_inb + hist_exp_outb
# hist_sim = hist_sim_inb + hist_sim_outb

# plt.hist(bins[:-1], bins, weights =  hist_exp_inb, histtype = 'step')
# plt.hist(bins[:-1], bins, weights =  hist_sim_inb, histtype = 'step')
# plt.hist(bins[:-1], bins, weights =  hist_exp_outb, histtype = 'step')
# plt.hist(bins[:-1], bins, weights =  hist_sim_outb, histtype = 'step')
# plt.hist(bins[:-1], bins, weights =  (hist_exp_inb ) / (hist_sim_inb), histtype = 'step')
# plt.hist(bins[:-1], bins, weights =  (hist_exp_outb ) / (hist_sim_outb), histtype = 'step')
# axs[1].hist(bins[:-1], bins, weights =  (hist_exp_inb + hist_exp_outb) / (hist_sim_inb + hist_sim_outb), histtype = 'stepfilled')

plt.subplots_adjust(hspace=0.05)
# plt.savefig("MM2_epg_distribution.pdf", bbox_inches = 'tight')
plt.savefig("addendum_v3/M2_X_distribution.new.2.pdf", bbox_inches = 'tight')
plt.close()

df_fall2018 = pd.concat([df_exp_epg_inbs, df_exp_epg_outbs])

fig, ax = plt.subplots(1, 1, figsize = (10, 6))
h = ax.hist2d(df_fall2018.xB, df_fall2018.Q2, norm = LogNorm(vmin = 0.9, vmax = 5000), cmap = parula_map, bins = [np.linspace(0.05, 0.65, 61), np.linspace(0, 6.5, 66)], rasterized = True)

cbar = plt.colorbar(h[3])
cbar.ax.set_yticks([1, 10, 100, 1000])
cbar.ax.set_yticklabels(['$1$', '$10$', '$10^2$', '$10^3$'])
cbar.set_label(r"$\mathrm{Events}/(0.01)/(0.1~\mathrm{GeV}^2/c^2)$")

x1 = 1/2/M/8.604
x2 = 1/(5-M**2)
x3 = (10.604/8.604-1)/M*10.604* (1-np.cos(np.radians(35)))
x4 = (1-(4-M**2)/2/10.604/M)/(1+(4-M**2)/2/10.604**2/(1-np.cos(np.radians(35))))
x5 = 1/ (2*10.604*M - M/10.604/(1-np.cos(np.radians(7.74))))
x6 = (2*10.604*M/(4-M**2) -1 )  / (2*10.604*M/(4-M**2) + M/10.604/(1-np.cos(np.radians(8))))

print(x1, x2, x3, x4, x5, x6)

l1 = np.linspace(x1, x3, 101)
plt.plot(l1, l1*2*M*(10.604-2), color = 'k', linewidth = 4, solid_capstyle='round')
l2 = np.linspace(x1, x5, 101)
plt.plot(l2, 1+l2*0, color = 'k', linewidth = 4, solid_capstyle='round')

l3 = np.linspace(x3, x4, 101)
plt.plot(l3, 2*10.604*M*l3/(1+M*l3/10.604/(1-np.cos(np.radians(35)))), color = 'k', linewidth = 4, solid_capstyle='round')
l4 = np.linspace(x6, x4, 101)
plt.plot(l4, (4 - M*M)*l4/(1 - l4), color = 'k', linewidth = 4, solid_capstyle='round')

l5 = np.linspace(x5, x6, 101)
plt.plot(l5, 2*10.604*M*l5/(1+M*l5/10.604/(1-np.cos(np.radians(7.74)))), color = 'k', linewidth = 4, solid_capstyle='round', rasterized = True)


for integrated_binnum in df_summary_table_rebinned.loc[df_summary_table_rebinned.active_bin==1, "integrated_binnum"].unique():
    xBmin, xBmax, Q2min, Q2max = np.unique(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max"]].to_numpy())
    draw_box(xBmin, xBmax, Q2min, Q2max, ax)

# ax.set_xlim([0, 0.6])
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.set_xticklabels(['${}$'.format(i) for i in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]])
ax.set_xticks(np.linspace(0, 0.65, 66), minor = True)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
ax.set_yticklabels(['${}$'.format(i) for i in [0, 1, 2, 3, 4, 5, 6]])
ax.set_yticks(np.linspace(0, 6.5, 66), minor = True)

ax.set_xlabel(r"$x_B$")
ax.set_ylabel(r"$Q^2\quad(\mathrm{GeV}^2/c^2)$")

ax.tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'major', length = 10)
ax.tick_params( top = True, left = True, bottom = True, right = True, direction = 'in', pad = 10, which = 'minor', length = 5)
# ax.annotate(xy = (0.3, 1), xytext = (0.3, 1), text = "CLAS12 (This work)")

#

x1 = np.linspace(1/2/0.9382721/(5.75-0.8),0.285, 101)
ax.plot(x1, x1*2*0.9382721*(5.75-0.8), color = 'tab:red', linewidth = 4, solid_capstyle='round')#, label = 'CLAS6 Phase Space')

x2 = np.linspace(1/2/0.9382721/(5.75-0.8),0.118, 101)
ax.plot(x2, 1+x2*0, color = 'tab:red', linewidth = 4, solid_capstyle='round')

x3 = np.linspace(0.12, 0.415, 101)
ax.plot(x3, 2*5.75*0.9382721*x3/(1+0.9382721*x3/5.75/(1-0.93358)), color = 'tab:red', linewidth = 4, solid_capstyle='round')

x4 = np.linspace(0.285, 0.614, 101)
ax.plot(x4, 2*5.75*0.9382721*x4/(1+0.9382721*x4/5.75/(1-0.707107)), color = 'tab:red', linewidth = 4, solid_capstyle='round')

x5 = np.linspace(0.415, 0.611, 101)
ax.plot(x5, (4 - 0.9382721*0.9382721)*x5/(1 - x5), color = 'tab:red', linewidth = 4, solid_capstyle='round')
# ax.annotate(xy = (0.45, 2), xytext = (0.45, 2), text = "CLAS", color = 'b')


# clas6_cond_1 = (df_rga_taken_data.Q2 < df_rga_taken_data.xB*2*0.9382721*(5.75-0.8))
# clas6_cond_2 = (df_rga_taken_data.Q2 > 1)
# clas6_cond_3 = (df_rga_taken_data.Q2 > 2*5.75*0.9382721*df_rga_taken_data.xB/(1+0.9382721*df_rga_taken_data.xB/5.75/(1-0.93358)))
# clas6_cond_4 = (df_rga_taken_data.Q2 < 2*5.75*0.9382721*df_rga_taken_data.xB/(1+0.9382721*df_rga_taken_data.xB/5.75/(1-0.707107)))


plt.savefig("addendum_v3/phasespace.new.pdf", bbox_inches = 'tight')
'''

'''

#final analysis
nominal_suffix        = schema_suffices[6]
narrow_cut_suffix     = schema_suffices[7]
wide_cut_suffix       = schema_suffices[8]
more_smearing_suffix  = schema_suffices[9]
less_smearing_suffix  = schema_suffices[10]
loose_fiducial_suffix = schema_suffices[11]
tight_fiducial_suffix = schema_suffices[12]
bkg_merging_suffix    = schema_suffices[13]

df_summary_table_rebinned_sig = pd.read_pickle("impact_study_dec2024/summary_table.sig.rebinned.pkl")
df_summary_table_rebinned_bkg = pd.read_pickle("impact_study_dec2024/summary_table.bkg.rebinned.pkl")
df_summary_table_rebinned_exp = pd.read_pickle("impact_study_dec2024/summary_table.exp.rebinned.pkl")
df_summary_table_rebinned_gen = pd.read_pickle("summary_table.gen.rebinned.pkl")
df_display = pd.read_pickle("summary_table_display.pkl")

epg_inb_exp_stat_err_debug = df_summary_table_rebinned.epg_inb_exp_stat_err
epg_outb_exp_stat_err_debug = df_summary_table_rebinned.epg_outb_exp_stat_err
epg_exp_stat_err_debug = df_summary_table_rebinned.epg_exp_stat_err

# df_summary_table_rebinned = pd.DataFrame([i for i in range(2693)])
for column in ['xBbin', 'Q2bin', 'tbin', 'phi_binnum', 'phi_width', 'xBmin', 'xBmax', 'Q2min', 'Q2max', 't1min', 't1max', 'integrated_binnum', 'active_bin_inb', 'active_bin_outb', 'active_bin']+["km15_0d005_ratio", "km15_0d0005_ratio", "pureBH_0d005_ratio", "pureBH_0d0005_ratio", "rc_factor", "rc_factor_stat_err_ratio", "rc_factor_syst_err_ratio", "rad_factor_syst_err_ratio", "fbin_factor_syst_err_ratio", "pureBH_cross_section_this_point_norad", "pureBH_km15_cross_section_this_point_norad", "km15_cross_section_this_point_norad", "vgg_cross_section_this_point_norad"]:
    df_summary_table_rebinned.loc[:, column] = pd.read_pickle("df_summary_table_rebinned.backup.pkl").loc[:, column]
df_summary_table_rebinned_sig.loc[:, "n_entry_err"] = np.sqrt(df_summary_table_rebinned_sig.n_entry)
df_summary_table_rebinned_bkg.loc[:, "n_entry_err"] = np.sqrt(df_summary_table_rebinned_bkg.n_entry)
df_summary_table_rebinned_exp.loc[:, "n_entry_err"] = np.sqrt(df_summary_table_rebinned_exp.n_entry)

for integrated_binnum in df_summary_table_rebinned.integrated_binnum:
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin"])

df = pd.read_pickle("summary_table.rad.suppressed.rebinned.pkl")

directory_norad = "km15gen/dvcs_km15/norad"
df_km15gen_pureBH_km15_norad_rebinned = df.loc[df.directory == directory_norad, :].reset_index()

columns = ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point", "phi_avg_this_point", "tmin_this_point", "tcol_this_point", "P1_this_point", "P2_this_point", "maxP1P2", "minP1P2"]
df_summary_table_rebinned.loc[:, columns] = df_km15gen_pureBH_km15_norad_rebinned.loc[:, columns]

df_volume = pd.read_csv("volume_list.csv")
bin_volume_bulk = df_volume.loc[df_volume.integrated_bin <=147, :].to_numpy()[:, 1]
bin_volume_fringe = df_volume.loc[df_volume.integrated_bin > 147, :].to_numpy()[:, 1]

df_summary_table_rebinned.loc[:, "this_bin_volume"] = 0
df_summary_table_rebinned.loc[:, "bin_volume"] = 0
for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "bin_volume"] = bin_volume_bulk[integrated_binnum-1]

df_summary_table_rebinned.loc[:, "this_bin_volume"] = df_summary_table_rebinned.bin_volume * df_summary_table_rebinned.phi_width / 24.

directory = "dvcs_km15/fall2018_inb3"
n_entry_sum = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "n_entry"].to_numpy()
weight_sum  = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_sum"].to_numpy()
weight_err  = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_err"].to_numpy()
df_summary_table_rebinned.loc[:, "gen_inb_sim"]          = weight_sum/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_inb
df_summary_table_rebinned.loc[:, "gen_inb_sim_stat_err"] = weight_err/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_inb
df_summary_table_rebinned.loc[:, "gen_inb_sim_stat_err_ratio"] = divideHist(df_summary_table_rebinned.gen_inb_sim_stat_err, df_summary_table_rebinned.gen_inb_sim)

directory = "dvcs_km15/fall2018_outb3"
n_entry_sum = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "n_entry"].to_numpy()
weight_sum  = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_sum"].to_numpy()
weight_err  = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_err"].to_numpy()
df_summary_table_rebinned.loc[:, "gen_outb_sim"]          = weight_sum/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_outb
df_summary_table_rebinned.loc[:, "gen_outb_sim_stat_err"] = weight_err/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_outb
df_summary_table_rebinned.loc[:, "gen_outb_sim_stat_err_ratio"] = divideHist(df_summary_table_rebinned.gen_outb_sim_stat_err, df_summary_table_rebinned.gen_outb_sim)


directory = "pureBH/fall2018_inb3"
n_entry_sum = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "n_entry"].to_numpy()
weight_sum  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_sum"].to_numpy()
weight_err  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_err"].to_numpy()
df_summary_table_rebinned.loc[:, "gen_inb_sim_pureBH"]          = weight_sum/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_inb
df_summary_table_rebinned.loc[:, "gen_inb_sim_pureBH_stat_err"] = weight_err/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_inb
df_summary_table_rebinned.loc[:, "gen_inb_sim_pureBH_stat_err_ratio"] = divideHist(df_summary_table_rebinned.gen_inb_sim_pureBH_stat_err, df_summary_table_rebinned.gen_inb_sim_pureBH)

directory = "pureBH/fall2018_outb"
n_entry_sum = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "n_entry"].to_numpy()
weight_sum  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_sum"].to_numpy()
weight_err  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_err"].to_numpy()
df_summary_table_rebinned.loc[:, "gen_outb_sim_pureBH"]          = weight_sum/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_outb
df_summary_table_rebinned.loc[:, "gen_outb_sim_pureBH_stat_err"] = weight_err/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_outb
df_summary_table_rebinned.loc[:, "gen_outb_sim_pureBH_stat_err_ratio"] = divideHist(df_summary_table_rebinned.gen_outb_sim_pureBH_stat_err, df_summary_table_rebinned.gen_outb_sim_pureBH)


directory = "dvcs_vgg/fall2018_inb"
n_entry_sum = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "n_entry"].to_numpy()
weight_sum  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_sum"].to_numpy()
weight_err  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_err"].to_numpy()
df_summary_table_rebinned.loc[:, "gen_inb_sim_vgg"]          = weight_sum/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_inb
df_summary_table_rebinned.loc[:, "gen_inb_sim_vgg_stat_err"] = weight_err/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_inb
df_summary_table_rebinned.loc[:, "gen_inb_sim_vgg_stat_err_ratio"] = divideHist(df_summary_table_rebinned.gen_inb_sim_vgg_stat_err, df_summary_table_rebinned.gen_inb_sim_vgg)

directory = "dvcs_vgg/fall2018_outb"
n_entry_sum = df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "n_entry"].to_numpy()
weight_sum  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_sum"].to_numpy()
weight_err  = 2e-3*np.pi*df_summary_table_rebinned_gen.loc[df_summary_table_rebinned_gen.directory == directory, "weight_err"].to_numpy()
df_summary_table_rebinned.loc[:, "gen_outb_sim_vgg"]          = weight_sum/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_outb
df_summary_table_rebinned.loc[:, "gen_outb_sim_vgg_stat_err"] = weight_err/n_entry_sum * df_summary_table_rebinned.this_bin_volume * luminosity_outb
df_summary_table_rebinned.loc[:, "gen_outb_sim_vgg_stat_err_ratio"] = divideHist(df_summary_table_rebinned.gen_outb_sim_vgg_stat_err, df_summary_table_rebinned.gen_outb_sim_vgg)

df_summary_table_rebinned.loc[:, "gen_sim"]                            = df_summary_table_rebinned.loc[:, "gen_inb_sim"] + df_summary_table_rebinned.loc[:, "gen_outb_sim"]
df_summary_table_rebinned.loc[:, "gen_sim_stat_err"]                   = np.sqrt(df_summary_table_rebinned.loc[:, "gen_inb_sim_stat_err"]**2 + df_summary_table_rebinned.loc[:, "gen_outb_sim_stat_err"]**2)
df_summary_table_rebinned.loc[:, "gen_sim_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.gen_sim_stat_err, df_summary_table_rebinned.gen_sim)

df_summary_table_rebinned.loc[:, "gen_sim_pureBH"]                     = df_summary_table_rebinned.loc[:, "gen_inb_sim_pureBH"] + df_summary_table_rebinned.loc[:, "gen_outb_sim_pureBH"]
df_summary_table_rebinned.loc[:, "gen_sim_pureBH_stat_err"]            = np.sqrt(df_summary_table_rebinned.loc[:, "gen_inb_sim_pureBH_stat_err"]**2 + df_summary_table_rebinned.loc[:, "gen_outb_sim_pureBH_stat_err"]**2)
df_summary_table_rebinned.loc[:, "gen_sim_pureBH_stat_err_ratio"]      = divideHist(df_summary_table_rebinned.gen_sim_pureBH_stat_err, df_summary_table_rebinned.gen_sim_pureBH)

df_summary_table_rebinned.loc[:, "gen_sim_vgg"]                        = df_summary_table_rebinned.loc[:, "gen_inb_sim_vgg"] + df_summary_table_rebinned.loc[:, "gen_outb_sim_vgg"]
df_summary_table_rebinned.loc[:, "gen_sim_vgg_stat_err"]               = np.sqrt(df_summary_table_rebinned.loc[:, "gen_inb_sim_vgg_stat_err"]**2 + df_summary_table_rebinned.loc[:, "gen_outb_sim_vgg_stat_err"]**2)
df_summary_table_rebinned.loc[:, "gen_sim_vgg_stat_err_ratio"]         = divideHist(df_summary_table_rebinned.gen_sim_vgg_stat_err, df_summary_table_rebinned.gen_sim_vgg)

chunks = chunks_inb

# df_summary_table_rebinned.loc[:, "active_bin_inb"] = 0

# 0. experimental efficiency is fixed.
df_summary_table_rebinned.loc[:, "epg_inb_exp"]                                           = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected
df_summary_table_rebinned.loc[:, "epg_inb_exp_stat_err"]                                  = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected_err
df_summary_table_rebinned.loc[:, "epg_inb_exp_stat_err_ratio"]                            = divideHist(df_summary_table_rebinned.epg_inb_exp_stat_err, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "pi0_inb_exp"]                                           = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/pi0") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected
df_summary_table_rebinned.loc[:, "pi0_inb_exp_stat_err"]                                  = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/pi0") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected_err
df_summary_table_rebinned.loc[:, "pi0_inb_exp_stat_err_ratio"]                            = divideHist(df_summary_table_rebinned.pi0_inb_exp_stat_err, df_summary_table_rebinned.pi0_inb_exp)

for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_exp_integrated"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_exp"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_exp_integrated_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_exp_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_inb_exp_integrated_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_inb_exp_integrated_stat_err, df_summary_table_rebinned.pi0_inb_exp_integrated)

# Using alternative modes
for suffix in schema_suffices[7:13]:
    if suffix in df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs")].variation.unique():
        df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs") & (df_summary_table_rebinned_exp.variation == suffix), :].reset_index().n_entry_eff_corrected
        df_summary_table_rebinned.loc[:, "pi0_inb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/pi0") & (df_summary_table_rebinned_exp.variation == suffix), :].reset_index().n_entry_eff_corrected
    else:
        df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected
        df_summary_table_rebinned.loc[:, "pi0_inb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/pi0") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_exp_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_exp_{}".format(suffix)])

chunks = chunks_outb

# df_summary_table_rebinned.loc[:, "active_bin_outb"] = 0

# 0. experimental efficiency is fixed.
df_summary_table_rebinned.loc[:, "epg_outb_exp"]                                           = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected
df_summary_table_rebinned.loc[:, "epg_outb_exp_stat_err"]                                  = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected_err
df_summary_table_rebinned.loc[:, "epg_outb_exp_stat_err_ratio"]                            = divideHist(df_summary_table_rebinned.epg_outb_exp_stat_err, df_summary_table_rebinned.epg_outb_exp)
df_summary_table_rebinned.loc[:, "pi0_outb_exp"]                                           = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/pi0") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected
df_summary_table_rebinned.loc[:, "pi0_outb_exp_stat_err"]                                  = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/pi0") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected_err
df_summary_table_rebinned.loc[:, "pi0_outb_exp_stat_err_ratio"]                            = divideHist(df_summary_table_rebinned.pi0_outb_exp_stat_err, df_summary_table_rebinned.pi0_outb_exp)

for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_exp_integrated"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_exp"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_exp_integrated_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_exp_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_outb_exp_integrated_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_outb_exp_integrated_stat_err, df_summary_table_rebinned.pi0_outb_exp_integrated)

# Using alternative modes
for suffix in schema_suffices[7:13]:
    if suffix in df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs")].variation.unique():
        df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs") & (df_summary_table_rebinned_exp.variation == suffix), :].reset_index().n_entry_eff_corrected
        df_summary_table_rebinned.loc[:, "pi0_outb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/pi0") & (df_summary_table_rebinned_exp.variation == suffix), :].reset_index().n_entry_eff_corrected
    else:
        df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected
        df_summary_table_rebinned.loc[:, "pi0_outb_exp_{}".format(suffix)]       = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/pi0") & (df_summary_table_rebinned_exp.variation == nominal_suffix), :].reset_index().n_entry_eff_corrected

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_exp_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_exp_{}".format(suffix)])

df_summary_table_rebinned.loc[:, "epg_exp"] = df_summary_table_rebinned.loc[:, "epg_inb_exp"] + df_summary_table_rebinned.loc[:, "epg_outb_exp"]
df_summary_table_rebinned.loc[:, "epg_exp_stat_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "epg_inb_exp_stat_err"]**2 + df_summary_table_rebinned.loc[:, "epg_outb_exp_stat_err"]**2)

df_summary_table_kinematics = pd.DataFrame(df_summary_table_rebinned.loc[(y(df_summary_table_rebinned.xB_avg_this_point, df_summary_table_rebinned.Q2_avg_this_point, 0, 0)>0.7), :].integrated_binnum.unique())
df_summary_table_kinematics = df_summary_table_kinematics.rename(columns = {df_summary_table_kinematics.columns[0]: "integrated_binnum"})


df_summary_table_rebinned.loc[:, "active_bin_inb_bkg_merging"] = 0
# Using nominal efficiency map
# background merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"]                            = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3_{}nA".format(sim_current_inb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err"]                   = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3_{}nA".format(sim_current_inb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_inb_sim_bkg_merging)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim"]                                        = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_stat_err"]                               = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_stat_err_ratio"]                         = divideHist(df_summary_table_rebinned.dvcs_inb_sim_stat_err, df_summary_table_rebinned.dvcs_inb_sim)

df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_bkg_merging"]                        = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"], df_summary_table_rebinned.loc[:, "dvcs_inb_sim"])
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_bkg_merging_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_inb_sim_stat_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_bkg_merging_stat_err"]               = df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging * df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging"]                     = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_max"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_min"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging - df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_stat_err"]            = 0.5*divideHist(np.abs(df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_max - df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_min), df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_stat_err_ratio"]      = divideHist(df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err, df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging)

df_summary_table_rebinned.loc[:, "bkg_inb_sim"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
df_summary_table_rebinned.loc[:, "bkg_inb_sim_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_inb_sim_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.bkg_inb_sim_stat_err, df_summary_table_rebinned.bkg_inb_sim)
        
df_summary_table_rebinned.loc[:, "pi0_inb_sim"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
df_summary_table_rebinned.loc[:, "pi0_inb_sim_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0))
df_summary_table_rebinned.loc[:, "pi0_inb_sim_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.pi0_inb_sim_stat_err, df_summary_table_rebinned.pi0_inb_sim)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"]                            = df_summary_table_rebinned.dvcs_inb_sim * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err"]                   = df_summary_table_rebinned.dvcs_inb_sim_bkg_merging * df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging"]                             = df_summary_table_rebinned.pi0_inb_sim * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_inb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_stat_err"]                    = df_summary_table_rebinned.pi0_inb_sim_bkg_merging * df_summary_table_rebinned.pi0_inb_sim_bkg_merging_stat_err_ratio
                                
df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging"]                             = df_summary_table_rebinned.bkg_inb_sim * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.bkg_inb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_inb_sim_bkg_merging * df_summary_table_rebinned.bkg_inb_sim_bkg_merging_stat_err_ratio
                                 
for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_integrated_bkg_merging"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_integrated_bkg_merging_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_inb_sim_integrated_bkg_merging_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging_stat_err, df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging)

df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging"]                      = divideHist(df_summary_table_rebinned.pi0_inb_exp_integrated, df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_stat_err_ratio"]       = np.sqrt(df_summary_table_rebinned.pi0_inb_exp_integrated_stat_err_ratio**2 + df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_stat_err"]             = df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging"]                             = df_summary_table_rebinned.bkg_inb_sim_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging"]                             = np.minimum(df_summary_table_rebinned.bkg_inb_exp_bkg_merging, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.bkg_inb_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_inb_exp_bkg_merging * df_summary_table_rebinned.bkg_inb_exp_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging"]                            = df_summary_table_rebinned.epg_inb_exp - df_summary_table_rebinned.bkg_inb_exp_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_stat_err"]                   = np.sqrt(df_summary_table_rebinned.epg_inb_exp_stat_err **2 + df_summary_table_rebinned.bkg_inb_exp_bkg_merging_stat_err**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_inb_exp_bkg_merging)

df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging"]                       = divideHist(df_summary_table_rebinned.bkg_inb_exp_bkg_merging, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging_stat_err_ratio"]        = np.sqrt(df_summary_table_rebinned.bkg_inb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.epg_inb_exp_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging_stat_err"]              = df_summary_table_rebinned.contamination_inb_bkg_merging * df_summary_table_rebinned.contamination_inb_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging"]                                  = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging, df_summary_table_rebinned.gen_inb_sim)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_stat_err_ratio"]                   = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.gen_inb_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_stat_err"]                         = df_summary_table_rebinned.acceptance_inb_sim_bkg_merging * df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging"]                              = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging"]                            = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity_inb)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_stat_err"]                   = df_summary_table_rebinned.xsec_inb_exp_bkg_merging * df_summary_table_rebinned.xsec_inb_exp_bkg_merging_stat_err_ratio

# Using alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "bkg_inb_sim_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
    df_summary_table_rebinned.loc[:, "pi0_inb_sim_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
    df_summary_table_rebinned.loc[:, "dvcs_inb_sim_{}".format(suffix)]      = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == suffix), :].reset_index().n_entry_corrected

    df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}".format(suffix)]      = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
    df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "bkg_inb_sim_{}".format(suffix)]   * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
    df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "pi0_inb_sim_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_{}".format(suffix)]    = divideHist(df_summary_table_rebinned.loc[:, "pi0_inb_exp_integrated_{}".format(suffix)], df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_integrated_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)]           = df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_{}".format(suffix)] * df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_{}".format(suffix)] # use standard π0 inb exp to sim
    df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)]           = np.minimum(df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)])
    df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_{}".format(suffix)]          = df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)] - df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)]

    df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging_{}".format(suffix)]     = divideHist(df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_{}".format(suffix)]                       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_inb_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_{}".format(suffix)]             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_{}".format(suffix)])


df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pureBH"]                                  = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pureBH_stat_err"]                         = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pureBH_stat_err_ratio"]                   = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pureBH_stat_err, df_summary_table_rebinned.dvcs_inb_sim_pureBH)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_vgg"]                                     = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_inb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_vgg_stat_err"]                            = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_inb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_vgg_stat_err_ratio"]                      = divideHist(df_summary_table_rebinned.dvcs_inb_sim_vgg_stat_err, df_summary_table_rebinned.dvcs_inb_sim_vgg)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_pureBH"]                     = df_summary_table_rebinned.dvcs_inb_sim_pureBH * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_pureBH_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_pureBH_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_pureBH_stat_err"]            = df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_pureBH * df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_pureBH_stat_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_vgg"]                        = df_summary_table_rebinned.dvcs_inb_sim_vgg * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_vgg_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_vgg_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_vgg_stat_err"]               = df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_vgg * df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_vgg_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_pureBH"]                = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_pureBH, df_summary_table_rebinned.gen_inb_sim_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_vgg"]                   = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_vgg, df_summary_table_rebinned.gen_inb_sim_vgg)

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_pureBH"]      = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_vgg"]         = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_vgg)

# Systematic uncertainty
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_9_4sigma), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_9_4sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_11_smearing90), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_11_smearing90 - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_12_loosefid), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid2"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_12_loosefid - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"]
for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio"]               = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid"]], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_syst_err"]        = df_summary_table_rebinned.bkg_inb_exp_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err"]       = df_summary_table_rebinned.bkg_inb_exp_bkg_merging_syst_err
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging_syst_err, df_summary_table_rebinned.dvcs_inb_exp_bkg_merging)

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing1"]    = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing2"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid2"]
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_model"]        = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_inb_exp_bkg_merging", "acceptance_corrected_yield_inb_exp_bkg_merging_vgg", "acceptance_corrected_yield_inb_exp_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err_ratio"]**2)

# df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + df_summary_table_rebinned.normalization_inb_syst_err_ratio**2+ 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + 0.0453**2)

df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_inb_exp_bkg_merging * df_summary_table_rebinned.xsec_inb_exp_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_inb_exp_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_inb_exp_bkg_merging_stat_err_ratio < 0.5), "active_bin_inb_bkg_merging"] = 1
# df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_inb_exp_bkg_merging_stat_err_ratio < 0.5, "active_bin_inb_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_inb_exp_bkg_merging == 0, "active_bin_inb_bkg_merging"] = 0

df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_inb_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_inb_bkg_merging"] = 0
df_summary_table_rebinned.loc[:, "active_bin_inb_bkg_merging"] = 0
# Using nominal efficiency map
# background merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"]                            = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3_{}nA".format(sim_current_inb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err"]                   = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3_{}nA".format(sim_current_inb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_inb_sim_bkg_merging)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim"]                                        = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_stat_err"]                               = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_stat_err_ratio"]                         = divideHist(df_summary_table_rebinned.dvcs_inb_sim_stat_err, df_summary_table_rebinned.dvcs_inb_sim)

df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_bkg_merging"]                        = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"], df_summary_table_rebinned.loc[:, "dvcs_inb_sim"])
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_bkg_merging_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_inb_sim_stat_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_bkg_merging_stat_err"]               = df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging * df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging"]                     = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_max"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_min"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging - df_summary_table_rebinned.bkg_to_nobkg_inb_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_stat_err"]            = 0.5*divideHist(np.abs(df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_max - df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_min), df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_bkg_merging_stat_err_ratio"]      = divideHist(df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err, df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging)

df_summary_table_rebinned.loc[:, "bkg_inb_sim"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
df_summary_table_rebinned.loc[:, "bkg_inb_sim_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_inb_sim_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.bkg_inb_sim_stat_err, df_summary_table_rebinned.bkg_inb_sim)
        
df_summary_table_rebinned.loc[:, "pi0_inb_sim"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
df_summary_table_rebinned.loc[:, "pi0_inb_sim_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0))
df_summary_table_rebinned.loc[:, "pi0_inb_sim_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.pi0_inb_sim_stat_err, df_summary_table_rebinned.pi0_inb_sim)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"]                            = df_summary_table_rebinned.dvcs_inb_sim * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err"]                   = df_summary_table_rebinned.dvcs_inb_sim_bkg_merging * df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging"]                             = df_summary_table_rebinned.pi0_inb_sim * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_inb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_stat_err"]                    = df_summary_table_rebinned.pi0_inb_sim_bkg_merging * df_summary_table_rebinned.pi0_inb_sim_bkg_merging_stat_err_ratio
                                
df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging"]                             = df_summary_table_rebinned.bkg_inb_sim * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.bkg_inb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_inb_sim_bkg_merging * df_summary_table_rebinned.bkg_inb_sim_bkg_merging_stat_err_ratio
                                 
for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_integrated_bkg_merging"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_integrated_bkg_merging_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_inb_sim_integrated_bkg_merging_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging_stat_err, df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging)

df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging"]                      = divideHist(df_summary_table_rebinned.pi0_inb_exp_integrated, df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_stat_err_ratio"]       = np.sqrt(df_summary_table_rebinned.pi0_inb_exp_integrated_stat_err_ratio**2 + df_summary_table_rebinned.pi0_inb_sim_integrated_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_stat_err"]             = df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging"]                             = df_summary_table_rebinned.bkg_inb_sim_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging"]                             = np.minimum(df_summary_table_rebinned.bkg_inb_exp_bkg_merging, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.bkg_inb_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_inb_exp_bkg_merging * df_summary_table_rebinned.bkg_inb_exp_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging"]                            = df_summary_table_rebinned.epg_inb_exp - df_summary_table_rebinned.bkg_inb_exp_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_stat_err"]                   = np.sqrt(df_summary_table_rebinned.epg_inb_exp_stat_err **2 + df_summary_table_rebinned.bkg_inb_exp_bkg_merging_stat_err**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_inb_exp_bkg_merging)

df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging"]                       = divideHist(df_summary_table_rebinned.bkg_inb_exp_bkg_merging, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging_stat_err_ratio"]        = np.sqrt(df_summary_table_rebinned.bkg_inb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.epg_inb_exp_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging_stat_err"]              = df_summary_table_rebinned.contamination_inb_bkg_merging * df_summary_table_rebinned.contamination_inb_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging"]                                  = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging, df_summary_table_rebinned.gen_inb_sim)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_stat_err_ratio"]                   = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.gen_inb_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_stat_err"]                         = df_summary_table_rebinned.acceptance_inb_sim_bkg_merging * df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging"]                              = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging"]                            = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity_inb)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_stat_err"]                   = df_summary_table_rebinned.xsec_inb_exp_bkg_merging * df_summary_table_rebinned.xsec_inb_exp_bkg_merging_stat_err_ratio

# Using alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "bkg_inb_sim_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
    df_summary_table_rebinned.loc[:, "pi0_inb_sim_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry for i in chunks_inb], axis = 0)
    df_summary_table_rebinned.loc[:, "dvcs_inb_sim_{}".format(suffix)]      = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == suffix), :].reset_index().n_entry_corrected

    df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}".format(suffix)]      = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
    df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "bkg_inb_sim_{}".format(suffix)]   * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
    df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "pi0_inb_sim_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_bkg_merging_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_{}".format(suffix)]    = divideHist(df_summary_table_rebinned.loc[:, "pi0_inb_exp_integrated_{}".format(suffix)], df_summary_table_rebinned.loc[:, "pi0_inb_sim_bkg_merging_integrated_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)]           = df_summary_table_rebinned.loc[:, "bkg_inb_sim_bkg_merging_{}".format(suffix)] * df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_{}".format(suffix)] # use standard π0 inb exp to sim
    df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)]           = np.minimum(df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)])
    df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_{}".format(suffix)]          = df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)] - df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)]

    df_summary_table_rebinned.loc[:, "contamination_inb_bkg_merging_{}".format(suffix)]     = divideHist(df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_{}".format(suffix)]                       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_inb_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_{}".format(suffix)]             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_{}".format(suffix)])


df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pureBH"]                                  = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pureBH_stat_err"]                         = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pureBH_stat_err_ratio"]                   = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pureBH_stat_err, df_summary_table_rebinned.dvcs_inb_sim_pureBH)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_vgg"]                                     = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_inb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_vgg_stat_err"]                            = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_inb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_vgg_stat_err_ratio"]                      = divideHist(df_summary_table_rebinned.dvcs_inb_sim_vgg_stat_err, df_summary_table_rebinned.dvcs_inb_sim_vgg)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_pureBH"]                     = df_summary_table_rebinned.dvcs_inb_sim_pureBH * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_pureBH_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_pureBH_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_pureBH_stat_err"]            = df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_pureBH * df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_pureBH_stat_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_vgg"]                        = df_summary_table_rebinned.dvcs_inb_sim_vgg * df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_vgg_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_vgg_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_vgg_stat_err"]               = df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_vgg * df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_vgg_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_pureBH"]                = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_pureBH, df_summary_table_rebinned.gen_inb_sim_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_vgg"]                   = divideHist(df_summary_table_rebinned.dvcs_inb_sim_bkg_merging_vgg, df_summary_table_rebinned.gen_inb_sim_vgg)

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_pureBH"]      = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_vgg"]         = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_vgg)

# Systematic uncertainty
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_9_4sigma), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_9_4sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_11_smearing90), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_11_smearing90 - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_12_loosefid), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid2"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_12_loosefid - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"]
for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio"]               = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid"]], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_syst_err"]        = df_summary_table_rebinned.bkg_inb_exp_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_bkg_merging_syst_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err"]       = df_summary_table_rebinned.bkg_inb_exp_bkg_merging_syst_err
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.dvcs_inb_exp_bkg_merging_syst_err, df_summary_table_rebinned.dvcs_inb_exp_bkg_merging)

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing1"]    = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing2"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_inb_sim_bkg_merging_syst_err_ratio_fid2"]
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_model"]        = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_inb_exp_bkg_merging", "acceptance_corrected_yield_inb_exp_bkg_merging_vgg", "acceptance_corrected_yield_inb_exp_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_inb_sim_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err_ratio"]**2)

# df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + df_summary_table_rebinned.normalization_inb_syst_err_ratio**2+ 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + 0.0453**2)

df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_inb_exp_bkg_merging * df_summary_table_rebinned.xsec_inb_exp_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_inb_exp_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_inb_exp_bkg_merging_stat_err_ratio < 0.5), "active_bin_inb_bkg_merging"] = 1
# df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_inb_exp_bkg_merging_stat_err_ratio < 0.5, "active_bin_inb_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_inb_exp_bkg_merging == 0, "active_bin_inb_bkg_merging"] = 0

df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_inb_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_inb_bkg_merging"] = 0


df_summary_table_rebinned.loc[:, "active_bin_outb_bkg_merging"] = 0
# Using nominal efficiency map
# background merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging"]                            = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3_{}nA".format(sim_current_outb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_stat_err"]                   = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3_{}nA".format(sim_current_outb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_outb_sim_bkg_merging)

df_summary_table_rebinned.loc[:, "dvcs_outb_sim"]                                        = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_stat_err"]                               = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_stat_err_ratio"]                         = divideHist(df_summary_table_rebinned.dvcs_outb_sim_stat_err, df_summary_table_rebinned.dvcs_outb_sim)

df_summary_table_rebinned.loc[:, "bkg_to_nobkg_outb_bkg_merging"]                        = divideHist(df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging"], df_summary_table_rebinned.loc[:, "dvcs_outb_sim"])
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_outb_bkg_merging_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_stat_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_outb_bkg_merging_stat_err"]               = df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging * df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_bkg_merging"]                     = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_outb/sim_current_outb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_bkg_merging_max"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_outb/sim_current_outb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging + df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_bkg_merging_min"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_outb/sim_current_outb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging - df_summary_table_rebinned.bkg_to_nobkg_outb_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_bkg_merging_stat_err"]            = 0.5*divideHist(np.abs(df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_max - df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_min), df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_bkg_merging_stat_err_ratio"]      = divideHist(df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_stat_err, df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging)

df_summary_table_rebinned.loc[:, "bkg_outb_sim"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_outb], axis = 0)
df_summary_table_rebinned.loc[:, "bkg_outb_sim_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_outb], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_outb_sim_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.bkg_outb_sim_stat_err, df_summary_table_rebinned.bkg_outb_sim)
        
df_summary_table_rebinned.loc[:, "pi0_outb_sim"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_outb], axis = 0)
df_summary_table_rebinned.loc[:, "pi0_outb_sim_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry for i in chunks_outb], axis = 0))
df_summary_table_rebinned.loc[:, "pi0_outb_sim_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.pi0_outb_sim_stat_err, df_summary_table_rebinned.pi0_outb_sim)

df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging"]                            = df_summary_table_rebinned.dvcs_outb_sim * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_stat_err"]                   = df_summary_table_rebinned.dvcs_outb_sim_bkg_merging * df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "pi0_outb_sim_bkg_merging"]                             = df_summary_table_rebinned.pi0_outb_sim * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
df_summary_table_rebinned.loc[:, "pi0_outb_sim_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_outb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_outb_sim_bkg_merging_stat_err"]                    = df_summary_table_rebinned.pi0_outb_sim_bkg_merging * df_summary_table_rebinned.pi0_outb_sim_bkg_merging_stat_err_ratio
                                
df_summary_table_rebinned.loc[:, "bkg_outb_sim_bkg_merging"]                             = df_summary_table_rebinned.bkg_outb_sim * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_outb_sim_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.bkg_outb_sim_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_outb_sim_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_outb_sim_bkg_merging * df_summary_table_rebinned.bkg_outb_sim_bkg_merging_stat_err_ratio
                                 
for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_integrated_bkg_merging"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_integrated_bkg_merging_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_bkg_merging_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_outb_sim_integrated_bkg_merging_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_outb_sim_integrated_bkg_merging_stat_err, df_summary_table_rebinned.pi0_outb_sim_integrated_bkg_merging)

df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging"]                      = divideHist(df_summary_table_rebinned.pi0_outb_exp_integrated, df_summary_table_rebinned.pi0_outb_sim_integrated_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_stat_err_ratio"]       = np.sqrt(df_summary_table_rebinned.pi0_outb_exp_integrated_stat_err_ratio**2 + df_summary_table_rebinned.pi0_outb_sim_integrated_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_stat_err"]             = df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging * df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging"]                             = df_summary_table_rebinned.bkg_outb_sim_bkg_merging * df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging"]                             = np.minimum(df_summary_table_rebinned.bkg_outb_exp_bkg_merging, df_summary_table_rebinned.epg_outb_exp)
df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.bkg_outb_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_outb_exp_bkg_merging * df_summary_table_rebinned.bkg_outb_exp_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging"]                            = df_summary_table_rebinned.epg_outb_exp - df_summary_table_rebinned.bkg_outb_exp_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_stat_err"]                   = np.sqrt(df_summary_table_rebinned.epg_outb_exp_stat_err **2 + df_summary_table_rebinned.bkg_outb_exp_bkg_merging_stat_err**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_outb_exp_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_outb_exp_bkg_merging)

df_summary_table_rebinned.loc[:, "contamination_outb_bkg_merging"]                       = divideHist(df_summary_table_rebinned.bkg_outb_exp_bkg_merging, df_summary_table_rebinned.epg_outb_exp)
df_summary_table_rebinned.loc[:, "contamination_outb_bkg_merging_stat_err_ratio"]        = np.sqrt(df_summary_table_rebinned.bkg_outb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.epg_outb_exp_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "contamination_outb_bkg_merging_stat_err"]              = df_summary_table_rebinned.contamination_outb_bkg_merging * df_summary_table_rebinned.contamination_outb_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging"]                                  = divideHist(df_summary_table_rebinned.dvcs_outb_sim_bkg_merging, df_summary_table_rebinned.gen_outb_sim)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_stat_err_ratio"]                   = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.gen_outb_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_stat_err"]                         = df_summary_table_rebinned.acceptance_outb_sim_bkg_merging * df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging"]                              = divideHist(df_summary_table_rebinned.dvcs_outb_exp_bkg_merging, df_summary_table_rebinned.acceptance_outb_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_outb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging"]                            = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity_outb)
df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging_stat_err"]                   = df_summary_table_rebinned.xsec_outb_exp_bkg_merging * df_summary_table_rebinned.xsec_outb_exp_bkg_merging_stat_err_ratio

# Using alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "bkg_outb_sim_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry for i in chunks_outb], axis = 0)
    df_summary_table_rebinned.loc[:, "pi0_outb_sim_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry for i in chunks_outb], axis = 0)
    df_summary_table_rebinned.loc[:, "dvcs_outb_sim_{}".format(suffix)]      = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3") & (df_summary_table_rebinned_sig.variation == suffix), :].reset_index().n_entry_corrected

    df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_{}".format(suffix)]      = df_summary_table_rebinned.loc[:, "dvcs_outb_sim_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
    df_summary_table_rebinned.loc[:, "bkg_outb_sim_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "bkg_outb_sim_{}".format(suffix)]   * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
    df_summary_table_rebinned.loc[:, "pi0_outb_sim_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "pi0_outb_sim_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_bkg_merging_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_bkg_merging_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_{}".format(suffix)]    = divideHist(df_summary_table_rebinned.loc[:, "pi0_outb_exp_integrated_{}".format(suffix)], df_summary_table_rebinned.loc[:, "pi0_outb_sim_bkg_merging_integrated_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_{}".format(suffix)]           = df_summary_table_rebinned.loc[:, "bkg_outb_sim_bkg_merging_{}".format(suffix)] * df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_{}".format(suffix)] # use standard π0 outb exp to sim
    df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_{}".format(suffix)]           = np.minimum(df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)])
    df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_{}".format(suffix)]          = df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)] - df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_{}".format(suffix)]

    df_summary_table_rebinned.loc[:, "contamination_outb_bkg_merging_{}".format(suffix)]     = divideHist(df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_{}".format(suffix)]                       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_outb_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging_{}".format(suffix)]             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_{}".format(suffix)])


df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pureBH"]                                  = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pureBH_stat_err"]                         = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pureBH_stat_err_ratio"]                   = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pureBH_stat_err, df_summary_table_rebinned.dvcs_outb_sim_pureBH)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_vgg"]                                     = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_vgg_stat_err"]                            = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_vgg_stat_err_ratio"]                      = divideHist(df_summary_table_rebinned.dvcs_outb_sim_vgg_stat_err, df_summary_table_rebinned.dvcs_outb_sim_vgg)

df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_pureBH"]                     = df_summary_table_rebinned.dvcs_outb_sim_pureBH * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_pureBH_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_pureBH_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_pureBH_stat_err"]            = df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_pureBH * df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_pureBH_stat_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_vgg"]                        = df_summary_table_rebinned.dvcs_outb_sim_vgg * df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_vgg_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_vgg_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_vgg_stat_err"]               = df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_vgg * df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_vgg_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_pureBH"]                = divideHist(df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_pureBH, df_summary_table_rebinned.gen_outb_sim_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_vgg"]                   = divideHist(df_summary_table_rebinned.dvcs_outb_sim_bkg_merging_vgg, df_summary_table_rebinned.gen_outb_sim_vgg)

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging_pureBH"]      = divideHist(df_summary_table_rebinned.dvcs_outb_exp_bkg_merging, df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging_vgg"]         = divideHist(df_summary_table_rebinned.dvcs_outb_exp_bkg_merging, df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_vgg)

# Systematic uncertainty
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_9_4sigma), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_9_4sigma - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_11_smearing90), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_11_smearing90 - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_fid1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_12_loosefid), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_fid2"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_12_loosefid - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_fid3"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_cut2"]
for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio"]               = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid"]], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_syst_err"]        = df_summary_table_rebinned.bkg_outb_exp_bkg_merging * df_summary_table_rebinned.pi0_outb_exp_to_sim_bkg_merging_syst_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_syst_err"]       = df_summary_table_rebinned.bkg_outb_exp_bkg_merging_syst_err
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.dvcs_outb_exp_bkg_merging_syst_err, df_summary_table_rebinned.dvcs_outb_exp_bkg_merging)

df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_cut1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_cut2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_cut3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing1"]    = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing2"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing3"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_fid1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_fid2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_fid3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_outb_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_outb_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_outb_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_outb_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_outb_sim_bkg_merging_syst_err_ratio_fid2"]
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_model"]        = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_outb_exp_bkg_merging", "acceptance_corrected_yield_outb_exp_bkg_merging_vgg", "acceptance_corrected_yield_outb_exp_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_outb_sim_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_syst_err_ratio"]**2)

# df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + df_summary_table_rebinned.normalization_outb_syst_err_ratio**2+ 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + 0.0453**2)

df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_outb_exp_bkg_merging * df_summary_table_rebinned.xsec_outb_exp_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_outb_exp_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_outb_exp_bkg_merging_stat_err_ratio < 0.5), "active_bin_outb_bkg_merging"] = 1
# df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_outb_exp_bkg_merging_stat_err_ratio < 0.5, "active_bin_outb_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_outb_exp_bkg_merging == 0, "active_bin_outb_bkg_merging"] = 0

df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_outb_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_outb_bkg_merging"] = 0


df_summary_table_rebinned.loc[:, "active_bin_bkg_merging"] = 0

df_summary_table_rebinned.loc[:, "bkg_exp_bkg_merging"] = df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging"] + df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging"]
df_summary_table_rebinned.loc[:, "bkg_exp_bkg_merging_stat_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_stat_err"]**2 + df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_stat_err"]**2)
df_summary_table_rebinned.loc[:, "bkg_exp_bkg_merging_syst_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "bkg_inb_exp_bkg_merging_syst_err"]**2 + df_summary_table_rebinned.loc[:, "bkg_outb_exp_bkg_merging_syst_err"]**2)

df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging"]                = df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging"] + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging"]
df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_stat_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_stat_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_stat_err"]**2)
df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_stat_err_ratio"] = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_stat_err"], df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging"])
df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_syst_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_syst_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_syst_err"]**2)
df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_syst_err"], df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging"])

#nominal
df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging"]                           = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging"] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging"]
df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_stat_err"]                  = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_stat_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_stat_err"]**2 + exp_to_BH_inb_mean_stat_err**2 + exp_to_BH_outb_mean_stat_err**2)
df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_stat_err_ratio"]            = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_stat_err"], df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging"])

df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging"]                     = divideHist(df_summary_table_rebinned.dvcs_sim_bkg_merging, df_summary_table_rebinned.gen_sim)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.gen_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_stat_err"]            = df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging"]*df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_stat_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging"]     = divideHist(df_summary_table_rebinned.dvcs_exp_bkg_merging, df_summary_table_rebinned.acceptance_sim_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_sim_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_stat_err_ratio

#alternative models
for model in ["pureBH", "vgg"]:
    df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}".format(model)]                             = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}".format(model)] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_{}".format(model)]
    df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}_stat_err".format(model)]                    = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}_stat_err".format(model)]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_{}_stat_err".format(model)]**2)
    df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}_stat_err_ratio".format(model)]              = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}_stat_err".format(model)], df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}".format(model)])

    df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}".format(model)]                       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}".format(model)], df_summary_table_rebinned.loc[:, "gen_sim_{}".format(model)])
    df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}_stat_err_ratio".format(model)]        = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}_stat_err_ratio".format(model)]**2 + df_summary_table_rebinned.loc[:, "gen_sim_{}_stat_err_ratio".format(model)])
    df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}_stat_err".format(model)]              = df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}".format(model)]* df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}_stat_err_ratio".format(model)]

    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_{}".format(model)]       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging"], df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}".format(model)])
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_stat_err_ratio_{}".format(model)]               = np.sqrt(df_summary_table_rebinned.dvcs_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_sim_bkg_merging_stat_err_ratio**2)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_stat_err_{}".format(model)]                     = df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_stat_err_ratio

#alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_{}".format(suffix)]                            = df_summary_table_rebinned.loc[:, "dvcs_inb_exp_bkg_merging_{}".format(suffix)] + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_bkg_merging_{}".format(suffix)]
    df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}".format(suffix)]                            = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_bkg_merging_{}".format(suffix)] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_bkg_merging_{}".format(suffix)]
    df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}".format(suffix)]                      = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_{}".format(suffix)]      = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_{}".format(suffix)])

df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_cut1"]              = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_cut2"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_cut3"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_smearing1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_smearing2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_smearing3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_fid1"]              = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_fid2"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_fid3"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)

df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_sim_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_sim_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_sim_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_sim_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_sim_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_sim_bkg_merging_syst_err_ratio_fid2"]

df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_model"]             = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_exp_bkg_merging", "acceptance_corrected_yield_exp_bkg_merging_vgg", "acceptance_corrected_yield_exp_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_sim_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_exp_bkg_merging_syst_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_bkg_merging_syst_err"] = df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_syst_err_ratio * df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging

df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging"] = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity)
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_stat_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_stat_err"] = df_summary_table_rebinned.xsec_exp_bkg_merging * df_summary_table_rebinned.xsec_exp_bkg_merging_stat_err_ratio
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err_ratio"]            = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_exp_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_exp_bkg_merging * df_summary_table_rebinned.xsec_exp_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_exp_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_exp_bkg_merging_stat_err_ratio < 0.5), "active_bin_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_exp_bkg_merging == 0, "active_bin_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_bkg_merging"] = 0

df_summary_table_rebinned.loc[(df_summary_table_rebinned.tbin == 0) & (df_summary_table_rebinned.phi_avg_this_point>90) & (df_summary_table_rebinned.phi_avg_this_point<270) & (df_summary_table_rebinned.xsec_exp_bkg_merging < 0.9* df_summary_table_rebinned.pureBH_cross_section_this_point_norad ), "active_bin_bkg_merging"]=0

print(df_summary_table_rebinned.active_bin_bkg_merging.sum())



df_summary_table_rebinned.loc[:, "active_bin_inb_pi0_eff_corrected_bkg_merging"] = 0
# Using nominal efficiency map
# background merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging"]                            = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3_{}nA".format(sim_current_inb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                   = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3_{}nA".format(sim_current_inb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected"]                                        = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_stat_err"]                               = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_stat_err_ratio"]                         = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_stat_err, df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected)

df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging"]                        = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging"], df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected"])
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_stat_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging_stat_err"]               = df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging"]                     = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_max"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging + df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_min"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_inb/sim_current_inb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging - df_summary_table_rebinned.bkg_to_nobkg_inb_pi0_eff_corrected_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err"]            = 0.5*divideHist(np.abs(df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_max - df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_min), df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio"]      = divideHist(df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_inb], axis = 0)
df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err**2 for i in chunks_inb], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected_stat_err, df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected)
        
df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_inb], axis = 0)
df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err**2 for i in chunks_inb], axis = 0))
df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.pi0_inb_sim_pi0_eff_corrected_stat_err, df_summary_table_rebinned.pi0_inb_sim_pi0_eff_corrected)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging"]                            = df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                   = df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_bkg_merging"]                             = df_summary_table_rebinned.pi0_inb_sim_pi0_eff_corrected * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_inb_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                    = df_summary_table_rebinned.pi0_inb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio
                                
df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_bkg_merging"]                             = df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio
                                 
for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_pi0_eff_corrected_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging"]                      = divideHist(df_summary_table_rebinned.pi0_inb_exp_integrated, df_summary_table_rebinned.pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]       = np.sqrt(df_summary_table_rebinned.pi0_inb_exp_integrated_stat_err_ratio**2 + df_summary_table_rebinned.pi0_inb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err"]             = df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging"]                             = df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging"]                             = np.minimum(df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.bkg_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging"]                            = df_summary_table_rebinned.epg_inb_exp - df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                   = np.sqrt(df_summary_table_rebinned.epg_inb_exp_stat_err **2 + df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "contamination_inb_pi0_eff_corrected_bkg_merging"]                       = divideHist(df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.epg_inb_exp)
df_summary_table_rebinned.loc[:, "contamination_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio"]        = np.sqrt(df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.epg_inb_exp_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "contamination_inb_pi0_eff_corrected_bkg_merging_stat_err"]              = df_summary_table_rebinned.contamination_inb_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.contamination_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging"]                                  = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.gen_inb_sim)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]                   = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.gen_inb_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                         = df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging"]                              = divideHist(df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "xsec_inb_exp_pi0_eff_corrected_bkg_merging"]                            = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity_inb)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                   = df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

# Using alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_inb], axis = 0)
    df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_inb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_inb], axis = 0)
    df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_{}".format(suffix)]      = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == suffix), :].reset_index().n_entry_eff_pi0_corrected

    df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]      = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
    df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_{}".format(suffix)]   * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
    df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]    = divideHist(df_summary_table_rebinned.loc[:, "pi0_inb_exp_integrated_{}".format(suffix)], df_summary_table_rebinned.loc[:, "pi0_inb_sim_pi0_eff_corrected_bkg_merging_integrated_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]           = df_summary_table_rebinned.loc[:, "bkg_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)] * df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)] # use standard π0 inb exp to sim
    df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]           = np.minimum(df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)])
    df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]          = df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)] - df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]

    df_summary_table_rebinned.loc[:, "contamination_inb_pi0_eff_corrected_bkg_merging_{}".format(suffix)]     = divideHist(df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_inb_exp_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]                             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_inb_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)])


df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_pureBH"]                                  = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_pureBH_stat_err"]                         = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_inb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_pureBH_stat_err_ratio"]                   = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_pureBH_stat_err, df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_pureBH)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_vgg"]                                     = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_inb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_vgg_stat_err"]                            = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_inb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_vgg_stat_err_ratio"]                      = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_vgg_stat_err, df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_vgg)

df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH"]                     = df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_pureBH * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_pureBH_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err"]            = df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH * df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg"]                        = df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_vgg * df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_vgg_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err"]               = df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg * df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_pureBH"]                = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH, df_summary_table_rebinned.gen_inb_sim_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_vgg"]                   = divideHist(df_summary_table_rebinned.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg, df_summary_table_rebinned.gen_inb_sim_vgg)

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_pureBH"]      = divideHist(df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_vgg"]         = divideHist(df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_vgg)

# Systematic uncertainty
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_9_4sigma), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_9_4sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_11_smearing90), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_11_smearing90 - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_12_loosefid), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_12_loosefid - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]               = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid"]], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err"]        = df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_inb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err"]       = df_summary_table_rebinned.bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err
df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err, df_summary_table_rebinned.dvcs_inb_exp_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing1"]    = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_model"]        = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging", "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_vgg", "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_inb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"]**2)

# df_summary_table_rebinned.loc[:, "xsec_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + df_summary_table_rebinned.normalization_inb_syst_err_ratio**2+ 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + 0.0453**2)

df_summary_table_rebinned.loc[:, "xsec_inb_exp_pi0_eff_corrected_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.5), "active_bin_inb_pi0_eff_corrected_bkg_merging"] = 1
# df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.5, "active_bin_inb_pi0_eff_corrected_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_inb_exp_pi0_eff_corrected_bkg_merging == 0, "active_bin_inb_pi0_eff_corrected_bkg_merging"] = 0

df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_inb_pi0_eff_corrected_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_inb_pi0_eff_corrected_bkg_merging"] = 0

df_summary_table_rebinned.loc[:, "active_bin_outb_pi0_eff_corrected_bkg_merging"] = 0
# Using nominal efficiency map
# background merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging"]                            = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3_{}nA".format(sim_current_outb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                   = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3_{}nA".format(sim_current_outb)) & (df_summary_table_rebinned_sig.variation == schema_suffices[-1]), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected"]                                        = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_stat_err"]                               = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_stat_err_ratio"]                         = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_stat_err, df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected)

df_summary_table_rebinned.loc[:, "bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging"]                        = divideHist(df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging"], df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected"])
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_stat_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging_stat_err"]               = df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging"]                     = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_outb/sim_current_outb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_max"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_outb/sim_current_outb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging + df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_min"]                 = np.maximum(np.zeros(len(df_summary_table_rebinned)), ( 1 + effective_current_outb/sim_current_outb * ( -1 + df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging - df_summary_table_rebinned.bkg_to_nobkg_outb_pi0_eff_corrected_bkg_merging_stat_err)))
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err"]            = 0.5*divideHist(np.abs(df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_max - df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_min), df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio"]      = divideHist(df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_outb], axis = 0)
df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err**2 for i in chunks_outb], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected_stat_err, df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected)
        
df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected"]                                          = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_outb], axis = 0)
df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_stat_err"]                                 = np.sqrt(np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err**2 for i in chunks_outb], axis = 0))
df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_stat_err_ratio"]                           = divideHist(df_summary_table_rebinned.pi0_outb_sim_pi0_eff_corrected_stat_err, df_summary_table_rebinned.pi0_outb_sim_pi0_eff_corrected)

df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging"]                            = df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                   = df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_bkg_merging"]                             = df_summary_table_rebinned.pi0_outb_sim_pi0_eff_corrected * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_outb_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                    = df_summary_table_rebinned.pi0_outb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio
                                
df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_bkg_merging"]                             = df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio
                                 
for integrated_binnum in range(1, 147+1):
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging"]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_pi0_eff_corrected_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err"]     = np.sqrt(np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]**2))

df_summary_table_rebinned.loc[:, "pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err_ratio"]   = divideHist(df_summary_table_rebinned.pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging"]                      = divideHist(df_summary_table_rebinned.pi0_outb_exp_integrated, df_summary_table_rebinned.pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]       = np.sqrt(df_summary_table_rebinned.pi0_outb_exp_integrated_stat_err_ratio**2 + df_summary_table_rebinned.pi0_outb_sim_integrated_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err"]             = df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging"]                             = df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging"]                             = np.minimum(df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.epg_outb_exp)
df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]              = np.sqrt(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.bkg_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                    = df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging"]                            = df_summary_table_rebinned.epg_outb_exp - df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                   = np.sqrt(df_summary_table_rebinned.epg_outb_exp_stat_err **2 + df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = divideHist(df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err, df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "contamination_outb_pi0_eff_corrected_bkg_merging"]                       = divideHist(df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.epg_outb_exp)
df_summary_table_rebinned.loc[:, "contamination_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio"]        = np.sqrt(df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.epg_outb_exp_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "contamination_outb_pi0_eff_corrected_bkg_merging_stat_err"]              = df_summary_table_rebinned.contamination_outb_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.contamination_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging"]                                  = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.gen_outb_sim)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]                   = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.gen_outb_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]                         = df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging"]                              = divideHist(df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

df_summary_table_rebinned.loc[:, "xsec_outb_exp_pi0_eff_corrected_bkg_merging"]                            = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity_outb)
df_summary_table_rebinned.loc[:, "xsec_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]             = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_outb_exp_pi0_eff_corrected_bkg_merging_stat_err"]                   = df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

# Using alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_1gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_outb], axis = 0)
    df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_{}".format(suffix)]       = np.sum([df_summary_table_rebinned_bkg.loc[(df_summary_table_rebinned_bkg.directory == "sim_rad_rec_fall2018_outb/pi0_2gamma/{}".format(i)) & (df_summary_table_rebinned_bkg.variation == suffix), :].reset_index().n_entry_eff_pi0_corrected for i in chunks_outb], axis = 0)
    df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_{}".format(suffix)]      = df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_km15/fall2018_outb3") & (df_summary_table_rebinned_sig.variation == suffix), :].reset_index().n_entry_eff_pi0_corrected

    df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]      = df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
    df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_{}".format(suffix)]   * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
    df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]       = df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_{}".format(suffix)]  * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging

    for integrated_binnum in range(1, 147+1):
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_integrated_{}".format(suffix)]              = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]    = divideHist(df_summary_table_rebinned.loc[:, "pi0_outb_exp_integrated_{}".format(suffix)], df_summary_table_rebinned.loc[:, "pi0_outb_sim_pi0_eff_corrected_bkg_merging_integrated_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]           = df_summary_table_rebinned.loc[:, "bkg_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)] * df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)] # use standard π0 outb exp to sim
    df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]           = np.minimum(df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)])
    df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]          = df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)] - df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]

    df_summary_table_rebinned.loc[:, "contamination_outb_pi0_eff_corrected_bkg_merging_{}".format(suffix)]     = divideHist(df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "epg_outb_exp_{}".format(suffix)])

    df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]                             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_outb_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]             = divideHist(df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)])


df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_pureBH"]                                  = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_pureBH_stat_err"]                         = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "pureBH/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_pureBH_stat_err_ratio"]                   = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_pureBH_stat_err, df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_pureBH)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_vgg"]                                     = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_vgg_stat_err"]                            = 2e-3*np.pi*df_summary_table_rebinned_sig.loc[(df_summary_table_rebinned_sig.directory == "dvcs_vgg/fall2018_outb") & (df_summary_table_rebinned_sig.variation == nominal_suffix), :].reset_index().n_entry_eff_pi0_corrected_err
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_vgg_stat_err_ratio"]                      = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_vgg_stat_err, df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_vgg)

df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH"]                     = df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_pureBH * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_pureBH_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err"]            = df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH * df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg"]                        = df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_vgg * df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err_ratio"]         = np.sqrt(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_vgg_stat_err_ratio**2 + df_summary_table_rebinned.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err"]               = df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg * df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err_ratio

df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_pureBH"]                = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH, df_summary_table_rebinned.gen_outb_sim_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_vgg"]                   = divideHist(df_summary_table_rebinned.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg, df_summary_table_rebinned.gen_outb_sim_vgg)

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_pureBH"]      = divideHist(df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_pureBH)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_vgg"]         = divideHist(df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_vgg)

# Systematic uncertainty
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_9_4sigma), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_9_4sigma - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_11_smearing90), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_11_smearing90 - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid1"]  = 0.5 * divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_12_loosefid), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_12_loosefid - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"]  = divideHist(np.abs(df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]               = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid"]], axis = 0))
df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err"]        = df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.pi0_outb_exp_to_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err"]       = df_summary_table_rebinned.bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err
df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err, df_summary_table_rebinned.dvcs_outb_exp_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing1"]    = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"]    = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_model"]        = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging", "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_vgg", "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_outb_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"]**2)

# df_summary_table_rebinned.loc[:, "xsec_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + df_summary_table_rebinned.normalization_outb_syst_err_ratio**2+ 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 + 0.0453**2)

df_summary_table_rebinned.loc[:, "xsec_outb_exp_pi0_eff_corrected_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.5), "active_bin_outb_pi0_eff_corrected_bkg_merging"] = 1
# df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.5, "active_bin_outb_pi0_eff_corrected_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_outb_exp_pi0_eff_corrected_bkg_merging == 0, "active_bin_outb_pi0_eff_corrected_bkg_merging"] = 0

df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_outb_pi0_eff_corrected_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_outb_pi0_eff_corrected_bkg_merging"] = 0

df_summary_table_rebinned.loc[:, "active_bin_pi0_eff_corrected_bkg_merging"] = 0

df_summary_table_rebinned.loc[:, "bkg_exp_pi0_eff_corrected_bkg_merging"]                = df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging"] + df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging"]
df_summary_table_rebinned.loc[:, "bkg_exp_pi0_eff_corrected_bkg_merging_stat_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err"]**2 + df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err"]**2)
df_summary_table_rebinned.loc[:, "bkg_exp_pi0_eff_corrected_bkg_merging_syst_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err"]**2 + df_summary_table_rebinned.loc[:, "bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err"]**2)

df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging"]                = df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging"] + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging"]
df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_stat_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err"]**2)
df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"] = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_stat_err"], df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging"])
df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err"]       = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err"]**2)
df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err"], df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging"])

#nominal
df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected"]                           = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected"] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected"]
df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_stat_err"]                  = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_stat_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_stat_err"]**2)
df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_stat_err_ratio"]            = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_stat_err"], df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected"])

df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging"]                           = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging"] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging"]
df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_stat_err"]                  = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err"]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err"]**2)
df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]            = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_stat_err"], df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging"])

df_summary_table_rebinned.loc[:, "eff_bkg_merging_pi0_eff_corrected_bkg_merging"]  = divideHist(df_summary_table_rebinned.dvcs_sim_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.dvcs_sim_pi0_eff_corrected)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err_ratio"]  = np.sqrt(df_summary_table_rebinned.dvcs_sim_pi0_eff_corrected_stat_err_ratio**2 + df_summary_table_rebinned.dvcs_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err"]  = df_summary_table_rebinned.loc[:, "eff_bkg_merging_pi0_eff_corrected_bkg_merging"] * df_summary_table_rebinned.loc[:, "eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err_ratio"] 

df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging"]                     = divideHist(df_summary_table_rebinned.dvcs_sim_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.gen_sim)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]      = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]**2 + df_summary_table_rebinned.gen_sim_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err"]            = df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging"]*df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging"]     = divideHist(df_summary_table_rebinned.dvcs_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"]               = np.sqrt(df_summary_table_rebinned.dvcs_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err"]                     = df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

#alternative models
for model in ["pureBH", "vgg"]:
    df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}".format(model)]                             = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(model)] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(model)]
    df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}_stat_err".format(model)]                    = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_{}_stat_err".format(model)]**2 + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_{}_stat_err".format(model)]**2)
    df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}_stat_err_ratio".format(model)]              = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}_stat_err".format(model)], df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}".format(model)])

    df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}".format(model)]                       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}".format(model)], df_summary_table_rebinned.loc[:, "gen_sim_{}".format(model)])
    df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}_stat_err_ratio".format(model)]        = np.sqrt(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}_stat_err_ratio".format(model)]**2 + df_summary_table_rebinned.loc[:, "gen_sim_{}_stat_err_ratio".format(model)])
    df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}_stat_err".format(model)]              = df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}".format(model)]* df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}_stat_err_ratio".format(model)]

    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_{}".format(model)]       = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging"], df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}".format(model)])
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio_{}".format(model)]               = np.sqrt(df_summary_table_rebinned.dvcs_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err_{}".format(model)]                     = df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio

#alternative modes
for suffix in schema_suffices[7:13]:
    df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]                            = df_summary_table_rebinned.loc[:, "dvcs_inb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)] + df_summary_table_rebinned.loc[:, "dvcs_outb_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]
    df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]                            = df_summary_table_rebinned.loc[:, "dvcs_inb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)] + df_summary_table_rebinned.loc[:, "dvcs_outb_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]
    df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)]                      = divideHist(df_summary_table_rebinned.loc[:, "dvcs_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.gen_sim)
    df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)]      = divideHist(df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_{}".format(suffix)], df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_{}".format(suffix)])

df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut1"]              = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_9_4sigma), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_9_4sigma - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_8_2sigma - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing1"]         = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_11_smearing90), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_11_smearing90 - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"]         = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_10_smearing110 - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid1"]              = 0.5 * divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_12_loosefid), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_13_tightfid - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"]              = divideHist(np.abs(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_12_loosefid - df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)

df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut3 > 0.5, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing3 > 0.5, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing2"]
df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3"] = df_summary_table_rebinned.loc[df_summary_table_rebinned.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid3 > 0.5, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid2"]

df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_model"]             = divideHist(np.std(df_summary_table_rebinned.loc[:, ["acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging", "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_vgg", "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_pureBH"]], axis = 1), df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging)

for syst_err_type in ["cut", "smearing", "fid"]:
    df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)] = 0.5* (df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}2".format(syst_err_type)] + df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}3".format(syst_err_type)])
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(np.sum([df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_{}".format(syst_err_type)]**2 for syst_err_type in ["cut", "smearing", "fid", "model"]], axis = 0))
df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err"]       = df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging"] * df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]

df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"] = np.sqrt(df_summary_table_rebinned.loc[:, "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"]**2 + df_summary_table_rebinned.loc[:, "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"]**2)
df_summary_table_rebinned.loc[:, "acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err"] = df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio * df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging

df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging"] = divideHist(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging, df_summary_table_rebinned.this_bin_volume * df_summary_table_rebinned.rc_factor * luminosity)
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio"] = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned.rc_factor_stat_err_ratio**2)
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_stat_err"] = df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"]            = np.sqrt(df_summary_table_rebinned.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + df_summary_table_rebinned.rc_factor_syst_err_ratio**2 )#+ 0.0453**2)
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err"] = df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio

df_summary_table_rebinned.loc[(df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio < 0.5) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.5), "active_bin_pi0_eff_corrected_bkg_merging"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging == 0, "active_bin_pi0_eff_corrected_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_pi0_eff_corrected_bkg_merging"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_pi0_eff_corrected_bkg_merging"] = 0

df_summary_table_rebinned.loc[(df_summary_table_rebinned.tbin == 0) & (df_summary_table_rebinned.phi_avg_this_point>90) & (df_summary_table_rebinned.phi_avg_this_point<270) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging < 0.9* df_summary_table_rebinned.pureBH_cross_section_this_point_norad ), "active_bin_pi0_eff_corrected_bkg_merging"]=0

print(df_summary_table_rebinned.active_bin_pi0_eff_corrected_bkg_merging.sum())

df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin"] = df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_bin_by_bin"]       = df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_global"]     = np.sqrt(0.3**2 + 0.0476**2)#0.3 is normalization 0.0476 is others
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_global"]           = df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_global

df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"]            = np.sqrt(df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin"]**2 + df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_global"]**2)
df_summary_table_rebinned.loc[:, "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err"]                  = df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging * df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio

print(df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin.mean())


df_summary_table_rebinned.loc[:, "active_bin_nominal"] = 0
# df_summary_table_rebinned.loc[(df_summary_table_rebinned.epg_exp > 10) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio < 0.4) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.4) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_up < 0.6) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_down < 0.6), "active_bin_nominal"] = 1
df_summary_table_rebinned.loc[(df_summary_table_rebinned.epg_exp > 10) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin < 0.4) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.4), "active_bin_nominal"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_exp == 0, "active_bin_nominal"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_nominal"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_nominal"] = 0
# df_summary_table_rebinned.loc[df_summary_table_rebinned.efficiency < 0.25, "active_bin_nominal"] = 0
df_summary_table_rebinned.loc[(df_summary_table_rebinned.tbin == 0) & (df_summary_table_rebinned.phi_avg_this_point>90) & (df_summary_table_rebinned.phi_avg_this_point<270) & (df_summary_table_rebinned.xsec_exp < 0.9* df_summary_table_rebinned.pureBH_cross_section_this_point_norad ), "active_bin_nominal"]=0
print(df_summary_table_rebinned.active_bin_nominal.sum())


models = ["pi0"]#, "km15", "bh", "vgg"]
for integrated_binnum in df_summary_table_rebinned.integrated_binnum:
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_nominal_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_nominal"])
    for model in models:
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb_{}_eff_corrected_bkg_merging_integrated".format(model)] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb_{}_eff_corrected_bkg_merging".format(model)])
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb_{}_eff_corrected_bkg_merging_integrated".format(model)] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb_{}_eff_corrected_bkg_merging".format(model)])
        df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_{}_eff_corrected_bkg_merging_integrated".format(model)] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_{}_eff_corrected_bkg_merging".format(model)])

    df_summary_table_rebinned.loc[:, "normalization_{}".format(model)] = divideHist(df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging"])     
    df_summary_table_rebinned.loc[:, "normalization_inb_{}".format(model)] = divideHist(df_summary_table_rebinned.loc[:, "xsec_inb_exp_{}_eff_corrected_bkg_merging".format(model)], df_summary_table_rebinned.loc[:, "xsec_inb_exp_bkg_merging"])         
    df_summary_table_rebinned.loc[:, "normalization_outb_{}".format(model)] = divideHist(df_summary_table_rebinned.loc[:, "xsec_outb_exp_{}_eff_corrected_bkg_merging".format(model)], df_summary_table_rebinned.loc[:, "xsec_outb_exp_bkg_merging"])                 


# report uncertainty

print("Background", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))#, "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"].std()))
print("Model", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_model"].mean()))
print("Cut", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut"].mean()))
print("Fid", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid"].mean()))
print("Smearing", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing"].mean()))
# print("{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))
# print("Background Merging", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "bkg_merging_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))
print("Rad Factor", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "rad_factor_syst_err_ratio"].mean()))
print("Fbin", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "fbin_factor_syst_err_ratio"].mean()))
print("Total Bin-by-bin ", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin"].mean()))
# print("Bin-by-bin", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin"].mean()))
# print("Efficiency", "{:.2f}".format(100*df_summary_table_rebinned.loc[((df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3)) & (df_summary_table_rebinned.efficiency>0), "efficiency_pi0_eff_corrected_bkg_merging_syst_err_ratio"].unique().mean()))
# print("Normalization", "{:.2f}".format(100*1/df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "normalization_pi0"].mean()))
print("Total", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.tbin == 3), "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))


print(df_summary_table_rebinned.active_bin_nominal.sum())

print(np.sqrt(0.3**2 + 0.0476**2))

# report uncertainty

print("Background", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))#, "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "dvcs_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"].std()))
print("Model", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_model"].mean()))
print("Cut", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_cut"].mean()))
print("Fid", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_fid"].mean()))
print("Smearing", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio_smearing"].mean()))
# print("{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))
# print("Background Merging", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "bkg_merging_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))
print("Rad Factor", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "rad_factor_syst_err_ratio"].mean()))
print("Fbin", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "fbin_factor_syst_err_ratio"].mean()))
print("Total Bin-by-bin ", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin"].mean()))
# print("Bin-by-bin", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin"].mean()))
# print("Efficiency", "{:.2f}".format(100*df_summary_table_rebinned.loc[((df_summary_table_rebinned.active_bin_nominal == 1)) & (df_summary_table_rebinned.efficiency>0), "efficiency_pi0_eff_corrected_bkg_merging_syst_err_ratio"].unique().mean()))
# print("Normalization", "{:.2f}".format(100*1/df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "normalization_pi0"].mean()))
print("Total", "{:.2f}".format(100*df_summary_table_rebinned.loc[(df_summary_table_rebinned.active_bin_nominal == 1), "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio"].mean()))


print(df_summary_table_rebinned.active_bin_nominal.sum())

df_summary_table_rebinned.to_pickle("addendum_v3/df_summary_table_rebinned.final_analysis.pkl")
# analysis done
'''

df_display                = pd.read_pickle("summary_table_display.pkl")
df_summary_table_rebinned = pd.read_pickle("addendum_v3/df_summary_table_rebinned.final_analysis.pkl")

# # just 5 plots
xB_panes = 8
Q2_panes = 7
models = ["pi0"]
label_scheme = ["$\mathrm{Data}$"]
color_scheme = ['k']

# xB_binnum = 3
# Q2_binnum = 3

# Q2_avg = []
# xB_avgs = []
# fig, axs = plt.subplots(5, 1, figsize = (10, 20))
# for t_binnum in range(1, 6):

#   df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
#   Q2_avg.append(df_this_bin.t_avg_this_point.unique()[0])
#   xB_avgs.append(df_this_bin.xB_avg_this_point.unique()[0])
#   integrated_binnum = df_this_bin.integrated_binnum.unique()[0]
#   axs[t_binnum-1].annotate("{}.~".format(t_binnum)+r"$\langle |t| \rangle = {:.2f}~\mathrm{{GeV}}^2$".format(df_this_bin.t_avg_this_point.unique()[0]), xy = (0.5, 0.83), xytext = (0.5, 0.83), xycoords = 'axes fraction', fontsize = 40, ha = 'center' )
#   axs[t_binnum-1].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=3, label = r"$\mathrm{KM15}$")
#   axs[t_binnum-1].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].phi_display), df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].xsec_VGG_display, color = 'tab:orange', lw=3, label = r"$\mathrm{VGG}$")
#   axs[t_binnum-1].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=3, label = r"$\mathrm{BH}$")

#   df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal"] == 1), :]
#   weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
#   weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
#   weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
#   weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err

#   for i, model in enumerate(models):
#       if i == 0:
#           dots = axs[t_binnum-1].scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
#           axs[t_binnum-1].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
#           axs[t_binnum-1].fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
#       else:
#           axs[t_binnum-1].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
#           pass
#       # print(df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)].max())
#   axs[t_binnum-1].set_yscale('log')

#   axs[t_binnum-1].set_xlim([0, 360])
#   axs[t_binnum-1].set_xticks(np.linspace(0, 360, 12+1), minor = True)
#   axs[t_binnum-1].set_xlabel("$\phi$ ($^{\circ}$)", fontsize = 60 )
#   if t_binnum - 1 == 4:
#       axs[t_binnum-1].set_xticks(np.linspace(0, 360, 4+1))
#   else:
#       axs[t_binnum-1].set_xticks(np.linspace(0, 360, 4+1), ['']*5)

  
#   axs[t_binnum-1].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 20, labelsize = 40)
#   axs[t_binnum-1].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 10)
#   axs[t_binnum-1].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
#   axs[t_binnum-1].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)

# axs[0].set_ylim([4.8e-2, 9.2e-1])
# axs[1].set_ylim([1.4e-2, 7.2e-1])
# axs[2].set_ylim([4.8e-3, 5.2e-1])
# axs[3].set_ylim([1.8e-3, 1.2e-1])
# axs[4].set_ylim([7.8e-4, 6.2e-2])

# handles, labels = axs[2].get_legend_handles_labels() 
# handles = [handles[-1], handles[0], handles[1], handles[2]]
# labels  = [labels [-1], labels [0], labels [1], labels [2]]

# # axs[0].set_title(r"$\langle x_B \rangle={:.2f},~\langle Q^2 \rangle={:.2f}~\mathrm{{GeV}}^2/c^2$".format(np.mean(xB_avgs), np.mean(Q2_avg)), fontsize= 40)
# fig.text(-0.05, 0.5, r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$", va='center', rotation = 'vertical', fontsize = 60)
 
# axs[0].legend(handles, labels, loc = 'lower left', bbox_to_anchor = (-.2, .85), title = r"$\langle x_B \rangle={:.2f},~\langle Q^2 \rangle={:.2f}~\mathrm{{GeV}}^2/c^2$".format(np.mean(xB_avgs), np.mean(Q2_avg)), title_fontsize = 40, fontsize = 40, alignment = 'center', ncol = 2, markerscale = 3, framealpha = 0, handlelength = 1, columnspacing = .5)
# # plt.tight_layout()
# plt.subplots_adjust(wspace = 0.2 , hspace = 0.0, left = .2 )
# plt.savefig("addendum_v3/xsec.pdf".format(t_binnum), bbox_inches = 'tight')
# plt.close()

df_display.loc[:, "weight_BH"] = weight_BH(df_display.xB_display, df_display.Q2_display, df_display.t_display, np.degrees(df_display.phi_display))

df_display.loc[:, "xsec_BH_display_w"]           = df_display.loc[:, "xsec_BH_display"]     *df_display.loc[:, "weight_BH"]
df_display.loc[:, "xsec_KM15_display_w"]         = df_display.loc[:, "xsec_KM15_display"]   *df_display.loc[:, "weight_BH"]
df_display.loc[:, "xsec_BH_KM15_display_w"]      = df_display.loc[:, "xsec_BH_KM15_display"]*df_display.loc[:, "weight_BH"]
df_display.loc[:, "xsec_VGG_display_w"]          = df_display.loc[:, "xsec_VGG_display"]    *df_display.loc[:, "weight_BH"]

model = 'pi0'

df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_w".format(model)]                       = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err_w".format(model)]              = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_w".format(model)]              = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_bin_by_bin_w".format(model)]   = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_bin_by_bin".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_w"]                                                      = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging"]*df_summary_table_rebinned.loc[:, "weight_BH"]
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_stat_err_w"]                                             = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_stat_err"]*df_summary_table_rebinned.loc[:, "weight_BH"]
df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err_w"]                                             = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err"]*df_summary_table_rebinned.loc[:, "weight_BH"]

n_sample    = 48
i           = 0
df_fittings = pd.DataFrame()
phi_dummy = np.linspace(0, 360, n_sample)
for integrated_binnum in np.sort(df_summary_table_rebinned.integrated_binnum.unique()):
  if integrated_binnum < 1:
    continue
  if integrated_binnum > 147:
    continue
  df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal"] == 1), :]
  df_display_this_bin = df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)]

  if len(df_summary_table_rebinned_this_bin)<4:
    continue
  xBs = df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0] * np.ones_like(phi_dummy)
  Q2s = df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0] * np.ones_like(phi_dummy)
  ts = df_summary_table_rebinned_this_bin.t_avg_this_point.unique()[0] * np.ones_like(phi_dummy)

  fig, ax = plt.subplots(1, 1, figsize = (10, 6))
  cosine_1_th_km15, cosine_1_th_km15_err, cosine_1_th_vgg, cosine_1_th_vgg_err, cosine_1_th_bh, cosine_1_th_bh_err, cosine_1_th_bh_km15, cosine_1_th_bh_km15_err, popt_1_th_km15, popt_1_th_vgg, popt_1_th_bh, popt_1_th_bh_km15 = cosine_1_fitting_one_bin_th(df_display_this_bin)
  cosine_2_th_km15, cosine_2_th_km15_err, cosine_2_th_vgg, cosine_2_th_vgg_err, cosine_2_th_bh, cosine_2_th_bh_err, cosine_2_th_bh_km15, cosine_2_th_bh_km15_err, popt_2_th_km15, popt_2_th_vgg, popt_2_th_bh, popt_2_th_bh_km15 = cosine_2_fitting_one_bin_th(df_display_this_bin)
  cosine_3_th_km15, cosine_3_th_km15_err, cosine_3_th_vgg, cosine_3_th_vgg_err, cosine_3_th_bh, cosine_3_th_bh_err, cosine_3_th_bh_km15, cosine_3_th_bh_km15_err, popt_3_th_km15, popt_3_th_vgg, popt_3_th_bh, popt_3_th_bh_km15 = cosine_3_fitting_one_bin_th(df_display_this_bin)
  plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display_w, color = 'cyan', lw=3, label = r"$\mathrm{KM15}$")
  plt.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_km15), color = 'cyan', ls = '--')
  plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].phi_display), df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].xsec_VGG_display_w, color = 'tab:orange', lw=3, label = r"$\mathrm{VGG}$")
  plt.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_vgg), color = 'tab:orange', ls = '--')
  plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display_w, color = 'tab:red', lw=3, label = r"$\mathrm{BH}$")
  plt.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_bh_km15), color = 'tab:red', ls = '--')

  plt.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_w".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err_w".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
  xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_bkgmerging_only, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_bkgmerging_only = cosine_1_fitting_one_bin_exp(df_summary_table_rebinned_this_bin, p0 = popt_1_th_km15)
  xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_bkgmerging_only, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_bkgmerging_only = cosine_2_fitting_one_bin_exp(df_summary_table_rebinned_this_bin, p0 = popt_2_th_km15)
  xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_bkgmerging_only, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_bkgmerging_only = cosine_3_fitting_one_bin_exp(df_summary_table_rebinned_this_bin, p0 = popt_3_th_km15)
  plt.plot(phi_dummy, cosine_fitting_1(np.radians(phi_dummy), *popt_1_pi0), color = 'k')
  # plt.plot(phi_dummy, cosine_fitting_2(np.radians(phi_dummy), *popt_2_pi0), color = 'k')
  chi2fit = np.sum((df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - cosine_fitting_2(np.radians(df_summary_table_rebinned_this_bin.phi_avg_this_point), *popt_2_pi0))**2/(df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w)**2)
  ndf = len(df_summary_table_rebinned_this_bin) - 3
  pvalue = (1-chi2.cdf(chi2fit, ndf))
  plt.savefig("addendum_v3/fitting/eachbin/modified_xsec_{}.pdf".format(integrated_binnum))
  plt.close()

  fig, ax = plt.subplots(1, 1, figsize = (10, 6))
  plt.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
  plt.plot(phi_dummy, cosine_fitting_1(np.radians(phi_dummy), *popt_1_pi0)/weight_BH(xBs, Q2s, ts, phi_dummy), color = 'k')
  plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=3, label = r"$\mathrm{KM15}$")
  plt.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_km15)/weight_BH(xBs, Q2s, ts, phi_dummy), color = 'tab:cyan')  
  plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].phi_display), df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].xsec_VGG_display, color = 'tab:orange', lw=3, label = r"$\mathrm{VGG}$")
  plt.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_vgg)/weight_BH(xBs, Q2s, ts, phi_dummy), color = 'tab:orange')  
  plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=3, label = r"$\mathrm{BH}$")
  plt.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_bh_km15)/weight_BH(xBs, Q2s, ts, phi_dummy), color = 'tab:brown')  

  plt.yscale('log')
  plt.savefig("addendum_v3/fitting/eachbin/xsec_{}.pdf".format(integrated_binnum))
  plt.close()

#   df_fitting = pd.DataFrame([{"integrated_binnum": integrated_binnum,
#                               "xB_avg": xB_avg,
#                               "Q2_avg": Q2_avg,
#                               "t_avg": t_avg,
#                               "cosine_1_exp_pi0": cosine_1_exp_pi0,
#                               "cosine_1_exp_pi0_stat_err": cosine_1_exp_pi0_stat_err,
#                               "cosine_1_exp_pi0_syst_err": cosine_1_exp_pi0_syst_err,
#                               "cosine_1_exp_pi0_norm_err_up": cosine_1_exp_pi0_norm_err_up,
#                               "cosine_1_exp_pi0_norm_err_down": cosine_1_exp_pi0_norm_err_down,
#                               "cosine_1_exp_bkgmerging_only": cosine_1_exp_bkgmerging_only,
#                               "popt_1_pi0": popt_1_pi0,
#                               "popt_1_pi0_min": popt_1_pi0_min,
#                               "popt_1_pi0_max": popt_1_pi0_max,
#                               "popt_1_bkgmerging_only": popt_1_bkgmerging_only,
#                               "cosine_2_exp_pi0": cosine_2_exp_pi0,
#                               "cosine_2_exp_pi0_stat_err": cosine_2_exp_pi0_stat_err,
#                               "cosine_2_exp_pi0_syst_err": cosine_2_exp_pi0_syst_err,
#                               "cosine_2_exp_pi0_norm_err_up": cosine_2_exp_pi0_norm_err_up,
#                               "cosine_2_exp_pi0_norm_err_down": cosine_2_exp_pi0_norm_err_down,
#                               "cosine_2_exp_bkgmerging_only": cosine_2_exp_bkgmerging_only,
#                               "popt_2_pi0": popt_2_pi0,
#                               "popt_2_pi0_min": popt_2_pi0_min,
#                               "popt_2_pi0_max": popt_2_pi0_max,
#                               "popt_2_bkgmerging_only": popt_2_bkgmerging_only,
#                               "cosine_3_exp_pi0": cosine_3_exp_pi0,
#                               "cosine_3_exp_pi0_stat_err": cosine_3_exp_pi0_stat_err,
#                               "cosine_3_exp_pi0_syst_err": cosine_3_exp_pi0_syst_err,
#                               "cosine_3_exp_pi0_norm_err_up": cosine_3_exp_pi0_norm_err_up,
#                               "cosine_3_exp_pi0_norm_err_down": cosine_3_exp_pi0_norm_err_down,
#                               "cosine_3_exp_bkgmerging_only": cosine_3_exp_bkgmerging_only,
#                               "popt_3_pi0": popt_3_pi0,
#                               "popt_3_pi0_min": popt_3_pi0_min,
#                               "popt_3_pi0_max": popt_3_pi0_max,
#                               "popt_3_bkgmerging_only": popt_3_bkgmerging_only,
#                               "pvalue": pvalue}], index = [integrated_binnum])
#   df_fittings = pd.concat([df_fittings, df_fitting])

# df_fittings.to_pickle("addendum_v3/fitting/eachbin/df_fittings.pkl")

# n_display = 2
# n_sample  = 5
# phi_dummy = np.linspace(0, 360, n_sample)


# df_t_dependence                          = pd.DataFrame()
# df_display_this_bin = df_display.loc[(df_display.integratedbin_display == args.integrated_binnum) & (df_display.phi_display == 0), :]

# xBbin = int(df_display_this_bin.xBbin_display)
# Q2bin = int(df_display_this_bin.Q2bin_display)
# tbin = int(df_display_this_bin.tbin_display)

# assert len(df_display_this_bin) == 1

# xBs  = float(df_display_this_bin.xB_display) * np.ones_like(phi_dummy)
# Q2s  = float(df_display_this_bin.Q2_display)  * np.ones_like(phi_dummy)
# ts   = np.linspace(-tmin(float(df_display_this_bin.xB_display.mean()), float(df_display_this_bin.Q2_display.mean()), 0, 0), 1, n_display)
# for index, t_running in enumerate(ts):
#   print("t dep", xBbin, Q2bin, index)
#   yd = y(float(df_display_this_bin.xB_display), float(df_display_this_bin.Q2_display), 0, 0)
#   if t_running < -tmin(float(df_display_this_bin.xB_display), float(df_display_this_bin.Q2_display), 0, 0):
#     continue
#   if (yd<=0.1) or (yd>0.9):
#     continue
#   Wd = W(float(df_display_this_bin.xB_display), float(df_display_this_bin.Q2_display), 0, 0)
#   if (Wd < 2):
#     continue
#   df_t_dependence_this_bin =  pd.DataFrame(data = {"integrated_binnum": args.integrated_binnum, "xBbin": xBbin * np.ones_like(phi_dummy), "Q2bin": Q2bin * np.ones_like(phi_dummy), 
#     "tbin": tbin * np.ones_like(phi_dummy), "tindex": index * np.ones_like(phi_dummy), "xB_display": xBs, "Q2_display": Q2s, "t_display": t_running * np.ones_like(phi_dummy), "phi_display": phi_dummy})
#   df_t_dependence_this_bin.loc[:, "xsec_KM15_display"]       = printKMarray (df_t_dependence_this_bin.xB_display, df_t_dependence_this_bin.Q2_display, df_t_dependence_this_bin.t_display, np.radians(df_t_dependence_this_bin.phi_display), mode = 5)
#   df_t_dependence_this_bin.loc[:, "xsec_BH_KM15_display"]    = printKMarray (df_t_dependence_this_bin.xB_display, df_t_dependence_this_bin.Q2_display, df_t_dependence_this_bin.t_display, np.radians(df_t_dependence_this_bin.phi_display), mode = 1)
#   df_t_dependence_this_bin.loc[:, "xsec_BH_display"]         = printBHarray (df_t_dependence_this_bin.xB_display, df_t_dependence_this_bin.Q2_display, df_t_dependence_this_bin.t_display, np.radians(df_t_dependence_this_bin.phi_display), local=True)
#   df_t_dependence_this_bin.loc[:, "xsec_VGG_display"]        = printVGGarray(df_t_dependence_this_bin.xB_display, df_t_dependence_this_bin.Q2_display, df_t_dependence_this_bin.t_display, np.radians(df_t_dependence_this_bin.phi_display), local=True)
#   df_t_dependence_this_bin.loc[:, "weight_BH"]               = weight_BH(df_t_dependence_this_bin.xB_display, df_t_dependence_this_bin.Q2_display, df_t_dependence_this_bin.t_display, df_t_dependence_this_bin.phi_display)
#   df_t_dependence_this_bin.loc[:, "xsec_BH_display_w"]       = df_t_dependence_this_bin.loc[:, "xsec_BH_display"]     *df_t_dependence_this_bin.loc[:, "weight_BH"]
#   df_t_dependence_this_bin.loc[:, "xsec_KM15_display_w"]     = df_t_dependence_this_bin.loc[:, "xsec_KM15_display"]   *df_t_dependence_this_bin.loc[:, "weight_BH"]
#   df_t_dependence_this_bin.loc[:, "xsec_BH_KM15_display_w"]  = df_t_dependence_this_bin.loc[:, "xsec_BH_KM15_display"]*df_t_dependence_this_bin.loc[:, "weight_BH"]
#   df_t_dependence_this_bin.loc[:, "xsec_VGG_display_w"]      = df_t_dependence_this_bin.loc[:, "xsec_VGG_display"]    *df_t_dependence_this_bin.loc[:, "weight_BH"]
#   df_t_dependence = pd.concat([df_t_dependence, df_t_dependence_this_bin])
# df_t_dependence.reset_index(inplace = True)
# df_t_dependence.rename(columns = {"index": "phiindex"}, inplace = True)
# df_t_dependence.to_pickle("addendum_v3/fitting/df_t_dependence_{}.pkl".format(args.integrated_binnum))

# df_Q2_dependence                          = pd.DataFrame()

# xBs   = float(df_display_this_bin.xB_display) * np.ones_like(phi_dummy)
# ts    = float(df_display_this_bin.t_display)  * np.ones_like(phi_dummy)
# Q2s   = np.linspace(1, 6, n_display)
# Q2index = []
# for index, Q2_running in enumerate(Q2s):
#   xBbin = int(df_display_this_bin.xBbin_display)
#   Q2bin = int(df_display_this_bin.Q2bin_display)
#   tbin  = int(df_display_this_bin.tbin_display)
#   print("Q2 dep", xBbin, tbin, index)
#   if float(df_display_this_bin.t_display) < -tmin(float(df_display_this_bin.xB_display), Q2_running, 0, 0):
#     continue
#   yd = y(float(df_display_this_bin.xB_display), Q2_running, 0, 0)
#   if (yd<=0.1) or (yd>0.9):
#     continue
#   Wd = W(float(df_display_this_bin.xB_display), Q2_running, 0, 0)
#   if (Wd < 2):
#     continue
#   df_Q2_dependence_this_bin =  pd.DataFrame(data = {"integrated_binnum": args.integrated_binnum, "xBbin": xBbin * np.ones_like(phi_dummy), "Q2bin": Q2bin * np.ones_like(phi_dummy),
#     "tbin": tbin * np.ones_like(phi_dummy), "Q2index": index * np.ones_like(phi_dummy), "xB_display": xBs, "Q2_display": Q2_running * np.ones_like(phi_dummy), "t_display": ts, "phi_display": phi_dummy})
#   df_Q2_dependence_this_bin.loc[:, "xsec_KM15_display"]    = printKMarray (df_Q2_dependence_this_bin.xB_display, df_Q2_dependence_this_bin.Q2_display, df_Q2_dependence_this_bin.t_display, np.radians(df_Q2_dependence_this_bin.phi_display), mode = 5)
#   df_Q2_dependence_this_bin.loc[:, "xsec_BH_KM15_display"] = printKMarray (df_Q2_dependence_this_bin.xB_display, df_Q2_dependence_this_bin.Q2_display, df_Q2_dependence_this_bin.t_display, np.radians(df_Q2_dependence_this_bin.phi_display), mode = 1)
#   df_Q2_dependence_this_bin.loc[:, "xsec_BH_display"]      = printBHarray (df_Q2_dependence_this_bin.xB_display, df_Q2_dependence_this_bin.Q2_display, df_Q2_dependence_this_bin.t_display, np.radians(df_Q2_dependence_this_bin.phi_display), local=True)
#   df_Q2_dependence_this_bin.loc[:, "xsec_VGG_display"]     = printVGGarray(df_Q2_dependence_this_bin.xB_display, df_Q2_dependence_this_bin.Q2_display, df_Q2_dependence_this_bin.t_display, np.radians(df_Q2_dependence_this_bin.phi_display), local=True)
#   df_Q2_dependence_this_bin.loc[:, "weight_BH"] = weight_BH(df_Q2_dependence_this_bin.xB_display, df_Q2_dependence_this_bin.Q2_display, df_Q2_dependence_this_bin.t_display, df_Q2_dependence_this_bin.phi_display)
#   df_Q2_dependence_this_bin.loc[:, "xsec_BH_display_w"]       = df_Q2_dependence_this_bin.loc[:, "xsec_BH_display"]     *df_Q2_dependence_this_bin.loc[:, "weight_BH"]
#   df_Q2_dependence_this_bin.loc[:, "xsec_KM15_display_w"]     = df_Q2_dependence_this_bin.loc[:, "xsec_KM15_display"]   *df_Q2_dependence_this_bin.loc[:, "weight_BH"]
#   df_Q2_dependence_this_bin.loc[:, "xsec_BH_KM15_display_w"]  = df_Q2_dependence_this_bin.loc[:, "xsec_BH_KM15_display"]*df_Q2_dependence_this_bin.loc[:, "weight_BH"]
#   df_Q2_dependence_this_bin.loc[:, "xsec_VGG_display_w"]      = df_Q2_dependence_this_bin.loc[:, "xsec_VGG_display"]    *df_Q2_dependence_this_bin.loc[:, "weight_BH"]
#   df_Q2_dependence = pd.concat([df_Q2_dependence, df_Q2_dependence_this_bin])
# df_Q2_dependence.reset_index(inplace = True)
# df_Q2_dependence.rename(columns = {"index": "phiindex"}, inplace = True)
# df_Q2_dependence.to_pickle("addendum_v3/fitting/df_Q2_dependence_{}.pkl".format(args.integrated_binnum))

# df_Q2_dependence = pd.DataFrame()
# n_sample    = 500

# #Q2 dependence
# for xBbin in df_summary_table_rebinned.xBbin.unique():
#   for Q2bin in df_summary_table_rebinned.Q2bin.unique():
#     df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbin) & (df_sumary_table_rebinned.Q2bin == Q2bin), :]
#     xB_avgs = np.unique(df_summary_table_rebinned_this_bin.xB_avg_this_point).mean() * np.ones(n_sample)
#     xB_stds = np.unique(df_summary_table_rebinned_this_bin.xB_avg_this_point).std()  * np.ones(n_sample)
#     Q2_avgs = np.unique(df_summary_table_rebinned_this_bin.Q2_avg_this_point).mean() * np.ones(n_sample)

#   if len(np.unique(df_summary_table_rebinned_this_bin.t_avg_this_point)) > 1:
#     continue
#   ts = np.unique(df_summary_table_rebinned_this_bin.t_avg_this_point)[0] * np.ones(n_sample)
#   phis = np.linspace(0, 2*np.pi, n_sample)
#   weight_BHs_theory = weight_BH(xB_theory, Q2_theory, t_theory, np.degrees(phi_theory))
#   xsec_KM15_theory = printKMarray(xB_theory, Q2_theory, t_theory, phi_theory, mode = 5)
#   xsec_BH_theory   = printBHarray(xB_theory, Q2_theory, t_theory, phi_theory, local=True)
#   xsec_VGG_theory  = printVGGarray(xB_theory, Q2_theory, t_theory, phi_theory, local=True)


# xBs = df_summary_table_rebinned.xB_avg_this_point
# Q2s = df_summary_table_rebinned.Q2_avg_this_point
# ts  = df_summary_table_rebinned.t_avg_this_point
# phis = df_summary_table_rebinned.phi_avg_this_point
# df_summary_table_rebinned.loc[:, "weight_BH"] = weight_BH(xBs, Q2s, ts, phis)

# df_summary_table_rebinned.loc[:, "km15_cross_section_this_point_norad_w"] = df_summary_table_rebinned.loc[:, "km15_cross_section_this_point_norad"]*df_summary_table_rebinned.loc[:, "weight_BH"]
# df_summary_table_rebinned.loc[:, "vgg_cross_section_this_point_norad_w"] = df_summary_table_rebinned.loc[:, "vgg_cross_section_this_point_norad"]*df_summary_table_rebinned.loc[:, "weight_BH"]
# df_summary_table_rebinned.loc[:, "pureBH_cross_section_this_point_norad_w"] = df_summary_table_rebinned.loc[:, "pureBH_cross_section_this_point_norad"]*df_summary_table_rebinned.loc[:, "weight_BH"]
# df_summary_table_rebinned.loc[:, "pureBH_km15_cross_section_this_point_norad_w"] = df_summary_table_rebinned.loc[:, "pureBH_km15_cross_section_this_point_norad"]*df_summary_table_rebinned.loc[:, "weight_BH"]

# df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_w".format(model)]                  = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
# df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_stat_err_w".format(model)]         = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_stat_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
# df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err_w".format(model)]         = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
# # df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err_bin_by_bin_w".format(model)]         = df_summary_table_rebinned.loc[:, "xsec_exp_bkg_merging_syst_err_bin_by_bin".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]


# for model in ['pi0']:
#     df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_w".format(model)]                  = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
#     df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err_w".format(model)]         = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
#     df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_w".format(model)]         = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
#     df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_bin_by_bin_w".format(model)]         = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_bin_by_bin".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
#     # df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_down_w".format(model)]    = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_down".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
#     # df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_up_w".format(model)]      = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_syst_err_up".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]
#     # df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_err_w".format(model)]              = df_summary_table_rebinned.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_err".format(model)]*df_summary_table_rebinned.loc[:, "weight_BH"]

# df_display.loc[:, "weight_BH"] = weight_BH(df_display.xB_display, df_display.Q2_display, df_display.t_display, np.degrees(df_display.phi_display))

# df_display.loc[:, "xsec_BH_display_w"]           = df_display.loc[:, "xsec_BH_display"]     *df_display.loc[:, "weight_BH"]
# df_display.loc[:, "xsec_KM15_display_w"]         = df_display.loc[:, "xsec_KM15_display"]   *df_display.loc[:, "weight_BH"]
# df_display.loc[:, "xsec_BH_KM15_display_w"]      = df_display.loc[:, "xsec_BH_KM15_display"]*df_display.loc[:, "weight_BH"]
# df_display.loc[:, "xsec_VGG_display_w"]          = df_display.loc[:, "xsec_VGG_display"]    *df_display.loc[:, "weight_BH"]


# def get_tdependece_th(xB_avg, Q2_avg, t1_min, t1_max, p0 = (1, 0), n_theory = 100, n_sample = 100):
#     cosine_1_ths_km15     = []
#     cosine_1_ths_km15_err = []
#     cosine_1_ths_vgg      = []
#     cosine_1_ths_vgg_err  = []
#     cosine_1_ths_bh       = []
#     cosine_1_ths_bh_err   = []

#     cosine_2_ths_km15     = []
#     cosine_2_ths_km15_err = []
#     cosine_2_ths_vgg      = []
#     cosine_2_ths_vgg_err  = []
#     cosine_2_ths_bh       = []
#     cosine_2_ths_bh_err   = []

#     cosine_3_ths_km15     = []
#     cosine_3_ths_km15_err = []
#     cosine_3_ths_vgg      = []
#     cosine_3_ths_vgg_err  = []
#     cosine_3_ths_bh       = []
#     cosine_3_ths_bh_err   = []

#     ts_theory = np.linspace(t1_min, t1_max, n_theory)
#     for t_theory in ts_theory:
#         xB_theory  = np.ones(n_sample) * xB_avg
#         Q2_theory  = np.ones(n_sample) * Q2_avg
#         t_theory   = np.ones(n_sample) * t_theory
#         phi_theory = np.linspace(0, 2*np.pi, n_sample)
#         weight_BHs_theory = weight_BH(xB_theory, Q2_theory, t_theory, np.degrees(phi_theory))
#         xsec_KM15_theory = printKMarray(xB_theory, Q2_theory, t_theory, phi_theory, mode = 5)
#         xsec_BH_theory   = printBHarray(xB_theory, Q2_theory, t_theory, phi_theory, local=True)
#         xsec_VGG_theory  = printVGGarray(xB_theory, Q2_theory, t_theory, phi_theory, local=True)

#         xsec_KM15_theory_w = weight_BHs_theory * xsec_KM15_theory
#         xsec_BH_theory_w   = weight_BHs_theory * xsec_BH_theory
#         xsec_VGG_theory_w  = weight_BHs_theory * xsec_VGG_theory

#         popt_1_th_km15, pcov = curve_fit(cosine_fitting_1, phi_theory, xsec_KM15_theory_w, p0 = p0)
#         cosine_1_th_km15 =  -popt_1_th_km15[1]
#         cosine_1_th_km15_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_1_ths_km15    .append(cosine_1_th_km15)
#         cosine_1_ths_km15_err.append(cosine_1_th_km15_err)

#         popt_2_th_km15, pcov = curve_fit(cosine_fitting_2, phi_theory, xsec_KM15_theory_w, p0 = (1, 0, 0))
#         cosine_2_th_km15 =  -popt_2_th_km15[1]
#         cosine_2_th_km15_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_2_ths_km15    .append(cosine_2_th_km15)
#         cosine_2_ths_km15_err.append(cosine_2_th_km15_err)

#         popt_3_th_km15, pcov = curve_fit(cosine_fitting_3, phi_theory, xsec_KM15_theory_w, p0 = (1, 0, 0, 0))
#         cosine_3_th_km15 =  -popt_3_th_km15[1]
#         cosine_3_th_km15_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_3_ths_km15    .append(cosine_3_th_km15)
#         cosine_3_ths_km15_err.append(cosine_3_th_km15_err)

#         popt_1_th_bh, pcov = curve_fit(cosine_fitting_1, phi_theory, xsec_BH_theory_w, p0 = p0)
#         cosine_1_th_bh =  -popt_1_th_bh[1]
#         cosine_1_th_bh_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_1_ths_bh    .append(cosine_1_th_bh)
#         cosine_1_ths_bh_err.append(cosine_1_th_bh_err)

#         popt_2_th_bh, pcov = curve_fit(cosine_fitting_2, phi_theory, xsec_BH_theory_w, p0 = (1, 0, 0))
#         cosine_2_th_bh =  -popt_2_th_bh[1]
#         cosine_2_th_bh_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_2_ths_bh    .append(cosine_2_th_bh)
#         cosine_2_ths_bh_err.append(cosine_2_th_bh_err)

#         popt_3_th_bh, pcov = curve_fit(cosine_fitting_3, phi_theory, xsec_BH_theory_w, p0 = (1, 0, 0, 0))
#         cosine_3_th_bh =  -popt_3_th_bh[1]
#         cosine_3_th_bh_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_3_ths_bh    .append(cosine_3_th_bh)
#         cosine_3_ths_bh_err.append(cosine_3_th_bh_err)
    
#         popt_1_th_vgg, pcov = curve_fit(cosine_fitting_1, phi_theory, xsec_VGG_theory_w, p0 = p0)
#         cosine_1_th_vgg =  -popt_1_th_vgg[1]
#         cosine_1_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_1_ths_vgg    .append(cosine_1_th_vgg)
#         cosine_1_ths_vgg_err.append(cosine_1_th_vgg_err)

#         popt_2_th_vgg, pcov = curve_fit(cosine_fitting_2, phi_theory, xsec_VGG_theory_w, p0 = (1, 0, 0))
#         cosine_2_th_vgg =  -popt_2_th_vgg[1]
#         cosine_2_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_2_ths_vgg    .append(cosine_2_th_vgg)
#         cosine_2_ths_vgg_err.append(cosine_2_th_vgg_err)

#         popt_3_th_vgg, pcov = curve_fit(cosine_fitting_3, phi_theory, xsec_VGG_theory_w, p0 = (1, 0, 0, 0))
#         cosine_3_th_vgg =  -popt_3_th_vgg[1]
#         cosine_3_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_3_ths_vgg    .append(cosine_3_th_vgg)
#         cosine_3_ths_vgg_err.append(cosine_3_th_vgg_err)

#     return ts_theory, cosine_1_ths_km15, cosine_1_ths_km15_err, cosine_1_ths_vgg, cosine_1_ths_vgg_err, cosine_1_ths_bh, cosine_1_ths_bh_err, cosine_2_ths_km15, cosine_2_ths_km15_err, cosine_2_ths_vgg, cosine_2_ths_vgg_err, cosine_2_ths_bh, cosine_2_ths_bh_err, cosine_3_ths_km15, cosine_3_ths_km15_err, cosine_3_ths_vgg, cosine_3_ths_vgg_err, cosine_3_ths_bh, cosine_3_ths_bh_err# popt_1_ths_km15, popt_1_ths_vgg, popt_1_ths_bh, popt_2_ths_km15, popt_2_ths_vgg, popt_2_ths_bh, popt_3_ths_km15, popt_3_ths_vgg, popt_3_ths_bh

# def get_Q2dependence_exp(xBbin, Q2bins, tbin):   
#     t_avgs = []
#     xB_avgs = []
#     Q2s = []
#     cosine_1_exps_pi0     = []
#     cosine_1_exps_pi0_stat_err = []
#     cosine_1_exps_pi0_syst_err_up = []
#     cosine_1_exps_pi0_syst_err_down = []
#     # cosine_1_exps_km15    = []
#     # cosine_1_exps_vgg     = []
#     # cosine_1_exps_bh      = []
#     # cosine_1_exps_global1 = []
#     cosine_1_exps_bkgmerging_only = []

#     cosine_2_exps_pi0     = []
#     cosine_2_exps_pi0_stat_err = []
#     cosine_2_exps_pi0_syst_err_up = []
#     cosine_2_exps_pi0_syst_err_down = []
#     # cosine_2_exps_km15    = []
#     # cosine_2_exps_vgg     = []
#     # cosine_2_exps_bh      = []
#     # cosine_2_exps_global1 = []
#     cosine_2_exps_bkgmerging_only = []

#     cosine_3_exps_pi0     = []
#     cosine_3_exps_pi0_stat_err = []
#     cosine_3_exps_pi0_syst_err_up = []
#     cosine_3_exps_pi0_syst_err_down = []
#     # cosine_3_exps_km15    = []
#     # cosine_3_exps_vgg     = []
#     # cosine_3_exps_bh      = []
#     # cosine_3_exps_global1 = []
#     cosine_3_exps_bkgmerging_only = []

#     for Q2bin in Q2bins:
#         df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbin) & (df_summary_table_rebinned.Q2bin == Q2bin) & (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), :]
#         if len(df_this_bin) < 5 :
#             continue
#         # xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_pi0_norm_err_up, cosine_1_exp_pi0_norm_err_down, cosine_1_exp_km15, cosine_1_exp_vgg, cosine_1_exp_bh, cosine_1_exp_global1, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_km15, popt_1_vgg, popt_1_bh, popt_1_global1 = cosine_1_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         # xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_pi0_norm_err_up, cosine_2_exp_pi0_norm_err_down, cosine_2_exp_km15, cosine_2_exp_vgg, cosine_2_exp_bh, cosine_2_exp_global1, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_km15, popt_2_vgg, popt_2_bh, popt_2_global1 = cosine_2_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         # xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_pi0_norm_err_up, cosine_3_exp_pi0_norm_err_down, cosine_3_exp_km15, cosine_3_exp_vgg, cosine_3_exp_bh, cosine_3_exp_global1, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_km15, popt_3_vgg, popt_3_bh, popt_3_global1 = cosine_3_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_pi0_norm_err_up, cosine_1_exp_pi0_norm_err_down, cosine_1_exp_bkgmerging_only, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_bkgmerging_only = cosine_1_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_pi0_norm_err_up, cosine_2_exp_pi0_norm_err_down, cosine_2_exp_bkgmerging_only, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_bkgmerging_only = cosine_2_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_pi0_norm_err_up, cosine_3_exp_pi0_norm_err_down, cosine_3_exp_bkgmerging_only, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_bkgmerging_only = cosine_3_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         t_avgs.append(t_avg)
#         xB_avgs.append(xB_avg)
#         Q2s.append(Q2_avg)

#         cosine_1_exps_pi0.append(cosine_1_exp_pi0)
#         cosine_1_exps_pi0_stat_err.append(cosine_1_exp_pi0_stat_err)
#         cosine_1_exps_pi0_syst_err_up.append(np.sqrt(cosine_1_exp_pi0_syst_err**2 + cosine_1_exp_pi0_norm_err_up**2))
#         cosine_1_exps_pi0_syst_err_down.append(np.sqrt(cosine_1_exp_pi0_syst_err**2 + cosine_1_exp_pi0_norm_err_down**2))
#         # cosine_1_exps_km15.append(cosine_1_exp_km15)
#         # cosine_1_exps_vgg.append(cosine_1_exp_vgg)
#         # cosine_1_exps_bh.append(cosine_1_exp_bh)
#         # cosine_1_exps_global1.append(cosine_1_exp_global1)
#         cosine_1_exps_bkgmerging_only.append(cosine_1_exp_bkgmerging_only)

#         cosine_2_exps_pi0.append(cosine_2_exp_pi0)
#         cosine_2_exps_pi0_stat_err.append(cosine_2_exp_pi0_stat_err)
#         cosine_2_exps_pi0_syst_err_up.append(np.sqrt(cosine_2_exp_pi0_syst_err**2 + cosine_2_exp_pi0_norm_err_up**2))
#         cosine_2_exps_pi0_syst_err_down.append(np.sqrt(cosine_2_exp_pi0_syst_err**2 + cosine_2_exp_pi0_norm_err_down**2))
#         # cosine_2_exps_km15.append(cosine_2_exp_km15)
#         # cosine_2_exps_vgg.append(cosine_2_exp_vgg)
#         # cosine_2_exps_bh.append(cosine_2_exp_bh)
#         # cosine_2_exps_global1.append(cosine_2_exp_global1)
#         cosine_2_exps_bkgmerging_only.append(cosine_2_exp_bkgmerging_only)

#         cosine_3_exps_pi0.append(cosine_3_exp_pi0)
#         cosine_3_exps_pi0_stat_err.append(cosine_3_exp_pi0_stat_err)
#         cosine_3_exps_pi0_syst_err_up.append(np.sqrt(cosine_3_exp_pi0_syst_err**2 + cosine_3_exp_pi0_norm_err_up**2))
#         cosine_3_exps_pi0_syst_err_down.append(np.sqrt(cosine_3_exp_pi0_syst_err**2 + cosine_3_exp_pi0_norm_err_down**2))
#         # cosine_3_exps_km15.append(cosine_3_exp_km15)
#         # cosine_3_exps_vgg.append(cosine_3_exp_vgg)
#         # cosine_3_exps_bh.append(cosine_3_exp_bh)
#         # cosine_3_exps_global1.append(cosine_3_exp_global1)
#         cosine_3_exps_bkgmerging_only.append(cosine_3_exp_bkgmerging_only)

#     return Q2s, xB_avgs, t_avgs, np.array(cosine_1_exps_pi0), np.array(cosine_1_exps_pi0_stat_err), np.array(cosine_1_exps_pi0_syst_err_up), np.array(cosine_1_exps_pi0_syst_err_down), np.array(cosine_1_exps_bkgmerging_only), np.array(cosine_2_exps_pi0), np.array(cosine_2_exps_pi0_stat_err), np.array(cosine_2_exps_pi0_syst_err_up), np.array(cosine_2_exps_pi0_syst_err_down), np.array(cosine_2_exps_bkgmerging_only), np.array(cosine_3_exps_pi0), np.array(cosine_3_exps_pi0_stat_err), np.array(cosine_3_exps_pi0_syst_err_up), np.array(cosine_3_exps_pi0_syst_err_down), np.array(cosine_3_exps_bkgmerging_only)

# def get_tdependence_exp(xBbin, Q2bin, tbins):   
#     ts = []
#     xB_avgs = []
#     Q2_avgs = []
#     cosine_1_exps_pi0     = []
#     cosine_1_exps_pi0_stat_err = []
#     cosine_1_exps_pi0_syst_err_up = []
#     cosine_1_exps_pi0_syst_err_down = []
#     # cosine_1_exps_km15    = []
#     # cosine_1_exps_vgg     = []
#     # cosine_1_exps_bh      = []
#     # cosine_1_exps_global1 = []
#     cosine_1_exps_bkgmerging_only = []

#     cosine_2_exps_pi0     = []
#     cosine_2_exps_pi0_stat_err = []
#     cosine_2_exps_pi0_syst_err_up = []
#     cosine_2_exps_pi0_syst_err_down = []
#     # cosine_2_exps_km15    = []
#     # cosine_2_exps_vgg     = []
#     # cosine_2_exps_bh      = []
#     # cosine_2_exps_global1 = []
#     cosine_2_exps_bkgmerging_only = []

#     cosine_3_exps_pi0     = []
#     cosine_3_exps_pi0_stat_err = []
#     cosine_3_exps_pi0_syst_err_up = []
#     cosine_3_exps_pi0_syst_err_down = []
#     # cosine_3_exps_km15    = []
#     # cosine_3_exps_vgg     = []
#     # cosine_3_exps_bh      = []
#     # cosine_3_exps_global1 = []
#     cosine_3_exps_bkgmerging_only = []

#     for tbin in tbins:
#         df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbin) & (df_summary_table_rebinned.Q2bin == Q2bin) & (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), :]
#         if len(df_this_bin) < 5 :
#             continue
#         # xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_pi0_norm_err_up, cosine_1_exp_pi0_norm_err_down, cosine_1_exp_km15, cosine_1_exp_vgg, cosine_1_exp_bh, cosine_1_exp_global1, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_km15, popt_1_vgg, popt_1_bh, popt_1_global1 = cosine_1_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         # xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_pi0_norm_err_up, cosine_2_exp_pi0_norm_err_down, cosine_2_exp_km15, cosine_2_exp_vgg, cosine_2_exp_bh, cosine_2_exp_global1, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_km15, popt_2_vgg, popt_2_bh, popt_2_global1 = cosine_2_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         # xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_pi0_norm_err_up, cosine_3_exp_pi0_norm_err_down, cosine_3_exp_km15, cosine_3_exp_vgg, cosine_3_exp_bh, cosine_3_exp_global1, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_km15, popt_3_vgg, popt_3_bh, popt_3_global1 = cosine_3_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_pi0_norm_err_up, cosine_1_exp_pi0_norm_err_down, cosine_1_exp_bkgmerging_only, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_bkgmerging_only = cosine_1_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_pi0_norm_err_up, cosine_2_exp_pi0_norm_err_down, cosine_2_exp_bkgmerging_only, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_bkgmerging_only = cosine_2_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_pi0_norm_err_up, cosine_3_exp_pi0_norm_err_down, cosine_3_exp_bkgmerging_only, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_bkgmerging_only = cosine_3_fitting_one_bin_exp(xBbin, Q2bin, tbin)
#         ts.append(t_avg)
#         xB_avgs.append(xB_avg)
#         Q2_avgs.append(Q2_avg)

#         cosine_1_exps_pi0.append(cosine_1_exp_pi0)
#         cosine_1_exps_pi0_stat_err.append(cosine_1_exp_pi0_stat_err)
#         cosine_1_exps_pi0_syst_err_up.append(np.sqrt(cosine_1_exp_pi0_syst_err**2 + cosine_1_exp_pi0_norm_err_up**2))
#         cosine_1_exps_pi0_syst_err_down.append(np.sqrt(cosine_1_exp_pi0_syst_err**2 + cosine_1_exp_pi0_norm_err_down**2))
#         # cosine_1_exps_km15.append(cosine_1_exp_km15)
#         # cosine_1_exps_vgg.append(cosine_1_exp_vgg)
#         # cosine_1_exps_bh.append(cosine_1_exp_bh)
#         # cosine_1_exps_global1.append(cosine_1_exp_global1)
#         cosine_1_exps_bkgmerging_only.append(cosine_1_exp_bkgmerging_only)

#         cosine_2_exps_pi0.append(cosine_2_exp_pi0)
#         cosine_2_exps_pi0_stat_err.append(cosine_2_exp_pi0_stat_err)
#         cosine_2_exps_pi0_syst_err_up.append(np.sqrt(cosine_2_exp_pi0_syst_err**2 + cosine_2_exp_pi0_norm_err_up**2))
#         cosine_2_exps_pi0_syst_err_down.append(np.sqrt(cosine_2_exp_pi0_syst_err**2 + cosine_2_exp_pi0_norm_err_down**2))
#         # cosine_2_exps_km15.append(cosine_2_exp_km15)
#         # cosine_2_exps_vgg.append(cosine_2_exp_vgg)
#         # cosine_2_exps_bh.append(cosine_2_exp_bh)
#         # cosine_2_exps_global1.append(cosine_2_exp_global1)
#         cosine_2_exps_bkgmerging_only.append(cosine_2_exp_bkgmerging_only)

#         cosine_3_exps_pi0.append(cosine_3_exp_pi0)
#         cosine_3_exps_pi0_stat_err.append(cosine_3_exp_pi0_stat_err)
#         cosine_3_exps_pi0_syst_err_up.append(np.sqrt(cosine_3_exp_pi0_syst_err**2 + cosine_3_exp_pi0_norm_err_up**2))
#         cosine_3_exps_pi0_syst_err_down.append(np.sqrt(cosine_3_exp_pi0_syst_err**2 + cosine_3_exp_pi0_norm_err_down**2))
#         # cosine_3_exps_km15.append(cosine_3_exp_km15)
#         # cosine_3_exps_vgg.append(cosine_3_exp_vgg)
#         # cosine_3_exps_bh.append(cosine_3_exp_bh)
#         # cosine_3_exps_global1.append(cosine_3_exp_global1)
#         cosine_3_exps_bkgmerging_only.append(cosine_3_exp_bkgmerging_only)

#     # return ts, xB_avgs, Q2_avgs, np.array(cosine_1_exps_pi0), np.array(cosine_1_exps_pi0_stat_err), np.array(cosine_1_exps_pi0_syst_err_up), np.array(cosine_1_exps_pi0_syst_err_down), np.array(cosine_1_exps_pi0_syst_err_down),,  np.array(cosine_1_exps_km15), np.array(cosine_1_exps_vgg), np.array(cosine_1_exps_bh), np.array(cosine_1_exps_global1), np.array(cosine_2_exps_pi0), np.array(cosine_2_exps_pi0_stat_err), np.array(cosine_2_exps_pi0_syst_err_up), np.array(cosine_2_exps_pi0_syst_err_down),,  np.array(cosine_2_exps_km15), np.array(cosine_2_exps_vgg), np.array(cosine_2_exps_bh), np.array(cosine_2_exps_global1), np.array(cosine_3_exps_pi0), np.array(cosine_3_exps_pi0_stat_err), np.array(cosine_3_exps_pi0_syst_err_up), np.array(cosine_3_exps_pi0_syst_err_down),,  np.array(cosine_3_exps_km15), np.array(cosine_3_exps_vgg), np.array(cosine_3_exps_bh), np.array(cosine_3_exps_global1)
#     # return ts, xB_avgs, Q2_avgs, np.array(cosine_1_exps_pi0), np.array(cosine_1_exps_pi0_stat_err), np.array(cosine_1_exps_pi0_syst_err_up), np.array(cosine_1_exps_pi0_syst_err_down),,  np.array(cosine_2_exps_pi0), np.array(cosine_2_exps_pi0_stat_err), np.array(cosine_2_exps_pi0_syst_err_up), np.array(cosine_2_exps_pi0_syst_err_down),,  np.array(cosine_3_exps_pi0), np.array(cosine_3_exps_pi0_stat_err), np.array(cosine_3_exps_pi0_syst_err_up), np.array(cosine_3_exps_pi0_syst_err_down), 
#     # return ts, xB_avgs, Q2_avgs, np.array(cosine_1_exps_pi0), np.array(cosine_1_exps_pi0_stat_err), np.array(cosine_1_exps_pi0_syst_err_up), np.array(cosine_1_exps_bkgmerging_only), np.array(cosine_2_exps_pi0), np.array(cosine_2_exps_pi0_stat_err), np.array(cosine_2_exps_pi0_syst_err_up), np.array(cosine_2_exps_pi0_syst_err_down), np.array(cosine_2_exps_bkgmerging_only), np.array(cosine_3_exps_pi0), np.array(cosine_3_exps_pi0_stat_err), np.array(cosine_3_exps_pi0_syst_err_up), np.array(cosine_3_exps_pi0_syst_err_down), np.array(cosine_3_exps_bkgmerging_only)
#     return ts, xB_avgs, Q2_avgs, np.array(cosine_1_exps_pi0), np.array(cosine_1_exps_pi0_stat_err), np.array(cosine_1_exps_pi0_syst_err_up), np.array(cosine_1_exps_pi0_syst_err_down), np.array(cosine_1_exps_bkgmerging_only), np.array(cosine_2_exps_pi0), np.array(cosine_2_exps_pi0_stat_err), np.array(cosine_2_exps_pi0_syst_err_up), np.array(cosine_2_exps_pi0_syst_err_down), np.array(cosine_2_exps_bkgmerging_only), np.array(cosine_3_exps_pi0), np.array(cosine_3_exps_pi0_stat_err), np.array(cosine_3_exps_pi0_syst_err_up), np.array(cosine_3_exps_pi0_syst_err_down), np.array(cosine_3_exps_bkgmerging_only)    

# def get_Q2dependece_th(xB_avg, Q2_min, Q2_max, t_avg, p0 = (1, 0), n_theory = 100, n_sample = 100):
#     cosine_1_ths_km15     = []
#     cosine_1_ths_km15_err = []
#     cosine_1_ths_vgg      = []
#     cosine_1_ths_vgg_err  = []
#     cosine_1_ths_bh       = []
#     cosine_1_ths_bh_err   = []

#     cosine_2_ths_km15     = []
#     cosine_2_ths_km15_err = []
#     cosine_2_ths_vgg      = []
#     cosine_2_ths_vgg_err  = []
#     cosine_2_ths_bh       = []
#     cosine_2_ths_bh_err   = []

#     cosine_3_ths_km15     = []
#     cosine_3_ths_km15_err = []
#     cosine_3_ths_vgg      = []
#     cosine_3_ths_vgg_err  = []
#     cosine_3_ths_bh       = []
#     cosine_3_ths_bh_err   = []

#     Q2s_theory = np.linspace(Q2_min, Q2_max, n_theory)
#     for Q2_theory in Q2s_theory:
#         xB_theory  = np.ones(n_sample) * xB_avg
#         Q2_theory  = np.ones(n_sample) * Q2_theory
#         t_theory   = np.ones(n_sample) * t_avg
#         phi_theory = np.linspace(0, 2*np.pi, n_sample)

#         weight_BHs_theory = weight_BH(xB_theory, Q2_theory, t_theory, np.degrees(phi_theory))
#         xsec_KM15_theory = printKMarray(xB_theory, Q2_theory, t_theory, phi_theory, mode = 5)
#         xsec_BH_theory   = printBHarray(xB_theory, Q2_theory, t_theory, phi_theory, local=True)
#         xsec_VGG_theory  = printVGGarray(xB_theory, Q2_theory, t_theory, phi_theory, local=True)

#         xsec_KM15_theory_w = weight_BHs_theory * xsec_KM15_theory
#         xsec_BH_theory_w   = weight_BHs_theory * xsec_BH_theory
#         xsec_VGG_theory_w  = weight_BHs_theory * xsec_VGG_theory

#         popt_1_th_km15, pcov = curve_fit(cosine_fitting_1, phi_theory, xsec_KM15_theory_w, p0 = p0)
#         cosine_1_th_km15 =  -popt_1_th_km15[1]
#         cosine_1_th_km15_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_1_ths_km15    .append(cosine_1_th_km15)
#         cosine_1_ths_km15_err.append(cosine_1_th_km15_err)

#         popt_2_th_km15, pcov = curve_fit(cosine_fitting_2, phi_theory, xsec_KM15_theory_w, p0 = (1, 0, 0))
#         cosine_2_th_km15 =  -popt_2_th_km15[1]
#         cosine_2_th_km15_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_2_ths_km15    .append(cosine_2_th_km15)
#         cosine_2_ths_km15_err.append(cosine_2_th_km15_err)

#         popt_3_th_km15, pcov = curve_fit(cosine_fitting_3, phi_theory, xsec_KM15_theory_w, p0 = (1, 0, 0, 0))
#         cosine_3_th_km15 =  -popt_3_th_km15[1]
#         cosine_3_th_km15_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_3_ths_km15    .append(cosine_3_th_km15)
#         cosine_3_ths_km15_err.append(cosine_3_th_km15_err)

#         popt_1_th_bh, pcov = curve_fit(cosine_fitting_1, phi_theory, xsec_BH_theory_w, p0 = p0)
#         cosine_1_th_bh =  -popt_1_th_bh[1]
#         cosine_1_th_bh_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_1_ths_bh    .append(cosine_1_th_bh)
#         cosine_1_ths_bh_err.append(cosine_1_th_bh_err)

#         popt_2_th_bh, pcov = curve_fit(cosine_fitting_2, phi_theory, xsec_BH_theory_w, p0 = (1, 0, 0))
#         cosine_2_th_bh =  -popt_2_th_bh[1]
#         cosine_2_th_bh_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_2_ths_bh    .append(cosine_2_th_bh)
#         cosine_2_ths_bh_err.append(cosine_2_th_bh_err)

#         popt_3_th_bh, pcov = curve_fit(cosine_fitting_3, phi_theory, xsec_BH_theory_w, p0 = (1, 0, 0, 0))
#         cosine_3_th_bh =  -popt_3_th_bh[1]
#         cosine_3_th_bh_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_3_ths_bh    .append(cosine_3_th_bh)
#         cosine_3_ths_bh_err.append(cosine_3_th_bh_err)
    
#         popt_1_th_vgg, pcov = curve_fit(cosine_fitting_1, phi_theory, xsec_VGG_theory_w, p0 = p0)
#         cosine_1_th_vgg =  -popt_1_th_vgg[1]
#         cosine_1_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_1_ths_vgg    .append(cosine_1_th_vgg)
#         cosine_1_ths_vgg_err.append(cosine_1_th_vgg_err)

#         popt_2_th_vgg, pcov = curve_fit(cosine_fitting_2, phi_theory, xsec_VGG_theory_w, p0 = (1, 0, 0))
#         cosine_2_th_vgg =  -popt_2_th_vgg[1]
#         cosine_2_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_2_ths_vgg    .append(cosine_2_th_vgg)
#         cosine_2_ths_vgg_err.append(cosine_2_th_vgg_err)

#         popt_3_th_vgg, pcov = curve_fit(cosine_fitting_3, phi_theory, xsec_VGG_theory_w, p0 = (1, 0, 0, 0))
#         cosine_3_th_vgg =  -popt_3_th_vgg[1]
#         cosine_3_th_vgg_err =  np.sqrt(np.diag(pcov))[1]
#         cosine_3_ths_vgg    .append(cosine_3_th_vgg)
#         cosine_3_ths_vgg_err.append(cosine_3_th_vgg_err)

#     return Q2s_theory, cosine_1_ths_km15, cosine_1_ths_km15_err, cosine_1_ths_vgg, cosine_1_ths_vgg_err, cosine_1_ths_bh, cosine_1_ths_bh_err, cosine_2_ths_km15, cosine_2_ths_km15_err, cosine_2_ths_vgg, cosine_2_ths_vgg_err, cosine_2_ths_bh, cosine_2_ths_bh_err, cosine_3_ths_km15, cosine_3_ths_km15_err, cosine_3_ths_vgg, cosine_3_ths_vgg_err, cosine_3_ths_bh, cosine_3_ths_bh_err# popt_1_ths_km15, popt_1_ths_vgg, popt_1_ths_bh, popt_2_ths_km15, popt_2_ths_vgg, popt_2_ths_bh, popt_3_ths_km15, popt_3_ths_vgg, popt_3_ths_bh    

# phi_dummy = np.linspace(0, 360, 101)
# fig, axs = plt.subplots(1, 1, figsize = (10, 6))
# integrated_binnum = 64

# df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal == 1), :] 
# xBbin = df_this_bin.xBbin.unique()[0]
# Q2bin = df_this_bin.Q2bin.unique()[0]
# tbin  = df_this_bin.tbin.unique()[0]

# # xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_pi0_norm_err_up, cosine_1_exp_pi0_norm_err_down, cosine_1_exp_km15, cosine_1_exp_vgg, cosine_1_exp_bh, cosine_1_exp_global1, cosine_1_th_km15, cosine_1_th_km15_err, cosine_1_th_vgg, cosine_1_th_vgg_err, cosine_1_th_bh, cosine_1_th_bh_err, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_km15, popt_1_vgg, popt_1_bh, popt_1_global1, popt_1_th_km15, popt_1_th_vgg, popt_1_th_bh = cosine_1_fitting_one_bin(xBbin, Q2bin, tbin)
# # xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_pi0_norm_err_up, cosine_2_exp_pi0_norm_err_down, cosine_2_exp_km15, cosine_2_exp_vgg, cosine_2_exp_bh, cosine_2_exp_global1, cosine_2_th_km15, cosine_2_th_km15_err, cosine_2_th_vgg, cosine_2_th_vgg_err, cosine_2_th_bh, cosine_2_th_bh_err, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_km15, popt_2_vgg, popt_2_bh, popt_2_global1, popt_2_th_km15, popt_2_th_vgg, popt_2_th_bh = cosine_2_fitting_one_bin(xBbin, Q2bin, tbin)
# # xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_pi0_norm_err_up, cosine_3_exp_pi0_norm_err_down, cosine_3_exp_km15, cosine_3_exp_vgg, cosine_3_exp_bh, cosine_3_exp_global1, cosine_3_th_km15, cosine_3_th_km15_err, cosine_3_th_vgg, cosine_3_th_vgg_err, cosine_3_th_bh, cosine_3_th_bh_err, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_km15, popt_3_vgg, popt_3_bh, popt_3_global1, popt_3_th_km15, popt_3_th_vgg, popt_3_th_bh = cosine_3_fitting_one_bin(xBbin, Q2bin, tbin)
# xB_avg, Q2_avg, t_avg, cosine_1_exp_pi0, cosine_1_exp_pi0_stat_err, cosine_1_exp_pi0_syst_err, cosine_1_exp_pi0_norm_err_up, cosine_1_exp_pi0_norm_err_down, cosine_1_exp_bkgmerging_only, cosine_1_th_km15, cosine_1_th_km15_err, cosine_1_th_vgg, cosine_1_th_vgg_err, cosine_1_th_bh, cosine_1_th_bh_err, popt_1_pi0, popt_1_pi0_min, popt_1_pi0_max, popt_1_bkgmerging_only, popt_1_th_km15, popt_1_th_vgg, popt_1_th_bh = cosine_1_fitting_one_bin(xBbin, Q2bin, tbin)
# xB_avg, Q2_avg, t_avg, cosine_2_exp_pi0, cosine_2_exp_pi0_stat_err, cosine_2_exp_pi0_syst_err, cosine_2_exp_pi0_norm_err_up, cosine_2_exp_pi0_norm_err_down, cosine_2_exp_bkgmerging_only, cosine_2_th_km15, cosine_2_th_km15_err, cosine_2_th_vgg, cosine_2_th_vgg_err, cosine_2_th_bh, cosine_2_th_bh_err, popt_2_pi0, popt_2_pi0_min, popt_2_pi0_max, popt_2_bkgmerging_only, popt_2_th_km15, popt_2_th_vgg, popt_2_th_bh = cosine_2_fitting_one_bin(xBbin, Q2bin, tbin)
# xB_avg, Q2_avg, t_avg, cosine_3_exp_pi0, cosine_3_exp_pi0_stat_err, cosine_3_exp_pi0_syst_err, cosine_3_exp_pi0_norm_err_up, cosine_3_exp_pi0_norm_err_down, cosine_3_exp_bkgmerging_only, cosine_3_th_km15, cosine_3_th_km15_err, cosine_3_th_vgg, cosine_3_th_vgg_err, cosine_3_th_bh, cosine_3_th_bh_err, popt_3_pi0, popt_3_pi0_min, popt_3_pi0_max, popt_3_bkgmerging_only, popt_3_th_km15, popt_3_th_vgg, popt_3_th_bh = cosine_3_fitting_one_bin(xBbin, Q2bin, tbin)
# plt.errorbar(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, yerr = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, ls = '', marker = 'o', color = 'k')
# plt.errorbar(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_bkg_merging_w, yerr = df_this_bin.xsec_exp_bkg_merging_stat_err_w, ls = '', marker = 'o', color = 'tab:purple')

# # plt.errorbar(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_bh_eff_corrected_bkg_merging_w, ls = '', marker = 'o', color = 'r')
# # plt.errorbar(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_vgg_eff_corrected_bkg_merging_w, ls = '', marker = 'o', color = 'tab:orange')
# plt.fill_between(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, alpha = .5)

# axs.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_pi0), color = 'k')
# # axs.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_pi0_min), color = 'r')
# # axs.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_pi0_max), color = 'k')
# # axs.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_km15), color = 'k', lw = 3)
# # axs.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_vgg), color = 'k', lw = 3)
# # axs.plot(phi_dummy, cosine_fitting_3(np.radians(phi_dummy), *popt_3_th_bh), color = 'k', lw = 3)

# df_display_this_bin = df_display.loc[df_display.integratedbin_display == integrated_binnum, :]

# axs.plot(np.degrees(df_display_this_bin.phi_display), df_display_this_bin.xsec_KM15_display_w, color = 'cyan')
# axs.plot(np.degrees(df_display_this_bin.phi_display), df_display_this_bin.xsec_BH_display_w, color = 'red')
# axs.plot(np.degrees(df_display_this_bin.phi_display), df_display_this_bin.xsec_VGG_display_w, color = 'tab:orange')

# plt.xlim([0, 360])
# plt.xticks([0, 90, 180, 270, 360])
# plt.savefig("addendum_v3/modified_one_bin_test.from_script.pdf")




# if "fitting_result" not in locals().keys():
#     fitting_result = {}

# with open('addendum_v3/fitting_result.npy', 'rb') as f:
#     for key in ["xBbins_panel00", "Q2bins_panel00", "tbins_panel00", "tmin_panel00", "tmax_panel00", "ts_panel00", "xB_avgs_panel00", "Q2_avgs_panel00", "cosine_1_exps_pi0_panel00", "cosine_1_exps_pi0_stat_err_panel00", "cosine_1_exps_pi0_syst_err_up_panel00", "cosine_1_exps_pi0_syst_err_down_panel00", "cosine_1_exps_km15_panel00", "cosine_1_exps_vgg_panel00", "cosine_1_exps_bh_panel00", "cosine_1_exps_global1_panel00", "cosine_2_exps_pi0_panel00", "cosine_2_exps_pi0_stat_err_panel00", "cosine_2_exps_pi0_syst_err_up_panel00", "cosine_2_exps_pi0_syst_err_down_panel00", "cosine_2_exps_km15_panel00", "cosine_2_exps_vgg_panel00", "cosine_2_exps_bh_panel00", "cosine_2_exps_global1_panel00", "cosine_3_exps_pi0_panel00", "cosine_3_exps_pi0_stat_err_panel00", "cosine_3_exps_pi0_syst_err_up_panel00", "cosine_3_exps_pi0_syst_err_down_panel00", "cosine_3_exps_km15_panel00", "cosine_3_exps_vgg_panel00", "cosine_3_exps_bh_panel00", "cosine_3_exps_global1_panel00", "ts_theory_panel00", "cosine_1_ths_km15_panel00", "cosine_1_ths_km15_err_panel00", "cosine_1_ths_vgg_panel00", "cosine_1_ths_vgg_err_panel00", "cosine_1_ths_bh_panel00", "cosine_1_ths_bh_err_panel00", "cosine_2_ths_km15_panel00", "cosine_2_ths_km15_err_panel00", "cosine_2_ths_vgg_panel00", "cosine_2_ths_vgg_err_panel00", "cosine_2_ths_bh_panel00", "cosine_2_ths_bh_err_panel00", "cosine_3_ths_km15_panel00", "cosine_3_ths_km15_err_panel00", "cosine_3_ths_vgg_panel00", "cosine_3_ths_vgg_err_panel00", "cosine_3_ths_bh_panel00", "cosine_3_ths_bh_err_panel00", "xBbins_panel10", "Q2bins_panel10", "tbins_panel10", "tmin_panel10", "tmax_panel10", "ts_panel10", "xB_avgs_panel10", "Q2_avgs_panel10", "cosine_1_exps_pi0_panel10", "cosine_1_exps_pi0_stat_err_panel10", "cosine_1_exps_pi0_syst_err_up_panel10", "cosine_1_exps_pi0_syst_err_down_panel10", "cosine_1_exps_km15_panel10", "cosine_1_exps_vgg_panel10", "cosine_1_exps_bh_panel10", "cosine_1_exps_global1_panel10", "cosine_2_exps_pi0_panel10", "cosine_2_exps_pi0_stat_err_panel10", "cosine_2_exps_pi0_syst_err_up_panel10", "cosine_2_exps_pi0_syst_err_down_panel10", "cosine_2_exps_km15_panel10", "cosine_2_exps_vgg_panel10", "cosine_2_exps_bh_panel10", "cosine_2_exps_global1_panel10", "cosine_3_exps_pi0_panel10", "cosine_3_exps_pi0_stat_err_panel10", "cosine_3_exps_pi0_syst_err_up_panel10", "cosine_3_exps_pi0_syst_err_down_panel10", "cosine_3_exps_km15_panel10", "cosine_3_exps_vgg_panel10", "cosine_3_exps_bh_panel10", "cosine_3_exps_global1_panel10", "ts_theory_panel10", "cosine_1_ths_km15_panel10", "cosine_1_ths_km15_err_panel10", "cosine_1_ths_vgg_panel10", "cosine_1_ths_vgg_err_panel10", "cosine_1_ths_bh_panel10", "cosine_1_ths_bh_err_panel10", "cosine_2_ths_km15_panel10", "cosine_2_ths_km15_err_panel10", "cosine_2_ths_vgg_panel10", "cosine_2_ths_vgg_err_panel10", "cosine_2_ths_bh_panel10", "cosine_2_ths_bh_err_panel10", "cosine_3_ths_km15_panel10", "cosine_3_ths_km15_err_panel10", "cosine_3_ths_vgg_panel10", "cosine_3_ths_vgg_err_panel10", "cosine_3_ths_bh_panel10", "cosine_3_ths_bh_err_panel10", "xBbins_panel20", "Q2bins_panel20", "tbins_panel20", "tmin_panel20", "tmax_panel20", "ts_panel20", "xB_avgs_panel20", "Q2_avgs_panel20", "cosine_1_exps_pi0_panel20", "cosine_1_exps_pi0_stat_err_panel20", "cosine_1_exps_pi0_syst_err_up_panel20", "cosine_1_exps_pi0_syst_err_down_panel20", "cosine_1_exps_km15_panel20", "cosine_1_exps_vgg_panel20", "cosine_1_exps_bh_panel20", "cosine_1_exps_global1_panel20", "cosine_2_exps_pi0_panel20", "cosine_2_exps_pi0_stat_err_panel20", "cosine_2_exps_pi0_syst_err_up_panel20", "cosine_2_exps_pi0_syst_err_down_panel20", "cosine_2_exps_km15_panel20", "cosine_2_exps_vgg_panel20", "cosine_2_exps_bh_panel20", "cosine_2_exps_global1_panel20", "cosine_3_exps_pi0_panel20", "cosine_3_exps_pi0_stat_err_panel20", "cosine_3_exps_pi0_syst_err_up_panel20", "cosine_3_exps_pi0_syst_err_down_panel20", "cosine_3_exps_km15_panel20", "cosine_3_exps_vgg_panel20", "cosine_3_exps_bh_panel20", "cosine_3_exps_global1_panel20", "ts_theory_panel20", "cosine_1_ths_km15_panel20", "cosine_1_ths_km15_err_panel20", "cosine_1_ths_vgg_panel20", "cosine_1_ths_vgg_err_panel20", "cosine_1_ths_bh_panel20", "cosine_1_ths_bh_err_panel20", "cosine_2_ths_km15_panel20", "cosine_2_ths_km15_err_panel20", "cosine_2_ths_vgg_panel20", "cosine_2_ths_vgg_err_panel20", "cosine_2_ths_bh_panel20", "cosine_2_ths_bh_err_panel20", "cosine_3_ths_km15_panel20", "cosine_3_ths_km15_err_panel20", "cosine_3_ths_vgg_panel20", "cosine_3_ths_vgg_err_panel20", "cosine_3_ths_bh_panel20", "cosine_3_ths_bh_err_panel20", "xBbins_panel01", "Q2bins_panel01", "tbins_panel01", "Q2min_panel01", "Q2max_panel01", "Q2s_panel01", "xB_avgs_panel01", "t_avgs_panel01", "cosine_1_exps_pi0_panel01", "cosine_1_exps_pi0_stat_err_panel01", "cosine_1_exps_pi0_syst_err_up_panel01", "cosine_1_exps_pi0_syst_err_down_panel01", "cosine_1_exps_km15_panel01", "cosine_1_exps_vgg_panel01", "cosine_1_exps_bh_panel01", "cosine_1_exps_global1_panel01", "cosine_2_exps_pi0_panel01", "cosine_2_exps_pi0_stat_err_panel01", "cosine_2_exps_pi0_syst_err_up_panel01", "cosine_2_exps_pi0_syst_err_down_panel01", "cosine_2_exps_km15_panel01", "cosine_2_exps_vgg_panel01", "cosine_2_exps_bh_panel01", "cosine_2_exps_global1_panel01", "cosine_3_exps_pi0_panel01", "cosine_3_exps_pi0_stat_err_panel01", "cosine_3_exps_pi0_syst_err_up_panel01", "cosine_3_exps_pi0_syst_err_down_panel01", "cosine_3_exps_km15_panel01", "cosine_3_exps_vgg_panel01", "cosine_3_exps_bh_panel01", "cosine_3_exps_global1_panel01", "Q2s_theory_panel01", "cosine_1_ths_km15_panel01", "cosine_1_ths_km15_err_panel01", "cosine_1_ths_vgg_panel01", "cosine_1_ths_vgg_err_panel01", "cosine_1_ths_bh_panel01", "cosine_1_ths_bh_err_panel01", "cosine_2_ths_km15_panel01", "cosine_2_ths_km15_err_panel01", "cosine_2_ths_vgg_panel01", "cosine_2_ths_vgg_err_panel01", "cosine_2_ths_bh_panel01", "cosine_2_ths_bh_err_panel01", "cosine_3_ths_km15_panel01", "cosine_3_ths_km15_err_panel01", "cosine_3_ths_vgg_panel01", "cosine_3_ths_vgg_err_panel01", "cosine_3_ths_bh_panel01", "cosine_3_ths_bh_err_panel01", "xBbins_panel11", "Q2bins_panel11", "tbins_panel11", "Q2min_panel11", "Q2max_panel11", "Q2s_panel11", "xB_avgs_panel11", "t_avgs_panel11", "cosine_1_exps_pi0_panel11", "cosine_1_exps_pi0_stat_err_panel11", "cosine_1_exps_pi0_syst_err_up_panel11", "cosine_1_exps_pi0_syst_err_down_panel11", "cosine_1_exps_km15_panel11", "cosine_1_exps_vgg_panel11", "cosine_1_exps_bh_panel11", "cosine_1_exps_global1_panel11", "cosine_2_exps_pi0_panel11", "cosine_2_exps_pi0_stat_err_panel11", "cosine_2_exps_pi0_syst_err_up_panel11", "cosine_2_exps_pi0_syst_err_down_panel11", "cosine_2_exps_km15_panel11", "cosine_2_exps_vgg_panel11", "cosine_2_exps_bh_panel11", "cosine_2_exps_global1_panel11", "cosine_3_exps_pi0_panel11", "cosine_3_exps_pi0_stat_err_panel11", "cosine_3_exps_pi0_syst_err_up_panel11", "cosine_3_exps_pi0_syst_err_down_panel11", "cosine_3_exps_km15_panel11", "cosine_3_exps_vgg_panel11", "cosine_3_exps_bh_panel11", "cosine_3_exps_global1_panel11", "Q2s_theory_panel11", "cosine_1_ths_km15_panel11", "cosine_1_ths_km15_err_panel11", "cosine_1_ths_vgg_panel11", "cosine_1_ths_vgg_err_panel11", "cosine_1_ths_bh_panel11", "cosine_1_ths_bh_err_panel11", "cosine_2_ths_km15_panel11", "cosine_2_ths_km15_err_panel11", "cosine_2_ths_vgg_panel11", "cosine_2_ths_vgg_err_panel11", "cosine_2_ths_bh_panel11", "cosine_2_ths_bh_err_panel11", "cosine_3_ths_km15_panel11", "cosine_3_ths_km15_err_panel11", "cosine_3_ths_vgg_panel11", "cosine_3_ths_vgg_err_panel11", "cosine_3_ths_bh_panel11", "cosine_3_ths_bh_err_panel11", "xBbins_panel21", "Q2bins_panel21", "tbins_panel21", "Q2min_panel21", "Q2max_panel21", "Q2s_panel21", "xB_avgs_panel21", "t_avgs_panel21", "cosine_1_exps_pi0_panel21", "cosine_1_exps_pi0_stat_err_panel21", "cosine_1_exps_pi0_syst_err_up_panel21", "cosine_1_exps_pi0_syst_err_down_panel21", "cosine_1_exps_km15_panel21", "cosine_1_exps_vgg_panel21", "cosine_1_exps_bh_panel21", "cosine_1_exps_global1_panel21", "cosine_2_exps_pi0_panel21", "cosine_2_exps_pi0_stat_err_panel21", "cosine_2_exps_pi0_syst_err_up_panel21", "cosine_2_exps_pi0_syst_err_down_panel21", "cosine_2_exps_km15_panel21", "cosine_2_exps_vgg_panel21", "cosine_2_exps_bh_panel21", "cosine_2_exps_global1_panel21", "cosine_3_exps_pi0_panel21", "cosine_3_exps_pi0_stat_err_panel21", "cosine_3_exps_pi0_syst_err_up_panel21", "cosine_3_exps_pi0_syst_err_down_panel21", "cosine_3_exps_km15_panel21", "cosine_3_exps_vgg_panel21", "cosine_3_exps_bh_panel21", "cosine_3_exps_global1_panel21", "Q2s_theory_panel21", "cosine_1_ths_km15_panel21", "cosine_1_ths_km15_err_panel21", "cosine_1_ths_vgg_panel21", "cosine_1_ths_vgg_err_panel21", "cosine_1_ths_bh_panel21", "cosine_1_ths_bh_err_panel21", "cosine_2_ths_km15_panel21", "cosine_2_ths_km15_err_panel21", "cosine_2_ths_vgg_panel21", "cosine_2_ths_vgg_err_panel21", "cosine_2_ths_bh_panel21", "cosine_2_ths_bh_err_panel21", "cosine_3_ths_km15_panel21", "cosine_3_ths_km15_err_panel21", "cosine_3_ths_vgg_panel21", "cosine_3_ths_vgg_err_panel21", "cosine_3_ths_bh_panel21", "cosine_3_ths_bh_err_panel21"]:
#         fitting_result[key] = np.load(f)
# ts_theory_panel00 = fitting_result["ts_theory_panel00"]
# cosine_1_ths_km15_panel00 = fitting_result["cosine_1_ths_km15_panel00"]
# cosine_1_ths_km15_err_panel00 = fitting_result["cosine_1_ths_km15_err_panel00"]
# cosine_1_ths_vgg_panel00 = fitting_result["cosine_1_ths_vgg_panel00"]
# cosine_1_ths_vgg_err_panel00 = fitting_result["cosine_1_ths_vgg_err_panel00"]
# cosine_1_ths_bh_panel00 = fitting_result["cosine_1_ths_bh_panel00"]
# cosine_1_ths_bh_err_panel00 = fitting_result["cosine_1_ths_bh_err_panel00"]
# cosine_2_ths_km15_panel00 = fitting_result["cosine_2_ths_km15_panel00"]
# cosine_2_ths_km15_err_panel00 = fitting_result["cosine_2_ths_km15_err_panel00"]
# cosine_2_ths_vgg_panel00 = fitting_result["cosine_2_ths_vgg_panel00"]
# cosine_2_ths_vgg_err_panel00 = fitting_result["cosine_2_ths_vgg_err_panel00"]
# cosine_2_ths_bh_panel00 = fitting_result["cosine_2_ths_bh_panel00"]
# cosine_2_ths_bh_err_panel00 = fitting_result["cosine_2_ths_bh_err_panel00"]
# cosine_3_ths_km15_panel00 = fitting_result["cosine_3_ths_km15_panel00"]
# cosine_3_ths_km15_err_panel00 = fitting_result["cosine_3_ths_km15_err_panel00"]
# cosine_3_ths_vgg_panel00 = fitting_result["cosine_3_ths_vgg_panel00"]
# cosine_3_ths_vgg_err_panel00 = fitting_result["cosine_3_ths_vgg_err_panel00"]
# cosine_3_ths_bh_panel00 = fitting_result["cosine_3_ths_bh_panel00"]
# cosine_3_ths_bh_err_panel00 = fitting_result["cosine_3_ths_bh_err_panel00"]
# ts_theory_panel10 = fitting_result["ts_theory_panel10"]
# cosine_1_ths_km15_panel10 = fitting_result["cosine_1_ths_km15_panel10"]
# cosine_1_ths_km15_err_panel10 = fitting_result["cosine_1_ths_km15_err_panel10"]
# cosine_1_ths_vgg_panel10 = fitting_result["cosine_1_ths_vgg_panel10"]
# cosine_1_ths_vgg_err_panel10 = fitting_result["cosine_1_ths_vgg_err_panel10"]
# cosine_1_ths_bh_panel10 = fitting_result["cosine_1_ths_bh_panel10"]
# cosine_1_ths_bh_err_panel10 = fitting_result["cosine_1_ths_bh_err_panel10"]
# cosine_2_ths_km15_panel10 = fitting_result["cosine_2_ths_km15_panel10"]
# cosine_2_ths_km15_err_panel10 = fitting_result["cosine_2_ths_km15_err_panel10"]
# cosine_2_ths_vgg_panel10 = fitting_result["cosine_2_ths_vgg_panel10"]
# cosine_2_ths_vgg_err_panel10 = fitting_result["cosine_2_ths_vgg_err_panel10"]
# cosine_2_ths_bh_panel10 = fitting_result["cosine_2_ths_bh_panel10"]
# cosine_2_ths_bh_err_panel10 = fitting_result["cosine_2_ths_bh_err_panel10"]
# cosine_3_ths_km15_panel10 = fitting_result["cosine_3_ths_km15_panel10"]
# cosine_3_ths_km15_err_panel10 = fitting_result["cosine_3_ths_km15_err_panel10"]
# cosine_3_ths_vgg_panel10 = fitting_result["cosine_3_ths_vgg_panel10"]
# cosine_3_ths_vgg_err_panel10 = fitting_result["cosine_3_ths_vgg_err_panel10"]
# cosine_3_ths_bh_panel10 = fitting_result["cosine_3_ths_bh_panel10"]
# cosine_3_ths_bh_err_panel10 = fitting_result["cosine_3_ths_bh_err_panel10"]
# ts_theory_panel20 = fitting_result["ts_theory_panel20"]
# cosine_1_ths_km15_panel20 = fitting_result["cosine_1_ths_km15_panel20"]
# cosine_1_ths_km15_err_panel20 = fitting_result["cosine_1_ths_km15_err_panel20"]
# cosine_1_ths_vgg_panel20 = fitting_result["cosine_1_ths_vgg_panel20"]
# cosine_1_ths_vgg_err_panel20 = fitting_result["cosine_1_ths_vgg_err_panel20"]
# cosine_1_ths_bh_panel20 = fitting_result["cosine_1_ths_bh_panel20"]
# cosine_1_ths_bh_err_panel20 = fitting_result["cosine_1_ths_bh_err_panel20"]
# cosine_2_ths_km15_panel20 = fitting_result["cosine_2_ths_km15_panel20"]
# cosine_2_ths_km15_err_panel20 = fitting_result["cosine_2_ths_km15_err_panel20"]
# cosine_2_ths_vgg_panel20 = fitting_result["cosine_2_ths_vgg_panel20"]
# cosine_2_ths_vgg_err_panel20 = fitting_result["cosine_2_ths_vgg_err_panel20"]
# cosine_2_ths_bh_panel20 = fitting_result["cosine_2_ths_bh_panel20"]
# cosine_2_ths_bh_err_panel20 = fitting_result["cosine_2_ths_bh_err_panel20"]
# cosine_3_ths_km15_panel20 = fitting_result["cosine_3_ths_km15_panel20"]
# cosine_3_ths_km15_err_panel20 = fitting_result["cosine_3_ths_km15_err_panel20"]
# cosine_3_ths_vgg_panel20 = fitting_result["cosine_3_ths_vgg_panel20"]
# cosine_3_ths_vgg_err_panel20 = fitting_result["cosine_3_ths_vgg_err_panel20"]
# cosine_3_ths_bh_panel20 = fitting_result["cosine_3_ths_bh_panel20"]
# cosine_3_ths_bh_err_panel20 = fitting_result["cosine_3_ths_bh_err_panel20"]
# Q2s_theory_panel01 = fitting_result["Q2s_theory_panel01"]
# cosine_1_ths_km15_panel01 = fitting_result["cosine_1_ths_km15_panel01"]
# cosine_1_ths_km15_err_panel01 = fitting_result["cosine_1_ths_km15_err_panel01"]
# cosine_1_ths_vgg_panel01 = fitting_result["cosine_1_ths_vgg_panel01"]
# cosine_1_ths_vgg_err_panel01 = fitting_result["cosine_1_ths_vgg_err_panel01"]
# cosine_1_ths_bh_panel01 = fitting_result["cosine_1_ths_bh_panel01"]
# cosine_1_ths_bh_err_panel01 = fitting_result["cosine_1_ths_bh_err_panel01"]
# cosine_2_ths_km15_panel01 = fitting_result["cosine_2_ths_km15_panel01"]
# cosine_2_ths_km15_err_panel01 = fitting_result["cosine_2_ths_km15_err_panel01"]
# cosine_2_ths_vgg_panel01 = fitting_result["cosine_2_ths_vgg_panel01"]
# cosine_2_ths_vgg_err_panel01 = fitting_result["cosine_2_ths_vgg_err_panel01"]
# cosine_2_ths_bh_panel01 = fitting_result["cosine_2_ths_bh_panel01"]
# cosine_2_ths_bh_err_panel01 = fitting_result["cosine_2_ths_bh_err_panel01"]
# cosine_3_ths_km15_panel01 = fitting_result["cosine_3_ths_km15_panel01"]
# cosine_3_ths_km15_err_panel01 = fitting_result["cosine_3_ths_km15_err_panel01"]
# cosine_3_ths_vgg_panel01 = fitting_result["cosine_3_ths_vgg_panel01"]
# cosine_3_ths_vgg_err_panel01 = fitting_result["cosine_3_ths_vgg_err_panel01"]
# cosine_3_ths_bh_panel01 = fitting_result["cosine_3_ths_bh_panel01"]
# cosine_3_ths_bh_err_panel01 = fitting_result["cosine_3_ths_bh_err_panel01"]
# Q2s_theory_panel11 = fitting_result["Q2s_theory_panel11"]
# cosine_1_ths_km15_panel11 = fitting_result["cosine_1_ths_km15_panel11"]
# cosine_1_ths_km15_err_panel11 = fitting_result["cosine_1_ths_km15_err_panel11"]
# cosine_1_ths_vgg_panel11 = fitting_result["cosine_1_ths_vgg_panel11"]
# cosine_1_ths_vgg_err_panel11 = fitting_result["cosine_1_ths_vgg_err_panel11"]
# cosine_1_ths_bh_panel11 = fitting_result["cosine_1_ths_bh_panel11"]
# cosine_1_ths_bh_err_panel11 = fitting_result["cosine_1_ths_bh_err_panel11"]
# cosine_2_ths_km15_panel11 = fitting_result["cosine_2_ths_km15_panel11"]
# cosine_2_ths_km15_err_panel11 = fitting_result["cosine_2_ths_km15_err_panel11"]
# cosine_2_ths_vgg_panel11 = fitting_result["cosine_2_ths_vgg_panel11"]
# cosine_2_ths_vgg_err_panel11 = fitting_result["cosine_2_ths_vgg_err_panel11"]
# cosine_2_ths_bh_panel11 = fitting_result["cosine_2_ths_bh_panel11"]
# cosine_2_ths_bh_err_panel11 = fitting_result["cosine_2_ths_bh_err_panel11"]
# cosine_3_ths_km15_panel11 = fitting_result["cosine_3_ths_km15_panel11"]
# cosine_3_ths_km15_err_panel11 = fitting_result["cosine_3_ths_km15_err_panel11"]
# cosine_3_ths_vgg_panel11 = fitting_result["cosine_3_ths_vgg_panel11"]
# cosine_3_ths_vgg_err_panel11 = fitting_result["cosine_3_ths_vgg_err_panel11"]
# cosine_3_ths_bh_panel11 = fitting_result["cosine_3_ths_bh_panel11"]
# cosine_3_ths_bh_err_panel11 = fitting_result["cosine_3_ths_bh_err_panel11"]
# Q2s_theory_panel21 = fitting_result["Q2s_theory_panel21"]
# cosine_1_ths_km15_panel21 = fitting_result["cosine_1_ths_km15_panel21"]
# cosine_1_ths_km15_err_panel21 = fitting_result["cosine_1_ths_km15_err_panel21"]
# cosine_1_ths_vgg_panel21 = fitting_result["cosine_1_ths_vgg_panel21"]
# cosine_1_ths_vgg_err_panel21 = fitting_result["cosine_1_ths_vgg_err_panel21"]
# cosine_1_ths_bh_panel21 = fitting_result["cosine_1_ths_bh_panel21"]
# cosine_1_ths_bh_err_panel21 = fitting_result["cosine_1_ths_bh_err_panel21"]
# cosine_2_ths_km15_panel21 = fitting_result["cosine_2_ths_km15_panel21"]
# cosine_2_ths_km15_err_panel21 = fitting_result["cosine_2_ths_km15_err_panel21"]
# cosine_2_ths_vgg_panel21 = fitting_result["cosine_2_ths_vgg_panel21"]
# cosine_2_ths_vgg_err_panel21 = fitting_result["cosine_2_ths_vgg_err_panel21"]
# cosine_2_ths_bh_panel21 = fitting_result["cosine_2_ths_bh_panel21"]
# cosine_2_ths_bh_err_panel21 = fitting_result["cosine_2_ths_bh_err_panel21"]
# cosine_3_ths_km15_panel21 = fitting_result["cosine_3_ths_km15_panel21"]
# cosine_3_ths_km15_err_panel21 = fitting_result["cosine_3_ths_km15_err_panel21"]
# cosine_3_ths_vgg_panel21 = fitting_result["cosine_3_ths_vgg_panel21"]
# cosine_3_ths_vgg_err_panel21 = fitting_result["cosine_3_ths_vgg_err_panel21"]
# cosine_3_ths_bh_panel21 = fitting_result["cosine_3_ths_bh_panel21"]
# cosine_3_ths_bh_err_panel21 = fitting_result["cosine_3_ths_bh_err_panel21"]


# cosine_2_ths_vgg_panel00 = np.array(cosine_2_ths_vgg_panel00)
# uvsp_res =UnivariateSpline(ts_theory_panel00[((ts_theory_panel00>0.49-0.05) & (ts_theory_panel00<0.49))| ((ts_theory_panel00>0.56) & (ts_theory_panel00<0.56+0.05))], cosine_2_ths_vgg_panel00[((ts_theory_panel00>0.49-0.05) & (ts_theory_panel00<0.49))| ((ts_theory_panel00>0.56) & (ts_theory_panel00<0.56+0.05))])
# cosine_2_ths_vgg_panel00_smooth = copy(cosine_2_ths_vgg_panel00)
# cosine_2_ths_vgg_panel00_smooth[((ts_theory_panel00>0.49) & (ts_theory_panel00<0.56))] = uvsp_res(ts_theory_panel00[((ts_theory_panel00>0.49) & (ts_theory_panel00<0.56))])
# cosine_2_ths_vgg_panel00_smooth = savgol_filter(cosine_2_ths_vgg_panel00_smooth, 51, 5)

# cosine_2_ths_bh_panel00 = np.array(cosine_2_ths_bh_panel00)
# uvsp_res =UnivariateSpline(ts_theory_panel00[((ts_theory_panel00>0.49-0.05) & (ts_theory_panel00<0.49))| ((ts_theory_panel00>0.56) & (ts_theory_panel00<0.56+0.05))], cosine_2_ths_bh_panel00[((ts_theory_panel00>0.49-0.05) & (ts_theory_panel00<0.49))| ((ts_theory_panel00>0.56) & (ts_theory_panel00<0.56+0.05))])
# cosine_2_ths_bh_panel00_smooth = copy(cosine_2_ths_bh_panel00)
# cosine_2_ths_bh_panel00_smooth[((ts_theory_panel00>0.49) & (ts_theory_panel00<0.56))] = uvsp_res(ts_theory_panel00[((ts_theory_panel00>0.49) & (ts_theory_panel00<0.56))])
# cosine_2_ths_bh_panel00_smooth = savgol_filter(cosine_2_ths_bh_panel00_smooth, 51, 5)


# cosine_2_ths_vgg_panel10 = np.array(cosine_2_ths_vgg_panel10)
# uvsp_res =UnivariateSpline(ts_theory_panel10[((ts_theory_panel10>0.83-0.05) & (ts_theory_panel10<0.83))| ((ts_theory_panel10>0.91) & (ts_theory_panel10<0.91+0.05))], cosine_2_ths_vgg_panel10[((ts_theory_panel10>0.83-0.05) & (ts_theory_panel10<0.83))| ((ts_theory_panel10>0.91) & (ts_theory_panel10<0.91+0.05))])
# cosine_2_ths_vgg_panel10_smooth = copy(cosine_2_ths_vgg_panel10)
# cosine_2_ths_vgg_panel10_smooth[((ts_theory_panel10>0.83) & (ts_theory_panel10<0.91))] = uvsp_res(ts_theory_panel10[((ts_theory_panel10>0.83) & (ts_theory_panel10<0.91))])
# cosine_2_ths_vgg_panel10_smooth = savgol_filter(cosine_2_ths_vgg_panel10_smooth, 51, 5)

# cosine_2_ths_bh_panel10 = np.array(cosine_2_ths_bh_panel10)
# uvsp_res =UnivariateSpline(ts_theory_panel10[((ts_theory_panel10>0.83-0.05) & (ts_theory_panel10<0.83))| ((ts_theory_panel10>0.91) & (ts_theory_panel10<0.91+0.05))], cosine_2_ths_bh_panel10[((ts_theory_panel10>0.83-0.05) & (ts_theory_panel10<0.83))| ((ts_theory_panel10>0.91) & (ts_theory_panel10<0.91+0.05))])
# cosine_2_ths_bh_panel10_smooth = copy(cosine_2_ths_bh_panel10)
# cosine_2_ths_bh_panel10_smooth[((ts_theory_panel10>0.83) & (ts_theory_panel10<0.91))] = uvsp_res(ts_theory_panel10[((ts_theory_panel10>0.83) & (ts_theory_panel10<0.91))])
# cosine_2_ths_bh_panel10_smooth = savgol_filter(cosine_2_ths_bh_panel10_smooth, 51, 5)


# cosine_2_ths_vgg_panel20 = np.array(cosine_2_ths_vgg_panel20)
# # uvsp_res =UnivariateSpline(ts_theory_panel20[((ts_theory_panel20>0.83-0.05) & (ts_theory_panel20<0.83))| ((ts_theory_panel20>0.91) & (ts_theory_panel20<0.91+0.05))], cosine_2_ths_vgg_panel20[((ts_theory_panel20>0.83-0.05) & (ts_theory_panel20<0.83))| ((ts_theory_panel20>0.91) & (ts_theory_panel20<0.91+0.05))])
# cosine_2_ths_vgg_panel20_smooth = copy(cosine_2_ths_vgg_panel20)
# # cosine_2_ths_vgg_panel20_smooth[((ts_theory_panel20>0.83) & (ts_theory_panel20<0.91))] = uvsp_res(ts_theory_panel20[((ts_theory_panel20>0.83) & (ts_theory_panel20<0.91))])
# cosine_2_ths_vgg_panel20_smooth = savgol_filter(cosine_2_ths_vgg_panel20_smooth, 51, 5)

# cosine_2_ths_bh_panel20 = np.array(cosine_2_ths_bh_panel20)
# # uvsp_res =UnivariateSpline(ts_theory_panel20[((ts_theory_panel20>0.83-0.05) & (ts_theory_panel20<0.83))| ((ts_theory_panel20>0.91) & (ts_theory_panel20<0.91+0.05))], cosine_2_ths_bh_panel20[((ts_theory_panel20>0.83-0.05) & (ts_theory_panel20<0.83))| ((ts_theory_panel20>0.91) & (ts_theory_panel20<0.91+0.05))])
# cosine_2_ths_bh_panel20_smooth = copy(cosine_2_ths_bh_panel20)
# # cosine_2_ths_bh_panel20_smooth[((ts_theory_panel20>0.83) & (ts_theory_panel20<0.91))] = uvsp_res(ts_theory_panel20[((ts_theory_panel20>0.83) & (ts_theory_panel20<0.91))])
# cosine_2_ths_bh_panel20_smooth = savgol_filter(cosine_2_ths_bh_panel20_smooth, 51, 5)

# cosine_2_ths_vgg_panel01 = np.array(cosine_2_ths_vgg_panel01)
# uvsp_res =UnivariateSpline(Q2s_theory_panel01[((Q2s_theory_panel01>1.95-0.3) & (Q2s_theory_panel01<1.95))| ((Q2s_theory_panel01>2.03) & (Q2s_theory_panel01<2.03+0.3))], cosine_2_ths_vgg_panel01[((Q2s_theory_panel01>1.95-0.3) & (Q2s_theory_panel01<1.95))| ((Q2s_theory_panel01>2.03) & (Q2s_theory_panel01<2.03+0.3))])
# cosine_2_ths_vgg_panel01_smooth = copy(cosine_2_ths_vgg_panel01)
# cosine_2_ths_vgg_panel01_smooth[((Q2s_theory_panel01>1.95) & (Q2s_theory_panel01<2.03))] = uvsp_res(Q2s_theory_panel01[((Q2s_theory_panel01>1.95) & (Q2s_theory_panel01<2.03))])
# cosine_2_ths_vgg_panel01_smooth = savgol_filter(cosine_2_ths_vgg_panel01_smooth, 51, 5)

# cosine_2_ths_bh_panel01 = np.array(cosine_2_ths_bh_panel01)
# uvsp_res =UnivariateSpline(Q2s_theory_panel01[((Q2s_theory_panel01>1.95-0.3) & (Q2s_theory_panel01<1.95))| ((Q2s_theory_panel01>2.03) & (Q2s_theory_panel01<2.03+0.3))], cosine_2_ths_bh_panel01[((Q2s_theory_panel01>1.95-0.3) & (Q2s_theory_panel01<1.95))| ((Q2s_theory_panel01>2.03) & (Q2s_theory_panel01<2.03+0.3))])
# cosine_2_ths_bh_panel01_smooth = copy(cosine_2_ths_bh_panel01)
# cosine_2_ths_bh_panel01_smooth[((Q2s_theory_panel01>1.95) & (Q2s_theory_panel01<2.03))] = uvsp_res(Q2s_theory_panel01[((Q2s_theory_panel01>1.95) & (Q2s_theory_panel01<2.03))])
# cosine_2_ths_bh_panel01_smooth = savgol_filter(cosine_2_ths_bh_panel01_smooth, 51, 5)


# cosine_2_ths_vgg_panel11 = np.array(cosine_2_ths_vgg_panel11)
# uvsp_res =UnivariateSpline(Q2s_theory_panel11[((Q2s_theory_panel11>2.89-0.3) & (Q2s_theory_panel11<2.89))| ((Q2s_theory_panel11>2.97) & (Q2s_theory_panel11<2.97+0.3))], cosine_2_ths_vgg_panel11[((Q2s_theory_panel11>2.89-0.3) & (Q2s_theory_panel11<2.89))| ((Q2s_theory_panel11>2.97) & (Q2s_theory_panel11<2.97+0.3))])
# cosine_2_ths_vgg_panel11_smooth = copy(cosine_2_ths_vgg_panel11)
# cosine_2_ths_vgg_panel11_smooth[((Q2s_theory_panel11>2.89) & (Q2s_theory_panel11<2.97))] = uvsp_res(Q2s_theory_panel11[((Q2s_theory_panel11>2.89) & (Q2s_theory_panel11<2.97))])
# cosine_2_ths_vgg_panel11_smooth = savgol_filter(cosine_2_ths_vgg_panel11_smooth, 51, 5)

# cosine_2_ths_bh_panel11 = np.array(cosine_2_ths_bh_panel11)
# uvsp_res =UnivariateSpline(Q2s_theory_panel11[((Q2s_theory_panel11>2.89-0.3) & (Q2s_theory_panel11<2.89))| ((Q2s_theory_panel11>2.97) & (Q2s_theory_panel11<2.97+0.3))], cosine_2_ths_bh_panel11[((Q2s_theory_panel11>2.89-0.3) & (Q2s_theory_panel11<2.89))| ((Q2s_theory_panel11>2.97) & (Q2s_theory_panel11<2.97+0.3))])
# cosine_2_ths_bh_panel11_smooth = copy(cosine_2_ths_bh_panel11)
# cosine_2_ths_bh_panel11_smooth[((Q2s_theory_panel11>2.89) & (Q2s_theory_panel11<2.97))] = uvsp_res(Q2s_theory_panel11[((Q2s_theory_panel11>2.89) & (Q2s_theory_panel11<2.97))])
# cosine_2_ths_bh_panel11_smooth = savgol_filter(cosine_2_ths_bh_panel11_smooth, 51, 5)

# cosine_2_ths_vgg_panel21 = np.array(cosine_2_ths_vgg_panel21)
# # uvsp_res =UnivariateSpline(Q2s_theory_panel21[((Q2s_theory_panel21>0.49-0.3) & (Q2s_theory_panel21<0.49))| ((Q2s_theory_panel21>0.56) & (Q2s_theory_panel21<0.56+0.3))], cosine_2_ths_vgg_panel21[((Q2s_theory_panel21>0.49-0.3) & (Q2s_theory_panel21<0.49))| ((Q2s_theory_panel21>0.56) & (Q2s_theory_panel21<0.56+0.3))])
# cosine_2_ths_vgg_panel21_smooth = copy(cosine_2_ths_vgg_panel21)
# # cosine_2_ths_vgg_panel21_smooth[((Q2s_theory_panel21>0.49) & (Q2s_theory_panel21<0.56))] = uvsp_res(Q2s_theory_panel21[((Q2s_theory_panel21>0.49) & (Q2s_theory_panel21<0.56))])
# cosine_2_ths_vgg_panel21_smooth = savgol_filter(cosine_2_ths_vgg_panel21_smooth, 51, 5)

# cosine_2_ths_bh_panel21 = np.array(cosine_2_ths_bh_panel21)
# # uvsp_res =UnivariateSpline(Q2s_theory_panel21[((Q2s_theory_panel21>0.49-0.3) & (Q2s_theory_panel21<0.49))| ((Q2s_theory_panel21>0.56) & (Q2s_theory_panel21<0.56+0.3))], cosine_2_ths_bh_panel21[((Q2s_theory_panel21>0.49-0.3) & (Q2s_theory_panel21<0.49))| ((Q2s_theory_panel21>0.56) & (Q2s_theory_panel21<0.56+0.3))])
# cosine_2_ths_bh_panel21_smooth = copy(cosine_2_ths_bh_panel21)
# # cosine_2_ths_bh_panel21_smooth[((Q2s_theory_panel21>0.49) & (Q2s_theory_panel21<0.56))] = uvsp_res(Q2s_theory_panel21[((Q2s_theory_panel21>0.49) & (Q2s_theory_panel21<0.56))])
# cosine_2_ths_bh_panel21_smooth = savgol_filter(cosine_2_ths_bh_panel21_smooth, 51, 5)


# n_theory = 500
# n_sample = 500

# print ("00", end = " ")
# xBbins_panel00, Q2bins_panel00, tbins_panel00, tmin_panel00, tmax_panel00 = 2, 3, [1, 2, 3, 4, 5, 6], 0.15, 1.0
# ts_panel00, xB_avgs_panel00, Q2_avgs_panel00, cosine_1_exps_pi0_panel00, cosine_1_exps_pi0_stat_err_panel00, cosine_1_exps_pi0_syst_err_up_panel00, cosine_1_exps_pi0_syst_err_down_panel00, cosine_1_exps_bkgmerging_only_panel00, cosine_2_exps_pi0_panel00, cosine_2_exps_pi0_stat_err_panel00, cosine_2_exps_pi0_syst_err_up_panel00, cosine_2_exps_pi0_syst_err_down_panel00, cosine_2_exps_bkgmerging_only_panel00, cosine_3_exps_pi0_panel00, cosine_3_exps_pi0_stat_err_panel00, cosine_3_exps_pi0_syst_err_up_panel00, cosine_3_exps_pi0_syst_err_down_panel00,cosine_3_exps_bkgmerging_only_panel00 =  get_tdependence_exp(xBbins_panel00, Q2bins_panel00, tbins_panel00)

# print ("10", end = " ")
# xBbins_panel10, Q2bins_panel10, tbins_panel10, tmin_panel10, tmax_panel10 = 3, 3, [1, 2, 3, 4, 5, 6], 0.15, 1.0
# ts_panel10, xB_avgs_panel10, Q2_avgs_panel10, cosine_1_exps_pi0_panel10, cosine_1_exps_pi0_stat_err_panel10, cosine_1_exps_pi0_syst_err_up_panel10, cosine_1_exps_pi0_syst_err_down_panel10, cosine_1_exps_bkgmerging_only_panel10, cosine_2_exps_pi0_panel10, cosine_2_exps_pi0_stat_err_panel10, cosine_2_exps_pi0_syst_err_up_panel10, cosine_2_exps_pi0_syst_err_down_panel10, cosine_2_exps_bkgmerging_only_panel10, cosine_3_exps_pi0_panel10, cosine_3_exps_pi0_stat_err_panel10, cosine_3_exps_pi0_syst_err_up_panel10, cosine_3_exps_pi0_syst_err_down_panel10,cosine_3_exps_bkgmerging_only_panel10 =  get_tdependence_exp(xBbins_panel10, Q2bins_panel10, tbins_panel10)

# print ("20", end = " ")
# xBbins_panel20, Q2bins_panel20, tbins_panel20, tmin_panel20, tmax_panel20 = 4, 3, [1, 2, 3, 4, 5, 6], 0.15, 1.0
# ts_panel20, xB_avgs_panel20, Q2_avgs_panel20, cosine_1_exps_pi0_panel20, cosine_1_exps_pi0_stat_err_panel20, cosine_1_exps_pi0_syst_err_up_panel20, cosine_1_exps_pi0_syst_err_down_panel20, cosine_1_exps_bkgmerging_only_panel20, cosine_2_exps_pi0_panel20, cosine_2_exps_pi0_stat_err_panel20, cosine_2_exps_pi0_syst_err_up_panel20, cosine_2_exps_pi0_syst_err_down_panel20, cosine_2_exps_bkgmerging_only_panel20, cosine_3_exps_pi0_panel20, cosine_3_exps_pi0_stat_err_panel20, cosine_3_exps_pi0_syst_err_up_panel20, cosine_3_exps_pi0_syst_err_down_panel20,cosine_3_exps_bkgmerging_only_panel20 =  get_tdependence_exp(xBbins_panel20, Q2bins_panel20, tbins_panel20)

# print ("01", end = " ")
# xBbins_panel01, Q2bins_panel01, tbins_panel01, Q2min_panel01, Q2max_panel01 = 2, [0, 1, 2, 3], 3, 1, 2.3
# Q2s_panel01, xB_avgs_panel01, t_avgs_panel01, cosine_1_exps_pi0_panel01, cosine_1_exps_pi0_stat_err_panel01, cosine_1_exps_pi0_syst_err_up_panel01, cosine_1_exps_pi0_syst_err_down_panel01, cosine_1_exps_bkgmerging_only_panel01, cosine_2_exps_pi0_panel01, cosine_2_exps_pi0_stat_err_panel01, cosine_2_exps_pi0_syst_err_up_panel01, cosine_2_exps_pi0_syst_err_down_panel01, cosine_2_exps_bkgmerging_only_panel01, cosine_3_exps_pi0_panel01, cosine_3_exps_pi0_stat_err_panel01, cosine_3_exps_pi0_syst_err_up_panel01, cosine_3_exps_pi0_syst_err_down_panel01,cosine_3_exps_bkgmerging_only_panel01 =  get_Q2dependence_exp(xBbins_panel01, Q2bins_panel01, tbins_panel01)

# print ("11", end = " ")
# xBbins_panel11, Q2bins_panel11, tbins_panel11, Q2min_panel11, Q2max_panel11 = 3, [0, 1, 2, 3, 4, 5, 6], 3, 1, 3
# Q2s_panel11, xB_avgs_panel11, t_avgs_panel11, cosine_1_exps_pi0_panel11, cosine_1_exps_pi0_stat_err_panel11, cosine_1_exps_pi0_syst_err_up_panel11, cosine_1_exps_pi0_syst_err_down_panel11, cosine_1_exps_bkgmerging_only_panel11, cosine_2_exps_pi0_panel11, cosine_2_exps_pi0_stat_err_panel11, cosine_2_exps_pi0_syst_err_up_panel11, cosine_2_exps_pi0_syst_err_down_panel11, cosine_2_exps_bkgmerging_only_panel11, cosine_3_exps_pi0_panel11, cosine_3_exps_pi0_stat_err_panel11, cosine_3_exps_pi0_syst_err_up_panel11, cosine_3_exps_pi0_syst_err_down_panel11,cosine_3_exps_bkgmerging_only_panel11 =  get_Q2dependence_exp(xBbins_panel11, Q2bins_panel11, tbins_panel11)

# print ("21", end = " ")
# xBbins_panel21, Q2bins_panel21, tbins_panel21, Q2min_panel21, Q2max_panel21 = 4, [0, 1, 2, 3, 4, 5, 6], 3, 1, 4
# Q2s_panel21, xB_avgs_panel21, t_avgs_panel21, cosine_1_exps_pi0_panel21, cosine_1_exps_pi0_stat_err_panel21, cosine_1_exps_pi0_syst_err_up_panel21, cosine_1_exps_pi0_syst_err_down_panel21, cosine_1_exps_bkgmerging_only_panel21, cosine_2_exps_pi0_panel21, cosine_2_exps_pi0_stat_err_panel21, cosine_2_exps_pi0_syst_err_up_panel21, cosine_2_exps_pi0_syst_err_down_panel21, cosine_2_exps_bkgmerging_only_panel21, cosine_3_exps_pi0_panel21, cosine_3_exps_pi0_stat_err_panel21, cosine_3_exps_pi0_syst_err_up_panel21, cosine_3_exps_pi0_syst_err_down_panel21,cosine_3_exps_bkgmerging_only_panel21 =  get_Q2dependence_exp(xBbins_panel21, Q2bins_panel21, tbins_panel21)

# integrated_binnums_modified = []
# for tbin in tbins_panel00:
#     if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel00) & (df_summary_table_rebinned.Q2bin == Q2bins_panel00) &  (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"]) > 4:
#         integrated_binnums_modified.append(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel00) & (df_summary_table_rebinned.Q2bin == Q2bins_panel00) &  (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"].unique())

# for tbin in tbins_panel10:
#     if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel10) & (df_summary_table_rebinned.Q2bin == Q2bins_panel10) &  (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"]) > 4:
#         integrated_binnums_modified.append(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel10) & (df_summary_table_rebinned.Q2bin == Q2bins_panel10) &  (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"].unique())

# for tbin in tbins_panel20:
#     if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel20) & (df_summary_table_rebinned.Q2bin == Q2bins_panel20) &  (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"]) > 4:
#         integrated_binnums_modified.append(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel20) & (df_summary_table_rebinned.Q2bin == Q2bins_panel20) &  (df_summary_table_rebinned.tbin == tbin) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"].unique())

# for Q2bin in Q2bins_panel01:
#     if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel01) & (df_summary_table_rebinned.Q2bin == Q2bin) &  (df_summary_table_rebinned.tbin == tbins_panel01) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"]) > 4:
#         integrated_binnums_modified.append(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel01) & (df_summary_table_rebinned.Q2bin == Q2bin) &  (df_summary_table_rebinned.tbin == tbins_panel01) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"].unique())

# for Q2bin in Q2bins_panel11:
#     if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel11) & (df_summary_table_rebinned.Q2bin == Q2bin) &  (df_summary_table_rebinned.tbin == tbins_panel11) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"]) > 4:
#         integrated_binnums_modified.append(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel11) & (df_summary_table_rebinned.Q2bin == Q2bin) &  (df_summary_table_rebinned.tbin == tbins_panel11) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"].unique())

# for Q2bin in Q2bins_panel21:
#     if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel21) & (df_summary_table_rebinned.Q2bin == Q2bin) &  (df_summary_table_rebinned.tbin == tbins_panel21) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"]) > 4:
#         integrated_binnums_modified.append(df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xBbins_panel21) & (df_summary_table_rebinned.Q2bin == Q2bin) &  (df_summary_table_rebinned.tbin == tbins_panel21) & (df_summary_table_rebinned.active_bin_nominal == 1), "integrated_binnum"].unique())

# integrated_binnums_modified = np.array(integrated_binnums_modified).flatten()


# for integrated_binnum in integrated_binnums_modified:#integrated_binnum = 88
#     fig, axs = plt.subplots(1, 1, figsize = (10, 6))
#     df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal == 1), :]
#     axs.errorbar(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, yerr = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, ls ='', marker = 'o', color = 'k', label = r'$\mathrm{With~normalization}$', zorder = 5, mfc ='None')
#     axs.errorbar(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_bkg_merging_w, yerr = df_this_bin.xsec_exp_bkg_merging_stat_err_w, ls ='', marker = 'o', color = 'tab:pink', label = r'$\mathrm{Without~normalization}$')
#     axs.fill_between(df_this_bin.phi_avg_this_point, df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w + df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_w, color = 'k', alpha = 0.5)
#     popt, pcov = curve_fit(cosine_fitting_2, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w, absolute_sigma = True,  p0 = (1, 0, 0))
#     print(popt)
#     chi2fit = np.sum((df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - cosine_fitting_2(np.radians(df_this_bin.phi_avg_this_point), *popt))**2/(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_w)**2)
#     ndf = len(df_this_bin) - 3
#     pvalue = (1-chi2.cdf(chi2fit, ndf))
#     axs.plot(phi_dummy, cosine_fitting_2(np.radians(phi_dummy), *popt), color = 'k', ls = ':', label = r"$\mathrm{Measurement~(This~Work)}$" + "\n" + "${}{}\cos\phi_{{\mathrm{{BMK}}}}{}\cos2\phi_{{\mathrm{{BMK}}}}$".format(engineering_to_latex2("{:.2e}".format(popt[0])), engineering_to_latex2("{:+.2e}".format(-popt[1])), engineering_to_latex2("{:+.2e}".format(popt[2]))) + "\n" + "$\chi^2 ={:.1f}, ndf = {:.0f}, p = {:.3f}$".format(chi2fit, ndf, pvalue) )
#     # popt, pcov = curve_fit(cosine_fitting, np.radians(df_this_bin.phi_avg_this_point), df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w, sigma = df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_err_w, absolute_sigma = True,  p0 = (1, 0))
#     # print(popt)
#     # chi2fit = np.sum((df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_w - cosine_fitting(np.radians(df_this_bin.phi_avg_this_point), *popt))**2/(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_err_w)**2)
#     # ndf = len(df_this_bin) - 3
#     # pvalue = (1-chi2.cdf(chi2fit, ndf))
#     # axs.plot(phi_dummy, cosine_fitting(np.radians(phi_dummy), *popt), color = 'k', ls = '--', label = "${}{}\cos\phi_{{\mathrm{{BMK}}}}$".format(engineering_to_latex2("{:.2e}".format(popt[0])), engineering_to_latex2("{:+.2e}".format(-popt[1]))) + "\n" + "$\chi^2 ={:.1f}, ndf = {:.0f}, p = {:.3f}$".format(chi2fit, ndf, pvalue) )

#     df_display_this_bin = df_display.loc[(df_display.integratedbin_display == integrated_binnum), :]


#     # axs.errorbar(np.degrees(df_display_this_bin.phi_display), df_display_this_bin.xsec_KM15_display_w, ls ='', marker = 'o', color = 'cyan')
#     popt, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_KM15_display_w, p0 = (1, 0, 0))
#     print(popt)

#     axs.plot(np.degrees(df_display_this_bin.phi_display), cosine_fitting_2(df_display_this_bin.phi_display, *popt), color = 'cyan', ls = '--', label = r"$\mathrm{Theory~(KM15)}$" + "\n" + "${}{}\cos\phi_{{\mathrm{{BMK}}}}{}\cos2\phi_{{\mathrm{{BMK}}}}$".format(engineering_to_latex2("{:.2e}".format(popt[0])), engineering_to_latex2("{:+.2e}".format(-popt[1])), engineering_to_latex2("{:+.2e}".format(popt[2]))))
#     # axs.errorbar(np.degrees(df_display_this_bin.phi_display), df_display_this_bin.xsec_VGG_display_w, ls ='', marker = 'o', color = 'tab:orange')
#     popt, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_VGG_display_w, p0 = (1, 0, 0))
#     print(popt)

#     axs.plot(np.degrees(df_display_this_bin.phi_display), cosine_fitting_2(df_display_this_bin.phi_display, *popt), color = 'tab:orange', ls = '--', label = r"$\mathrm{Theory~(VGG)}$" + "\n" + "${}{}\cos\phi_{{\mathrm{{BMK}}}}{}\cos2\phi_{{\mathrm{{BMK}}}}$".format(engineering_to_latex2("{:.2e}".format(popt[0])), engineering_to_latex2("{:+.2e}".format(-popt[1])), engineering_to_latex2("{:+.2e}".format(popt[2]))))
#     # axs.errorbar(np.degrees(df_display_this_bin.phi_display), df_display_this_bin.xsec_BH_display_w, ls ='', marker = 'o', color = 'tab:red')
#     popt, pcov = curve_fit(cosine_fitting_2, df_display_this_bin.phi_display, df_display_this_bin.xsec_BH_display_w, p0 = (1, 0, 0))
#     print(popt)

#     axs.plot(np.degrees(df_display_this_bin.phi_display), cosine_fitting_2(df_display_this_bin.phi_display, *popt), color = 'tab:red', ls = '--', label = r"$\mathrm{Theory~(BH)}$" + "\n" + "${}{}\cos\phi_{{\mathrm{{BMK}}}}{}\cos2\phi_{{\mathrm{{BMK}}}}$".format(engineering_to_latex2("{:.2e}".format(popt[0])), engineering_to_latex2("{:+.2e}".format(-popt[1])), engineering_to_latex2("{:+.2e}".format(popt[2]))))
#     axs.set_ylabel(r"$2\pi \frac{P_1 P_2}{\int d\phi P_1 P_2} \times \frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~(\mathrm{nb/GeV}^4)$")
#     axs.set_xlabel("$\phi$~($^\circ$)")
#     axs.set_xticks([0, 90, 180, 270, 360])
#     axs.set_xlim([0, 360])

#     xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
#     print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)
#     xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]

#     # xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
#     # Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
#     # t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
#     # plt.figlegend(loc = 'upper left', title = "${}.$\n".format(integrated_binnum)+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
#     xBheader = "$\\langle x_B\\rangle  = {:.3f}$\n".format(xBmean)#, xBmax)
#     Q2header = "$\\langle Q^2\\rangle  = {:.3f}~\mathrm{{GeV}}^2/c^2$\n".format(Q2mean)#, Q2max)
#     t1header = "$\\langle |t|\\rangle~ = {:.3f}~\mathrm{{GeV}}^2$".format(t1mean)#, t1max)
#     handles, labels = axs.get_legend_handles_labels() 
#     orders  = [4, 5, 0, 1, 2, 3]
#     handles = [handles[order] for order in orders]
#     labels  = [ labels[order] for order in orders]
#     plt.figlegend(handles, labels, loc = 'upper left', title = "${}.$\n".format(integrated_binnum)+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 1.0))

#     plt.savefig("addendum_v3/modified_xsec_{}.from_script.pdf".format(integrated_binnum), bbox_inches = 'tight')
#     plt.close()

# label_scheme = ["$\pi^0~\mathrm{Norm.~(Nominal)}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]

# fig, axs = plt.subplots(3, 2, figsize = (10, 8))
# #t dependence_1
# axs[0, 0].errorbar(ts_panel00, cosine_3_exps_pi0_panel00, yerr = cosine_3_exps_pi0_stat_err_panel00, marker = 'o', color = 'k', mfc = 'None', zorder = 5, label = "$\pi^0~\mathrm{Norm.~(Nominal)}$")
# # axs[0, 0].scatter(ts_panel00, cosine_3_exps_km15_panel00, marker = 'o', color = 'cyan', zorder = 5, label = r"$\mathrm{KM15~Norm.}$")
# # axs[0, 0].scatter(ts_panel00, cosine_3_exps_vgg_panel00, marker = 'o', color = 'tab:orange', zorder = 5, label = r"$\mathrm{VGG~Norm.}$")
# # axs[0, 0].scatter(ts_panel00, cosine_3_exps_bh_panel00, marker = 'o', color = 'tab:red', zorder = 5, label = r"$\mathrm{BH~Norm.}$")
# # axs[0, 0].scatter(ts_panel00, cosine_3_exps_global1_panel00, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Global~Norm.}$")
# axs[0, 0].scatter(ts_panel00, cosine_3_exps_bkgmerging_only_panel00, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Without~Normalization}$")
# axs[0, 0].fill_between(ts_panel00, cosine_3_exps_pi0_panel00 - cosine_3_exps_pi0_syst_err_down_panel00, cosine_3_exps_pi0_panel00 + cosine_3_exps_pi0_syst_err_up_panel00, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[0, 0].errorbar(ts_theory_panel00, cosine_2_ths_km15_panel00       , color = 'cyan', label = r"$\mathrm{Theory~(KM15)}$")
# axs[0, 0].errorbar(ts_theory_panel00, cosine_2_ths_vgg_panel00_smooth , color = 'tab:orange', label = r"$\mathrm{Theory~(VGG)}$")
# axs[0, 0].errorbar(ts_theory_panel00, cosine_2_ths_bh_panel00_smooth  , color = 'tab:red', label = r"$\mathrm{Theory~(BH)}$")

# axs[1, 0].errorbar(ts_panel10, cosine_3_exps_pi0_panel10, yerr = cosine_3_exps_pi0_stat_err_panel10, marker = 'o', color = 'k', mfc = 'None', zorder = 5, label = "$\pi^0~\mathrm{Norm.~(Nominal)}$")
# # axs[1, 0].scatter(ts_panel10, cosine_3_exps_km15_panel10, marker = 'o', color = 'cyan', zorder = 5, label = r"$\mathrm{KM15~Norm.}$")
# # axs[1, 0].scatter(ts_panel10, cosine_3_exps_vgg_panel10, marker = 'o', color = 'tab:orange', zorder = 5, label = r"$\mathrm{VGG~Norm.}$")
# # axs[1, 0].scatter(ts_panel10, cosine_3_exps_bh_panel10, marker = 'o', color = 'tab:red', zorder = 5, label = r"$\mathrm{BH~Norm.}$")
# # axs[1, 0].scatter(ts_panel10, cosine_3_exps_global1_panel10, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Global~Norm.}$")
# axs[1, 0].scatter(ts_panel10, cosine_3_exps_bkgmerging_only_panel10, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Without~Normalization}$")
# axs[1, 0].fill_between(ts_panel10, cosine_3_exps_pi0_panel10 - cosine_3_exps_pi0_syst_err_down_panel10, cosine_3_exps_pi0_panel10 + cosine_3_exps_pi0_syst_err_up_panel10, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[1, 0].errorbar(ts_theory_panel10, cosine_2_ths_km15_panel10       , color = 'cyan', label = r"$\mathrm{Theory~(KM15)}$")
# axs[1, 0].errorbar(ts_theory_panel10, cosine_2_ths_vgg_panel10_smooth , color = 'tab:orange', label = r"$\mathrm{Theory~(VGG)}$")
# axs[1, 0].errorbar(ts_theory_panel10, cosine_2_ths_bh_panel10_smooth  , color = 'tab:red', label = r"$\mathrm{Theory~(BH)}$")

# axs[2, 0].errorbar(ts_panel20, cosine_3_exps_pi0_panel20, yerr = cosine_3_exps_pi0_stat_err_panel20, marker = 'o', color = 'k', mfc = 'None', zorder = 5, label = "$\pi^0~\mathrm{Norm.~(Nominal)}$")
# # axs[2, 0].scatter(ts_panel20, cosine_3_exps_km15_panel20, marker = 'o', color = 'cyan', zorder = 5, label = r"$\mathrm{KM15~Norm.}$")
# # axs[2, 0].scatter(ts_panel20, cosine_3_exps_vgg_panel20, marker = 'o', color = 'tab:orange', zorder = 5, label = r"$\mathrm{VGG~Norm.}$")
# # axs[2, 0].scatter(ts_panel20, cosine_3_exps_bh_panel20, marker = 'o', color = 'tab:red', zorder = 5, label = r"$\mathrm{BH~Norm.}$")
# # axs[2, 0].scatter(ts_panel20, cosine_3_exps_global1_panel20, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Global~Norm.}$")
# axs[2, 0].scatter(ts_panel20, cosine_3_exps_bkgmerging_only_panel20, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Without~Normalization}$")
# axs[2, 0].fill_between(ts_panel20, cosine_3_exps_pi0_panel20 - cosine_3_exps_pi0_syst_err_down_panel20, cosine_3_exps_pi0_panel20 + cosine_3_exps_pi0_syst_err_up_panel20, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[2, 0].errorbar(ts_theory_panel20, cosine_2_ths_km15_panel20       , color = 'cyan', label = r"$\mathrm{Theory~(KM15)}$")
# axs[2, 0].errorbar(ts_theory_panel20, cosine_2_ths_vgg_panel20_smooth , color = 'tab:orange', label = r"$\mathrm{Theory~(VGG)}$")
# axs[2, 0].errorbar(ts_theory_panel20, cosine_2_ths_bh_panel20_smooth  , color = 'tab:red', label = r"$\mathrm{Theory~(BH)}$")

# axs[0, 1].errorbar(Q2s_panel01, cosine_3_exps_pi0_panel01, yerr = cosine_3_exps_pi0_stat_err_panel01, marker = 'o', color = 'k', mfc = 'None', zorder = 5, label = "$\pi^0~\mathrm{Norm.~(Nominal)}$")
# # axs[0, 1].scatter(Q2s_panel01, cosine_3_exps_km15_panel01, marker = 'o', color = 'cyan', zorder = 5, label = r"$\mathrm{KM15~Norm.}$")
# # axs[0, 1].scatter(Q2s_panel01, cosine_3_exps_vgg_panel01, marker = 'o', color = 'tab:orange', zorder = 5, label = r"$\mathrm{VGG~Norm.}$")
# # axs[0, 1].scatter(Q2s_panel01, cosine_3_exps_bh_panel01, marker = 'o', color = 'tab:red', zorder = 5, label = r"$\mathrm{BH~Norm.}$")
# # axs[0, 1].scatter(Q2s_panel01, cosine_3_exps_global1_panel01, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Global~Norm.}$")
# axs[0, 1].scatter(Q2s_panel01, cosine_3_exps_bkgmerging_only_panel01, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Without~Normalization}$")
# axs[0, 1].fill_between(Q2s_panel01, cosine_3_exps_pi0_panel01 - cosine_3_exps_pi0_syst_err_down_panel01, cosine_3_exps_pi0_panel01 + cosine_3_exps_pi0_syst_err_up_panel01, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[0, 1].errorbar(Q2s_theory_panel01, cosine_2_ths_km15_panel01       , color = 'cyan', label = r"$\mathrm{Theory~(KM15)}$")
# axs[0, 1].errorbar(Q2s_theory_panel01, cosine_2_ths_vgg_panel01_smooth , color = 'tab:orange', label = r"$\mathrm{Theory~(VGG)}$")
# axs[0, 1].errorbar(Q2s_theory_panel01, cosine_2_ths_bh_panel01_smooth  , color = 'tab:red', label = r"$\mathrm{Theory~(BH)}$")

# axs[1, 1].errorbar(Q2s_panel11, cosine_3_exps_pi0_panel11, yerr = cosine_3_exps_pi0_stat_err_panel11, marker = 'o', color = 'k', mfc = 'None', zorder = 5, label = "$\pi^0~\mathrm{Norm.~(Nominal)}$")
# # axs[1, 1].scatter(Q2s_panel11, cosine_3_exps_km15_panel11, marker = 'o', color = 'cyan', zorder = 5, label = r"$\mathrm{KM15~Norm.}$")
# # axs[1, 1].scatter(Q2s_panel11, cosine_3_exps_vgg_panel11, marker = 'o', color = 'tab:orange', zorder = 5, label = r"$\mathrm{VGG~Norm.}$")
# # axs[1, 1].scatter(Q2s_panel11, cosine_3_exps_bh_panel11, marker = 'o', color = 'tab:red', zorder = 5, label = r"$\mathrm{BH~Norm.}$")
# # axs[1, 1].scatter(Q2s_panel11, cosine_3_exps_global1_panel11, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Global~Norm.}$")
# axs[1, 1].scatter(Q2s_panel11, cosine_3_exps_bkgmerging_only_panel11, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Without~Normalization}$")
# axs[1, 1].fill_between(Q2s_panel11, cosine_3_exps_pi0_panel11 - cosine_3_exps_pi0_syst_err_down_panel11, cosine_3_exps_pi0_panel11 + cosine_3_exps_pi0_syst_err_up_panel11, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[1, 1].errorbar(Q2s_theory_panel11, cosine_2_ths_km15_panel11       , color = 'cyan', label = r"$\mathrm{Theory~(KM15)}$")
# axs[1, 1].errorbar(Q2s_theory_panel11, cosine_2_ths_vgg_panel11_smooth , color = 'tab:orange', label = r"$\mathrm{Theory~(VGG)}$")
# axs[1, 1].errorbar(Q2s_theory_panel11, cosine_2_ths_bh_panel11_smooth  , color = 'tab:red', label = r"$\mathrm{Theory~(BH)}$")

# axs[2, 1].errorbar(Q2s_panel21, cosine_3_exps_pi0_panel21, yerr = cosine_3_exps_pi0_stat_err_panel21, marker = 'o', color = 'k', mfc = 'None', zorder = 5, label = "$\pi^0~\mathrm{Norm.~(Nominal)}$")
# # axs[2, 1].scatter(Q2s_panel21, cosine_3_exps_km15_panel21, marker = 'o', color = 'cyan', zorder = 5, label = r"$\mathrm{KM15~Norm.}$")
# # axs[2, 1].scatter(Q2s_panel21, cosine_3_exps_vgg_panel21, marker = 'o', color = 'tab:orange', zorder = 5, label = r"$\mathrm{VGG~Norm.}$")
# # axs[2, 1].scatter(Q2s_panel21, cosine_3_exps_bh_panel21, marker = 'o', color = 'tab:red', zorder = 5, label = r"$\mathrm{BH~Norm.}$")
# # axs[2, 1].scatter(Q2s_panel21, cosine_3_exps_global1_panel21, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Global~Norm.}$")
# axs[2, 1].scatter(Q2s_panel21, cosine_3_exps_bkgmerging_only_panel21, marker = 'o', color = 'tab:pink', zorder = 5, label = r"$\mathrm{Without~Normalization}$")
# axs[2, 1].fill_between(Q2s_panel21, cosine_3_exps_pi0_panel21 - cosine_3_exps_pi0_syst_err_down_panel21, cosine_3_exps_pi0_panel21 + cosine_3_exps_pi0_syst_err_up_panel21, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[2, 1].errorbar(Q2s_theory_panel21, cosine_2_ths_km15_panel21       , color = 'cyan', label = r"$\mathrm{Theory~(KM15)}$")
# axs[2, 1].errorbar(Q2s_theory_panel21, cosine_2_ths_vgg_panel21_smooth , color = 'tab:orange', label = r"$\mathrm{Theory~(VGG)}$")
# axs[2, 1].errorbar(Q2s_theory_panel21, cosine_2_ths_bh_panel21_smooth  , color = 'tab:red', label = r"$\mathrm{Theory~(BH)}$")

# axs[0, 0].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel00)), xycoords = 'axes fraction')
# axs[1, 0].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel10)), xycoords = 'axes fraction')
# axs[2, 0].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel20)), xycoords = 'axes fraction')
# axs[0, 1].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel01)), xycoords = 'axes fraction')
# axs[1, 1].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel11)), xycoords = 'axes fraction')
# axs[2, 1].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel21)), xycoords = 'axes fraction')

# for ax in axs[:, 0]:
#     ax.set_xlim([0.15, 0.95])
#     ax.set_xticks(np.linspace(0.16, 1.0, 43), minor = True)
#     ax.set_xticks(np.linspace(0.2, 1.0, 9), ['$0.2$', '', '$0.4$', '', '$0.6$', '', '$0.8$', '', '$1.0$'])
#     ax.tick_params(right = True, left = True, top = True, bottom = True, direction = 'in', which = 'both')

# for ax in axs[:, 1]:
#     ax.set_xlim([1, 4])
#     ax.set_xticks(np.linspace(1, 4, 31), minor = True)
#     ax.set_xticks(np.linspace(1, 4, 4))
#     ax.tick_params(right = True, left = True, top = True, bottom = True, direction = 'in', which = 'both')

# for ax in axs[:2, 0]:
#     ax.set_xticks(np.linspace(0.2, 1.0, 9), ['']*9)
#     # ax.set_xticks(np.linspace(0.15, 0.95, 81), minor = True)

# for ax in axs[:2, 1]:
#     ax.set_xticks(np.linspace(1, 4, 4), ['']*4)
#     # ax.set_xticks(np.linspace(1, 4, 31), minor = True)

# handles, labels = axs[2, 1].get_legend_handles_labels()
# orders = [2, 1, 0, 3, 4, 5]
# handles = [handles[i] for i in orders]
# labels = [labels[i] for i in orders]
# # handles = [handles[5], handles[4], handles[2], handles[0], handles[1], handles[3], handles[8], handles[6], handles[7]]
# # labels = [labels[5], labels[4], labels[2], labels[0], labels[1], labels[3], labels[8], labels[6], labels[7]]
# plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (.9, 0.75), title = r"$d\sigma_{ep\rightarrow e'p'\gamma}^{\cos \phi,w}~(\mathrm{nb/GeV}^4)$", fontsize = 20, title_fontsize = 20, markerscale = 1.3)

# axs[-1, 0].set_xlabel(r"$|t|~(\mathrm{GeV}^2)$")
# axs[-1, 1].set_xlabel(r"$Q^2~(\mathrm{GeV}^2/c^2)$")
# plt.subplots_adjust(hspace = 0., wspace = 0.23)


# axs[0, 0].set_ylim([-0.06, 0.11])
# axs[0, 0].set_yticks([-0.05, 0.00, 0.05, 0.1])
# axs[0, 0].set_yticks(np.linspace(-0.06, 0.11, 18), minor = True)


# axs[0, 1].set_ylim([-0.03, 0.15])
# axs[0, 1].set_yticks([0, 0.05, 0.1, 0.15])
# axs[0, 1].set_yticks(np.linspace(-0.03, 0.15, 19), minor = True)


# axs[1, 0].set_ylim([-0.03, 0.07])
# axs[1, 0].set_yticks([-0.02, 0, 0.02, 0.04, 0.06])
# axs[1, 0].set_yticks(np.linspace(-0.03, 0.07, 11), minor = True)

# axs[1, 1].set_ylim([-0.01, 0.07])
# axs[1, 1].set_yticks([0, 0.02, 0.04, 0.06])
# axs[1, 1].set_yticks(np.linspace(-0.01, 0.07, 9), minor = True)

# axs[2, 0].set_ylim([-0.01, 0.055])
# axs[2, 0].set_yticks([0, 0.02, 0.04])
# axs[2, 0].set_yticks(np.linspace(-0.01, 0.05, 7), minor = True)

# axs[2, 1].set_ylim([-0.01, 0.07])
# axs[2, 1].set_yticks([0, 0.02, 0.04, 0.06])
# axs[2, 1].set_yticks(np.linspace(-0.01, 0.07, 9), minor = True)

# axs[0, 0].set_title("$\langle Q^2 \\rangle = {:.3f}~\mathrm{{GeV}}/c^2$".format(np.mean([*Q2_avgs_panel00, *Q2_avgs_panel10, *Q2_avgs_panel20])))
# axs[0, 1].set_title("$\langle |t| \\rangle = {:.3f}~\mathrm{{GeV}}^2$".format(np.mean([*t_avgs_panel01, *t_avgs_panel11, *t_avgs_panel21])))
# # plt.savefig("addendum_v3/modified_xsec_entire_w_different_normalizations.from_script.pdf", bbox_inches = 'tight')
# plt.savefig("addendum_v3/modified_xsec_entire_w_different_normalizations.from_script.pdf", bbox_inches = 'tight')

# label_scheme = ["$\pi^0~\mathrm{Norm.~(Nominal)}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]

# fig, axs = plt.subplots(3, 2, figsize = (10, 8))
# #t dependence_1
# axs[0, 0].errorbar(ts_panel00, cosine_3_exps_pi0_panel00, yerr = cosine_3_exps_pi0_stat_err_panel00, marker = 'o', color = 'k', zorder = 5, label = r"$\mathrm{Measurement}$" + "\n"+ r"$\mathrm{(This~work)}$")
# axs[0, 0].fill_between(ts_panel00, cosine_3_exps_pi0_panel00 - cosine_3_exps_pi0_syst_err_down_panel00, cosine_3_exps_pi0_panel00 + cosine_3_exps_pi0_syst_err_up_panel00, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[0, 0].errorbar(ts_theory_panel00, cosine_2_ths_km15_panel00, color = 'cyan', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(KM15)}$")
# axs[0, 0].errorbar(ts_theory_panel00, cosine_2_ths_vgg_panel00_smooth , color = 'tab:orange', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(VGG)}$")
# axs[0, 0].errorbar(ts_theory_panel00, cosine_2_ths_bh_panel00_smooth  , color = 'tab:red', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(BH)}$")

# axs[1, 0].errorbar(ts_panel10, cosine_3_exps_pi0_panel10, yerr = cosine_3_exps_pi0_stat_err_panel10, marker = 'o', color = 'k', zorder = 5, label = r"$\mathrm{Measurement}$" + "\n"+ r"$\mathrm{(This~work)}$")
# axs[1, 0].fill_between(ts_panel10, cosine_3_exps_pi0_panel10 - cosine_3_exps_pi0_syst_err_down_panel10, cosine_3_exps_pi0_panel10 + cosine_3_exps_pi0_syst_err_up_panel10, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[1, 0].errorbar(ts_theory_panel10, cosine_2_ths_km15_panel10, color = 'cyan', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(KM15)}$")
# axs[1, 0].errorbar(ts_theory_panel10, cosine_2_ths_vgg_panel10_smooth , color = 'tab:orange', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(VGG)}$")
# axs[1, 0].errorbar(ts_theory_panel10, cosine_2_ths_bh_panel10_smooth  , color = 'tab:red', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(BH)}$")

# axs[2, 0].errorbar(ts_panel20, cosine_3_exps_pi0_panel20, yerr = cosine_3_exps_pi0_stat_err_panel20, marker = 'o', color = 'k', zorder = 5, label = r"$\mathrm{Measurement}$" + "\n"+ r"$\mathrm{(This~work)}$")
# axs[2, 0].fill_between(ts_panel20, cosine_3_exps_pi0_panel20 - cosine_3_exps_pi0_syst_err_down_panel20, cosine_3_exps_pi0_panel20 + cosine_3_exps_pi0_syst_err_up_panel20, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[2, 0].errorbar(ts_theory_panel20, cosine_2_ths_km15_panel20, color = 'cyan', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(KM15)}$")
# axs[2, 0].errorbar(ts_theory_panel20, cosine_2_ths_vgg_panel20_smooth , color = 'tab:orange', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(VGG)}$")
# axs[2, 0].errorbar(ts_theory_panel20, cosine_2_ths_bh_panel20_smooth  , color = 'tab:red', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(BH)}$")

# axs[0, 1].errorbar(Q2s_panel01, cosine_3_exps_pi0_panel01, yerr = cosine_3_exps_pi0_stat_err_panel01, marker = 'o', color = 'k', zorder = 5, label = r"$\mathrm{Measurement}$" + "\n"+ r"$\mathrm{(This~work)}$")
# axs[0, 1].fill_between(Q2s_panel01, cosine_3_exps_pi0_panel01 - cosine_3_exps_pi0_syst_err_down_panel01, cosine_3_exps_pi0_panel01 + cosine_3_exps_pi0_syst_err_up_panel01, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[0, 1].errorbar(Q2s_theory_panel01, cosine_2_ths_km15_panel01, color = 'cyan', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(KM15)}$")
# axs[0, 1].errorbar(Q2s_theory_panel01, cosine_2_ths_vgg_panel01_smooth , color = 'tab:orange', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(VGG)}$")
# axs[0, 1].errorbar(Q2s_theory_panel01, cosine_2_ths_bh_panel01_smooth  , color = 'tab:red', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(BH)}$")

# axs[1, 1].errorbar(Q2s_panel11, cosine_3_exps_pi0_panel11, yerr = cosine_3_exps_pi0_stat_err_panel11, marker = 'o', color = 'k', zorder = 5, label = r"$\mathrm{Measurement}$" + "\n"+ r"$\mathrm{(This~work)}$")
# axs[1, 1].fill_between(Q2s_panel11, cosine_3_exps_pi0_panel11 - cosine_3_exps_pi0_syst_err_down_panel11, cosine_3_exps_pi0_panel11 + cosine_3_exps_pi0_syst_err_up_panel11, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[1, 1].errorbar(Q2s_theory_panel11, cosine_2_ths_km15_panel11, color = 'cyan', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(KM15)}$")
# axs[1, 1].errorbar(Q2s_theory_panel11, cosine_2_ths_vgg_panel11_smooth , color = 'tab:orange', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(VGG)}$")
# axs[1, 1].errorbar(Q2s_theory_panel11, cosine_2_ths_bh_panel11_smooth  , color = 'tab:red', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(BH)}$")

# axs[2, 1].errorbar(Q2s_panel21, cosine_3_exps_pi0_panel21, yerr = cosine_3_exps_pi0_stat_err_panel21, marker = 'o', color = 'k', zorder = 5, label = r"$\mathrm{Measurement}$" + "\n"+ r"$\mathrm{(This~work)}$")
# axs[2, 1].fill_between(Q2s_panel21, cosine_3_exps_pi0_panel21 - cosine_3_exps_pi0_syst_err_down_panel21, cosine_3_exps_pi0_panel21 + cosine_3_exps_pi0_syst_err_up_panel21, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
# axs[2, 1].errorbar(Q2s_theory_panel21, cosine_2_ths_km15_panel21, color = 'cyan', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(KM15)}$")
# axs[2, 1].errorbar(Q2s_theory_panel21, cosine_2_ths_vgg_panel21_smooth , color = 'tab:orange', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(VGG)}$")
# axs[2, 1].errorbar(Q2s_theory_panel21, cosine_2_ths_bh_panel21_smooth  , color = 'tab:red', label = r"$\mathrm{Theory}$" + "\n" + r"$\mathrm{(BH)}$")

# axs[0, 0].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel00)), xycoords = 'axes fraction')
# axs[1, 0].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel10)), xycoords = 'axes fraction')
# axs[2, 0].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel20)), xycoords = 'axes fraction')
# axs[0, 1].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel01)), xycoords = 'axes fraction')
# axs[1, 1].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel11)), xycoords = 'axes fraction')
# axs[2, 1].annotate(xy = (0.45, 0.8), xytext = (0.45, 0.8), text = r"$\langle x_B \rangle = {:.3f}$".format(np.mean(xB_avgs_panel21)), xycoords = 'axes fraction')

# for ax in axs[:, 0]:
#     ax.set_xlim([0.15, 0.95])
#     ax.set_xticks(np.linspace(0.16, 1.0, 43), minor = True)
#     ax.set_xticks(np.linspace(0.2, 1.0, 9), ['$0.2$', '', '$0.4$', '', '$0.6$', '', '$0.8$', '', '$1.0$'])
#     ax.tick_params(right = True, left = True, top = True, bottom = True, direction = 'in', which = 'both')

# for ax in axs[:, 1]:
#     ax.set_xlim([1, 4])
#     ax.set_xticks(np.linspace(1, 4, 31), minor = True)
#     ax.set_xticks(np.linspace(1, 4, 4))
#     ax.tick_params(right = True, left = True, top = True, bottom = True, direction = 'in', which = 'both')

# for ax in axs[:2, 0]:
#     ax.set_xticks(np.linspace(0.2, 1.0, 9), ['']*9)
#     # ax.set_xticks(np.linspace(0.15, 0.95, 81), minor = True)

# for ax in axs[:2, 1]:
#     ax.set_xticks(np.linspace(1, 4, 4), ['']*4)
#     # ax.set_xticks(np.linspace(1, 4, 31), minor = True)

# handles, labels = axs[2, 1].get_legend_handles_labels()
# # handles = []
# # labels = [labels[5], labels[4], labels[2], labels[0], labels[1], labels[3], labels[8], labels[6], labels[7]]
# plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (.9, 0.75), title = r"$d\sigma_{ep\rightarrow e'p'\gamma}^{\cos \phi,w}~(\mathrm{nb/GeV}^4)$", fontsize = 20, title_fontsize = 20, markerscale = 1.3)

# axs[-1, 0].set_xlabel(r"$|t|~(\mathrm{GeV}^2)$")
# axs[-1, 1].set_xlabel(r"$Q^2~(\mathrm{GeV}^2/c^2)$")
# plt.subplots_adjust(hspace = 0., wspace = 0.23)


# axs[0, 0].set_ylim([-0.06, 0.11])
# axs[0, 0].set_yticks([-0.05, 0.00, 0.05, 0.1])
# axs[0, 0].set_yticks(np.linspace(-0.06, 0.11, 18), minor = True)


# axs[0, 1].set_ylim([-0.03, 0.15])
# axs[0, 1].set_yticks([0, 0.05, 0.1, 0.15])
# axs[0, 1].set_yticks(np.linspace(-0.03, 0.15, 19), minor = True)


# axs[1, 0].set_ylim([-0.03, 0.07])
# axs[1, 0].set_yticks([-0.02, 0, 0.02, 0.04, 0.06])
# axs[1, 0].set_yticks(np.linspace(-0.03, 0.07, 11), minor = True)

# axs[1, 1].set_ylim([-0.01, 0.07])
# axs[1, 1].set_yticks([0, 0.02, 0.04, 0.06])
# axs[1, 1].set_yticks(np.linspace(-0.01, 0.07, 9), minor = True)

# axs[2, 0].set_ylim([-0.025, 0.045])
# axs[2, 0].set_yticks([-0.02, 0, 0.02, 0.04])
# axs[2, 0].set_yticks(np.linspace(-0.02, 0.04, 4), minor = True)

# axs[2, 1].set_ylim([-0.005, 0.045])
# axs[2, 1].set_yticks([0, 0.02, 0.04])
# axs[2, 1].set_yticks(np.linspace(-0.00, 0.04, 5), minor = True)

# axs[0, 0].set_title("$\langle Q^2 \\rangle = {:.3f}~\mathrm{{GeV}}/c^2$".format(np.mean([*Q2_avgs_panel00, *Q2_avgs_panel10, *Q2_avgs_panel20])))
# axs[0, 1].set_title("$\langle |t| \\rangle = {:.3f}~\mathrm{{GeV}}^2$".format(np.mean([*t_avgs_panel01, *t_avgs_panel11, *t_avgs_panel21])))
# plt.savefig("addendum_v3/modified_xsec_entire.from_script.pdf", bbox_inches = 'tight')

'''
# global normalization
fig, axs = plt.subplots(7, 4, figsize = (30, 40))

labeled = 0

for i, integrated_binnum in enumerate(df_summary_table_rebinned.loc[(y(df_summary_table_rebinned.xB_avg_this_point, df_summary_table_rebinned.Q2_avg_this_point, 0, 0)>0.7), :].integrated_binnum.unique()):
    df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.km15_0d005_ratio>0.9) & (df_summary_table_rebinned.pureBH_0d005_ratio>0.9), :]
    df_display_this_bin = df_display.loc[df_display.integratedbin_display == integrated_binnum, :]
    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9), :]
    if len(df_summary_table_rebinned_this_bin)>0:
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].xsec_KM15_display, color = 'cyan', label = "KM15")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].xsec_VGG_display, color = 'tab:orange', label = "VGG")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].xsec_BH_KM15_display, color = 'red', label = "BH")
        axs[i%7, i//7].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging, yerr = df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging_stat_err, ls = '', marker = 'o', color = 'k', label = "Exp. (Bkg merging only)")
        axs[i%7, i//7].fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging -  df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging_syst_err, df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging +  df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging_syst_err,  color = 'k', alpha = .5)
        if not labeled:
            handles, labels = axs[i%7, i//7].get_legend_handles_labels()
            plt.figlegend(handles, labels, bbox_to_anchor = (.95, 0.5), loc = 'upper left')
        if len(df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), :])>0:
            xB_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "xB_avg_this_point"].to_numpy()
            Q2_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "Q2_avg_this_point"].to_numpy()
            t_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "t_avg_this_point"].to_numpy()
            phi_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "phi_avg_this_point"].to_numpy()
            # integrated_binnum = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "integrated_binnum"].to_numpy()
            xsec_BH_model   = printBHarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), local = True)
            xsec_VGG_model   = printVGGarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), local = True)
            xsec_KM15_model = printKMarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 5)
            xsec_BH_KM15_model = printKMarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 1)
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_BH_model"] = xsec_BH_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_VGG_model"] = xsec_VGG_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_KM15_model"] = xsec_KM15_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_BH_KM15_model"] = xsec_BH_KM15_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "exp_to_BH_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_inb_exp_bkg_merging"]/ xsec_BH_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "exp_to_BH_inb_stat_err"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_inb_exp_bkg_merging_stat_err"]/ xsec_BH_model

        
        ymax1 = np.max(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
        ymax2 = np.max(df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging + df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging_stat_err)
        ymax = 1.1*np.max([ymax1, ymax2])
        ymin1 = np.min(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
        ymin2 = np.min(df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging - df_summary_table_rebinned_this_bin.xsec_inb_exp_bkg_merging_stat_err)
        ymin = 0.9*np.min([ymin1, ymin2])
    else:
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].xsec_KM15_display, color = 'cyan', label = "KM15")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].xsec_VGG_display, color = 'tab:orange', label = "VGG")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].xsec_BH_KM15_display, color = 'red', label = "BH")
        ymax = 1.1*np.max(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
        ymin = 0.9*np.min(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
    # axs[i%7, i//7].fill_between(np.linspace(0, 45, 101), 0.8*ymax*np.ones(101), ymax*np.ones(101), color = 'k', alpha = 0.3)
    # axs[i%7, i//7].fill_between(np.linspace(315, 360, 101), 0.8*ymax*np.ones(101), ymax*np.ones(101), color = 'k', alpha = 0.3)
    axs[i%7, i//7].set_xlim([0, 360])
    axs[i%7, i//7].set_yscale('log')
    axs[i%7, i//7].set_ylim([ymin, ymax])
    axs[i%7, i//7].set_xticks([0, 90, 180, 270, 360], ['']*5)
    axs[i%7, i//7].set_xticks(np.linspace(0, 360, 9), minor = True)
    axs[i%7, i//7].tick_params(axis = 'both', direction = 'inout', which = 'major', right = True, top = True, length = 10)
    axs[i%7, i//7].tick_params(axis = 'both', direction = 'inout', which = 'minor', right = True, top = True, length = 5)
    df_summary_table_rebinned_this_bin_dummy = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), :]
    axs[i%7, i//7].annotate(r"${}.$".format(integrated_binnum),  xy = (0.12, 0.9), xytext = (0.12, 0.9), xycoords = 'axes fraction')
    axs[i%7, i//7].annotate(r"$x_B \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned_this_bin_dummy.xBmin.unique()[0], df_summary_table_rebinned_this_bin_dummy.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_summary_table_rebinned_this_bin_dummy.Q2min.unique()[0], df_summary_table_rebinned_this_bin_dummy.Q2max.unique()[0]) + "\n" + r"$|t|~ \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned_this_bin_dummy.t1min.unique()[0], df_summary_table_rebinned_this_bin_dummy.t1max.unique()[0]), xy = (0.21, 0.73), xytext = (0.21, 0.73), xycoords = 'axes fraction')
for ax in axs[-1, :].flatten():
    ax.set_xticks([0, 90, 180, 270, 360], ['${}$'.format(i) for i in [0, 90, 180, 270, 360]])
    ax.set_xticks(np.linspace(0, 360, 9), minor = True)
    ax.set_xlabel("$\phi$" + r" ($^{\circ}$)")
for ax in axs[:, 0].flatten():
    ax.set_ylabel(r"$\frac{d\sigma}{d x_B dQ^2 d |t| d\phi}$" + " " + r"$(\mathrm{{nb/GeV}}^4)$")


fig.subplots_adjust(hspace = 0)
plt.savefig("addendum_v3/inbending_global_normalization.from_script.pdf", bbox_inches = 'tight')
plt.close()

fig, axs = plt.subplots(7, 4, figsize = (30, 40))

labeled = 0

for i, integrated_binnum in enumerate(df_summary_table_rebinned.loc[(y(df_summary_table_rebinned.xB_avg_this_point, df_summary_table_rebinned.Q2_avg_this_point, 0, 0)>0.7), :].integrated_binnum.unique()):
    df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.km15_0d005_ratio>0.9) & (df_summary_table_rebinned.pureBH_0d005_ratio>0.9), :]
    df_display_this_bin = df_display.loc[df_display.integratedbin_display == integrated_binnum, :]
    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9), :]
    if len(df_summary_table_rebinned_this_bin)>0:
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].xsec_KM15_display, color = 'cyan', label = "KM15")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].xsec_VGG_display, color = 'tab:orange', label = "VGG")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].xsec_BH_KM15_display, color = 'red', label = "BH")
        axs[i%7, i//7].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging, yerr = df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging_stat_err, ls = '', marker = 'o', color = 'k', label = "Exp. (Bkg merging only)")
        axs[i%7, i//7].fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging -  df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging_syst_err, df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging +  df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging_syst_err,  color = 'k', alpha = .5)
        if not labeled:
            handles, labels = axs[i%7, i//7].get_legend_handles_labels()
            plt.figlegend(handles, labels, bbox_to_anchor = (.95, 0.5), loc = 'upper left')
        if len(df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), :])>0:
            xB_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "xB_avg_this_point"].to_numpy()
            Q2_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "Q2_avg_this_point"].to_numpy()
            t_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "t_avg_this_point"].to_numpy()
            phi_avg_this_point = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "phi_avg_this_point"].to_numpy()
            # integrated_binnum = df_summary_table_rebinned_this_bin.loc[(df_summary_table_rebinned_this_bin.phi_avg_this_point<=45) | (df_summary_table_rebinned_this_bin.phi_avg_this_point>315), "integrated_binnum"].to_numpy()
            xsec_BH_model   = printBHarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), local = True)
            xsec_VGG_model   = printVGGarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), local = True)
            xsec_KM15_model = printKMarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 5)
            xsec_BH_KM15_model = printKMarray(xB_avg_this_point, Q2_avg_this_point, t_avg_this_point, np.radians(phi_avg_this_point), mode = 1)
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_BH_model"] = xsec_BH_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_VGG_model"] = xsec_VGG_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_KM15_model"] = xsec_KM15_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_BH_KM15_model"] = xsec_BH_KM15_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "exp_to_BH_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_outb_exp_bkg_merging"]/ xsec_BH_model
            df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "exp_to_BH_outb_stat_err"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_bkg_merging == 1) & (df_summary_table_rebinned.km15_0d005_ratio > 0.9)  & (df_summary_table_rebinned.pureBH_0d005_ratio > 0.9) & ((df_summary_table_rebinned.phi_avg_this_point<=45) | (df_summary_table_rebinned.phi_avg_this_point>315)), "xsec_outb_exp_bkg_merging_stat_err"]/ xsec_BH_model

        
        ymax1 = np.max(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
        ymax2 = np.max(df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging + df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging_stat_err)
        ymax = 1.1*np.max([ymax1, ymax2])
        ymin1 = np.min(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
        ymin2 = np.min(df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging - df_summary_table_rebinned_this_bin.xsec_outb_exp_bkg_merging_stat_err)
        ymin = 0.9*np.min([ymin1, ymin2])
    else:
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_KM15_display>0].xsec_KM15_display, color = 'cyan', label = "KM15")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0].xsec_VGG_display, color = 'tab:orange', label = "VGG")
        axs[i%7, i//7].plot(df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].phi_display*180./np.pi, df_display_this_bin.loc[df_display_this_bin.xsec_BH_KM15_display>0].xsec_BH_KM15_display, color = 'red', label = "BH")
        ymax = 1.1*np.max(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
        ymin = 0.9*np.min(df_display_this_bin.loc[df_display_this_bin.xsec_VGG_display>0, ["xsec_KM15_display", "xsec_VGG_display", "xsec_BH_KM15_display"]])
    # axs[i%7, i//7].fill_between(np.linspace(0, 45, 101), 0.8*ymax*np.ones(101), ymax*np.ones(101), color = 'k', alpha = 0.3)
    # axs[i%7, i//7].fill_between(np.linspace(315, 360, 101), 0.8*ymax*np.ones(101), ymax*np.ones(101), color = 'k', alpha = 0.3)
    axs[i%7, i//7].set_xlim([0, 360])
    axs[i%7, i//7].set_yscale('log')
    axs[i%7, i//7].set_ylim([ymin, ymax])
    axs[i%7, i//7].set_xticks([0, 90, 180, 270, 360], ['']*5)
    axs[i%7, i//7].set_xticks(np.linspace(0, 360, 9), minor = True)
    axs[i%7, i//7].tick_params(axis = 'both', direction = 'inout', which = 'major', right = True, top = True, length = 10)
    axs[i%7, i//7].tick_params(axis = 'both', direction = 'inout', which = 'minor', right = True, top = True, length = 5)
    df_summary_table_rebinned_this_bin_dummy = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum), :]
    axs[i%7, i//7].annotate(r"${}.$".format(integrated_binnum),  xy = (0.12, 0.9), xytext = (0.12, 0.9), xycoords = 'axes fraction')
    axs[i%7, i//7].annotate(r"$x_B \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned_this_bin_dummy.xBmin.unique()[0], df_summary_table_rebinned_this_bin_dummy.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_summary_table_rebinned_this_bin_dummy.Q2min.unique()[0], df_summary_table_rebinned_this_bin_dummy.Q2max.unique()[0]) + "\n" + r"$|t|~ \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned_this_bin_dummy.t1min.unique()[0], df_summary_table_rebinned_this_bin_dummy.t1max.unique()[0]), xy = (0.21, 0.73), xytext = (0.21, 0.73), xycoords = 'axes fraction')
for ax in axs[-1, :].flatten():
    ax.set_xticks([0, 90, 180, 270, 360], ['${}$'.format(i) for i in [0, 90, 180, 270, 360]])
    ax.set_xticks(np.linspace(0, 360, 9), minor = True)
    ax.set_xlabel("$\phi$" + r" ($^{\circ}$)")
for ax in axs[:, 0].flatten():
    ax.set_ylabel(r"$\frac{d\sigma}{d x_B dQ^2 d |t| d\phi}$" + " " + r"$(\mathrm{{nb/GeV}}^4)$")


fig.subplots_adjust(hspace = 0)
plt.savefig("addendum_v3/outbending_global_normalization.from_script.pdf", bbox_inches = 'tight')
plt.close()


for integrated_binnum in  df_summary_table_kinematics.integrated_binnum.unique():
    df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "xB_avg_this_point"] = float(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum].xB_avg_this_point.unique()[0])
    df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "Q2_avg_this_point"] = float(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum].Q2_avg_this_point.unique()[0])
    df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "t_avg_this_point"]  = float(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum].t_avg_this_point.unique()[0])
    df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "xBbin"]             =   int(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum].xBbin.unique()[0])
    df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "Q2bin"]             =   int(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum].Q2bin.unique()[0])
    df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "tbin"]              =   int(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum].tbin.unique()[0])
    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[ (df_summary_table_rebinned.integrated_binnum == integrated_binnum)  & (df_summary_table_rebinned.exp_to_BH_inb>0), :]
    if len(df_summary_table_rebinned_this_bin)> 0:
    # df_summary_table_kinematics.loc[:, "exp_to_BH_mean"] = df_summary_table_rebinned_this_bin.exp
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_inb_mean"]          = np.sum(df_summary_table_rebinned_this_bin.exp_to_BH_inb / df_summary_table_rebinned_this_bin.exp_to_BH_inb_stat_err/ df_summary_table_rebinned_this_bin.exp_to_BH_inb_stat_err)/ np.sum(1/df_summary_table_rebinned_this_bin.exp_to_BH_inb_stat_err/df_summary_table_rebinned_this_bin.exp_to_BH_inb_stat_err)
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_inb_mean_stat_err"] = np.sqrt(1/ np.sum(1/df_summary_table_rebinned_this_bin.exp_to_BH_inb_stat_err/df_summary_table_rebinned_this_bin.exp_to_BH_inb_stat_err) )
    else:
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_inb_mean"]          = 0
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_inb_mean_stat_err"] = 0

    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[ (df_summary_table_rebinned.integrated_binnum == integrated_binnum)  & (df_summary_table_rebinned.exp_to_BH_outb>0), :]
    if len(df_summary_table_rebinned_this_bin)> 0:
    # df_summary_table_kinematics.loc[:, "exp_to_BH_mean"] = df_summary_table_rebinned_this_bin.exp
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_outb_mean"]          = np.sum(df_summary_table_rebinned_this_bin.exp_to_BH_outb / df_summary_table_rebinned_this_bin.exp_to_BH_outb_stat_err/ df_summary_table_rebinned_this_bin.exp_to_BH_outb_stat_err)/ np.sum(1/df_summary_table_rebinned_this_bin.exp_to_BH_outb_stat_err/df_summary_table_rebinned_this_bin.exp_to_BH_outb_stat_err)
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_outb_mean_stat_err"] = np.sqrt(1/ np.sum(1/df_summary_table_rebinned_this_bin.exp_to_BH_outb_stat_err/df_summary_table_rebinned_this_bin.exp_to_BH_outb_stat_err) )
    else:
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_outb_mean"]          = 0
        df_summary_table_kinematics.loc[df_summary_table_kinematics.integrated_binnum == integrated_binnum, "exp_to_BH_outb_mean_stat_err"] = 0


fig, axs = plt.subplots(1, 1, figsize = (12, 8))
for tbin in range(6):
    # for Q2bin in range(7):
    #     axs[tbin%3, tbin//3].errorbar(df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].xB_avg_this_point, df_summary_table_kinematics.loc[(df_summary_table_kinematics.Q2bin == Q2bin) & (df_summary_table_kinematics.tbin == tbin)].exp_to_BH_inb_mean, df_summary_table_kinematics.loc[(df_summary_table_kinematics.Q2bin == Q2bin) & (df_summary_table_kinematics.tbin == tbin)].exp_to_BH_inb_mean_stat_err, ls = '', marker = 'o')
    t1min = df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin == tbin, "t1min"].unique()[0]
    t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin == tbin, "t1max"].unique()[0]
    axs.errorbar(df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].Q2_avg_this_point, df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].exp_to_BH_inb_mean, df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].exp_to_BH_inb_mean_stat_err, ls = '', marker = 'o', label = r"${}<|t|<{}$".format(t1min,t1max))
    axs.set_ylim([0.3, 1.2])

exp_to_BH_inb_mean_mean     = df_summary_table_kinematics.loc[(df_summary_table_kinematics.exp_to_BH_inb_mean>0)].exp_to_BH_inb_mean.mean()
exp_to_BH_inb_mean_stat_err = df_summary_table_kinematics.loc[(df_summary_table_kinematics.exp_to_BH_inb_mean>0)].exp_to_BH_inb_mean.std()
axs.fill_between(np.linspace(1, 5.5, 101), np.ones(101)*(exp_to_BH_inb_mean_mean - exp_to_BH_inb_mean_stat_err), np.ones(101)*(exp_to_BH_inb_mean_mean + exp_to_BH_inb_mean_stat_err), alpha = 0.5, label = r'$\mathrm{average}$')
axs.legend(loc = 'upper left', bbox_to_anchor = (0.4, 0.45), ncol = 2, title = r"$\mathrm{Inbending}$")
axs.set_xlim([1, 5.5])
axs.set_xlabel(r"$\langle Q^2 \rangle ~(\mathrm{GeV}^2/c^2)$")
axs.set_ylabel(r"$\mathrm{Exp.~to~BH}$")
plt.savefig("addendum_v3/inbending_global_norm.from_script.pdf", bbox_inches = 'tight')

fig, axs = plt.subplots(1, 1, figsize = (12, 8))
for tbin in range(6):
    # for Q2bin in range(7):
    #     axs[tbin%3, tbin//3].errorbar(df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].xB_avg_this_point, df_summary_table_kinematics.loc[(df_summary_table_kinematics.Q2bin == Q2bin) & (df_summary_table_kinematics.tbin == tbin)].exp_to_BH_outb_mean, df_summary_table_kinematics.loc[(df_summary_table_kinematics.Q2bin == Q2bin) & (df_summary_table_kinematics.tbin == tbin)].exp_to_BH_outb_mean_stat_err, ls = '', marker = 'o')
    t1min = df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin == tbin, "t1min"].unique()[0]
    t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin == tbin, "t1max"].unique()[0]
    axs.errorbar(df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].Q2_avg_this_point, df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].exp_to_BH_outb_mean, df_summary_table_kinematics.loc[(df_summary_table_kinematics.tbin == tbin)].exp_to_BH_outb_mean_stat_err, ls = '', marker = 'o', label = r"${}<|t|<{}$".format(t1min,t1max))
    axs.set_ylim([0.3, 1.2])

exp_to_BH_outb_mean_mean     = df_summary_table_kinematics.loc[(df_summary_table_kinematics.exp_to_BH_outb_mean>0)].exp_to_BH_outb_mean.mean()
exp_to_BH_outb_mean_stat_err = df_summary_table_kinematics.loc[(df_summary_table_kinematics.exp_to_BH_outb_mean>0)].exp_to_BH_outb_mean.std()
axs.fill_between(np.linspace(1, 5.5, 101), np.ones(101)*(exp_to_BH_outb_mean_mean - exp_to_BH_outb_mean_stat_err), np.ones(101)*(exp_to_BH_outb_mean_mean + exp_to_BH_outb_mean_stat_err), alpha = 0.5, label = r'$\mathrm{average}$')
axs.legend(loc = 'upper left', bbox_to_anchor = (0.5, 0.98), ncol = 2, title = r"$\mathrm{Outbending}$")
axs.set_xlim([1, 5.5])
axs.set_xlabel(r"$\langle Q^2 \rangle ~(\mathrm{GeV}^2/c^2)$")
axs.set_ylabel(r"$\mathrm{Exp.~to~BH}$")
plt.savefig("addendum_v3/outbending_global_norm.from_script.pdf", bbox_inches = 'tight')
'''

'''
#stamp plots for the cross sections

# stamp plots for normalization
df_summary_table_rebinned.loc[:, "active_bin_nominal"] = 0
# df_summary_table_rebinned.loc[(df_summary_table_rebinned.epg_exp > 10) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio < 0.4) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.4), "active_bin_nominal"] = 1
# df_summary_table_rebinned.loc[(df_summary_table_rebinned.epg_exp > 10) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio_bin_by_bin < 0.6) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.4), "active_bin_nominal"] = 1
df_summary_table_rebinned.loc[(df_summary_table_rebinned.epg_exp > 10) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio < 0.6) & (df_summary_table_rebinned.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio < 0.4), "active_bin_nominal"] = 1
df_summary_table_rebinned.loc[df_summary_table_rebinned.xsec_exp == 0, "active_bin_nominal"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.pureBH_0d0005_ratio < 0.9, "active_bin_nominal"] = 0
df_summary_table_rebinned.loc[df_summary_table_rebinned.km15_0d0005_ratio < 0.9, "active_bin_nominal"] = 0
# df_summary_table_rebinned.loc[df_summary_table_rebinned.efficiency < 0.25, "active_bin_nominal"] = 0
# df_summary_table_rebinned.loc[(df_summary_table_rebinned.tbin == 0) & (df_summary_table_rebinned.phi_avg_this_point>90) & (df_summary_table_rebinned.phi_avg_this_point<270) & (df_summary_table_rebinned.xsec_exp < 0.9* df_summary_table_rebinned.pureBH_cross_section_this_point_norad ), "active_bin_nominal"]=0
df_summary_table_rebinned.loc[(df_summary_table_rebinned.tbin == 0) & (df_summary_table_rebinned.phi_avg_this_point>90) & (df_summary_table_rebinned.phi_avg_this_point<270) & (df_summary_table_rebinned.xsec_exp < 0.9* df_summary_table_rebinned.pureBH_cross_section_this_point_norad ), "active_bin_nominal"]=0
print(df_summary_table_rebinned.active_bin_nominal.sum())


fig, axs = plt.subplots(1, 1, figsize = (10, 6))
integrated_binnum = 88

# models         = ["pi0", "bh", "km15", "vgg", "global1"]
# # label_scheme = [r"$\mathrm{Global1~Norm.}$", r"$\mathrm{Global2~Norm.}$", "$\pi^0~\mathrm{Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{VGG~Norm.}$"]
# label_scheme = ["$\pi^0~\mathrm{Norm.~(Nominal)}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
# color_scheme = ['k', 'red', 'cyan', 'tab:orange', 'tab:blue']

models = ["", "pi0_eff_corrected_"]#, "bh_eff_corrected_", "km15_eff_corrected_", "vgg_eff_corrected_", "global1_eff_corrected_"]
# label_scheme = [r"$\mathrm{Global1~Norm.}$", r"$\mathrm{Global2~Norm.}$", "$\pi^0~\mathrm{Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{VGG~Norm.}$"]
label_scheme = ["$\mathrm{Without}$" + "\n" + "$\pi^0~\mathrm{normalization}$", "$\mathrm{With}$" + "\n" + "$\pi^0~\mathrm{normalization}$"]#, r"$\mathrm{BH-Local}$", r"$\mathrm{KM15-Local}$", r"$\mathrm{VGG-Local}$", r"$\mathrm{BH-Global}$"]
color_scheme = ['tab:pink', 'k']#, 'red', 'cyan', 'tab:orange', 'tab:blue']


df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal==1), :]

weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
# weights_err_up        = np.sqrt(weights_stat_err**2 + weights_syst_err_up**2)
# weights_err_down      = np.sqrt(weights_stat_err**2 + weights_syst_err_down**2)

# axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Exp.}$", histtype = 'step', color = 'k')
# axs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = "$\mathrm{Measurement~(This~work)}$")

for i, model in enumerate(models):
    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal".format(model)] == 1), :]
    if i == 0:
        dots = plt.scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
        plt.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
        axs.fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst. Uncertainty}$')
    else:
        plt.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)

plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=2, ls = '--', label = r"$\mathrm{Theory~(BH)}$")
plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=2, ls = '--', label = r"$\mathrm{Theory~(KM)}$")
plt.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_VGG_display, color = 'tab:orange', lw=2, ls = '--', label = r"$\mathrm{Theory~(VGG)}$")


xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]
print(xBmean, Q2mean, t1mean)
axs.annotate("$(x_B, Q^2) = ({:.2f}, {:.2f})$".format(df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0], df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.17, 0.88), xytext = (0.17, 0.88), xycoords = 'axes fraction', fontsize = 40)


axs.set_ylabel(r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$")
axs.set_xlabel("$\phi$~($^\circ$)")
axs.set_xlim([0, 360])
axs.set_xticks([0, 90, 180, 270, 360])
handles, labels = axs.get_legend_handles_labels()
# xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
# Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
# t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
orders = [5, 4, 0, 1, 2, 3, ]
handles = [handles[order] for order in orders]
labels = [labels[order] for order in orders]

xBheader = "$\\langle x_B \\rangle = {:.3f}$ \n".format(xBmean)
Q2header = "$\\langle Q^2 \\rangle = {:.3f}~\mathrm{{GeV}}^2/c^2$ \n".format(Q2mean)
t1header = "$\\langle |t| \\rangle~= {:.3f}~\mathrm{{GeV}}^2$".format(t1mean)

plt.figlegend(handles, labels, loc = 'upper left', title = r"${}.$".format(integrated_binnum)+"\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, .9), ncol = 1, alignment='left', columnspacing=0)

# axs.annotate(xy = (0.01, 0.7), xytext = (0.01, 0.7), text = xBheader + Q2header + t1header, xycoords = 'axes fraction', fontsize = 20)
plt.savefig("addendum_v3/xsec_{}_normalization.new.from_script.pdf".format(integrated_binnum), bbox_inches = 'tight')


xB_panes = 8
Q2_panes = 7
display_factor = 1

# models        = ["", "pi0_eff_corrected_", "bh_eff_corrected_", "km15_eff_corrected_", "vgg_eff_corrected_", "global1_eff_corrected_"]
# # label_scheme = [r"$\mathrm{Global1~Norm.}$", r"$\mathrm{Global2~Norm.}$", "$\pi^0~\mathrm{Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{VGG~Norm.}$"]
# label_scheme = ["$\mathrm{Without}$" + "\n" + "$\mathrm{additional~normalization}$", "$\mathrm{Nominal}~(\pi^0-\mathrm{Local})$", r"$\mathrm{BH-Local}$", r"$\mathrm{KM15-Local}$", r"$\mathrm{VGG-Local}$", r"$\mathrm{BH-Global}$"]
# color_scheme = ['tab:pink', 'k', 'red', 'cyan', 'tab:orange', 'tab:blue']
with_binnumber = 1

models = ["", "pi0_eff_corrected_"]#, "bh_eff_corrected_", "km15_eff_corrected_", "vgg_eff_corrected_", "global1_eff_corrected_"]
# label_scheme = [r"$\mathrm{Global1~Norm.}$", r"$\mathrm{Global2~Norm.}$", "$\pi^0~\mathrm{Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{VGG~Norm.}$"]
label_scheme = ["$\mathrm{Without}$" + "\n" + "$\pi^0~\mathrm{normalization}$", "$\mathrm{With}~\pi^0~\mathrm{normalization}$", "test"]#, r"$\mathrm{BH-Local}$", r"$\mathrm{KM15-Local}$", r"$\mathrm{VGG-Local}$", r"$\mathrm{BH-Global}$"]
color_scheme = ['tab:pink', 'k']#, 'red', 'cyan', 'tab:orange', 'tab:blue']

for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (65, 20))
    gs = GridSpec(Q2_panes, xB_panes, figure=fig)
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_nominal.sum() < 4:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            
            t_avgs.append(df_this_bin.t_avg_this_point.unique()[0])
            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]
            if integrated_binnum in [54, 106, 126]:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue
            if with_binnumber:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}.~({:.2f}, {:.2f})$".format(integrated_binnum, df_this_bin.xB_avg_this_point.unique()[0], df_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.23, 0.82), xytext = (0.23, 0.82), xycoords = 'axes fraction', fontsize = 40)
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("$({:.2f}, {:.2f})$".format(df_this_bin.xB_avg_this_point.unique()[0], df_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.3, 0.8), xytext = (0.3, 0.8), xycoords = 'axes fraction', fontsize = 40)                    
                    
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=3, label = r"$\mathrm{Theory~(KM15)}$")
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].phi_display), df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].xsec_VGG_display, color = 'tab:orange', lw=3, label = r"$\mathrm{Theory~(VGG)}$")
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=3, label = r"$\mathrm{Theory~(BH)}$")

            df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal"] == 1), :]
            weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
            weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
            weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
            for i, model in enumerate(models):
                if i == 0:
                    dots = axs[Q2_panes - Q2_binnum - 1, xB_binnum].scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
                else:
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)

            if integrated_binnum == {0: 63, 1: 87, 2: 88, 3: 89, 4: 90, 5: 91}[t_binnum]:
                # print(Q2_panes, xB_binnum)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].patch.set_linewidth(12)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].patch.set_edgecolor('tab:purple')
            
            if not labeled:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                # orders = [4, 5, 6, 7, 8, 3, 2, 0, 1 ]
                # handles = [handles[order] for order in orders]
                # labels = [labels[order] for order in orders]
                
                labeled = 1
            
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_yscale('log')
            ymin = np.min(df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display)
            ymin = np.log10(ymin)
            ymax = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.phi_avg_this_point > 30)  & (df_summary_table_rebinned.phi_avg_this_point < 330), ["xsec_exp_pi0_eff_corrected_bkg_merging", "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err", "xsec_exp_pi0_eff_corrected_bkg_merging_stat_err"]].sum(axis = 1).max() 
            ymax = np.log10(ymax)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 12+1), minor = True)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)", fontsize = 40 )
            if integrated_binnum in [4, 12, 32, 60, 84, 111, 130]:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 4+1))
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 4+1), ['']*5)

            adjustment_min = 0.2
            adjustment_max = 0.2
            if np.abs(ymin)%1 < 0.2:
                adjustment_min = 0.5
            if np.abs(ymax)%1 < 0.2:
                adjustment_max = 0.5
            if np.floor(ymax) - np.ceil(ymin) < 1  :
                if np.ceil(ymax) - ymax > ymin - np.floor(ymin):
                    ymin = np.floor(ymin)- 0.005
                else:
                    ymax = np.ceil(ymax) + 0.005
            if (xB_binnum == 6) & (Q2_binnum == 4):
                ymax = ymax + .2
            if (xB_binnum == 5) & (Q2_binnum == 3):
                ymax = ymax + .4
            if (xB_binnum == 5) & (Q2_binnum == 5):
                ymin = ymin + .1
            if (xB_binnum == 4) & (Q2_binnum == 3):
                ymax = ymax + .4
            if (xB_binnum == 6) & (Q2_binnum == 6):
                ymin = ymin + .2
            if (xB_binnum == 1) & (Q2_binnum == 0):
                ymin = ymin - .4
                ymax = ymax + .3
            if (xB_binnum == 4) & (Q2_binnum == 4):
                ymax = ymax - .3
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([10**(ymin-adjustment_min), 10**(ymax+adjustment_max)])
            ticks = np.logspace(np.ceil(ymin-adjustment_min), np.floor(ymax+adjustment_max), int( np.floor(ymax+adjustment_max)-np.ceil(ymin-adjustment_min)+1))
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_yticks(ticks, log_ticker(ticks))
            
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 40)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
    plt.figlegend(handles, labels, loc = 'lower right', bbox_to_anchor = (0.910, 0.050), fontsize = 40, title = r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$" + "\n" + r"$\langle |t| \rangle={:.3f}$".format(np.mean(t_avgs))+ r"$~\mathrm{GeV}^2$" + "\n" + r"$(\langle x_B \rangle, \langle Q^2 \rangle)$" + "\n" + r"$\mathrm{is~annotated~in~each~panel}$", title_fontsize = 40, alignment = 'left', ncol = 1, markerscale = 3, framealpha = 1, edgecolor = 'k')
    # axs[0, 0].annotate(xy = (0.0, 0.4), xytext = (0.0, 0.4), text = r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$" , xycoords = "axes fraction", fontsize = 60)
    plt.subplots_adjust(wspace = 0.2 , hspace = 0.0 )

    # if t_binnum == 3:

    ax_gs = plt.subplot(gs.new_subplotspec((0, 0), rowspan=3, colspan=2))
    integrated_binnum = {0: 63, 1: 87, 2: 88, 3: 89, 4: 90, 5: 91}[t_binnum]
    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal==1), :]
    
    weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
    weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
    weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
    weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
    
    # ax_gs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Exp.}$", histtype = 'step', color = 'k')
    # ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = "$\mathrm{Measurement~(This~work)}$")
    
    for i, model in enumerate(models):
        df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal".format(model)] == 1), :]
        if i == 0:
            dots = plt.scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
            ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
            ax_gs.fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
        else:
            ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
    
    ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=4, ls = '--', label = r"$\mathrm{Theory~(BH)}$")
    ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=4, ls = '--', label = r"$\mathrm{Theory~(KM)}$")
    ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_VGG_display, color = 'tab:orange', lw=4, ls = '--', label = r"$\mathrm{Theory~(VGG)}$")
    
    
    xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
    # print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)
    
    xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]
    # print(xBmean, Q2mean, t1mean)
    
    
    # ax_gs.set_ylabel(r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$")
    if with_binnumber:
        ax_gs.annotate("${}.~(x_B, Q^2) = ({:.2f}, {:.2f})$".format(integrated_binnum, df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0], df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.27, 0.9), xytext = (0.27, 0.9), xycoords = 'axes fraction', fontsize = 40)
    else:
        ax_gs.annotate("$(x_B, Q^2) = ({:.2f}, {:.2f})$".format(df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0], df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.30, 0.88), xytext = (0.30, 0.88), xycoords = 'axes fraction', fontsize = 40)
    ax_gs.set_xlabel("$\phi$~($^\circ$)", fontsize = 40)
    ax_gs.set_xlim([0, 360])
    ax_gs.set_xticks([0, 90, 180, 270, 360])
    ax_gs.set_xticklabels(['${}$'.format(n) for n  in [0, 90, 180, 270, 360]], fontsize = 40)
    if t_binnum > 0:
        ax_gs.set_ylim([0, 0.14])
        ax_gs.set_yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])
        ax_gs.set_yticklabels(['${:.2f}$'.format(n) for n  in [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14]])
    else:
        ax_gs.set_ylim([0, 1])
        ax_gs.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax_gs.set_yticklabels(['${:.2f}$'.format(n) for n  in [0, 2, 0.4, 0.6, 0.8, 1.0]])
    ax_gs.set_frame_on(True)
    ax_gs.patch.set_linewidth(12)
    ax_gs.patch.set_edgecolor('tab:purple')
    ax_gs.tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 40)
    ax_gs.tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
    ax_gs.tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
    ax_gs.tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
    
    # plt.savefig("addendum_v3/xsec_t{}_paper_w_different_normalizations.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.savefig("addendum_v3/xsec_t{}_an_w_different_normalizations.new.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    # plt.close()

xB_panes = 7
Q2_panes = 7
display_factor = 1

models = ["pi0", "bh", "km15", "vgg", "global1"]
# label_scheme = [r"$\mathrm{Global1~Norm.}$", r"$\mathrm{Global2~Norm.}$", "$\pi^0~\mathrm{Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{VGG~Norm.}$"]
label_scheme = ["$\mathrm{Measurement~(This~work)}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
color_scheme = ['k', 'red', 'cyan', 'tab:orange', 'tab:blue']



# stamp plot for the cross section
with_binnumber = 0
for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (65, 20))
    gs = GridSpec(Q2_panes, xB_panes, figure=fig)
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_nominal.sum() < 4:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            
            t_avgs.append(df_this_bin.t_avg_this_point.unique()[0])
            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]
            if integrated_binnum in [54, 106, 126]:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue
            if with_binnumber:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}.~({:.2f}, {:.2f})$".format(integrated_binnum, df_this_bin.xB_avg_this_point.unique()[0], df_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.3, 0.8), xytext = (0.3, 0.8), xycoords = 'axes fraction', fontsize = 40)
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("$({:.2f}, {:.2f})$".format(df_this_bin.xB_avg_this_point.unique()[0], df_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.3, 0.8), xytext = (0.3, 0.8), xycoords = 'axes fraction', fontsize = 40)                    
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=3, label = r"$\mathrm{Theory~(KM15)}$")
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].phi_display), df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].xsec_VGG_display, color = 'tab:orange', lw=3, label = r"$\mathrm{Theory~(VGG)}$")
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=3, label = r"$\mathrm{Theory~(BH)}$")

            df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal"] == 1), :]
            weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
            weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
            weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err

            for i, model in enumerate(models):
                if i == 0:
                    dots = axs[Q2_panes - Q2_binnum - 1, xB_binnum].scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
                else:
                    # axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                    pass

            if integrated_binnum == 89:
                # print(Q2_panes, xB_binnum)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].patch.set_linewidth(12)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].patch.set_edgecolor('tab:purple')
            
            if not labeled:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()                
                neworders = [3, 0, 1, 2]
                handles = [handles[i] for i in neworders]
                labels = [labels[i] for i in neworders]
                labeled = 1
            
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_yscale('log')
            ymin = np.min(df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display)
            ymin = np.log10(ymin)
            ymax = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.phi_avg_this_point > 30)  & (df_summary_table_rebinned.phi_avg_this_point < 330), ["xsec_exp_pi0_eff_corrected_bkg_merging", "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err", "xsec_exp_pi0_eff_corrected_bkg_merging_stat_err"]].sum(axis = 1).max() 
            ymax = np.log10(ymax)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 12+1), minor = True)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)", fontsize = 40 )
            if integrated_binnum in [4, 12, 32, 60, 84, 111, 130]:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 4+1))
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 4+1), ['']*5)

            adjustment_min = 0.2
            adjustment_max = 0.2
            if np.abs(ymin)%1 < 0.2:
                adjustment_min = 0.5
            if np.abs(ymax)%1 < 0.2:
                adjustment_max = 0.5
            if np.floor(ymax) - np.ceil(ymin) < 1  :
                if np.ceil(ymax) - ymax > ymin - np.floor(ymin):
                    ymin = np.floor(ymin)- 0.005
                else:
                    ymax = np.ceil(ymax) + 0.005
            if (xB_binnum == 6) & (Q2_binnum == 4):
                ymax = ymax + .2
            if (xB_binnum == 5) & (Q2_binnum == 3):
                ymax = ymax + .4
            if (xB_binnum == 5) & (Q2_binnum == 5):
                ymin = ymin + .1
            if (xB_binnum == 4) & (Q2_binnum == 3):
                ymax = ymax + .4
            if (xB_binnum == 6) & (Q2_binnum == 6):
                ymin = ymin + .2
            if (xB_binnum == 1) & (Q2_binnum == 0):
                ymin = ymin - .4
                ymax = ymax + .3
            if (xB_binnum == 4) & (Q2_binnum == 4):
                ymax = ymax - .3
            # print(ymin, ymax, df_this_bin.active_bin_nominal.sum())
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([10**(ymin-adjustment_min), 10**(ymax+adjustment_max)])
            ticks = np.logspace(np.ceil(ymin-adjustment_min), np.floor(ymax+adjustment_max), int( np.floor(ymax+adjustment_max)-np.ceil(ymin-adjustment_min)+1))
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_yticks(ticks, log_ticker(ticks))
            
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 40)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.760, 0.385), fontsize = 40, title = r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$" + "\n" + r"$\langle |t| \rangle={:.3f}$".format(np.mean(t_avgs))+ r"$~\mathrm{GeV}^2$"+ "\n" + r"$(\langle x_B \rangle, \langle Q^2 \rangle)~\mathrm{is~annotated~in~each~panel}$", title_fontsize = 40, alignment = 'left', ncol = 1, markerscale = 3, framealpha = 1, edgecolor = 'k')
    plt.subplots_adjust(wspace = 0.2 , hspace = 0.0 )

    if t_binnum == 3:
        ax_gs = plt.subplot(gs.new_subplotspec((0, 0), rowspan=3, colspan=2))
        integrated_binnum = 89
        df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal==1), :]
        
        weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
        weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
        weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
        weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
        
        # ax_gs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Exp.}$", histtype = 'step', color = 'k')
        # ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = "$\mathrm{Measurement~(This~work)}$")
        
        for i, model in enumerate(models):
            df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal".format(model)] == 1), :]
            if i == 0:
                dots = plt.scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
                ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                ax_gs.fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
            else:
                # ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                pass
        
        ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=4, ls = '--', label = r"$\mathrm{Theory~(BH)}$")
        ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=4, ls = '--', label = r"$\mathrm{Theory~(KM)}$")
        ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_VGG_display, color = 'tab:orange', lw=4, ls = '--', label = r"$\mathrm{Theory~(VGG)}$")
        
        
        xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
        # print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)
        
        xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]
        # print(xBmean, Q2mean, t1mean)
        
        
        # ax_gs.set_ylabel(r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$")
        if with_binnumber:
            ax_gs.annotate("${}.~(x_B, Q^2) = ({:.2f}, {:.2f})$".format(integrated_binnum, df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0], df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.30, 0.88), xytext = (0.30, 0.88), xycoords = 'axes fraction', fontsize = 40)
        else:
            ax_gs.annotate("$(x_B, Q^2) = ({:.2f}, {:.2f})$".format(df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0], df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.30, 0.88), xytext = (0.30, 0.88), xycoords = 'axes fraction', fontsize = 40)
        ax_gs.set_xlabel("$\phi$~($^\circ$)", fontsize = 40)
        ax_gs.set_xlim([0, 360])
        ax_gs.set_ylim([0, 0.14])
        ax_gs.set_xticks([0, 90, 180, 270, 360])
        ax_gs.set_xticklabels(['${}$'.format(n) for n  in [0, 90, 180, 270, 360]], fontsize = 40)
        ax_gs.set_yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])
        ax_gs.set_yticklabels(['${:.2f}$'.format(n) for n  in [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14]])
        ax_gs.set_frame_on(True)
        ax_gs.patch.set_linewidth(12)
        ax_gs.patch.set_edgecolor('tab:purple')
        ax_gs.tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 40)
        ax_gs.tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
        ax_gs.tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
        ax_gs.tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
    
    plt.savefig("addendum_v3/xsec_t{}_paper.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    # plt.close()

# stamp plot for a talk
xB_panes = 7
Q2_panes = 7
display_factor = 1

models = ["pi0", "bh", "km15", "vgg", "global1"]
# label_scheme = [r"$\mathrm{Global1~Norm.}$", r"$\mathrm{Global2~Norm.}$", "$\pi^0~\mathrm{Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{VGG~Norm.}$"]
label_scheme = ["$\mathrm{Measurement~(This~work)}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
color_scheme = ['k', 'red', 'cyan', 'tab:orange', 'tab:blue']

for t_binnum in [3]:
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (60, 30))
    gs = GridSpec(Q2_panes, xB_panes, figure=fig)
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if not df_this_bin.active_bin_nominal.sum():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            
            t_avgs.append(df_this_bin.t_avg_this_point.unique()[0])
            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]
            if integrated_binnum in [54, 106, 126]:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("$({:.2f}, {:.2f})$".format(df_this_bin.xB_avg_this_point.unique()[0], df_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.3, 0.85), xytext = (0.3, 0.85), xycoords = 'axes fraction', fontsize = 40)
                    
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=3, label = r"$\mathrm{Theory~(KM15)}$")
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].phi_display), df_display.loc[(df_display.integratedbin_display == integrated_binnum) & (df_display.phi_display >0.05) & (df_display.phi_display<2*np.pi-0.05)].xsec_VGG_display, color = 'tab:orange', lw=3, label = r"$\mathrm{Theory~(VGG)}$")
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=3, label = r"$\mathrm{Theory~(BH)}$")

            df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal"] == 1), :]
            weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
            weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
            weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
            for i, model in enumerate(models):
                if i == 0:
                    dots = axs[Q2_panes - Q2_binnum - 1, xB_binnum].scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                    axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3)#, label = r'$\mathrm{Syst.~Uncertainty}$')
                else:
                    # axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
                    pass

            if integrated_binnum == 89:
                # print(Q2_panes, xB_binnum)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].patch.set_linewidth(12)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].patch.set_edgecolor('tab:purple')
            
            if not labeled:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels() 
                neworders = [3, 0, 1, 2]
                handles = [handles[i] for i in neworders]
                labels = [labels[i] for i in neworders]
                labeled = 1
            
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_yscale('log')
            ymin = np.min(df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display)
            ymin = np.log10(ymin)
            ymax = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal == 1) & (df_summary_table_rebinned.phi_avg_this_point > 30)  & (df_summary_table_rebinned.phi_avg_this_point < 330), ["xsec_exp_pi0_eff_corrected_bkg_merging", "xsec_exp_pi0_eff_corrected_bkg_merging_syst_err", "xsec_exp_pi0_eff_corrected_bkg_merging_stat_err"]].sum(axis = 1).max() 
            ymax = np.log10(ymax)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 12+1), minor = True)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)", fontsize = 40 )
            if integrated_binnum in [4, 12, 32, 60, 84, 111, 130]:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 4+1))
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0, 360, 4+1), ['']*5)

            adjustment_min = 0.2
            adjustment_max = 0.2
            if np.abs(ymin)%1 < 0.2:
                adjustment_min = 0.5
            if np.abs(ymax)%1 < 0.2:
                adjustment_max = 0.5
            if np.floor(ymax) - np.ceil(ymin) < 1  :
                if np.ceil(ymax) - ymax > ymin - np.floor(ymin):
                    ymin = np.floor(ymin)- 0.005
                else:
                    ymax = np.ceil(ymax) + 0.005
            if (xB_binnum == 6) & (Q2_binnum == 4):
                ymax = ymax + .2
            if (xB_binnum == 5) & (Q2_binnum == 3):
                ymax = ymax + .4
            if (xB_binnum == 5) & (Q2_binnum == 5):
                ymin = ymin + .1
            if (xB_binnum == 4) & (Q2_binnum == 3):
                ymax = ymax + .4
            if (xB_binnum == 6) & (Q2_binnum == 6):
                ymin = ymin + .2
            if (xB_binnum == 1) & (Q2_binnum == 0):
                ymin = ymin - .4
                ymax = ymax + .3
            if (xB_binnum == 4) & (Q2_binnum == 4):
                ymax = ymax - .3
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([10**(ymin-adjustment_min), 10**(ymax+adjustment_max)])
            ticks = np.logspace(np.ceil(ymin-adjustment_min), np.floor(ymax+adjustment_max), int( np.floor(ymax+adjustment_max)-np.ceil(ymin-adjustment_min)+1))
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_yticks(ticks, log_ticker(ticks))
            
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 40)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.73, 0.35), fontsize = 40, title = r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$" + "\n" + r"$\langle |t| \rangle={:.3f}$".format(np.mean(t_avgs))+ r"$~\mathrm{GeV}^2$"+ "\n" + r"$(\langle x_B \rangle, \langle Q^2 \rangle)~\mathrm{is~annotated~in~each~panel}$", title_fontsize = 40, alignment = 'left', ncol = 1, markerscale = 3, framealpha = 1, edgecolor = 'k')
    # axs[0, 0].annotate(xy = (0.0, 0.4), xytext = (0.0, 0.4), text = r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$" , xycoords = "axes fraction", fontsize = 60)
    plt.subplots_adjust(wspace = 0.2 , hspace = 0.0 )

    ax_gs = plt.subplot(gs.new_subplotspec((0, 0), rowspan=3, colspan=2))
    integrated_binnum = 89
    df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_nominal==1), :]
    
    weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
    weights_stat_err      = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err
    weights               = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging
    weights_syst_err   = df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_syst_err
    
    # ax_gs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Exp.}$", histtype = 'step', color = 'k')
    # ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = "$\mathrm{Measurement~(This~work)}$")
    
    for i, model in enumerate(models):
        df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.loc[:, "active_bin_nominal".format(model)] == 1), :]
        if i == 0:
            dots = plt.scatter(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], color = color_scheme[i], zorder = 10)
            ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], yerr =df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging_stat_err".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
            ax_gs.fill_between(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights+weights_syst_err, weights-weights_syst_err, color = 'k', alpha = 0.3, label = r'$\mathrm{Syst.~Uncertainty}$')
        else:
            # ax_gs.errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, df_summary_table_rebinned_this_bin.loc[:, "xsec_exp_{}_eff_corrected_bkg_merging".format(model)], ls = '', marker = 'o', label = label_scheme[i], color = color_scheme[i], zorder = 5)
            pass
    
    ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_BH_KM15_display, color = 'tab:red', lw=4, ls = '--', label = r"$\mathrm{Theory~(BH)}$")
    ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_KM15_display, color = 'cyan', lw=4, ls = '--', label = r"$\mathrm{Theory~(KM)}$")
    ax_gs.plot(np.degrees(df_display.loc[(df_display.integratedbin_display == integrated_binnum)].phi_display), df_display.loc[df_display.integratedbin_display == integrated_binnum].xsec_VGG_display, color = 'tab:orange', lw=4, ls = '--', label = r"$\mathrm{Theory~(VGG)}$")
    
    
    xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
    # print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)
    
    xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]
    # print(xBmean, Q2mean, t1mean)
    
    
    # ax_gs.set_ylabel(r"$\frac{d\sigma_{ep\rightarrow e'p'\gamma}}{dx_B dQ^2 d|t| d\phi}~~\mathrm{[nb/GeV^4]}\\$")
    ax_gs.annotate("$(x_B, Q^2) = ({:.2f}, {:.2f})$".format(df_summary_table_rebinned_this_bin.xB_avg_this_point.unique()[0], df_summary_table_rebinned_this_bin.Q2_avg_this_point.unique()[0]), xy = (0.30, 0.88), xytext = (0.30, 0.88), xycoords = 'axes fraction', fontsize = 40)
    ax_gs.set_xlabel("$\phi$~($^\circ$)", fontsize = 40)
    ax_gs.set_xlim([0, 360])
    ax_gs.set_ylim([0, 0.14])
    ax_gs.set_xticks([0, 90, 180, 270, 360])
    ax_gs.set_xticklabels(['${}$'.format(n) for n  in [0, 90, 180, 270, 360]], fontsize = 40)
    ax_gs.set_yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])
    ax_gs.set_yticklabels(['${:.2f}$'.format(n) for n  in [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14]])
    ax_gs.set_frame_on(True)
    ax_gs.patch.set_linewidth(12)
    ax_gs.patch.set_edgecolor('tab:purple')
    ax_gs.tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 40)
    ax_gs.tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
    ax_gs.tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 40)
    ax_gs.tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
    
    plt.savefig("addendum_v3/xsec_t{}_presentation.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    # plt.close()
'''

'''
# AN
exp_pi0_fall2018_inb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("exp_pi0_fall2018_inb"))
sim_pi0_fall2018_inb_1 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_pi0_fall2018_inb_1"))
sim_pi0_fall2018_inb_2 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_pi0_fall2018_inb_2"))
sim_pi0_fall2018_inb_3 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_pi0_fall2018_inb_3"))
sim_pi0_fall2018_inb   = pd.concat([sim_pi0_fall2018_inb_1, sim_pi0_fall2018_inb_2])

sim_pi0_fall2018_inb.loc[:, "weights"]  = pi0_sigma_inb_in_nb * luminosity_inb * survival_rate_inb/10**8/2

exp_pi0_fall2018_inb.loc[:, "eff_bkg_merging_inb"] = 0
sim_pi0_fall2018_inb.loc[:, "eff_bkg_merging_inb"] = 0
for integrated_binnum in df_summary_table_rebinned.integrated_binnum.unique():
    for phi_binnum in range(24):
        if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy()):
            exp_pi0_fall2018_inb.loc[(exp_pi0_fall2018_inb.integrated_binnum == integrated_binnum) & (exp_pi0_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_pi0_fall2018_inb.loc[(sim_pi0_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_pi0_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
        else:
            exp_pi0_fall2018_inb.loc[(exp_pi0_fall2018_inb.integrated_binnum == integrated_binnum) & (exp_pi0_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0.7#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_pi0_fall2018_inb.loc[(sim_pi0_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_pi0_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0.7#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()


exp_pi0_fall2018_inb = exp_pi0_fall2018_inb.loc[exp_pi0_fall2018_inb.eff_bkg_merging_inb>0, :]
sim_pi0_fall2018_inb = sim_pi0_fall2018_inb.loc[sim_pi0_fall2018_inb.eff_bkg_merging_inb>0, :]


exp_pi0_fall2018_outb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("exp_pi0_fall2018_outb"))
sim_pi0_fall2018_outb_1 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_pi0_fall2018_outb_1"))
sim_pi0_fall2018_outb_2 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_pi0_fall2018_outb_2"))
sim_pi0_fall2018_outb_3 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_pi0_fall2018_outb_3"))
sim_pi0_fall2018_outb   = pd.concat([sim_pi0_fall2018_outb_2, sim_pi0_fall2018_outb_3])

sim_pi0_fall2018_outb.loc[:, "weights"]  = pi0_sigma_outb_in_nb * luminosity_outb * survival_rate_outb/10**8/2

exp_pi0_fall2018_outb.loc[:, "eff_bkg_merging_outb"] = 0
sim_pi0_fall2018_outb.loc[:, "eff_bkg_merging_outb"] = 0
for integrated_binnum in df_summary_table_rebinned.integrated_binnum.unique():
    for phi_binnum in range(24):
        if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy()):
            exp_pi0_fall2018_outb.loc[(exp_pi0_fall2018_outb.integrated_binnum == integrated_binnum) & (exp_pi0_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_pi0_fall2018_outb.loc[(sim_pi0_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_pi0_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
        else:
            exp_pi0_fall2018_outb.loc[(exp_pi0_fall2018_outb.integrated_binnum == integrated_binnum) & (exp_pi0_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0.7#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_pi0_fall2018_outb.loc[(sim_pi0_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_pi0_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0.7#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()


exp_pi0_fall2018_outb = exp_pi0_fall2018_outb.loc[exp_pi0_fall2018_outb.eff_bkg_merging_outb>0, :]
sim_pi0_fall2018_outb = sim_pi0_fall2018_outb.loc[sim_pi0_fall2018_outb.eff_bkg_merging_outb>0, :]

exp_dvcs_fall2018_inb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("exp_dvcs_fall2018_inb"))
sim_km15_fall2018_inb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_km15_fall2018_inb"))
sim_vgg_fall2018_inb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_vgg_fall2018_inb"))
sim_bh_fall2018_inb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bh_fall2018_inb"))
sim_bkg_fall2018_inb_1 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bkg_fall2018_inb_1"))
sim_bkg_fall2018_inb_2 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bkg_fall2018_inb_2"))
sim_bkg_fall2018_inb_3 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bkg_fall2018_inb_3"))
sim_bkg_fall2018_inb   = pd.concat([sim_bkg_fall2018_inb_1, sim_bkg_fall2018_inb_2])


exp_dvcs_fall2018_inb.loc[:, "eff_bkg_merging_inb"] = 0
sim_km15_fall2018_inb.loc[:, "eff_bkg_merging_inb"] = 0
sim_vgg_fall2018_inb.loc[:, "eff_vgg_merging_inb"] = 0
sim_bh_fall2018_inb.loc[:, "eff_bh_merging_inb"] = 0
sim_bkg_fall2018_inb.loc[:, "eff_bkg_merging_inb"] = 0
for integrated_binnum in df_summary_table_rebinned.integrated_binnum.unique():
    for phi_binnum in range(24):
        if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy()):
            exp_dvcs_fall2018_inb.loc[(exp_dvcs_fall2018_inb.integrated_binnum == integrated_binnum) & (exp_dvcs_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_km15_fall2018_inb.loc[(sim_km15_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_km15_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_vgg_fall2018_inb.loc[(sim_vgg_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_vgg_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"]    = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_bh_fall2018_inb.loc[(sim_bh_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_bh_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"]    = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_bkg_fall2018_inb .loc[(sim_bkg_fall2018_inb .integrated_binnum == integrated_binnum) & (sim_bkg_fall2018_inb .phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
        else:
            exp_dvcs_fall2018_inb.loc[(exp_dvcs_fall2018_inb.integrated_binnum == integrated_binnum) & (exp_dvcs_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_km15_fall2018_inb.loc[(sim_km15_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_km15_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_bkg_fall2018_inb .loc[(sim_bkg_fall2018_inb .integrated_binnum == integrated_binnum) & (sim_bkg_fall2018_inb .phi_binnum == phi_binnum), "eff_bkg_merging_inb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_vgg_fall2018_inb.loc[(sim_vgg_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_vgg_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"]    = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()
            sim_bh_fall2018_inb.loc[(sim_bh_fall2018_inb.integrated_binnum == integrated_binnum) & (sim_bh_fall2018_inb.phi_binnum == phi_binnum), "eff_bkg_merging_inb"]    = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_inb"].to_numpy().mean()


exp_dvcs_fall2018_inb = exp_dvcs_fall2018_inb.loc[exp_dvcs_fall2018_inb.eff_bkg_merging_inb>0, :]
sim_km15_fall2018_inb = sim_km15_fall2018_inb.loc[sim_km15_fall2018_inb.eff_bkg_merging_inb>0, :]
sim_vgg_fall2018_inb = sim_vgg_fall2018_inb.loc[sim_vgg_fall2018_inb.eff_bkg_merging_inb>0, :]
sim_bh_fall2018_inb = sim_bh_fall2018_inb.loc[sim_bh_fall2018_inb.eff_bkg_merging_inb>0, :]
sim_bkg_fall2018_inb  = sim_bkg_fall2018_inb .loc[sim_bkg_fall2018_inb. eff_bkg_merging_inb>0, :]


exp_dvcs_fall2018_outb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("exp_dvcs_fall2018_outb"))
sim_km15_fall2018_outb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_km15_fall2018_outb"))
sim_vgg_fall2018_outb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_vgg_fall2018_outb"))
sim_bh_fall2018_outb = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bh_fall2018_outb"))
# sim_bkg_fall2018_outb_1 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bkg_fall2018_outb_1"))
sim_bkg_fall2018_outb_2 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bkg_fall2018_outb_2"))
sim_bkg_fall2018_outb_3 = pd.read_pickle("impact_study_dec2024/{}.pkl".format("sim_bkg_fall2018_outb_3"))
sim_bkg_fall2018_outb   = pd.concat([sim_bkg_fall2018_outb_2, sim_bkg_fall2018_outb_3])


exp_dvcs_fall2018_outb.loc[:, "eff_bkg_merging_outb"] = 0
sim_km15_fall2018_outb.loc[:, "eff_bkg_merging_outb"] = 0
sim_vgg_fall2018_outb.loc[:, "eff_vgg_merging_outb"] = 0
sim_bh_fall2018_outb.loc[:, "eff_bh_merging_outb"] = 0
sim_bkg_fall2018_outb.loc[:, "eff_bkg_merging_outb"] = 0
for integrated_binnum in df_summary_table_rebinned.integrated_binnum.unique():
    for phi_binnum in range(24):
        if len(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy()):
            exp_dvcs_fall2018_outb.loc[(exp_dvcs_fall2018_outb.integrated_binnum == integrated_binnum) & (exp_dvcs_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_km15_fall2018_outb.loc[(sim_km15_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_km15_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_vgg_fall2018_outb.loc[(sim_vgg_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_vgg_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"]    = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_bh_fall2018_outb.loc[(sim_bh_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_bh_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"]    = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_bkg_fall2018_outb .loc[(sim_bkg_fall2018_outb .integrated_binnum == integrated_binnum) & (sim_bkg_fall2018_outb .phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
        else:
            exp_dvcs_fall2018_outb.loc[(exp_dvcs_fall2018_outb.integrated_binnum == integrated_binnum) & (exp_dvcs_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_km15_fall2018_outb.loc[(sim_km15_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_km15_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_bkg_fall2018_outb .loc[(sim_bkg_fall2018_outb .integrated_binnum == integrated_binnum) & (sim_bkg_fall2018_outb .phi_binnum == phi_binnum), "eff_bkg_merging_outb"] = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_vgg_fall2018_outb.loc[(sim_vgg_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_vgg_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"]    = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()
            sim_bh_fall2018_outb.loc[(sim_bh_fall2018_outb.integrated_binnum == integrated_binnum) & (sim_bh_fall2018_outb.phi_binnum == phi_binnum), "eff_bkg_merging_outb"]    = 0#df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.phi_binnum == phi_binnum), "eff_bkg_merging_outb"].to_numpy().mean()


exp_dvcs_fall2018_outb = exp_dvcs_fall2018_outb.loc[exp_dvcs_fall2018_outb.eff_bkg_merging_outb>0, :]
sim_km15_fall2018_outb = sim_km15_fall2018_outb.loc[sim_km15_fall2018_outb.eff_bkg_merging_outb>0, :]
sim_vgg_fall2018_outb = sim_vgg_fall2018_outb.loc[sim_vgg_fall2018_outb.eff_bkg_merging_outb>0, :]
sim_bh_fall2018_outb = sim_bh_fall2018_outb.loc[sim_bh_fall2018_outb.eff_bkg_merging_outb>0, :]
sim_bkg_fall2018_outb  = sim_bkg_fall2018_outb .loc[sim_bkg_fall2018_outb. eff_bkg_merging_outb>0, :]

fig, axs = plt.subplots(4, 3, figsize = (15, 12))

variables = ["Ep", "Etheta", "Ephi", "Pp", "Ptheta", "Pphi", "Gp", "Gtheta", "Gphi", "Gp2", "Gtheta2", "Gphi2"]
labels    = ["$p_{e'}$", r"$\theta_{e'}$", r"$\phi_{e'}$", "$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", "$p_{\gamma1}$", r"$\theta_{\gamma1}$", r"$\phi_{\gamma1}$", "$p_{\gamma2}$", r"$\theta_{\gamma2}$", r"$\phi_{\gamma2}$"]
units     = ["$\mathrm{GeV}/c$", "$^{\circ}$", "$^{\circ}$"]*4
xmins      = [0 ,  5, -180,   0, 15, -180,  0,  0, -180, 0,  0, -180]
xmaxs      = [12, 35,  180, 1.2, 65,  180, 12, 35,  180, 4, 35,  180]
xticks     = [7 ,  7,    5,   7,  6,    5,  7,  8,    5, 5,  8,    5]


for i, j in itertools.product(range(4), range(3)):
    var   = variables[3*i + j]
    label = labels[3*i + j]
    unit  = units[3*i + j]
    xmin  = xmins[3*i + j]
    xmax  = xmaxs[3*i + j]
    xtick = xticks[3*i + j]
    hist_exp, bins    = np.histogram(exp_pi0_fall2018_inb.loc[:, var], bins = 100)
    hist_sim, _    = np.histogram(sim_pi0_fall2018_inb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_inb.weights)#, histtype = 'step')
    hist_sim2, _   = np.histogram(sim_pi0_fall2018_inb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_inb.weights * sim_pi0_fall2018_inb.eff_bkg_merging_inb)#, histtype = 'step')
    hist_sim3, _   = np.histogram(sim_pi0_fall2018_inb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_inb.weights * sim_pi0_fall2018_inb.eff_bkg_merging_inb * sim_pi0_fall2018_inb.efficiency_pi0)#, histtype = 'step')
    # axs[i, j].hist(bins[:-1], bins, weights = hist_pi0, histtype = 'step', color = '#004D40', label = "Experimental Data")
    axs[i, j].hist(bins[:-1], bins, weights = hist_exp, histtype = 'step', color = 'k', label = "Signal Yield (Exp.)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim, histtype = 'step', color = 'tab:red', label = "Raw Yield (Sim.,\n Efficiency not corrected)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim2, histtype = 'step', color = 'tab:purple', label = "Signal Yield (Sim.,\n Background merging only)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim3, histtype = 'step', color = 'cyan', label = "Signal Yield (Sim.)")

    axs[i, j].set_xlabel(label+ " [" + unit + "]")
    axs[i, j].set_xlim(xmin, xmax)
    axs[i, j].set_xticks(np.linspace(xmin, xmax, xtick))
    # axs[i, j].set_yscale('log')
axs[1, 1].set_ylim(bottom = 1)
legends = axs[i, j].get_legend_handles_labels()
plt.figlegend(*legends, loc = 'upper left', bbox_to_anchor = (1.0, 0.6), fontsize = 20)

plt.suptitle("DV$\pi^0$P RG-A Fall2018 Inbending", fontsize = 30)
plt.tight_layout()
plt.savefig("addendum_v3/kinematics_inbending.from_script.pdf", bbox_inches = 'tight')
plt.close()
fig, axs = plt.subplots(2, 2, figsize = (15, 10))

variables = ["xB", "Q2", "t1", "phi1"]
labels    = ["$x_B$", "$Q^2$", "$|t|$", "$\phi$"]
units     = ["", "$\mathrm{GeV^2}/c^2$", "$\mathrm{GeV}^2$", "$^{\circ}$"]
xmins      = [0,   1, 0,   0]
xmaxs      = [0.6, 6, 1, 360]
xticks     = [7,   6, 6,   5]

for i, j in itertools.product(range(2), range(2)):
    var   = variables[2*i + j]
    label = labels[2*i + j]
    unit  = units[2*i + j]
    xmin  = xmins[2*i + j]
    xmax  = xmaxs[2*i + j]
    xtick = xticks[2*i + j]
    hist_exp, bins    = np.histogram(exp_pi0_fall2018_inb.loc[:, var], bins = 100)
    hist_sim, _    = np.histogram(sim_pi0_fall2018_inb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_inb.weights)#, histtype = 'step')
    hist_sim2, _   = np.histogram(sim_pi0_fall2018_inb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_inb.weights * sim_pi0_fall2018_inb.eff_bkg_merging_inb)#, histtype = 'step')
    hist_sim3, _   = np.histogram(sim_pi0_fall2018_inb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_inb.weights * sim_pi0_fall2018_inb.eff_bkg_merging_inb * sim_pi0_fall2018_inb.efficiency_pi0)#, histtype = 'step')
    # axs[i, j].hist(bins[:-1], bins, weights = hist_pi0, histtype = 'step', color = '#004D40', label = "Experimental Data")
    axs[i, j].hist(bins[:-1], bins, weights = hist_exp, histtype = 'step', color = 'k', label = "Signal Yield (Exp.)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim, histtype = 'step', color = 'tab:red', label = "Raw Yield (Sim.,\n Efficiency not corrected)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim2, histtype = 'step', color = 'tab:purple', label = "Signal Yield (Sim.,\n Background merging only)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim3, histtype = 'step', color = 'cyan', label = "Signal Yield (Sim.)")
    if unit:
        axs[i, j].set_xlabel(label+ " [" + unit + "]")
    else:
        axs[i, j].set_xlabel(label)
    axs[i, j].set_xlim(xmin, xmax)
    axs[i, j].set_xticks(np.linspace(xmin, xmax, xtick))
legends = axs[i, j].get_legend_handles_labels()
plt.figlegend(*legends, loc = 'upper left', bbox_to_anchor = (1.0, 0.6), fontsize = 20)
    
plt.suptitle("DV$\pi^0$P RG-A Fall2018 Inbending", fontsize = 30)
plt.tight_layout()
plt.savefig("addendum_v3/binning_vars_inbending.from_script.pdf", bbox_inches = 'tight')
plt.close()

fig, axs = plt.subplots(4, 3, figsize = (15, 12))

variables = ["Ep", "Etheta", "Ephi", "Pp", "Ptheta", "Pphi", "Gp", "Gtheta", "Gphi", "Gp2", "Gtheta2", "Gphi2"]
labels    = ["$p_{e'}$", r"$\theta_{e'}$", r"$\phi_{e'}$", "$p_{p'}$", r"$\theta_{p'}$", r"$\phi_{p'}$", "$p_{\gamma1}$", r"$\theta_{\gamma1}$", r"$\phi_{\gamma1}$", "$p_{\gamma2}$", r"$\theta_{\gamma2}$", r"$\phi_{\gamma2}$"]
units     = ["$\mathrm{GeV}/c$", "$^{\circ}$", "$^{\circ}$"]*4
xmins      = [0 ,  5, -180,   0, 15, -180,  0,  0, -180, 0,  0, -180]
xmaxs      = [12, 35,  180, 1.2, 65,  180, 12, 35,  180, 4, 35,  180]
xticks     = [7 ,  7,    5,   7,  6,    5,  7,  8,    5, 5,  8,    5]


for i, j in itertools.product(range(4), range(3)):
    var   = variables[3*i + j]
    label = labels[3*i + j]
    unit  = units[3*i + j]
    xmin  = xmins[3*i + j]
    xmax  = xmaxs[3*i + j]
    xtick = xticks[3*i + j]
    hist_exp, bins    = np.histogram(exp_pi0_fall2018_outb.loc[:, var], bins = 100)
    hist_sim, _    = np.histogram(sim_pi0_fall2018_outb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_outb.weights)#, histtype = 'step')
    hist_sim2, _   = np.histogram(sim_pi0_fall2018_outb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_outb.weights * sim_pi0_fall2018_outb.eff_bkg_merging_outb)#, histtype = 'step')
    hist_sim3, _   = np.histogram(sim_pi0_fall2018_outb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_outb.weights * sim_pi0_fall2018_outb.eff_bkg_merging_outb * sim_pi0_fall2018_outb.efficiency_pi0)#, histtype = 'step')
    # axs[i, j].hist(bins[:-1], bins, weights = hist_pi0, histtype = 'step', color = '#004D40', label = "Experimental Data")
    axs[i, j].hist(bins[:-1], bins, weights = hist_exp, histtype = 'step', color = 'k', label = "Signal Yield (Exp.)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim, histtype = 'step', color = 'tab:red', label = "Raw Yield (Sim.,\n Efficiency not corrected)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim2, histtype = 'step', color = 'tab:purple', label = "Signal Yield (Sim.,\n Background merging only)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim3, histtype = 'step', color = 'cyan', label = "Signal Yield (Sim.)")

    axs[i, j].set_xlabel(label+ " [" + unit + "]")
    axs[i, j].set_xlim(xmin, xmax)
    axs[i, j].set_xticks(np.linspace(xmin, xmax, xtick))
    # axs[i, j].set_yscale('log')
axs[1, 1].set_ylim(bottom = 1)
legends = axs[i, j].get_legend_handles_labels()
plt.figlegend(*legends, loc = 'upper left', bbox_to_anchor = (1.0, 0.6), fontsize = 20)

plt.suptitle("DV$\pi^0$P RG-A Fall2018 Outbending", fontsize = 30)
plt.tight_layout()
plt.savefig("addendum_v3/kinematics_outbending.from_script.pdf", bbox_inches = 'tight')
plt.close()
fig, axs = plt.subplots(2, 2, figsize = (15, 10))

variables = ["xB", "Q2", "t1", "phi1"]
labels    = ["$x_B$", "$Q^2$", "$|t|$", "$\phi$"]
units     = ["", "$\mathrm{GeV^2}/c^2$", "$\mathrm{GeV}^2$", "$^{\circ}$"]
xmins      = [0,   1, 0,   0]
xmaxs      = [0.6, 6, 1, 360]
xticks     = [7,   6, 6,   5]

for i, j in itertools.product(range(2), range(2)):
    var   = variables[2*i + j]
    label = labels[2*i + j]
    unit  = units[2*i + j]
    xmin  = xmins[2*i + j]
    xmax  = xmaxs[2*i + j]
    xtick = xticks[2*i + j]
    hist_exp, bins    = np.histogram(exp_pi0_fall2018_outb.loc[:, var], bins = 100)
    hist_sim, _    = np.histogram(sim_pi0_fall2018_outb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_outb.weights)#, histtype = 'step')
    hist_sim2, _   = np.histogram(sim_pi0_fall2018_outb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_outb.weights * sim_pi0_fall2018_outb.eff_bkg_merging_outb)#, histtype = 'step')
    hist_sim3, _   = np.histogram(sim_pi0_fall2018_outb.loc[:, var], bins = bins, weights = sim_pi0_fall2018_outb.weights * sim_pi0_fall2018_outb.eff_bkg_merging_outb * sim_pi0_fall2018_outb.efficiency_pi0)#, histtype = 'step')
    # axs[i, j].hist(bins[:-1], bins, weights = hist_pi0, histtype = 'step', color = '#004D40', label = "Experimental Data")
    axs[i, j].hist(bins[:-1], bins, weights = hist_exp, histtype = 'step', color = 'k', label = "Signal Yield (Exp.)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim, histtype = 'step', color = 'tab:red', label = "Raw Yield (Sim.,\n Efficiency not corrected)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim2, histtype = 'step', color = 'tab:purple', label = "Signal Yield (Sim.,\n Background merging only)")
    axs[i, j].hist(bins[:-1], bins, weights = hist_sim3, histtype = 'step', color = 'cyan', label = "Signal Yield (Sim.)")
    if unit:
        axs[i, j].set_xlabel(label+ " [" + unit + "]")
    else:
        axs[i, j].set_xlabel(label)
    axs[i, j].set_xlim(xmin, xmax)
    axs[i, j].set_xticks(np.linspace(xmin, xmax, xtick))
legends = axs[i, j].get_legend_handles_labels()
plt.figlegend(*legends, loc = 'upper left', bbox_to_anchor = (1.0, 0.6), fontsize = 20)
    
plt.suptitle("DV$\pi^0$P RG-A Fall2018 Outbending", fontsize = 30)
plt.tight_layout()
plt.savefig("addendum_v3/binning_vars_outbending.from_script.pdf", bbox_inches = 'tight')
plt.close()
'''


# df_summary_table_rebinned.to_pickle("addendum_v3/df_summary_table_rebinned_approved.pkl")
'''
# # AN

fig, axs = plt.subplots(1, 1, figsize = (10, 6))
integrated_binnum = 88

# models = ["contamination_inb_pi0_eff_corrected_bkg_merging", "contamination_inb_bh_eff_corrected_bkg_merging", "contamination_inb_km15_eff_corrected_bkg_merging", "contamination_inb_vgg_eff_corrected_bkg_merging", "contamination_inb_global1_eff_corrected_bkg_merging", ]
# labels = [r"$\pi^0~\mathrm{Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
# colors = ['k', 'r', 'cyan', 'tab:orange', 'tab:blue', ]

models = ["contamination_inb_pi0_eff_corrected_bkg_merging", "contamination_inb_bkg_merging"]
labels = [r"$\mathrm{With}~\pi^0~\mathrm{normalization}$", r"$\mathrm{Without~normalization}$"]
colors = ['k', 'tab:pink', 'cyan', 'tab:orange', 'tab:blue', ]

for i, model in enumerate(models):
    cont = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), "{}".format(model)]
    cont_stat_err = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), "{}_stat_err".format(model)]
    phi_avg = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), "phi_avg_this_point"]
    plt.errorbar(phi_avg, cont, cont_stat_err, ls = '--', marker = 'o', label = labels[i], color = colors[i])


xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]
print(xBmean, Q2mean, t1mean)


axs.set_ylabel(r"$\mathrm{Contamination~Ratio}$")
axs.set_xlabel("$\phi$~($^\circ$)")
axs.set_xlim([0, 360])
axs.set_xticks([0, 90, 180, 270, 360])
# handles, labels = axs.get_legend_handles_labels()
# xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
# Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
# t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
# handles = [handles[-1], handles[0], handles[1], handles[2]]
# labels = [labels[-1], labels[0], labels[1], labels[2]]

xBheader = "$\\langle x_B \\rangle = {:.3f}$ \n".format(xBmean)
Q2header = "$\\langle Q^2 \\rangle = {:.3f}~\mathrm{{GeV}}^2/c^2$ \n".format(Q2mean)
t1header = "$\\langle |t| \\rangle~= {:.3f}~\mathrm{{GeV}}^2$".format(t1mean)
# axs.annotate(xy = (0.01, 0.7), xytext = (0.01, 0.7), text = xBheader + Q2header + t1header, xycoords = 'axes fraction', fontsize = 20)
plt.figlegend(loc = 'upper left', title = r"${}.~\mathrm{{Inbending}}$".format(integrated_binnum)+"\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))

plt.savefig("addendum_v3/contamination_inb_{}.new.from_script.pdf".format(integrated_binnum), bbox_inches = 'tight')

fig, axs = plt.subplots(1, 1, figsize = (10, 6))
integrated_binnum = 88

# models = ["contamination_outb_pi0_eff_corrected_bkg_merging", "contamination_outb_bh_eff_corrected_bkg_merging", "contamination_outb_km15_eff_corrected_bkg_merging", "contamination_outb_vgg_eff_corrected_bkg_merging", "contamination_outb_global1_eff_corrected_bkg_merging", ]
# labels = [r"$\pi^0~\mathrm{Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
# colors = ['k', 'r', 'cyan', 'tab:orange', 'tab:blue', ]

models = ["contamination_outb_pi0_eff_corrected_bkg_merging", "contamination_outb_bkg_merging"]
labels = [r"$\mathrm{With}~\pi^0~\mathrm{normalization}$", r"$\mathrm{Without~normalization}$"]
colors = ['k', 'tab:pink', 'cyan', 'tab:orange', 'tab:blue', ]

for i, model in enumerate(models):
    cont = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), "{}".format(model)]
    cont_stat_err = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), "{}_stat_err".format(model)]
    phi_avg = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), "phi_avg_this_point"]
    plt.errorbar(phi_avg, cont, cont_stat_err, ls = '--', marker = 'o', label = labels[i], color = colors[i])


xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

xBmean, Q2mean, t1mean = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xB_avg_this_point", "Q2_avg_this_point", "t_avg_this_point"]].to_numpy().T[:, 0]
print(xBmean, Q2mean, t1mean)


axs.set_ylabel(r"$\mathrm{Contamination~Ratio}$")
axs.set_xlabel("$\phi$~($^\circ$)")
axs.set_xlim([0, 360])
axs.set_xticks([0, 90, 180, 270, 360])
# handles, labels = axs.get_legend_handles_labels()
# xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
# Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
# t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
# handles = [handles[-1], handles[0], handles[1], handles[2]]
# labels = [labels[-1], labels[0], labels[1], labels[2]]

xBheader = "$\\langle x_B \\rangle = {:.3f}$ \n".format(xBmean)
Q2header = "$\\langle Q^2 \\rangle = {:.3f}~\mathrm{{GeV}}^2/c^2$ \n".format(Q2mean)
t1header = "$\\langle |t| \\rangle~= {:.3f}~\mathrm{{GeV}}^2$".format(t1mean)
# axs.annotate(xy = (0.01, 0.7), xytext = (0.01, 0.7), text = xBheader + Q2header + t1header, xycoords = 'axes fraction', fontsize = 20)
plt.figlegend(loc = 'upper left', title = r"${}.~\mathrm{{Outbending}}$".format(integrated_binnum)+"\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))


plt.savefig("addendum_v3/contamination_outb_{}.new.from_script.pdf".format(integrated_binnum), bbox_inches = 'tight')


fig, axs = plt.subplots(1, 2, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix) & (df_summary_table_rebinned_exp.integrated_binnum == integrated_binnum), :]

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.n_entry_FD
weights_stat_err    = np.sqrt(df_summary_table_rebinned_this_bin.n_entry_FD)
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~FD,~FD)}$", histtype = 'step', color = 'k')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
weights             = df_summary_table_rebinned_this_bin.n_entry_CD
weights_stat_err    = np.sqrt(df_summary_table_rebinned_this_bin.n_entry_CD)
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FD)}$", histtype = 'step', color = 'r')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
weights             = df_summary_table_rebinned_this_bin.n_entry_CDFT
weights_stat_err    = np.sqrt(df_summary_table_rebinned_this_bin.n_entry_CDFT)
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FT)}$", histtype = 'step', color = 'tab:blue')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:blue')
axs[0].annotate(xy = (180, 230), xytext = (200, 230), text = '$\mathrm{Inbending}$')

df_summary_table_rebinned_this_bin = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix) & (df_summary_table_rebinned_exp.integrated_binnum == integrated_binnum), :]


phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.n_entry_FD
weights_stat_err    = np.sqrt(df_summary_table_rebinned_this_bin.n_entry_FD)
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~FD,~FD)}$", histtype = 'step', color = 'k')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
weights             = df_summary_table_rebinned_this_bin.n_entry_CD
weights_stat_err    = np.sqrt(df_summary_table_rebinned_this_bin.n_entry_CD)
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FD)}$", histtype = 'step', color = 'r')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
weights             = df_summary_table_rebinned_this_bin.n_entry_CDFT
weights_stat_err    = np.sqrt(df_summary_table_rebinned_this_bin.n_entry_CDFT)
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FT)}$", histtype = 'step', color = 'tab:blue')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:blue')
axs[1].annotate(xy = (180, 230), xytext = (180, 230), text = '$\mathrm{Outbending}$')

# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (KM15)", histtype = 'step')
# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pureBH
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (BH)", histtype = 'step')
# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_vgg
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (VGG)", histtype = 'step')
# axs[0].set_yscale('log')
xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

for ax in axs:
    ax.set_xlabel("$\phi$~($^\circ$)")
    ax.set_xlim([0, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim([0, 250])
handles, labels = axs[0].get_legend_handles_labels()
xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(handles, labels, loc = 'upper left', title = "$88.$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/Raw_yield_topology.from_script.pdf", bbox_inches = 'tight')
# plt.close()

fig, axs = plt.subplots(1, 2, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, :]

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.epg_inb_exp
weights_stat_err    = df_summary_table_rebinned_this_bin.epg_inb_exp_stat_err
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Yield~(Exp.)}$", histtype = 'step', color = 'k')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')

weights             = df_summary_table_rebinned_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Background~Yield~(Exp.)}$", histtype = 'step', color = 'r')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
weights_min         = df_summary_table_rebinned_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging - df_summary_table_rebinned_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err
weights_max         = df_summary_table_rebinned_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging + df_summary_table_rebinned_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[0].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'r', alpha = 0.3)

axs[0].annotate(xy = (180, 230), xytext = (200, 230), text = '$\mathrm{Inbending}$')

df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, :]

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.epg_outb_exp
weights_stat_err    = df_summary_table_rebinned_this_bin.epg_outb_exp_stat_err
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Yield~(Exp.)}$", histtype = 'step', color = 'k')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')

weights             = df_summary_table_rebinned_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Background~Yield~(Exp.)}$", histtype = 'step', color = 'r')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
weights_min         = df_summary_table_rebinned_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging - df_summary_table_rebinned_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err
weights_max         = df_summary_table_rebinned_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging + df_summary_table_rebinned_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'r', alpha = 0.3)



axs[1].annotate(xy = (180, 230), xytext = (180, 230), text = '$\mathrm{Outbending}$')

# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (KM15)", histtype = 'step')
# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pureBH
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (BH)", histtype = 'step')
# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_vgg
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (VGG)", histtype = 'step')
# axs[0].set_yscale('log')
xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

for ax in axs:
    ax.set_xlabel("$\phi$~($^\circ$)")
    ax.set_xlim([0, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim([0, 250])
handles, labels = axs[0].get_legend_handles_labels()
xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(handles, labels, loc = 'upper left', title = "$88.$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/Background_yield.from_script.pdf", bbox_inches = 'tight')


fig, axs = plt.subplots(1, 2, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, :]

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Exp.)}$", histtype = 'step', color = 'k')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging * (1 - np.sqrt(df_summary_table_rebinned_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio **2 + 0.3**2 + 0.0476**2))
weights_max         = df_summary_table_rebinned_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging * (1 + np.sqrt(df_summary_table_rebinned_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio **2 + 0.3**2 + 0.0476**2))
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[0].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)
        
weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~KM15)}$", histtype = 'step', color = 'cyan')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging *0.7#- df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_syst_err
weights_max         = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging *1.3#+ df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

# axs[0].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~VGG)}$", histtype = 'step', color = 'tab:orange')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg *0.7#- df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err
weights_max         = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg *1.3#+ df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

# axs[0].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:orange', alpha = 0.5)

weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err
axs[0].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~BH)}$", histtype = 'step', color = 'tab:red')
axs[0].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH *0.7#- df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err
weights_max         = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH *1.3#+ df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

# axs[0].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)



axs[0].annotate(xy = (180, 230), xytext = (200, 230), text = '$\mathrm{Inbending}$')

df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, :]


phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Exp.)}$", histtype = 'step', color = 'k')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging * (1 - np.sqrt(df_summary_table_rebinned_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio **2 + 0.3**2 + 0.0476**2))
weights_max         = df_summary_table_rebinned_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging * (1 + np.sqrt(df_summary_table_rebinned_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio **2 + 0.3**2 + 0.0476**2))
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

weights             = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.)}$", histtype = 'step', color = 'cyan')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging *0.7#- df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_syst_err
weights_max         = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging *1.3#+ df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

# axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)

weights             = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~VGG)}$", histtype = 'step', color = 'tab:orange')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg *0.7#- df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err
weights_max         = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg *1.3#+ df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

# axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:orange', alpha = 0.5)

weights             = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH
weights_stat_err    = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err
axs[1].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~BH)}$", histtype = 'step', color = 'tab:red')
axs[1].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')
weights_min         = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH *0.7#- df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err
weights_max         = df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH *1.3#+ df_summary_table_rebinned_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

# axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)



axs[1].annotate(xy = (180, 230), xytext = (180, 230), text = '$\mathrm{Outbending}$')

# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (KM15)", histtype = 'step')
# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (BH)", histtype = 'step')
# weights             = df_summary_table_rebinned_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg
# axs[0].hist(phi_binnum_rebinned[:-1], phi_binnum_rebinned, weights = weights, label = "Simlation (VGG)", histtype = 'step')
# axs[0].set_yscale('log')
xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

for ax in axs:
    ax.set_xlabel("$\phi$~($^\circ$)")
    ax.set_xlim([0, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim([0, 250])
handles, labels = axs[0].get_legend_handles_labels()
xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(handles, labels, loc = 'upper left', title = "$88.$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/Signal_yield.from_script.pdf", bbox_inches = 'tight')

fig, axs = plt.subplots(1, 1, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, :]

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Exp.}$", histtype = 'step', color = 'k')
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
weights_min         = df_summary_table_rebinned_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging * ( 1 - np.sqrt(0.3**2 + 0.0476**2 + df_summary_table_rebinned_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2))
weights_max         = df_summary_table_rebinned_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging * ( 1 + np.sqrt(0.3**2 + 0.0476**2 + df_summary_table_rebinned_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2))
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs.fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)
        
weights             = df_summary_table_rebinned_this_bin.gen_sim
weights_stat_err    = df_summary_table_rebinned_this_bin.gen_sim_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Sim.,~KM15}$", histtype = 'step', color = 'cyan', lw = 3, alpha = 0.8)
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')

weights             = df_summary_table_rebinned_this_bin.gen_sim_vgg
weights_stat_err    = df_summary_table_rebinned_this_bin.gen_sim_vgg_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Sim.,~VGG}$", histtype = 'step', color = 'tab:orange', lw = 3, alpha = 0.8)
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')

weights             = df_summary_table_rebinned_this_bin.gen_sim_pureBH
weights_stat_err    = df_summary_table_rebinned_this_bin.gen_sim_pureBH_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Sim.,~BH}$", histtype = 'step', color = 'tab:red', lw = 3, alpha = 0.8)
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')

xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

axs.set_xlabel("$\phi$~($^\circ$)")
axs.set_xlim([0, 360])
axs.set_xticks([0, 90, 180, 270, 360])
handles, labels = axs.get_legend_handles_labels()
xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(handles, labels, loc = 'upper left', title = "$88.$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.ylabel(r"$\mathrm{Acceptance~Corrected}$"+"\n"+r"$\mathrm{Signal~Yield}$")

plt.savefig("addendum_v3/Acceptance_corrected_yields.from_script.pdf", bbox_inches = 'tight')

fig, axs = plt.subplots(1, 1, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, :]

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{KM15}$", histtype = 'step', color = 'cyan', zorder = 10)
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', zorder = 5, marker = 'o', mfc = 'cyan')
weights_min         = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging * (1 - np.sqrt(df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio ** 2 + 0.3**2 + 0.0476**2))
weights_max         = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging * (1 + np.sqrt(df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio ** 2 + 0.3**2 + 0.0476**2))
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs.fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = .3, zorder = 1)


phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg
weights_stat_err    = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{VGG}$", histtype = 'step', color = 'tab:orange', zorder = 10)
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange', zorder = 5)
weights_min         = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg #- df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err
weights_max         = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg #+ df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs.fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:orange', alpha = 0.5)

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH
weights_stat_err    = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err
axs.hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{BH}$", histtype = 'step', color = 'tab:red', zorder = 110)
axs.errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red', zorder = 10)
weights_min         = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH #- df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err
weights_max         = df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH #+ df_summary_table_rebinned_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs.fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)



axs.set_xlabel("$\phi$~($^\circ$)")
axs.set_xlim([0, 360])
axs.set_xticks([0, 90, 180, 270, 360])
handles, labels = axs.get_legend_handles_labels()
xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.ylabel("$\mathrm{Acceptance}$")
plt.figlegend(handles, labels, loc = 'upper left', title = "$88.$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/Acceptance.from_script.pdf", bbox_inches = 'tight')


fig, axs = plt.subplots(1, 3, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb == 1), :]

# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

weights             = df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err
axs[0].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}}$')


# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0)
weights_stat_err    = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) * df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio
axs[1].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'tab:red', marker = 'o', label = r'$\epsilon_{\mathrm{additional}}$')

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

weights_min         = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) *0.7#- inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0_syst_err)
weights_max         = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) *1.3#+ inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0_syst_err)
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)

# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
weights_stat_err    = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging * np.sqrt(df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
axs[2].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'tab:cyan', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}} \times \epsilon_{\mathrm{additional}}$')#r'$\epsilon_{\mathrm{additional~merging}}} \times \epsilon_{\mathrm{additional~merging}}$')

weights_min         = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging *0.7#- df_summary_table_rebinned_this_bin.efficiency_inb_syst_err
weights_max         = inverseHist(df_summary_table_rebinned_this_bin.normalization_inb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging *1.3#+ df_summary_table_rebinned_this_bin.efficiency_inb_syst_err
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[2].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


for ax in axs:
    ax.set_xlabel("$\phi$~($^\circ$)")
    ax.set_xlim([0, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim([0, 1.4])
    ax.tick_params(left = True, right = True, axis = 'y', which = 'both', direction = 'inout', length = 10)
axs[0].set_xticks([0, 90, 180, 270, 360], ['$0$', '$90$', '$180$', '$270$', '$360~0$'])
axs[1].set_xticks([0, 90, 180, 270, 360], ['', '$90$', '$180$', '$270$', '$360~0$'])
axs[2].set_xticks([0, 90, 180, 270, 360], ['', '$90$', '$180$', '$270$', '$360$'])

axs[1].set_yticks(np.linspace(0, 1.4, 8), ['']*8)
axs[2].set_yticks(np.linspace(0, 1.4, 8), ['']*8)
# axs[2].axhline(0.25, color = 'k', ls = '--')
ax.tick_params(left = True, right = True, axis = 'y', which = 'both', direction = 'inout', length = 10)
plt.subplots_adjust(wspace = 0)
    
xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(loc = 'upper left', title = "$88.~\mathrm{Inbending~Data~Set}$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/efficiency_inbending.from_script.pdf", bbox_inches = 'tight')

fig, axs = plt.subplots(1, 3, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb == 1), :]

# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

weights             = df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err
axs[0].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label =r'$\epsilon_{\mathrm{bkg.~merging}}$')# \times \epsilon_{\mathrm{additional}}$')


# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0)
weights_stat_err    = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) * df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio
axs[1].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'tab:red', marker = 'o', label = r'$\epsilon_{\mathrm{additional}}$')

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

weights_min         = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) *0.7#- inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0_syst_err)
weights_max         = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) *1.3#+ inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0_syst_err)
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)

# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
weights_stat_err    = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging * np.sqrt(df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
axs[2].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'tab:cyan', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}}\times\epsilon_{\mathrm{additional}}$')

weights_min         = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging *0.7#- inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0_syst_err)
weights_max         = inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging *1.3#+ inverseHist(df_summary_table_rebinned_this_bin.normalization_outb_pi0_syst_err)
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[2].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


for ax in axs:
    ax.set_xlabel("$\phi$~($^\circ$)")
    ax.set_xlim([0, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim([0, 1.4])
    ax.tick_params(left = True, right = True, axis = 'y', which = 'both', direction = 'inout', length = 10)
axs[0].set_xticks([0, 90, 180, 270, 360], ['$0$', '$90$', '$180$', '$270$', '$360~0$'])
axs[1].set_xticks([0, 90, 180, 270, 360], ['', '$90$', '$180$', '$270$', '$360~0$'])
axs[2].set_xticks([0, 90, 180, 270, 360], ['', '$90$', '$180$', '$270$', '$360$'])

axs[1].set_yticks(np.linspace(0, 1.4, 8), ['']*8)
axs[2].set_yticks(np.linspace(0, 1.4, 8), ['']*8)
# axs[2].axhline(0.25, color = 'k', ls = '--')
ax.tick_params(left = True, right = True, axis = 'y', which = 'both', direction = 'inout', length = 10)
plt.subplots_adjust(wspace = 0)
    
xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(loc = 'upper left', title = "$88.~\mathrm{Outbending~Data~Set}$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/efficiency_outbending.from_script.pdf", bbox_inches = 'tight')

fig, axs = plt.subplots(1, 3, figsize = (10, 6))
integrated_binnum = 88
df_summary_table_rebinned_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), :]

# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

weights             = df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging
weights_stat_err    = df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err
axs[0].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}}$')


# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0)
weights_stat_err    = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) * df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio
axs[1].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'tab:red', marker = 'o', label = r'$\epsilon_{\mathrm{additional}}$')

phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
phi_bin_rebinned    = 15*phi_binnum_rebinned
phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

weights_min         = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) *0.7#- inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0_syst_err)
weights_max         = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) *1.3#+ inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0_syst_err)
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[1].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)

# phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_this_bin.phi_binnum.to_numpy(), [24]])
# phi_bin_rebinned    = 15*phi_binnum_rebinned
# phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
weights             = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging
weights_stat_err    = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging * np.sqrt(df_summary_table_rebinned_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
axs[2].errorbar(df_summary_table_rebinned_this_bin.phi_avg_this_point, weights, yerr = weights_stat_err, ls = '', color = 'tab:cyan', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}} \times \epsilon_{\mathrm{additional}}$')

weights_min         = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging *0.7#- inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0_syst_err)
weights_max         = inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0) * df_summary_table_rebinned_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging *1.3#+ inverseHist(df_summary_table_rebinned_this_bin.normalization_pi0_syst_err)
phi_bin_rebinned_fill_between = []
weights_min_fill_between      = []
weights_max_fill_between      = []
for i in range(len(phi_bin_rebinned)):
    if (i > 0):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i-1])
        weights_max_fill_between.append(weights_max.to_numpy()[i-1])
    if (i<len(phi_bin_rebinned)-1):
        phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
        weights_min_fill_between.append(weights_min.to_numpy()[i])
        weights_max_fill_between.append(weights_max.to_numpy()[i])

axs[2].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


for ax in axs:
    ax.set_xlabel("$\phi$~($^\circ$)")
    ax.set_xlim([0, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim([0, 1.4])
    ax.tick_params(left = True, right = True, axis = 'y', which = 'both', direction = 'inout', length = 10)
axs[0].set_xticks([0, 90, 180, 270, 360], ['$0$', '$90$', '$180$', '$270$', '$360~0$'])
axs[1].set_xticks([0, 90, 180, 270, 360], ['', '$90$', '$180$', '$270$', '$360~0$'])
axs[2].set_xticks([0, 90, 180, 270, 360], ['', '$90$', '$180$', '$270$', '$360$'])

axs[1].set_yticks(np.linspace(0, 1.4, 8), ['']*8)
axs[2].set_yticks(np.linspace(0, 1.4, 8), ['']*8)
# axs[2].axhline(0.25, color = 'k', ls = '--')
ax.tick_params(left = True, right = True, axis = 'y', which = 'both', direction = 'inout', length = 10)
plt.subplots_adjust(wspace = 0)
    
xBmin, xBmax, Q2min, Q2max, t1min, t1max = df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, ["xBmin", "xBmax", "Q2min", "Q2max", "t1min", "t1max"]].to_numpy().T[:, 0]
print(xBmin, xBmax, Q2min, Q2max, t1min, t1max)

xBheader = "$x_B \\in [{:.3f}, {:.3f}]$\n".format(xBmin, xBmax)
Q2header = "$Q^2 \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2/c^2$\n".format(Q2min, Q2max)
t1header = "$|t|~ \\in [{:.3f}, {:.3f}]~\mathrm{{GeV}}^2$".format(t1min, t1max)
plt.figlegend(loc = 'upper left', title = "$88.~\mathrm{Combined~Data~Set}$\n"+xBheader +Q2header + t1header,  bbox_to_anchor = (0.9, 0.8))
plt.savefig("addendum_v3/efficiency.from_script.pdf", bbox_inches = 'tight')

for integrated_binnum in df_summary_table_rebinned.integrated_binnum:
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb_pi0_eff_corrected_bkg_merging_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_inb_pi0_eff_corrected_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb_pi0_eff_corrected_bkg_merging_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_outb_pi0_eff_corrected_bkg_merging"])
    df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_nominal_integrated"] = np.sum(df_summary_table_rebinned.loc[df_summary_table_rebinned.integrated_binnum == integrated_binnum, "active_bin_nominal"])

xB_panes = 8
Q2_panes = 7

labeled = 0

# models = ["contamination_inb_pi0_eff_corrected_bkg_merging", "contamination_inb_bh_eff_corrected_bkg_merging", "contamination_inb_km15_eff_corrected_bkg_merging", "contamination_inb_vgg_eff_corrected_bkg_merging", "contamination_inb_global1_eff_corrected_bkg_merging", ]
# plot_labels = [r"$\pi^0~\mathrm{Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
# colors = ['k', 'r', 'cyan', 'tab:orange', 'tab:blue', ]
models = ["contamination_inb_pi0_eff_corrected_bkg_merging", "contamination_inb_bkg_merging"]
plot_labels = [r"$\mathrm{With}~\pi^0~\mathrm{normalization}$", r"$\mathrm{Without~normalization}$"]
colors = ['k', 'tab:pink', 'cyan', 'tab:orange', 'tab:blue', ]

for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()<=4:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            for i, model in enumerate(models):
                cont = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), "{}".format(model)]
                cont_stat_err = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), "{}_stat_err".format(model)]
                phi_avg = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), "phi_avg_this_point"]
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_avg, cont, cont_stat_err, ls = '--', marker = 'o', label = plot_labels[i], color = colors[i])
            if not labeled:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
                
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>4), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim(bottom = 0)
    # handles = [handles[-1], handles[0], handles[1]]
    # labels = [labels[-1], labels[0], labels[1]]
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.25, hspace = 0.25)
    
    plt.savefig("addendum_v3/contamination_all_bin_inb_{}.new.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()    

xB_panes = 8
Q2_panes = 7

labeled = 0

# models = ["contamination_outb_pi0_eff_corrected_bkg_merging", "contamination_outb_bh_eff_corrected_bkg_merging", "contamination_outb_km15_eff_corrected_bkg_merging", "contamination_outb_vgg_eff_corrected_bkg_merging", "contamination_outb_global1_eff_corrected_bkg_merging", ]
# plot_labels = [r"$\pi^0~\mathrm{Norm.}$", r"$\mathrm{BH~Norm.}$", r"$\mathrm{KM15~Norm.}$", r"$\mathrm{VGG~Norm.}$", r"$\mathrm{Global~Norm.}$"]
# colors = ['k', 'r', 'cyan', 'tab:orange', 'tab:blue', ]
models = ["contamination_outb_pi0_eff_corrected_bkg_merging", "contamination_outb_bkg_merging"]
plot_labels = [r"$\mathrm{With}~\pi^0~\mathrm{normalization}$", r"$\mathrm{Without~normalization}$"]
colors = ['k', 'tab:pink', 'cyan', 'tab:orange', 'tab:blue', ]

for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()<=4:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            for i, model in enumerate(models):
                cont = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), "{}".format(model)]
                cont_stat_err = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), "{}_stat_err".format(model)]
                phi_avg = df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), "phi_avg_this_point"]
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_avg, cont, cont_stat_err, ls = '--', marker = 'o', label = plot_labels[i], color = colors[i])
            if not labeled:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
                
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>4), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim(bottom = 0)
    # handles = [handles[-1], handles[0], handles[1]]
    # labels = [labels[-1], labels[0], labels[1]]
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.25, hspace = 0.25)
    
    plt.savefig("addendum_v3/contamination_all_bin_outb_{}.new.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]
            df_summary_table_rebinned_exp_this_bin = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_inb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix) & (df_summary_table_rebinned_exp.integrated_binnum == integrated_binnum), :]


            phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_exp_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_summary_table_rebinned_exp_this_bin.n_entry_FD
            weights_stat_err    = np.sqrt(df_summary_table_rebinned_exp_this_bin.n_entry_FD)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~FD,~FD)}$", histtype = 'step', color = 'k')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            weights             = df_summary_table_rebinned_exp_this_bin.n_entry_CD
            weights_stat_err    = np.sqrt(df_summary_table_rebinned_exp_this_bin.n_entry_CD)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FD)}$", histtype = 'step', color = 'r')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
            weights             = df_summary_table_rebinned_exp_this_bin.n_entry_CDFT
            weights_stat_err    = np.sqrt(df_summary_table_rebinned_exp_this_bin.n_entry_CDFT)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FT)}$", histtype = 'step', color = 'tab:blue')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:blue')

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            
            ymax1 = np.max(df_summary_table_rebinned_exp_this_bin.loc[:, ["n_entry_FD", "n_entry_CD", "n_entry_CDFT"]])*1.2
            ymax2 = np.max(df_summary_table_rebinned_exp_this_bin.loc[(df_summary_table_rebinned_exp_this_bin.phi_binnum>6) & (df_summary_table_rebinned_exp_this_bin.phi_binnum<18), ["n_entry_FD", "n_entry_CD", "n_entry_CDFT"]])*1.4
            ymax  = np.max([ymax1, ymax2])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, ymax])
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 20, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 10)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
            # ymin = np.min([df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum)].km15_cross_section_this_point_norad.min(), np.min(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] - df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_syst_err"]), np.min(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] - df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_stat_err"])])
            # ymin = np.floor(np.log10(ymin)) 
            # ymax = np.max([df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum)].km15_cross_section_this_point_norad.max(), np.max(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] + df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_syst_err"]), np.max(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] + df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_stat_err"])])
            # ymax = np.ceil(np.log10(ymax))
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([10**(ymin), 10**ymax])
    # handles = [handles[-1], handles[0], handles[1]]
    # labels = [labels[-1], labels[0], labels[1]]
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/raw_yield_inbending_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]
            df_summary_table_rebinned_exp_this_bin = df_summary_table_rebinned_exp.loc[(df_summary_table_rebinned_exp.directory == "exp_fall2018_outb/dvcs") & (df_summary_table_rebinned_exp.variation == nominal_suffix) & (df_summary_table_rebinned_exp.integrated_binnum == integrated_binnum), :]


            phi_binnum_rebinned = np.concatenate([df_summary_table_rebinned_exp_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_summary_table_rebinned_exp_this_bin.n_entry_FD
            weights_stat_err    = np.sqrt(df_summary_table_rebinned_exp_this_bin.n_entry_FD)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~FD,~FD)}$", histtype = 'step', color = 'k')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            weights             = df_summary_table_rebinned_exp_this_bin.n_entry_CD
            weights_stat_err    = np.sqrt(df_summary_table_rebinned_exp_this_bin.n_entry_CD)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FD)}$", histtype = 'step', color = 'r')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
            weights             = df_summary_table_rebinned_exp_this_bin.n_entry_CDFT
            weights_stat_err    = np.sqrt(df_summary_table_rebinned_exp_this_bin.n_entry_CDFT)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Raw~Yield~(Exp.,~CD,~FT)}$", histtype = 'step', color = 'tab:blue')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:blue')

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            
            ymax1 = np.max(df_summary_table_rebinned_exp_this_bin.loc[:, ["n_entry_FD", "n_entry_CD", "n_entry_CDFT"]])*1.2
            ymax2 = np.max(df_summary_table_rebinned_exp_this_bin.loc[(df_summary_table_rebinned_exp_this_bin.phi_binnum>6) & (df_summary_table_rebinned_exp_this_bin.phi_binnum<18), ["n_entry_FD", "n_entry_CD", "n_entry_CDFT"]])*1.4
            ymax  = np.max([ymax1, ymax2])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, ymax])
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
            # ymin = np.min([df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum)].km15_cross_section_this_point_norad.min(), np.min(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] - df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_syst_err"]), np.min(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] - df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_stat_err"])])
            # ymin = np.floor(np.log10(ymin)) 
            # ymax = np.max([df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum)].km15_cross_section_this_point_norad.max(), np.max(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] + df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_syst_err"]), np.max(df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp"] + df_summary_table_rebinned.loc[(df_summary_table_rebinned.integrated_binnum == integrated_binnum) & (df_summary_table_rebinned.active_bin == 1), "xsec_exp_stat_err"])])
            # ymax = np.ceil(np.log10(ymax))
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([10**(ymin), 10**ymax])
    # handles = [handles[-1], handles[0], handles[1]]
    # labels = [labels[-1], labels[0], labels[1]]
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/raw_yield_outbending_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.epg_inb_exp
            weights_stat_err    = df_this_bin.epg_inb_exp_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Yield~(Exp.)}$", histtype = 'step', color = 'k')
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            ymax1 = np.max(weights + weights_stat_err)
            ymax2 = np.max( 1.2 * (weights + weights_stat_err)[(phi_bin_centers>90) & (phi_bin_centers<270)])
            ymax  = np.max([ymax1, ymax2])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, ymax])

            weights             = df_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Background~Yield~(Exp.)}$", histtype = 'step', color = 'r')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
            weights_min         = df_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging - df_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err
            weights_max         = df_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging + df_this_bin.bkg_inb_exp_pi0_eff_corrected_bkg_merging_syst_err
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'r', alpha = 0.3)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/background_yield_inbending_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.epg_outb_exp
            weights_stat_err    = df_this_bin.epg_outb_exp_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Yield~(Exp.)}$", histtype = 'step', color = 'k')
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            ymax1 = np.max(weights + weights_stat_err)
            ymax2 = np.max( 1.2 * (weights + weights_stat_err)[(phi_bin_centers>90) & (phi_bin_centers<270)])
            ymax  = np.max([ymax1, ymax2])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, ymax])

            weights             = df_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Background~Yield~(Exp.)}$", histtype = 'step', color = 'r')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'r')
            weights_min         = df_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging - df_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err
            weights_max         = df_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging + df_this_bin.bkg_outb_exp_pi0_eff_corrected_bkg_merging_syst_err
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'r', alpha = 0.3)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/background_yield_outbending_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Exp.)}$", histtype = 'step', color = 'k')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            weights_min         = df_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging * (1 - np.sqrt(df_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + 0.3**2 + 0.0476**2))
            weights_max         = df_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging * (1 + np.sqrt(df_this_bin.dvcs_inb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + 0.3**2 + 0.0476**2))
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

            weights             = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~KM15)}$", histtype = 'step', color = 'cyan')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')
            weights_min         = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging #* 0.7#- df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_syst_err
            weights_max         = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging #* 1.3#+ df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_syst_err
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


            weights             = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg
            weights_stat_err    = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~VGG)}$", histtype = 'step', color = 'tab:orange')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')
            weights_min         = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg #- df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err_down
            weights_max         = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg #+ df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err_up
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:orange', alpha = 0.5)

            weights             = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH
            weights_stat_err    = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~BH)}$", histtype = 'step', color = 'tab:red')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')
            weights_min         = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH #- df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err_down
            weights_max         = df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH #+ df_this_bin.dvcs_inb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err_up
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim(bottom = 0)
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/signal_yield_inbending_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Exp.)}$", histtype = 'step', color = 'k')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            weights_min         = df_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging * (1 - np.sqrt(df_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + 0.3**2 + 0.0476**2))
            weights_max         = df_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging * (1 + np.sqrt(df_this_bin.dvcs_outb_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + 0.3**2 + 0.0476**2))
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

            weights             = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~KM15)}$", histtype = 'step', color = 'cyan')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')
            weights_min         = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging #* 0.7#- df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_syst_err
            weights_max         = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging #* 1.3#+ df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_syst_err
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


            weights             = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg
            weights_stat_err    = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~VGG)}$", histtype = 'step', color = 'tab:orange')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')
            weights_min         = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg #- df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err_down
            weights_max         = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg #+ df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err_up
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:orange', alpha = 0.5)

            weights             = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH
            weights_stat_err    = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Signal~Yield~(Sim.,~BH)}$", histtype = 'step', color = 'tab:red')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')
            weights_min         = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH #- df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err_down
            weights_max         = df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH #+ df_this_bin.dvcs_outb_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err_up
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.13, 1.04), xytext = (0.13, 1.04), xycoords = 'axes fraction', horizontalalignment = 'right')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.15, 0.9), xytext = (0.15, 0.9), xycoords = 'axes fraction')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim(bottom = 0)
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = True, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/signal_yield_outbending_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_nominal == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_nominal.sum()<5:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Exp.}$", histtype = 'step', color = 'k')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k')
            weights_min         = df_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging * ( 1 - np.sqrt(0.3**2 + 0.0476**2 + df_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2))
            weights_max         = df_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging * ( 1 + np.sqrt(0.3**2 + 0.0476**2 + df_this_bin.acceptance_corrected_yield_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio**2))
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

            weights             = df_this_bin.gen_sim
            weights_stat_err    = df_this_bin.gen_sim_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Sim.,~KM15}$", histtype = 'step', color = 'cyan', lw = 3, alpha = 0.8)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')

            weights             = df_this_bin.gen_sim_vgg
            weights_stat_err    = df_this_bin.gen_sim_vgg_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Sim.,~VGG}$", histtype = 'step', color = 'tab:orange', lw = 3, alpha = 0.8)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')

            weights             = df_this_bin.gen_sim_pureBH
            weights_stat_err    = df_this_bin.gen_sim_pureBH_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{Sim.,~BH}$", histtype = 'step', color = 'tab:red', lw = 3, alpha = 0.8)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            
            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_nominal_integrated>4), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Combined}$" +"\n" + r"$\mathrm{Acceptance~corrected~yield}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/acceptance_corrected_yield_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_nominal == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_nominal.sum()<5:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{KM15}$", histtype = 'step', color = 'cyan')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'cyan')
            weights_min         = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging * ( 1 - np.sqrt(df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + 0.3**2 + 0.0476**2))
            weights_max         = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging * ( 1 + np.sqrt(df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_syst_err_ratio**2 + 0.3**2 + 0.0476**2))
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'cyan', alpha = 0.5)


            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg
            weights_stat_err    = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{VGG}$", histtype = 'step', color = 'tab:orange')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:orange')
            weights_min         = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg# - df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err_normalization
            weights_max         = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg# + df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_vgg_syst_err_normalization
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:orange', alpha = 0.5)

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH
            weights_stat_err    = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].hist(phi_bin_rebinned[:-1], phi_bin_rebinned, weights = weights, label = "$\mathrm{BH}$", histtype = 'step', color = 'tab:red')
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'tab:red')
            weights_min         = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH# - df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err_normalization
            weights_max         = df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH# + df_this_bin.acceptance_sim_pi0_eff_corrected_bkg_merging_pureBH_syst_err_normalization
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'tab:red', alpha = 0.5)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_nominal_integrated>4), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Combined}$" +"\n" + r"$\mathrm{Acceptance}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/acceptance_all_bin_{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

            weights             = df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r"$\epsilon_{\mathrm{bkg.~merging}}$")

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.2])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/bkg_merging_inbending_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

            weights             = df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r"$\epsilon_{\mathrm{bkg.~merging}}$")

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.2])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/bkg_merging_outbending_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_nominal == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_nominal.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)

            weights             = df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging
            weights_stat_err    = df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r"$\epsilon_{\mathrm{bkg.~merging}}$")

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.2])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_nominal_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Combined~Data~Set}$" + "\n"+ r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/bkg_merging_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = inverseHist(df_this_bin.normalization_pi0)#df_this_bin.weights_inb
            weights_stat_err    = inverseHist(df_this_bin.normalization_pi0) *  df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio#df_this_bin.weights_inb_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r"$\epsilon_{\mathrm{additional}}$")


            weights_min         = inverseHist(df_this_bin.normalization_pi0) * 0.7
            weights_max         = inverseHist(df_this_bin.normalization_pi0) * 1.3
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.6])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/weight_inbending_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = inverseHist(df_this_bin.normalization_pi0)#df_this_bin.weights_outb
            weights_stat_err    = inverseHist(df_this_bin.normalization_pi0) *  df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio#df_this_bin.weights_outb_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o',  label = r"$\epsilon_{\mathrm{additional}}$")


            weights_min         = inverseHist(df_this_bin.normalization_pi0) * 0.7
            weights_max         = inverseHist(df_this_bin.normalization_pi0) * 1.3
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.6])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/weight_outbending_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_nominal == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_nominal.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = inverseHist(df_this_bin.normalization_pi0)#df_this_bin.weights
            weights_stat_err    = inverseHist(df_this_bin.normalization_pi0) *  df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio#df_this_bin.weights_stat_err
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o',  label = r"$\epsilon_{\mathrm{additional}}$")


            weights_min         = inverseHist(df_this_bin.normalization_pi0) * 0.7
            weights_max         = inverseHist(df_this_bin.normalization_pi0) * 1.3
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.6])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_nominal_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Combined~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/weight_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_inb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = inverseHist(df_this_bin.normalization_inb_pi0) * df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging
            weights_stat_err    = inverseHist(df_this_bin.normalization_inb_pi0) * df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging * np.sqrt(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}} \times \epsilon_{\mathrm{additional}}$')

            weights_min         = inverseHist(df_this_bin.normalization_inb_pi0) * df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging * 0.7
            weights_max         = inverseHist(df_this_bin.normalization_inb_pi0) * df_this_bin.eff_bkg_merging_inb_pi0_eff_corrected_bkg_merging * 1.3
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].axhline(0.25, color = 'k', ls = '--')

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.2])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_inb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Inbending~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/efficiency_inbending_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_outb_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = inverseHist(df_this_bin.normalization_outb_pi0) * df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging
            weights_stat_err    = inverseHist(df_this_bin.normalization_outb_pi0) * df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging * np.sqrt(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}} \times \epsilon_{\mathrm{additional}}$')

            weights_min         = inverseHist(df_this_bin.normalization_outb_pi0) * df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging * 0.7
            weights_max         = inverseHist(df_this_bin.normalization_outb_pi0) * df_this_bin.eff_bkg_merging_outb_pi0_eff_corrected_bkg_merging * 1.3
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].axhline(0.25, color = 'k', ls = '--')

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.2])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_outb_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Outbending~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/efficiency_outbending_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()

xB_panes = 8
Q2_panes = 7


for t_binnum in range(6):
    fig, axs = plt.subplots(Q2_panes, xB_panes, figsize = (45, 32.5))
    labeled = 0
    t_avgs = []
    for xB_binnum in range(xB_panes):
        for Q2_binnum in range(Q2_panes):
    
            df_this_bin = df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.Q2bin == Q2_binnum) & (df_summary_table_rebinned.tbin == t_binnum)  & (df_summary_table_rebinned.active_bin_pi0_eff_corrected_bkg_merging == 1), :]
            # if not len(df_this_bin)>0:
            if df_this_bin.active_bin_pi0_eff_corrected_bkg_merging.sum()==0:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].axis('off')
                continue

            integrated_binnum = df_this_bin.integrated_binnum.unique()[0]

            phi_binnum_rebinned = np.concatenate([df_this_bin.phi_binnum.to_numpy(), [24]])
            phi_bin_rebinned    = 15*phi_binnum_rebinned
            phi_bin_centers     = np.mean([phi_bin_rebinned[:-1], phi_bin_rebinned[1:]], axis = 0)
            weights             = inverseHist(df_this_bin.normalization_pi0) * df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging
            weights_stat_err    = inverseHist(df_this_bin.normalization_pi0) * df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging * np.sqrt(df_this_bin.xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio**2 + df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging_stat_err_ratio**2)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].errorbar(phi_bin_centers, weights, yerr = weights_stat_err, ls = '', color = 'k', marker = 'o', label = r'$\epsilon_{\mathrm{bkg.~merging}} \times \epsilon_{\mathrm{additional}}$')

            weights_min         = inverseHist(df_this_bin.normalization_pi0) * df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging * 0.7
            weights_max         = inverseHist(df_this_bin.normalization_pi0) * df_this_bin.eff_bkg_merging_pi0_eff_corrected_bkg_merging * 1.3
            phi_bin_rebinned_fill_between = []
            weights_min_fill_between      = []
            weights_max_fill_between      = []
            for i in range(len(phi_bin_rebinned)):
                if (i > 0):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i-1])
                    weights_max_fill_between.append(weights_max.to_numpy()[i-1])
                if (i<len(phi_bin_rebinned)-1):
                    phi_bin_rebinned_fill_between.append(phi_bin_rebinned[i])
                    weights_min_fill_between.append(weights_min.to_numpy()[i])
                    weights_max_fill_between.append(weights_max.to_numpy()[i])

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].fill_between(phi_bin_rebinned_fill_between, weights_min_fill_between, weights_max_fill_between, color = 'k', alpha = 0.5)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].axhline(0.25, color = 'k', ls = '--')

            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate("${}$.".format(integrated_binnum) , xy = (0.33, 1.03), xytext = (0.33, 1.03), xycoords = 'axes fraction', horizontalalignment = 'right', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].annotate(r"$ x_B \in[{:.3f}, {:.3f}]$".format(df_this_bin.xBmin.unique()[0], df_this_bin.xBmax.unique()[0]) +"\n"+ r"$Q^2 \in [{:.3f}, {:.3f}]$".format(df_this_bin.Q2min.unique()[0], df_this_bin.Q2max.unique()[0]), xy = (0.35, 0.92), xytext = (0.35, 0.92), xycoords = 'axes fraction', fontsize = 20)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlim([0, 360])
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_ylim([0, 1.2])
            # locator0 = MaxNLocator()
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].yaxis.set_major_locator(locator0)
            # axs[Q2_panes - Q2_binnum - 1, xB_binnum].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if Q2_binnum == df_summary_table_rebinned.loc[(df_summary_table_rebinned.xBbin == xB_binnum) & (df_summary_table_rebinned.tbin == t_binnum) & (df_summary_table_rebinned.active_bin_pi0_eff_corrected_bkg_merging_integrated>0), "Q2bin"].min():
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1))
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 12+1), minor = True)
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xlabel("$\phi$ ($^{\circ}$)" )
            else:
                axs[Q2_panes - Q2_binnum - 1, xB_binnum].set_xticks(np.linspace(0,360, 4+1), ['']*5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'major', direction = 'in', width = 4, length = 10, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(left = True, right = False, axis = 'y', which = 'minor', direction = 'in', width = 2, length = 5)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'major', direction = 'in', width = 1, length = 5, labelsize = 25)
            axs[Q2_panes - Q2_binnum - 1, xB_binnum].tick_params(axis = 'x', which = 'minor', direction = 'in', width = 1, length = 3)
            if labeled:
                pass
            else:
                handles, labels = axs[Q2_panes - Q2_binnum - 1, xB_binnum].get_legend_handles_labels()
                labeled = 1
    plt.figlegend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.2, 0.8), fontsize = 40, title = r"$\mathrm{Combined~Data~Set}$" +"\n" + r"$|t| \in[{:.3f}, {:.3f}]$".format(df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1min.unique()[0], df_summary_table_rebinned.loc[df_summary_table_rebinned.tbin== t_binnum].t1max.unique()[0]) + r" $\mathrm{GeV}^2$", title_fontsize = 40)
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    
    plt.savefig("addendum_v3/efficiency_all_bin{}.from_script.pdf".format(t_binnum), bbox_inches = 'tight')
    plt.close()
'''


'''
# #To Melany
# df_summary_table = pd.read_pickle("addendum_v3/df_summary_table_rebinned_approved.pkl")
# df_summary_table_report = df_summary_table.loc[df_summary_table.active_bin_nominal == 1, :]

# xB_report            = df_summary_table_report.xB_avg_this_point
# Q2_report            = df_summary_table_report.Q2_avg_this_point
# t_report             = df_summary_table_report.t_avg_this_point
# phi_report           = df_summary_table_report.phi_avg_this_point
# xsec_report          = df_summary_table_report.xsec_exp_pi0_eff_corrected_bkg_merging
# xsec_stat_err_ratio_report = xsec_exp_pi0_eff_corrected_bkg_merging_stat_err_ratio
# xsec_syst_err_ratio_report = xsec_exp_pi0_eff_corrected_bkg_merging_syst_err_ratio

# header_string = "xB, Q2, -t, phi, XUU, stat, sysm, sysp\n , GeV2/c2, GeV2, deg, %, %, %"
# data_for_Melany = np.array([xB_report, Q2_report, t_report, phi_report, xsec_report, xsec_stat_err_report, xsec_syst_err_report, xsec_syst_err_report]).T
# np.savetxt("addendum_v3/data_for_Melany.csv", data_for_Melany, delimiter=",", header=header_string, fmt = '%.3f, %.3f, %.3f, %.3f, %.3e, %.3f, %.3f, %.3f')
'''