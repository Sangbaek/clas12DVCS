#!/usr/bin/env python3
"""
Convert to np histogram for the binning
"""
import pandas as pd
import numpy as np
from utils.const import *
from utils.physics import *
import argparse
import os

def divideHist(df1, df2):
	return np.divide(df1, df2, where = df2!=0, out = np.zeros(df2.shape, dtype = float))

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-j","--jobnum", help="jobnum", default="3987")
	parser.add_argument("-u","--uri", help="uri", default = "/volatile/clas12/sangbaek/nov2021/convPkl_full/")
	parser.add_argument("-o","--output", help="output", default = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/")
	parser.add_argument("-kstart","--kstart", help="binning scheme number starting", default=None)

	args = parser.parse_args()
	jobnum = int(args.jobnum)
	output = args.output

	if args.kstart:
		kstart = 0
	else:
		kstart = len(collection_xBbins) - 1

	if jobnum in [*runs_inb_vgg45nA, *runs_inb_vgg50nA, *runs_inb_vgg55nA]:
		pol = "inb"
		mode = "dvcs"

	if jobnum in [*runs_outb_vgg50nA, *runs_outb_vgg40nA, *runs_outb_vgg40nAT]:
		pol = "outb"
		mode = "dvcs"

	if jobnum in [*runs_inb_bh45nA, *runs_inb_bh50nA, *runs_inb_bh55nA]:
		pol = "inb"
		mode = "bh"

	if jobnum in [*runs_outb_bh50nA, *runs_outb_bh40nA, *runs_outb_bh40nAT]:
		pol = "outb"
		mode = "bh"

	if "Gen" in args.uri:
		recgen = "Gen"
	else:
		recgen = "Rec"

	infile = "{}/{}/{}/{}.pkl".format(args.uri, pol, mode, jobnum)
	print("reading {}".format(infile))
	df = pd.read_pickle(infile)

	nu  = df.Q2/(2*M*df.xB)
	yb= nu/10.604
	W2 = M*M+2.0*M*10.604*yb-df.Q2
	tmin1 = np.abs(tmin(df.xB, df.Q2, df.t1, df.phi1))
	cond = (nu>2)&(nu<10.604 - 2) & (W2>4) & (df.t1>tmin1) & (Kfac2(df.xB, df.Q2, df.t1, df.phi1)>0)
	df = df.loc[cond, :]

	os.makedirs("{}".format(output), exist_ok = True)
	for k in range(kstart, len(collection_xBbins)):

		os.makedirs("{}binscheme{}".format(output, k), exist_ok = True)
		xBbins  = collection_xBbins[k]
		Q2bins  = collection_Q2bins[k]
		tbins   = collection_tbins [k]
		phibins = collection_phibins[k]

		#all
		out = "{}binscheme{}/{}{}.npz".format(output, k, jobnum, recgen)
		hist_all, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		np.savez(out, hist = hist_all)
		# helicity plus
		out = "{}binscheme{}/{}{}plus.npz".format(output, k, jobnum, recgen)
		hist_plus, _ = np.histogramdd(df.loc[df.helicity == 1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		np.savez(out, hist = hist_plus)
		# helicity minus
		out = "{}binscheme{}/{}{}minus.npz".format(output, k, jobnum, recgen)
		hist_minus, _ = np.histogramdd(df.loc[df.helicity == -1 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		np.savez(out, hist = hist_minus)

		for config in df.config.unique():
			#all
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, config)
			hist, _ = np.histogramdd(df.loc[df.config == config , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			np.savez(out, hist = hist)
			# helicity plus
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, config)
			hist, _ = np.histogramdd(df.loc[(df.config == config) & (df.helicity == 1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			np.savez(out, hist = hist)
			# helicity minus
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, config)
			hist, _ = np.histogramdd(df.loc[(df.config == config) & (df.helicity == -1) , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			np.savez(out, hist = hist)

		if recgen == "Gen":
			#all
			#phi
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "phi")
			hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = df.phi1)
			np.savez(out, hist = divideHist(hist, hist_all))
			#BornWeight
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "born")
			hist0, _ = np.histogramdd(df.loc[df.BornWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			hist, _ = np.histogramdd(df.loc[df.BornWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[df.BornWeight>0, "BornWeight"])
			np.savez(out, hist = divideHist(hist0, hist)*0.001*(2*np.pi))
			#GenWeight
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "rad")
			hist0, _ = np.histogramdd(df.loc[df.GenWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			hist, _ = np.histogramdd(df.loc[df.GenWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[df.GenWeight>0, "GenWeight"])
			np.savez(out, hist = divideHist(hist0, hist)*0.001*(2*np.pi))

			#Integrated
			# out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "Int")
			hist0, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins])
			# np.savez(out, hist = hist)
			#xB
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "xB")
			hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.xB)
			np.savez(out, hist = divideHist(hist, hist0))
			#Q2
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "Q2")
			hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.Q2)
			np.savez(out, hist = divideHist(hist, hist0))
			#t
			out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "t1")
			hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.t1)
			np.savez(out, hist = divideHist(hist, hist0))
			#helicity plus
			#phi
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "phi")
			hist, _ = np.histogramdd(df.loc[df.helicity == 1, ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = df.loc[df.helicity == 1, "phi1"])
			np.savez(out, hist = divideHist(hist, hist_plus))
			#BornWeight
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "born")
			hist0, _ = np.histogramdd(df.loc[(df.BornWeight>0) & (df.helicity == 1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			hist, _ = np.histogramdd(df.loc[(df.BornWeight>0) & (df.helicity == 1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[(df.BornWeight>0) & (df.helicity == 1), "BornWeight"])
			np.savez(out, hist = divideHist(hist0, hist)*0.001*(2*np.pi))
			#GenWeight
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "rad")
			hist0, _ = np.histogramdd(df.loc[(df.GenWeight>0) & (df.helicity == 1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			hist, _ = np.histogramdd(df.loc[(df.GenWeight>0) & (df.helicity == 1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[(df.GenWeight>0) & (df.helicity == 1), "GenWeight"])
			np.savez(out, hist = divideHist(hist0, hist)*0.001*(2*np.pi))
			#Integrated
			# out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "Int")
			hist0, _ = np.histogramdd(df.loc[df.helicity == 1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins])
			# np.savez(out, hist = divideHist(hist, hist0))
			#xB
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "xB")
			hist, _ = np.histogramdd(df.loc[df.helicity == 1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.loc[df.helicity == 1, "xB"])
			np.savez(out, hist = divideHist(hist, hist0))
			#Q2
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "Q2")
			hist, _ = np.histogramdd(df.loc[df.helicity == 1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.loc[df.helicity == 1, "Q2"])
			np.savez(out, hist = divideHist(hist, hist0))
			#t
			out = "{}binscheme{}/{}{}{}plus.npz".format(output, k, jobnum, recgen, "t1")
			hist, _ = np.histogramdd(df.loc[df.helicity == 1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.loc[df.helicity == 1, "t1"])
			np.savez(out, hist = divideHist(hist, hist0))

			#helicity minus
			#phi
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "phi")
			hist, _ = np.histogramdd(df.loc[df.helicity == -1, ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = df.loc[df.helicity == -1, "phi1"])
			np.savez(out, hist = divideHist(hist, hist_minus))
			#BornWeight
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "born")
			hist0, _ = np.histogramdd(df.loc[(df.BornWeight>0) & (df.helicity == -1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			hist, _ = np.histogramdd(df.loc[(df.BornWeight>0) & (df.helicity == -1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[(df.BornWeight>0) & (df.helicity == -1), "BornWeight"])
			np.savez(out, hist = divideHist(hist0, hist)*0.001*(2*np.pi))
			#GenWeight
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "rad")
			hist0, _ = np.histogramdd(df.loc[(df.GenWeight>0) & (df.helicity == -1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
			hist, _ = np.histogramdd(df.loc[(df.GenWeight>0) & (df.helicity == -1), ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[(df.GenWeight>0) & (df.helicity == -1), "GenWeight"])
			np.savez(out, hist = divideHist(hist0, hist)*0.001*(2*np.pi))
			#Integrated
			# out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "Int")
			hist0, _ = np.histogramdd(df.loc[df.helicity == -1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins])
			# np.savez(out, hist = divideHist(hist, hist0))
			#xB
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "xB")
			hist, _ = np.histogramdd(df.loc[df.helicity == -1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.loc[df.helicity == -1, "xB"])
			np.savez(out, hist = divideHist(hist, hist0))
			#Q2
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "Q2")
			hist, _ = np.histogramdd(df.loc[df.helicity == -1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.loc[df.helicity == -1, "Q2"])
			np.savez(out, hist = divideHist(hist, hist0))
			#t
			out = "{}binscheme{}/{}{}{}minus.npz".format(output, k, jobnum, recgen, "t1")
			hist, _ = np.histogramdd(df.loc[df.helicity == -1, ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.loc[df.helicity == -1, "t1"])
			np.savez(out, hist = divideHist(hist, hist0))

			# #volume count
			# out = "{}binscheme{}/{}{}{}.npz".format(output, k, jobnum, recgen, "binVolume")
			# finexBbins, fineQ2bins, finetbins, finephibins = [], [], [], []
			# for xBind in range(len(xBbins)-1):
			# 	finexBbins.append(np.linspace(xBbins[xBind], xBbins[xBind+1], 6+1))
			# for Q2ind in range(len(Q2bins)-1):
			# 	fineQ2bins.append(np.linspace(Q2bins[Q2ind], Q2bins[Q2ind+1], 6+1))
			# for tind in range(len(tbins)-1):
			# 	finetbins.append(np.linspace(tbins[tind], tbins[tind+1], 6+1))
			# for phiind in range(len(phibins)-1):
			# 	finephibins.append(np.linspace(phibins[phiind], phibins[phiind+1], 6+1))

			# finexBbins = np.unique(np.concatenate(finexBbins))
			# fineQ2bins = np.unique(np.concatenate(fineQ2bins))
			# finetbins = np.unique(np.concatenate(finetbins))
			# finephibins = np.unique(np.concatenate(finephibins))
			# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [finexBbins, fineQ2bins, finetbins, finephibins])
			# hist = np.divide(hist,hist, where = hist>0, out = np.zeros_like(hist))
			# np.savez(out, hist=hist)