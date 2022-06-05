#!/usr/bin/env python3
"""
Convert to np histogram for the binning
"""
import pandas as pd
import numpy as np
from utils.const import *
import argparse

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-j","--jobnum", help="jobnum", default="3987")
	parser.add_argument("-u","--uri", help="uri", default = "/volatile/clas12/sangbaek/nov2021/convPkl_full/")
	parser.add_argument("-i","--binscheme", help="binning scheme number", default="0")

	args = parser.parse_args()

	i = int(args.binscheme)

	xBbins  = collection_xBbins[i]
	Q2bins  = collection_Q2bins[i]
	tbins   = collection_tbins [i]
	phibins = collection_phibins[i]

	jobnum = int(args.jobnum)


	if jobnum in [*runs_inb_vgg50nA, *runs_inb_vgg55nA, *runs_inb_vgg45nA, *runs_inb_vgg0nA]:
		pol = "inb"
		mode = "dvcs"

	if jobnum in [*runs_outb_vgg50nA, *runs_outb_vgg40nA, *runs_outb_vgg0nA, *runs_outb_vgg40nAT]:
		pol = "outb"
		mode = "dvcs"

	if jobnum in [*runs_inb_bh50nA, *runs_inb_bh45nA]:
		pol = "inb"
		mode = "bh"

	if jobnum in runs_outb_bh50nA:
		pol = "outb"
		mode = "bh"

	if "Gen" in args.uri:
		recgen = "Gen"
	else:
		recgen = "Rec"

	infile = "{}/{}/{}/{}.pkl".format(args.uri, pol, mode, jobnum)
	print("reading {}".format(infile))
	df = pd.read_pickle(infile)

	# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}.npz".format(jobnum, recgen)
	# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	# np.savez(out, hist = hist)

	# for config in df.config.unique():
	# 	out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, config)
	# 	hist, _ = np.histogramdd(df.loc[df.config == config , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	# 	np.savez(out, hist = hist)

	if recgen == "Gen":
		# #phi
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "phi")
		# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = df.phi1)
		# np.savez(out, hist = hist)
		# #BornWeight
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "born")
		# hist, _ = np.histogramdd(df.loc[df.BornWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[df.BornWeight>0, "BornWeight"])
		# np.savez(out, hist = hist)
		# #GenWeight
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "rad")
		# hist, _ = np.histogramdd(df.loc[df.GenWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[df.GenWeight>0, "GenWeight"])
		# np.savez(out, hist = hist)

		# #Integrated
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "Int")
		# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins])
		# np.savez(out, hist = hist)
		# #xB
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "xB")
		# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.xB)
		# np.savez(out, hist = hist)
		# #Q2
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "Q2")
		# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.Q2)
		# np.savez(out, hist = hist)
		# #t
		# out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "t1")
		# hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.t1)
		# np.savez(out, hist = hist)

		#volume count
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "binVolume")
		finexBbins, fineQ2bins, finetbins, finephibins = [], [], [], []
		for xBind in range(len(xBbins)-1):
			finexBbins.append(np.linspace(xBbins[xBind], xBbins[xBind+1], 6))
		for Q2ind in range(len(Q2bins)-1):
			fineQ2bins.append(np.linspace(Q2bins[Q2ind], Q2bins[Q2ind+1], 6))
		for tind in range(len(tbins)-1):
			finetbins.append(np.linspace(tbins[tind], tbins[tind+1], 6))
		for phiind in range(len(phibins)-1):
			finephibins.append(np.linspace(phibins[phiind], phibins[phiind+1], 6))

		finexBbins = np.unique(np.concatenate(finexBbins))
		fineQ2bins = np.unique(np.concatenate(fineQ2bins))
		finetbins = np.unique(np.concatenate(finetbins))
		finephibins = np.unique(np.concatenate(finephibins))
		hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [finexBbins, fineQ2bins, finetbins, finephibins])
		hist = np.divide(hist,hist, where = hist>0, out = np.zeros(hist))
		# print("saving xB{}Q2{}t{}phi{} bins occupied {}".format(xBind, Q2ind, tind, phiind, np.sum(hist)/6**4))
		np.savez(out, hist=hist)