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

	out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}.npz".format(jobnum, recgen)
	hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	np.savez(out, hist = hist)

	for config in df.config.unique():
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, config)
		hist, _ = np.histogramdd(df.loc[df.config == config , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		np.savez(out, hist = hist)

	if recgen == "Gen":
		#phi
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "phi")
		hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = df.phi1)
		np.savez(out, hist = hist)
		#BornWeight
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "born")
		hist, _ = np.histogramdd(df.loc[df.BornWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[df.BornWeight>0, "BornWeight"])
		np.savez(out, hist = hist)
		#GenWeight
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "rad")
		hist, _ = np.histogramdd(df.loc[df.GenWeight>0 , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins], weights = 1/df.loc[df.GenWeight>0, "GenWeight"])
		np.savez(out, hist = hist)
		#xB
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "xB")
		hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.xB)
		hist = np.stack([hist]*(len(phibins)-1), axis = 3)
		np.savez(out, hist = hist)
		#Q2
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "Q2")
		hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.Q2)
		hist = np.stack([hist]*(len(phibins)-1), axis = 3)
		np.savez(out, hist = hist)
		#t
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, "t1")
		hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1"]].to_numpy(), bins = [xBbins, Q2bins, tbins], weights = df.t1)
		hist = np.stack([hist]*(len(phibins)-1), axis = 3)
		np.savez(out, hist = hist)