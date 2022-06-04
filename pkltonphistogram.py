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

	runs_inb_vgg50nA = [3987, 4124, 4139, 4181, 4182, 4397, 4528, 4529, 4535, 4539]
	runs_inb_vgg55nA = [4186, 4545]
	runs_inb_vgg45nA = [4188, 4547]
	runs_inb_vgg0nA = [4192, 4561]
	runs_inb_bh50nA = [4238, 4542]
	runs_inb_bh45nA = [4740, 4742, 4745, 4751, 4760]

	runs_outb_vgg50nA = [4240, 4250, 4251, 4252, 4255, 4398, 4532, 4534, 4540, 4541, 4717]
	runs_outb_vgg40nA = [4263, 4546]
	runs_outb_vgg0nA = [4262, 4554]
	runs_outb_vgg40nAT = [4266, 4562]
	runs_outb_bh50nA = [4249, 4544, 4780, 4808, 4812, 4818]


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

	df = pd.read_pickle("{}/{}/{}/{}.pkl".format(args.uri, pol, mode, jobnum))


	out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}.npz".format(jobnum, recgen)
	hist, _ = np.histogramdd(df.loc[: , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
	np.savez(out, hist = hist)

	for config in df.configs.unique():
		out = "/volatile/clas12/sangbaek/clas12DVCS/nphistograms/{}{}{}.npz".format(jobnum, recgen, config)
		hists, _ = np.histogramdd(df.loc[df.config == config , ["xB", "Q2", "t1", "phi1"]].to_numpy(), bins = [xBbins, Q2bins, tbins, phibins])
		np.savez(out, hist = hist)