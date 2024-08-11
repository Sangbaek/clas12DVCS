import numpy as np
import pandas as pd
from utils.const import *
from utils.physics import *
import os

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


def main(mode):

	suffix = schema_suffices[mode]

	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl_{}".format(suffix)):
		return

	print(suffix)

	df_summary = pd.DataFrame()
	binnum, bin_volume = np.loadtxt('volume_list.csv', skiprows = 1, delimiter = ',').T
	bin_volume = {int(binnum[i]): bin_volume[i] for i in range(len(binnum))}

	df_merged = {}

	for integrated_binnum in range(1, 147+1):
		df_merged[integrated_binnum] = []
	for integrated_binnum_gen in range(1, 147+1):
		print(integrated_binnum_gen)
		try:
			df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/pkl_{}/{}".format(suffix, integrated_binnum_gen))
		except:
			continue
		df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
		df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
		df.loc[:, "n_gen"] = 0
		df.loc[:, "bin_volume_gen"] = 0
		for phi_binnum_gen in range(24):
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_inb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
			 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
		for integrated_binnum in range(1,147+1):
			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
			df_merged[integrated_binnum].append(df_this_bin)
	for integrated_binnum_gen in range(148, 147+159+1):
		print(integrated_binnum_gen)
		try:
			df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15_fringe/excl_level_1/pkl_{}/{}".format(suffix, integrated_binnum_gen - 147))
		except:
			continue
		df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
		df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
		df.loc[:, "n_gen"] = 0
		df.loc[:, "bin_volume_gen"] = 0
		for phi_binnum_gen in range(24):
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_inb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
			 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
		for integrated_binnum in range(1,147+1):
			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
			df_merged[integrated_binnum].append(df_this_bin)

	for integrated_binnum in range(1, 147+1):
		df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_inb/dvcs_km15/excl_level_1/restructured_{}/{}.pkl".format(suffix, integrated_binnum))


	df_merged = {}

	for integrated_binnum in range(1, 147+1):
		df_merged[integrated_binnum] = []
	for integrated_binnum_gen in range(1, 147+1):
		print(integrated_binnum_gen)
		try:
			df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/pkl_{}/{}".format(suffix, integrated_binnum_gen))
		except:
			continue
		df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
		df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
		df.loc[:, "n_gen"] = 0
		df.loc[:, "bin_volume_gen"] = 0
		for phi_binnum_gen in range(24):
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_outb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
			 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_outb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
		for integrated_binnum in range(1,147+1):
			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
			df_merged[integrated_binnum].append(df_this_bin)
	for integrated_binnum_gen in range(148, 147+159+1):
		print(integrated_binnum_gen)
		try:
			df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15_fringe/excl_level_1/pkl_{}/{}".format(suffix, integrated_binnum_gen - 147))
		except:
			continue
		df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
		df = df.astype({"integrated_binnum_gen": int, "phi_binnum_gen": int})
		df.loc[:, "n_gen"] = 0
		df.loc[:, "bin_volume_gen"] = 0
		for phi_binnum_gen in range(24):
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == "dvcs_km15/fall2018_outb3") & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
			 * df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_outb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
		for integrated_binnum in range(1,147+1):
			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
			df_merged[integrated_binnum].append(df_this_bin)

	for integrated_binnum in range(1, 147+1):
		df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/sim_rad_rec_fall2018_outb/dvcs_km15/excl_level_1/restructured_{}/{}.pkl".format(suffix, integrated_binnum))

	return



if __name__ == "__main__":

	for mode in range(6, 14):
		main(mode)


