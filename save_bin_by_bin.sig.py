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

# df_volume = pd.read_csv("volume_list.csv")
# bin_volume_bulk = df_volume.loc[df_volume.integrated_bin <=147, :].to_numpy()[:, 1]
# bin_volume_fringe = df_volume.loc[df_volume.integrated_bin > 147, :].to_numpy()[:, 1]
binnum, bin_volume = np.loadtxt('volume_list.csv', skiprows = 1, delimiter = ',').T
bin_volume = {int(binnum[i]): bin_volume[i] for i in range(len(binnum))}


df_gen = pd.read_pickle("summary_table.gen.pkl")


def main_memory_efficient(mode, sig_directory, gen_directory):

	suffix = schema_suffices[mode]

	bulk_directory = "/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}".format(sig_directory, suffix)
	fringe_directory = "/volatile/clas12/sangbaek/dvcs_related/{}_fringe/excl_level_1/pkl_{}".format(sig_directory, suffix)

	if not os.path.exists(bulk_directory):
		return

	print(suffix)

	for integrated_binnum in range(1, 147+1):
		print("Integrated bin number {}".format(integrated_binnum))
		df_merged = []
		for integrated_binnum_gen in range(1, 147+159+1):
			if integrated_binnum_gen <= 147:
				print ("Reading bulk {}".format(integrated_binnum_gen))
				try:
					df = pd.read_pickle("{}/{}".format(bulk_directory, integrated_binnum_gen))
				except:
					print("{}/{} does not exist.".format(bulk_directory, integrated_binnum_gen))
					continue
			else:
				print ("Reading fringe {}".format(integrated_binnum_gen-147))
				try:
					df = pd.read_pickle("{}/{}".format(fringe_directory, integrated_binnum_gen-147))
				except:
					print("{}/{} does not exist.".format(fringe_directory, integrated_binnum_gen-147))
					continue

			df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
			for phi_binnum_gen in range(24):
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == gen_directory) & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
				df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
					* df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen

			df_merged.append(df.loc[df.integrated_binnum == integrated_binnum])

		df_merged = pd.concat(df_merged)
		df_merged = df_merged.reset_index()
		df_merged = df_merged.loc[:, df_merged.columns[1:]]
		print("Saving /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}".format(sig_directory, suffix ))
		df_merged.to_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}/{}.pkl".format(sig_directory, suffix, integrated_binnum))

	return

def main_fast(mode, sig_directory, gen_directory):

	suffix = schema_suffices[mode]

	bulk_directory = "/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}".format(sig_directory, suffix)
	fringe_directory = "/volatile/clas12/sangbaek/dvcs_related/{}_fringe/excl_level_1/pkl_{}".format(sig_directory, suffix)

	if not os.path.exists(bulk_directory):
		return

	print(suffix)

	df_merged = {}
	for integrated_binnum in range(1, 147+1):
		df_merged[integrated_binnum] = []
	for integrated_binnum_gen in range(1, 147+159+1):
		if integrated_binnum_gen <= 147:
			print ("Reading bulk {}".format(integrated_binnum_gen))
			try:
				df = pd.read_pickle("{}/{}".format(bulk_directory, integrated_binnum_gen))
			except:
				print("{}/{} does not exist.".format(bulk_directory, integrated_binnum_gen))
				continue
		else:
			print ("Reading fringe {}".format(integrated_binnum_gen-147))
			try:
				df = pd.read_pickle("{}/{}".format(fringe_directory, integrated_binnum_gen-147))
			except:
				print("{}/{} does not exist.".format(fringe_directory, integrated_binnum_gen-147))
				continue

		df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
		for phi_binnum_gen in range(24):
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == gen_directory) & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
				* df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen
		for integrated_binnum in range(1,147+1):
			df_this_bin  = df.loc[(df.integrated_binnum == integrated_binnum), :]
			df_merged[integrated_binnum].append(df_this_bin)
	print("Saving /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}".format(sig_directory, suffix ))
	for integrated_binnum in range(1, 147+1):
		print("Integrated bin number {}".format(integrated_binnum))
		df_merged[integrated_binnum] = pd.concat(df_merged[integrated_binnum])
		df_merged[integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}/{}.pkl".format(sig_directory, suffix, integrated_binnum))
		df_merged[integrated_binnum] = 0
	return

def main_fast_3(mode, sig_directory, gen_directory):

	suffix = schema_suffices[mode]

	bulk_directory = "/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_3_{}".format(sig_directory, suffix)
	# fringe_directory = "/volatile/clas12/sangbaek/dvcs_related/{}_fringe/excl_level_1/pkl_3_{}".format(sig_directory, suffix)
	fringe_directory = "/volatile/clas12/sangbaek/dvcs_related/{}_fringe/excl_level_1/pkl_{}".format(sig_directory, suffix)

	if not os.path.exists(bulk_directory):
		return

	print(suffix)

	df_merged = []
	for integrated_binnum_gen in range(1, 147+159+1):
		if integrated_binnum_gen <= 147:
			print ("Reading bulk {}".format(integrated_binnum_gen))
			try:
				df = pd.read_pickle("{}/{}".format(bulk_directory, integrated_binnum_gen))
			except:
				print("{}/{} does not exist.".format(bulk_directory, integrated_binnum_gen))
				continue
		else:
			print ("Reading fringe {}".format(integrated_binnum_gen-147))
			gen_directory = "pureBH/fall2018_inb"
			try:
				df = pd.read_pickle("{}/{}".format(fringe_directory, integrated_binnum_gen-147))
			except:
				print("{}/{} does not exist.".format(fringe_directory, integrated_binnum_gen-147))
				continue

		df.loc[:, "integrated_binnum_gen"] = integrated_binnum_gen
		for phi_binnum_gen in range(24):
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "n_gen"] = df_gen.loc[(df_gen.directory == gen_directory) & (df_gen.integrated_binnum_gen == integrated_binnum_gen) & (df_gen.phi_binnum_gen == phi_binnum_gen)].n_entry.to_numpy()[0]
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] = bin_volume[integrated_binnum_gen]/24.
			df.loc[df.phi_binnum_gen == phi_binnum_gen, "weights"] = df.loc[df.phi_binnum_gen == phi_binnum_gen].GenWeight \
				* df.loc[df.phi_binnum_gen == phi_binnum_gen, "bin_volume_gen"] * luminosity_inb / df.loc[df.phi_binnum_gen == phi_binnum_gen].n_gen

		df_merged.append(df)
	df_merged = pd.concat(df_merged)
	df_merged = df_merged.reset_index()
	df_merged = df_merged.loc[:, df_merged.columns[1:]]
	print("Saving /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_3_{}".format(sig_directory, suffix ))
	for integrated_binnum in range(1, 147+1):
		print("Integrated bin number {}".format(integrated_binnum))
		df_merged.loc[df_merged.integrated_binnum == integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_3_{}/{}.pkl".format(sig_directory, suffix, integrated_binnum))

	return

def main(mode, sig_directory, gen_directory):
	# if mode == 0:
	# 	main_memory_efficient(mode, sig_directory, gen_directory)
	# else:
	main_fast(mode, sig_directory, gen_directory)


if __name__ == "__main__":

	for mode in range(1, 14):
		main(mode, "sim_rad_rec_fall2018_inb/dvcs_km15", "dvcs_km15/fall2018_inb3")
		main(mode, "sim_rad_rec_fall2018_outb/dvcs_km15", "dvcs_km15/fall2018_outb3")
		main(mode, "sim_rad_rec_fall2018_inb/dvcs_vgg", "dvcs_vgg/fall2018_inb")
		main(mode, "sim_rad_rec_fall2018_outb/dvcs_vgg", "dvcs_vgg/fall2018_outb")
		# main(mode, "sim_rad_rec_fall2018_inb/pureBH", "pureBH/fall2018_inb")
		main_fast_3(mode, "sim_rad_rec_fall2018_inb/pureBH", "pureBH/fall2018_inb3")
		main(mode, "sim_rad_rec_fall2018_outb/pureBH", "pureBH/fall2018_outb")
		# break

	# main(13, "sim_rad_rec_fall2018_inb/dvcs_km15", "dvcs_km15/fall2018_inb3_45nA")
	# main(13, "sim_rad_rec_fall2018_outb/dvcs_km15", "dvcs_km15/fall2018_outb3_50nA")

