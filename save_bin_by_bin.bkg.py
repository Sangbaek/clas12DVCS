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


def main_memory_efficient(mode, pi0_directory):

	suffix = schema_suffices[mode]

	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}".format(pi0_directory, suffix)):
		return

	print(suffix)

	for chunk in [1, 2, 3]:
		print("Reading /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}/{}".format(pi0_directory, suffix, chunk))
		for integrated_binnum in range(1, 147+1):
			print("Integrated bin number {}".format(integrated_binnum))
			df_merged = []
			for filenum in range(20):
				print(chunk, filenum+1)
				df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}/{}/{}".format(pi0_directory, suffix, chunk, filenum+1))
				df_merged.append(df.loc[df.integrated_binnum == integrated_binnum])
			df_merged = pd.concat(df_merged)
			df_merged = df_merged.reset_index()
			df_merged = df_merged.loc[:, df_merged.columns[1:]]
			print("Saving /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}/{}".format(pi0_directory, suffix, chunk))
			df_merged.to_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}/{}/{}.pkl".format(pi0_directory, suffix, chunk, integrated_binnum))


	return

def main_fast(mode, pi0_directory):

	suffix = schema_suffices[mode]

	if not os.path.exists("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}".format(pi0_directory, suffix)):
		return

	print(suffix)

	for chunk in [1, 2, 3]:
		print("Reading /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}/{}".format(pi0_directory, suffix, chunk))
		df_merged = []
		for filenum in range(20):
			print(chunk, filenum+1)
			df = pd.read_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/pkl_{}/{}/{}".format(pi0_directory, suffix, chunk, filenum+1))
			df_merged.append(df)
		df_merged = pd.concat(df_merged)
		df_merged = df_merged.reset_index()
		df_merged = df_merged.loc[:, df_merged.columns[1:]]
		print("Saving /volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}/{}".format(pi0_directory, suffix, chunk))
		for integrated_binnum in range(1, 147+1):
			print("Integrated bin number {}".format(integrated_binnum))
			df_merged.loc[df_merged.integrated_binnum == integrated_binnum].to_pickle("/volatile/clas12/sangbaek/dvcs_related/{}/excl_level_1/restructured_{}/{}/{}.pkl".format(pi0_directory, suffix, chunk, integrated_binnum))


	return

def main(mode, pi0_directory):
	if mode < 2:
		main_memory_efficient(mode, pi0_directory)
	else:
		main_fast(mode, pi0_directory)


if __name__ == "__main__":

	for mode in [11, 12]:#range(1, 14):
		main(mode, "sim_rad_rec_fall2018_inb/pi0_1gamma")
		main(mode, "sim_rad_rec_fall2018_inb/pi0_2gamma")
		main(mode, "sim_rad_rec_fall2018_outb/pi0_1gamma")
		main(mode, "sim_rad_rec_fall2018_outb/pi0_2gamma")
